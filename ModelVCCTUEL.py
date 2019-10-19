#! /usr/python2
# _*_ coding: utf-8 _*_

"""
本程序用于修改使用Seam建模后的层合板模型的inp文件
在裂纹尖端插入参照以下引文的VCCT UEL
>Xie D, Biggers Jr S B. Progressive crack growth analysis using interface element based on the virtual crack closure technique[J]. Finite Elements in Analysis and Design, 2006, 42(11): 977-984.
具体做法为
将裂纹尖端节点分为两个节点，然后创建UEL单元

实际实现的时候我可以这么来干:
  把Seam的crack 尺寸设置的比真实的裂纹尺寸大一个单元宽度
  那么输出的inp文件就会将我们想的裂尖节点分离成两个无关的节点
  这时候我们找到这两个节点 建立UEL单元
"""
import numpy as np
from abaqus import *
from abaqusConstants import *
import regionToolset
import mesh
import visualization
import xyPlot
from odbAccess import *

if __name__ == '__main__' and __package__ is None:
    """

    """
    import inspect
    from os import sys, path
    #curpath=path.dirname(path.dirname(path.abspath(__file__)))
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    curpath = os.path.dirname(os.path.abspath(filename))
    pth=os.path.dirname(curpath)
    if pth not in sys.path:
        sys.path.append(pth)

from Abaqus.lipeng import edge2vector,feature2datum,labels2sequence,Pipes_Pagano_INP

ACISEPS=1e-6

def modelVCCTUEL(model1,ls,ws,plies,halfStructure,EquationOrDisplacement,ex=0.1,cracks=[]):
    """
    ls: tuple (层合板长度,长度方向的单元数目,double seed 的长度比例)
    ws: tuple (层合板宽度,宽度方向上的单元数目,double seed 的长度比例)
    """
    l,xdiv,xratio=ls
    w,ydiv,yratio=ws
    
    ## 根据半结构设置生成自定义参数 不要修改
    lam_len,lam_wid=l/(2.0 if halfStructure[0]==1 else 1),w/(2.0 if halfStructure[1]==1 else 1)
    
    plynum=len(plies)//(2 if halfStructure[2]==1 else 1)
    plyAngles=[plies[i][0] for i in xrange(plynum)]
    plyThicks=[plies[i][1] for i in xrange(plynum)]
    plyProps=[plies[i][2] for i in xrange(plynum)]

    lam_hgt=sum(plyThicks)
    plyprefix=",".join(["%02.0f"%x for x in plyAngles])
    plyLowZs=np.cumsum(plyThicks,dtype='f')
    plyLowZs=plyLowZs[-1]-plyLowZs

    print(plyLowZs)
    ## Part
    part1=model1.Part(name="Part-1",dimensionality=THREE_D,type=DEFORMABLE_BODY)
    s=model1.ConstrainedSketch(name="Sketch-1",sheetSize=200.0)
    s.rectangle(point1=(0,0),point2=(lam_len,lam_wid))
    part1.BaseSolidExtrude(sketch=s,depth=lam_hgt)

    for i in xrange(plynum-1):
        dt=part1.DatumPlaneByPrincipalPlane(offset=plyLowZs[i], principalPlane=XYPLANE)
        part1.PartitionCellByDatumPlane(cells=part1.cells,datumPlane=part1.datums[dt.id])
    part1.regenerate()                  

    ## Composite Laminate 
    Dt=part1.DatumCsysByThreePoints(name='Datum csys-1', coordSysType=CARTESIAN, origin=(
      0.0, 0.0, 0.0), line1=(1.0, 0.0, 0.0), line2=(0.0, 1.0, 0.0))
    layupOrientation = part1.datums[Dt.id]
    compositeLayup = part1.CompositeLayup(
      name="Laminate", description='/'.join(( str(x) for x in plyAngles)), elementType=SOLID, 
      symmetric=False, thicknessAssignment=FROM_SECTION)
    compositeLayup.ReferenceOrientation(orientationType=SYSTEM, 
    localCsys=layupOrientation, fieldName='', 
    additionalRotationType=ROTATION_NONE, angle=0.0, 
    additionalRotationField='', axis=AXIS_3, stackDirection=STACK_3)

    print("part1.cells:",len(part1.cells))
    plysetkeys=[]
    for i in xrange(plynum):
        print("ply ",i," start")
        z=plyLowZs[i]+plyThicks[i]/2.0
        r=part1.Set(name="Ply-%1d"%(i),cells=part1.cells.findAt(((0,0,z),)) )
        #r=regionToolset.Region(cells=part1.cells.findAt(((0,0,z),)))
        plyname='Ply-%1d-%3.1f'%(i,plyAngles[i])
        plysetkeys.append(plyname)
        print("ply ",i," 3","region:",r)
        compositeLayup.CompositePly(suppressed=False, plyName=plyname, region=r, 
          material=plyProps[i], thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
          orientationType=SPECIFY_ORIENT, orientationValue=plyAngles[i], 
          additionalRotationType=ROTATION_NONE, additionalRotationField='', 
          axis=AXIS_3, angle=0.0, numIntPoints=1) # must use 1 for integral point 
        print("ply ",i," end")

    ## Step
    model1.StaticStep(name='Step-1',noStop=OFF, previous='Initial', timeIncrementationMethod=AUTOMATIC, initialInc=1.0,maxNumInc=1000,minInc=1e-05)
    model1.fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'UR','RF', 'CF', 'CSTRESS', 'CDISP', 
    'SDV','NFORC'))

    ## Assembly
    rasm=model1.rootAssembly
    model1.rootAssembly.DatumCsysByDefault(CARTESIAN)
    myInstance1=model1.rootAssembly.Instance(dependent=OFF, name='Part-1-1', 
      part=part1)

    ## 创建Seam
    for crack_ind,crack in enumerate(cracks):
        position=crack[0]
        cracklen=crack[1]
        cracktype=crack[2]
        real_cracklen=cracklen+w/ydiv
        if cracktype==0:
            "crack 从y=0开始"
            offset=real_cracklen
            boundingBox={
                "yMin":0,"yMax":offset,
                "zMin":plyLowZs[position]-ACISEPS,"zMax":plyLowZs[position]+ACISEPS}
        elif cracktype==1:
            "crack 到y=lam_w结束"
            offset=w-real_cracklen
            boundingBox={
                "yMin":offset,"yMax":w,
                "zMin":plyLowZs[position]-ACISEPS,"zMax":plyLowZs[position]+ACISEPS}
        else:
            raise NotImplementedError("not support crack type , only can be 0 & 1 now!")
        dp_feature=rasm.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=offset)
        rasm.PartitionCellByDatumPlane(datumPlane=feature2datum(dp_feature,rasm), cells=myInstance1.cells)
        f1 = myInstance1.faces.getByBoundingBox(**boundingBox) 
        seam = rasm.Set(faces=f1, name='Seam_%d'%crack_ind)
        rasm.engineeringFeatures.assignSeam(regions=seam)

    ## Mesh
    rasm.setElementType(elemTypes=(mesh.ElemType(elemCode=C3D8R, elemLibrary=STANDARD, secondOrderAccuracy=OFF, kinematicSplit=AVERAGE_STRAIN, hourglassControl=ENHANCED, distortionControl=DEFAULT), 
                                    mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD), 
                                    mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)),
                        regions=(myInstance1.cells,))
    rasm.assignStackDirection(cells=myInstance1.cells, referenceRegion=myInstance1.faces.getByBoundingBox(zMin=lam_hgt)[0])

    xedges=myInstance1.edges.findAt(coordinates=map(lambda eg:eg.pointOn[0],
            filter(lambda eg:abs(abs(edge2vector(myInstance1.vertices,eg)[0])-lam_len)<ACISEPS,myInstance1.edges)))
    rasm.seedEdgeByBias(biasMethod=DOUBLE, centerEdges=xedges, constraint=FINER, number=xdiv, ratio=xratio)

    yedges=labels2sequence([e.index for e in filter(lambda eg:sum([abs(x) for x in edge2vector(myInstance1.vertices,eg)[0:3:2]])<ACISEPS,myInstance1.edges)],myInstance1.edges)
    #myInstance1.edges.findAt(coordinates=map(lambda eg:eg.pointOn[0],filter(lambda eg:sum([abs(x) for x in edge2vector(myInstance1.vertices,eg)[0:3:2]])<ACISEPS,myInstance1.edges)))
    #rasm.seedEdgeByBias(biasMethod=DOUBLE, constraint=FINER, endEdges=yedges, number=ydiv, ratio=yratio)
    rasm.seedEdgeBySize(edges=yedges, size=w/ydiv, deviationFactor=0.1, constraint=FINER)
    
    for i in range(plynum):
        zdiv=plies[i][3]
        zratio=plies[i][4]
        zedges=myInstance1.edges.findAt(coordinates=map(lambda eg:eg.pointOn[0],
            filter(lambda eg:sum([abs(x) for x in edge2vector(myInstance1.vertices,eg)[0:2]])<ACISEPS,
                myInstance1.edges.getByBoundingBox(zMin=plyLowZs[i]-ACISEPS,zMax=(plyLowZs[i-1] if i>0 else lam_hgt)+ACISEPS))))
        rasm.seedEdgeByBias(biasMethod=DOUBLE, constraint=FINER, endEdges=zedges, number=zdiv,ratio=zratio)
    rasm.generateMesh(regions=(myInstance1,))
    
    ## BC
    Origin = rasm.Set(vertices=myInstance1.vertices.findAt(((0,0,0),)), name='Origin')
    X0=rasm.Set(faces=myInstance1.faces.getByBoundingBox(xMax=0.0), name='X0')
    X1=rasm.Set(faces=myInstance1.faces.getByBoundingBox(xMin=lam_len), name='X1')
    Y0=rasm.Set(faces=myInstance1.faces.getByBoundingBox(yMax=0.0), name='Y0')
    Y1=rasm.Set(faces=myInstance1.faces.getByBoundingBox(yMin=lam_wid), name='Y1')
    Z0=rasm.Set(faces=myInstance1.faces.getByBoundingBox(zMax=0.0), name='Z0')
    Z1=rasm.Set(faces=myInstance1.faces.getByBoundingBox(zMin=lam_hgt), name='Z1')
    model1.rootAssembly.regenerate()

    notEqNodeLabel=set()
    model1.PinnedBC(name='X0Y0Z0-Pin', createStepName='Initial', region=Origin, localCsys=None)
    notEqNodeLabel.update([n.label for n in Origin.nodes])
    if halfStructure[0]==1:
        model1.XsymmBC(name='X=0-XSYMM', createStepName='Initial', region=X0, localCsys=None)    
        notEqNodeLabel.update([n.label for n in X0.nodes])
        print(len(notEqNodeLabel))
    else:
        pass
        #model1.EncastreBC(name='X=0-ENCAST',createStepName='Initial',region=X0, localCsys=None)
    if halfStructure[1]==1:
        #model1.YsymmBC(name='Y=0-YSYMM', createStepName='Initial', region=Y0, localCsys=None)
        model1.DisplacementBC(name='Y=0-LAMINATE_SYMM', createStepName='Initial', 
            region=Y0, u1=UNSET, u2=SET, u3=UNSET, ur1=SET, ur2=SET, ur3=UNSET, 
            amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
        notEqNodeLabel.update([n.label for n in Y0.nodes])
        print(len(notEqNodeLabel))
    
    if halfStructure[2]==1:
        model1.ZsymmBC(name='Z=0-ZSYMM', createStepName='Initial', region=Z0, localCsys=None)
        notEqNodeLabel.update([n.label for n in Z0.nodes])
        print(len(notEqNodeLabel))
    
    if EquationOrDisplacement==0:
        part2 = model1.Part(name='Part-2', dimensionality=THREE_D,
            type=DEFORMABLE_BODY)
        part2.ReferencePoint(point=(lam_len+1, 0.0, 0.0))
        myInstance2=model1.rootAssembly.Instance(dependent=ON, name='Part-2-1', 
          part=part2)
        faces=myInstance1.faces.getByBoundingBox(xMin=lam_len-ACISEPS,xMax=lam_len+ACISEPS)
        rasm.Set(faces=faces,name='x-load')
        r1 = myInstance2.referencePoints.findAt((lam_len+1.0,0,0))
        refPoints1=(r1,)
        xref=rasm.Set(referencePoints = refPoints1, name='x-ref')
        model1.Equation(name='Constraint-1', terms=((1.0, 'x-load',
            1), (-1.0, 'x-ref', 1)))
        model1.TabularAmplitude(name='disp-amp', timeSpan=TOTAL,
            smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (1.0, 1.0)))
        model1.DisplacementBC(name='X=1-Tension', createStepName='Step-1',
            region=xref, u1=ex*lam_len, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET,
            ur3=UNSET, amplitude='disp-amp', fixed=OFF, distributionType=UNIFORM,
            localCsys=None)
        model1.historyOutputRequests['H-Output-1'].setValues(
            frequency=10)
        model1.HistoryOutputRequest(name='H-Output-Ref',
            createStepName='Step-1', variables=('U', 'RF'), region=xref,
            sectionPoints=DEFAULT, rebar=EXCLUDE)
    else:
        model1.DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='X=1-TENSION',
            region=X1, u1=ex*lam_len, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

    #model1.boundaryConditions['X=0-ENCAST'].suppress()
    model1.constraints['Constraint-1'].suppress()
    
    ## Job    
    descinfo='%s\r\n %s\r\n %s\r\n'%(str((lam_len,lam_wid,lam_hgt)),",".join([str(x) for x in (xdiv,ydiv)]),str(plyAngles))
    suffix="".join(("XYZ"[i] if halfStructure[i]==1 else "" for i in xrange(3)))
    suffix=suffix+("" if len(suffix)==0 else "SYMM")
    #jobname="Job-%s-%s-%sX2"%(plyprefix,suffix,'Eq' if EquationOrDisplacement==0 else 'Disp')
    jobname="Job-%s"%model1.name
    job=mdb.Job(name=jobname,
        description=descinfo,model=model1,numCpus=1,numGPUs=0)

    job.writeInput(consistencyChecking=OFF)
    del job
    
    ## 插入Equation 约束
    import time
    #print(time.time())
    SetINP,EquationINP=Pipes_Pagano_INP(model1,notInclude=notEqNodeLabel)
    #print(time.time())

    ## 插入UEL
    import operator
    UELNodeOrders=[]
    for crack_ind,crack in enumerate(cracks):
        position=crack[0]
        cracklen=crack[1]
        cracktype=crack[2]
        real_cracklen=cracklen+w/ydiv
        zmin,zmax=plyLowZs[position]-ACISEPS,plyLowZs[position]+ACISEPS
        if cracktype==0:
            "crack 从y=0开始"
            ymin,ymax=cracklen-w/ydiv-ACISEPS,cracklen+ACISEPS
            ns=myInstance1.nodes.getByBoundingBox(zMin=zmin,zMax=zmax,yMin=ymin,yMax=ymax)
            assert len(ns)==8,"crack tip bounding box error"
            ymid=(ymin+ymax)/2.0
            xmid=lam_len/2.0
            nslabel=set((n.label for n in ns))
            e1=myInstance1.elements.getByBoundingBox(zMin=zmin,zMax=zmin+plyThicks[position]+2*ACISEPS,yMin=ymin,yMax=ymax)
            blankList=[]
            for e in e1:
                blankList.extend([myInstance1.nodes[nodeind].label for nodeind in e.connectivity])
            ns5678=[myInstance1.nodes[l-1] for l in (nslabel & set(blankList))]

            e2=myInstance1.elements.getByBoundingBox(zMax=zmax,zMin=plyLowZs[position+1]-ACISEPS,yMin=ymin,yMax=ymax)
            blankList2=[]
            for e in e2:
                blankList2.extend([myInstance1.nodes[nodeind].label for nodeind in e.connectivity])
            ns1234=[myInstance1.nodes[l-1] for l in (nslabel & set(blankList2))]
            
            nodeorder=[0 for _ in range(8)]
            i=0
            for xf in [operator.lt,operator.gt]:
                for yf in [operator.lt,operator.gt]:
                    nodeorder[i]=filter(lambda n:xf(n.coordinates[0],xmid) and  yf(n.coordinates[1],ymid),ns1234)[0].label
                    nodeorder[i+4]=filter(lambda n:xf(n.coordinates[0],xmid) and  yf(n.coordinates[1],ymid),ns5678)[0].label
                    i=i+1
            UELNodeOrders.append(nodeorder)
        elif cracktype==1:
            "crack 到y=lam_w结束"
            ymin,ymax=w-cracklen-ACISEPS,w-cracklen+w/ydiv+ACISEPS
            ns=myInstance1.nodes.getByBoundingBox(zMin=zmin,zMax=zmax,yMin=ymin,yMax=ymax)
            assert len(ns)==8,"crack tip bounding box error" 
            ymid=(ymin+ymax)/2.0
            xmid=lam_len/2.0
            nslabel=set((n.label for n in ns))
            e1=myInstance1.elements.getByBoundingBox(zMin=zmin,zMax=zmin+plyThicks[position]+2*ACISEPS,yMin=ymin,yMax=ymax)
            blankList=[]
            for e in e1:
                blankList.extend([myInstance1.nodes[nodeind].label for nodeind in e.connectivity])
            ns5678=[myInstance1.nodes[l-1] for l in (nslabel & set(blankList))]

            e2=myInstance1.elements.getByBoundingBox(zMax=zmax,zMin=plyLowZs[position+1]-ACISEPS,yMin=ymin,yMax=ymax)
            blankList2=[]
            for e in e2:
                blankList2.extend([myInstance1.nodes[nodeind].label for nodeind in e.connectivity])
            ns1234=[myInstance1.nodes[l-1] for l in (nslabel & set(blankList2))]

            #print(xmid,ymid)
            #print(blankList,blankList2)
            #print([(n.label,n.coordinates) for n in ns1234],[(n.label,n.coordinates) for n in ns5678])
            nodeorder=[0 for _ in range(8)]
            i=0
            for xf in [operator.lt,operator.gt]:
                for yf in [operator.gt,operator.lt]:
                    nodeorder[i]=filter(lambda n:xf(n.coordinates[0],xmid) and  yf(n.coordinates[1],ymid),ns1234)[0].label
                    nodeorder[i+4]=filter(lambda n:xf(n.coordinates[0],xmid) and  yf(n.coordinates[1],ymid),ns5678)[0].label
                    i=i+1
            UELNodeOrders.append(nodeorder)
        else:
            raise NotImplementedError("not support crack type , only can be 0 & 1 now!")
        
    #print(UELNodeOrders)
    UelINP=""
    if len(cracks)>0:
        UelINP+="** UEL of Interface element for VCCT\n*USER ELEMENT, TYPE=U1, NODES=8, COORD=3, PROPERTIES=1, VARIABLES=3\n1,2,3\n"
        UelINP+="*ELEMENT, TYPE=U1, ELSET=Cracks\n"
        for i,_ in enumerate(cracks):
            UelINP+="%d, %s\n"%(len(myInstance1.elements)+1+i,",".join([str(label) for label in UELNodeOrders[i]]))
        UelINP+="*UEL PROPERTY, ELSET=Cracks\n1.E10\n"
    
    #print(UelINP)
    fnmae="./%s_Modified.inp"%jobname
    src=open("./%s.inp"%jobname,"r")
    dst=open(fnmae,"w+")

    isUelInserted=False
    isEquationInserted=False
    isSetInserted=False
    instanceCnt=0
    for line in src.readlines():
        if not isUelInserted and line.upper().startswith("*END INSTANCE"):
            dst.write(UelINP)
            isUelInserted=True
        if not isEquationInserted and line.upper().startswith("*END ASSEMBLY"):
            dst.write(EquationINP)
            isEquationInserted=True
        dst.write(line)
        if not isSetInserted and line.upper().startswith("*END INSTANCE"):
            if instanceCnt==1:
                dst.write(SetINP)
                isSetInserted=True
            instanceCnt=instanceCnt+1
    src.close()
    dst.close()

    job=mdb.JobFromInputFile(name='%s'%jobname, 
        inputFileName=fnmae, type=ANALYSIS, atTime=None, 
        waitMinutes=0, waitHours=0, queue=None, memory=90, memoryUnits=PERCENTAGE, 
        getMemoryFromAnalysis=True, explicitPrecision=SINGLE, 
        nodalOutputPrecision=SINGLE, 
        userSubroutine='E:/1730895/Work/NumricalMethods/VCCT/VCCT_UEL.for', scratch='', 
        resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, numDomains=1, 
        activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=1)
    return job

materials=[
    ## 吴庆欣
    ('T300-7901',(137.78e3,8.91e3,8.91e3,0.3,0.3,0.48,4.41e3,4.41e3,3.01e3),None,None), #
    ('Pipes-Pagano',(20.0, 2.1, 2.1, 0.21, 0.21, 0.21, 0.85, 0.85, 0.85),None,None), #MPsi
    ## Uguen A, Zubillaga L, Turon A, et al. Comparison of cohesive zone models used to predict delamination initiated from free-edges: validation against experimental results[C]//ECCM-16TH European Conference on Composite Materials. 2014.
    ### GIc,GIIc,GIIIc,sigma_zz0,tau_xz0,tau_yz0=0.24,0.74,0.74,46,75,75
    ('T800-M21',(130.0e3,8.0e3,8.0e3,0.31,0.31,0.45,4.0e3,4.0e3,4.0e3),None,None), #MPa
    ('S2/SP250 GlasdEpoxy',None,None,None),
    ## Initiation of free-edge delamination in composite laminates
    ('G947/M18',(97.6e3, 8.0e3 ,8.0e3,0.37,0.37,0.5, 3.1e3, 3.1e3, 2.7e3)),
    ('Material_A',(10e6,0.3),None,None,None),
    ('Material_B',(19.231e6,0.0),None,None,None),
    ('Material_B1',(30e6,0.0),None,None,None), 
    ## 
    ('Expoxy_7901',(3.17e3, 0.355),None,None)

]

if __name__=="__main__":
    import os
    #os.chdir("E:/UEL")

    for i in range(1):
        crack_len=2.0
        model1=mdb.Model(name="Resin-Layer-2-")
        
        ## Material
        for mat in materials:
            if mat[1]:
                material1=model1.Material(name=mat[0])
                if len(mat[1])==9:
                    material1.Elastic(type=ENGINEERING_CONSTANTS, table=(mat[1], ))
                elif len(mat[1])==2:
                    material1.Elastic(type=ISOTROPIC, table=(mat[1], ))
        
        ## 输入参数
        ls=(0.05,1,1.0) #  层合板长度(x方向) , 单元数目，以及double seed ratio (end/center)
        ws=(8.0,400,1.0)  #  层合板长度(y方向) , 单元数目，以及double seed ratio (center/end)
        plies=[ (45,    0.125,    "T300-7901",  8,    1.0),
                (0,     0.01,   "Expoxy_7901",  1,    1.0),
                (0,     0.01,   "Expoxy_7901",  1,    1.0),
                (-45,   0.125,   "T300-7901",   8,    1.0),
                (-45,   0.125,    "T300-7901",  8,    1.0),
                (0,     0.01,   "Expoxy_7901",  1,    1.0),
                (0,     0.01,   "Expoxy_7901",  1,    1.0),
                (45,    0.125,    "T300-7901",  8,    1.0),] # 从顶部到底部的每层的角度、厚度、材料以及厚度
        
        halfStructure=[0,0,1] # 每个tuple分别包含是否在x,y,z方向上是否使用半结构 0否 1是
        # 注意宽度方向上的对称条件不知道为什么一直错误 所以宽度方向尽量不要使用半结构
        EquationOrDisplacement=0 #在拉伸边界上使用 参考点Equation约束(0) 或者 位移边界(1)
        cracks=[(0,crack_len,0),(0,crack_len,1)]

        print 'Start Model VCCT'
        job=modelVCCTUEL(model1,ls,ws,plies,halfStructure,EquationOrDisplacement,ex=0.01,cracks=cracks)

    #job.submit()
    #job.waitForCompletion()

    #import glob
    #Gs=[]
    #dst=open("E:/result.txt","w")
    #for fn in glob.glob("./*.log"):
    #    with open(fn,"r") as fp:
    #        for line in fp.readlines():
    #            if line.find("GI,GII,GIII") !=-1:
    #                G=line
    #            if line.find("TIP")!=-1:
    #                coord=float(line.split()[3])
    #        dst.write("%s %f : %s"%(fn,coord,G))
    #        Gs.append((coord,[float(x) for x in G.split()[1:]]))
    #dst.close()