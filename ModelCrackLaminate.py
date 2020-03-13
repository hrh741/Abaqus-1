#! /usr/python2
# _*_ coding: utf-8 _*_

"""
本程序用于生成带有端部裂纹(edge delaminated)的层合板结构
"""
from __future__ import print_function

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
    添加当前模块的路径
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
from Abaqus.Materials import Materials,UMATMaterial

ACISEPS=1e-6

def modelCrackLaminate(model1,ls,ws,plies,halfStructure,EquationOrDisplacement,ex=0.1,cracks=[]):
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
        position,crackLen,cracktype=crack
        if position>=plynum-1 or crackLen<=0 :
            # 裂纹面在中间 或者 不包含裂纹
            continue
        if cracktype==0:
            "crack 从y=0开始"
            offset=crackLen
            boundingBox={
                "yMin":0,"yMax":offset,
                "zMin":plyLowZs[position]-ACISEPS,"zMax":plyLowZs[position]+ACISEPS}
        elif cracktype==1:
            "crack 到y=lam_w结束"
            offset=lam_wid-crackLen
            boundingBox={
                "yMin":offset,"yMax":lam_wid,
                "zMin":plyLowZs[position]-ACISEPS,"zMax":plyLowZs[position]+ACISEPS}
        else:
            raise NotImplementedError("not support crack type , only can be 0 & 1 now!")
        print(crack,offset)
        dp_feature=rasm.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=offset)
        try:
            rasm.PartitionCellByDatumPlane(datumPlane=feature2datum(dp_feature,rasm), cells=myInstance1.cells)
        except:
            pass
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
    
    crackDict={crack[0]:crack[1] for crack in cracks}
    if (plynum-1) in crackDict:
        # 包含中面裂纹
        crackLen=crackDict[plynum-1]
        boundingBox={"zMax":0.0,
            "yMin":(0.0 if halfStructure[1]==1 else crackLen)-ACISEPS,
            "yMax":(lam_wid-crackLen if halfStructure[2]==1 else lam_wid)+ACISEPS,}
        Z0=rasm.Set(faces=myInstance1.faces.getByBoundingBox(**boundingBox), name='Z0')
    else :
        Z0=rasm.Set(faces=myInstance1.faces.getByBoundingBox(zMax=0.0), name='Z0')
    Z1=rasm.Set(faces=myInstance1.faces.getByBoundingBox(zMin=lam_hgt), name='Z1')
    X0Y0=rasm.Set(edges=myInstance1.edges.getByBoundingBox(xMax=0.0,yMax=0.0), name='X0Y0')
    X1Y0=rasm.Set(edges=myInstance1.edges.getByBoundingBox(xMin=lam_len,yMax=0.0), name='X1Y0')
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
    if halfStructure[1]==1:
        #model1.YsymmBC(name='Y=0-YSYMM', createStepName='Initial', region=Y0, localCsys=None)
        
        ### 这种有问题，应力在中间有波动
        #model1.DisplacementBC(name='X1Y0-Symm', 
        #    createStepName='Step-1', region=X1Y0, u1=0.0, u2=0.0, u3=UNSET, ur1=0.0, 
        #    ur2=UNSET, ur3=UNSET, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
        #    fieldName='', localCsys=None)
        
        ### 完美
        model1.DisplacementBC(name='Y0-Symm', 
            createStepName='Step-1', region=Y0, u1=UNSET, u2=0.0, u3=UNSET, ur1=0.0, 
            ur2=UNSET, ur3=UNSET, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
            fieldName='', localCsys=None)
        model1.DisplacementBC(name='X0Y0-Symm', 
            createStepName='Step-1', region=X0Y0, u1=0.0, u2=UNSET, u3=UNSET, ur1=UNSET, 
            ur2=UNSET, ur3=UNSET, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
            fieldName='', localCsys=None)
        model1.DisplacementBC(name='X1Y0-Symm', 
            createStepName='Step-1', region=X1Y0, u1=ex*lam_len, u2=UNSET, u3=UNSET, ur1=UNSET, 
            ur2=UNSET, ur3=UNSET, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
            fieldName='', localCsys=None)
        notEqNodeLabel.update([n.label for n in Y0.nodes])
        print(len(notEqNodeLabel))
    
    if halfStructure[2]==1:
        model1.ZsymmBC(name='Z=0-ZSYMM', createStepName='Initial', region=Z0, localCsys=None)
        ## 不要不要将Z0放入不添加约束的集合里，U1,U2的约束没有的话计算结果会出现错误
        #notEqNodeLabel.update([n.label for n in Z0.nodes])
        #print(len(notEqNodeLabel))

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
    model1.HistoryOutputRequest(name='H-Output-Ref',
        createStepName='Step-1', variables=('U', 'RF'), region=xref,
        sectionPoints=DEFAULT, rebar=EXCLUDE)

    if EquationOrDisplacement==0:
        #model1.Equation(name='Constraint-1', terms=((1.0, 'x-load',1), (-1.0, 'x-ref', 1)))
        model1.TabularAmplitude(name='disp-amp', timeSpan=TOTAL,
            smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (1.0, 1.0)))
        model1.DisplacementBC(name='X=1-Tension', createStepName='Step-1',
            region=xref, u1=ex*lam_len, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET,
            ur3=UNSET, amplitude='disp-amp', fixed=OFF, distributionType=UNIFORM,
            localCsys=None)
        model1.historyOutputRequests['H-Output-1'].setValues(
            frequency=1)
    else:
        ## 注意不要将边界条件都约束住 U2,U3 否则 会发现中间区域的应力与经典层合板理论不符合
        #model1.EncastreBC(name='X=0-ENCAST',createStepName='Initial',region=X0, localCsys=None)
        #model1.DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        #    distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='X=1-TENSION',
        #    region=X1, u1=ex*lam_len, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

        """
        model1.DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='X=0-TENSION',
            region=X0, u1=0.0, u2=0, u3=0, ur1=UNSET, ur2=UNSET, ur3=UNSET)
        model1.DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='X=1-TENSION',
            region=X1, u1=ex*lam_len, u2=0, u3=0, ur1=UNSET, ur2=UNSET, ur3=UNSET)
        """
        model1.DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='X=0-TENSION',
            region=X0, u1=0.0, u2=0, u3=0, ur1=UNSET, ur2=UNSET, ur3=UNSET)
        model1.DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='X=1-TENSION',
            region=xref, u1=ex*lam_len, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
        model1.Equation(name='Constraint-1', terms=((1.0, 'X1',1), (-1.0, 'x-ref', 1)))
    
    ## Job    
    descinfo='%s\r\n %s\r\n %s\r\n'%(str((lam_len,lam_wid,lam_hgt)),",".join([str(x) for x in (xdiv,ydiv)]),str(plyAngles))
    suffix="".join(("XYZ"[i] if halfStructure[i]==1 else "" for i in xrange(3)))
    suffix=suffix+("" if len(suffix)==0 else "SYMM")
    #jobname="Job-%s-%s-%sX2"%(plyprefix,suffix,'Eq' if EquationOrDisplacement==0 else 'Disp')
    jobname="Job-%s"%model1.name
    job=mdb.Job(name=jobname,
        description=descinfo,model=model1,numCpus=1,numGPUs=0)
    
    if EquationOrDisplacement==0:
        job.writeInput(consistencyChecking=OFF)
        del job
        

        ## 插入Equation 约束
        import time
        #print(time.time())
        SetINP,EquationINP=Pipes_Pagano_INP(model1,notInclude=notEqNodeLabel)
        #print(time.time())
        
        fname="./%s_PBC.inp"%jobname
        src=open("./%s.inp"%jobname,"r")
        dst=open(fname,"w+")

        isEquationInserted=False
        isSetInserted=False
        instanceCnt=0
        for line in src.readlines():
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
            inputFileName=fname, type=ANALYSIS, atTime=None, 
            waitMinutes=0, waitHours=0, queue=None, memory=90, memoryUnits=PERCENTAGE, 
            getMemoryFromAnalysis=True, explicitPrecision=SINGLE, 
            nodalOutputPrecision=SINGLE, scratch='', 
            resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, numDomains=1, 
            activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=1)
    return job

if __name__=="__main__":
    
    ## 输入参数
    ls=(180,20,5.0) #  层合板长度(x方向) , 单元数目，以及double seed ratio (end/center)
    ws=(25.0,50,1.0)  #  层合板宽度(y方向) , 单元数目，以及double seed ratio (center/end)
    plies=[ (45,    0.125,    "T300-7901",  8,    1.0),
            (0,     0.001,   "Epoxy7901",  1,    1.0),
            (0,     0.001,   "Epoxy7901",  1,    1.0),
            (-45,   0.125,   "T300-7901",   8,    1.0),
            (-45,   0.125,    "T300-7901",  8,    1.0),
            (0,     0.001,   "Epoxy7901",  1,    1.0),
            (0,     0.001,   "Epoxy7901",  1,    1.0),
            (45,    0.125,    "T300-7901",  8,    1.0),] # 从顶部到底部的每层的角度、厚度、材料以及厚度
    
    plies=[ (0,     0.135,    "T300-7901",  8,    1.0),
            (0,     0.001,   "Epoxy7901",  1,    1.0),
            (0,     0.001,   "Epoxy7901",  1,    1.0),
            (30,    0.135,    "T300-7901",  8,    1.0),
            (0,     0.001,   "Epoxy7901",  1,    1.0),
            (0,     0.001,   "Epoxy7901",  1,    1.0),
            (-30,   0.135,    "T300-7901",  8,    1.0),
            (0,     0.001,   "Epoxy7901",  1,    1.0),
            (0,     0.001,   "Epoxy7901",  1,    1.0),
            (60,    0.135,    "T300-7901",  8,    1.0),
            (0,     0.001,   "Epoxy7901",  1,    1.0),
            (0,     0.001,   "Epoxy7901",  1,    1.0),
            (-60,   0.135,    "T300-7901",  8,    1.0),
            (0,     0.001,   "Epoxy7901",  1,    1.0),
            (0,     0.001,   "Epoxy7901",  1,    1.0),
            (90,    0.135,    "T300-7901",  8,    1.0),
            (90,    0.135,    "T300-7901",  8,    1.0),
            (0,     0.001,   "Epoxy7901",  1,    1.0),
            (0,     0.001,   "Epoxy7901",  1,    1.0),
            (-60,    0.135,    "T300-7901",  8,    1.0),
            (0,     0.001,   "Epoxy7901",  1,    1.0),
            (0,     0.001,   "Epoxy7901",  1,    1.0),
            (60,    0.135,    "T300-7901",  8,    1.0),
            (0,     0.001,   "Epoxy7901",  1,    1.0),
            (0,     0.001,   "Epoxy7901",  1,    1.0),
            (-30,    0.135,    "T300-7901",  8,    1.0),
            (0,     0.001,   "Epoxy7901",  1,    1.0),
            (0,     0.001,   "Epoxy7901",  1,    1.0),
            (30,    0.135,    "T300-7901",  8,    1.0),
            (0,     0.001,   "Epoxy7901",  1,    1.0),
            (0,     0.001,   "Epoxy7901",  1,    1.0),
            (0,    0.135,    "T300-7901",  8,    1.0),] # 从顶部到底部的每层的角度、厚度、材料以及厚度

    plies=[ (0,     0.135,    "T300-7901",  8,    1.0),
            (30,    0.135,    "T300-7901",  8,    1.0),
            (-30,   0.135,    "T300-7901",  8,    1.0),
            (60,    0.135,    "T300-7901",  8,    1.0),
            (-60,   0.135,    "T300-7901",  8,    1.0),
            (90,    0.135,    "T300-7901",  8,    1.0),
            (90,    0.135,    "T300-7901",  8,    1.0),
            (-60,    0.135,    "T300-7901",  8,    1.0),
            (60,    0.135,    "T300-7901",  8,    1.0),
            (-30,    0.135,    "T300-7901",  8,    1.0),
            (30,    0.135,    "T300-7901",  8,    1.0),
            (0,    0.135,    "T300-7901",  8,    1.0),] # 从顶部到底部的每层的角度、厚度、材料以及厚度

    plies=[
        (30,    0.135,  "T300-648", 5, 1.0),
        (-30,   0.135,  "T300-648", 5, 1.0),
        (30,    0.135,  "T300-648", 5, 1.0),
        (-30,   0.135,  "T300-648", 5, 1.0),
        (30,    0.135,  "T300-648", 5, 1.0),
        (-30,   0.135,  "T300-648", 5, 1.0),
        (90,    0.135,  "T300-648", 5, 1.0),
        (90,    0.135,  "T300-648", 5, 1.0),
        (90,    0.135,  "T300-648", 5, 1.0),
        (90,    0.135,  "T300-648", 5, 1.0),
        (-30,   0.135,  "T300-648", 5, 1.0),
        (30,    0.135,  "T300-648", 5, 1.0),
        (-30,   0.135,  "T300-648", 5, 1.0),
        (30,    0.135,  "T300-648", 5, 1.0),
        (-30,   0.135,  "T300-648", 5, 1.0),
        (30,    0.135,  "T300-648", 5, 1.0),
    ]
    halfStructure=[0,0,0] # 每个tuple分别包含是否在x,y,z方向上是否使用半结构 0否 1是
    EquationOrDisplacement=1 #在拉伸边界上使用 参考点Equation约束(0) 或者 位移边界(1)

    ind=10
    modelNameFormat='YeLin-%d'
    for plyAngles in [[30,-30,30,-30,30,-30,90,90,90,90,-30,30,-30,30,-30,30],
                [45,-45,45,-45,0,0,90,90,90,90,0,0,-45,45,-45,45]]:
        plyAngles=np.array(plyAngles)
        plyThicks=0.135*np.ones_like(plyAngles)
        for mat in ['T300-648','T300-634']:
            modelName=modelNameFormat%ind
            ind=ind+1

            plies=[(plyAngles[i],plyThicks[i],mat,5,1.0) for i in range(len(plyAngles))]

            cracks=[(5,5.0,1),(5,5.0,0),(9,5.0,1),(9,5.0,0),]

            model1=mdb.Model(name=modelName)
            ## Material
            Materials(model1)
            modelCrackLaminate(model1,ls,ws,plies,halfStructure,EquationOrDisplacement,ex=0.01,cracks=cracks)

    """
    for ind,crack_len in enumerate([0.5*x for x in range(22,23)]):
        cracks=[(len(plies)/2-1,crack_len,1)] #中面裂纹
        modelName=modelnameFormat%ind

        model1=mdb.Model(name=modelName)
        ## Material
        Materials(model1)

        job=modelCrackLaminate(model1,ls,ws,plies,halfStructure,EquationOrDisplacement,ex=0.01,cracks=cracks)
    """
    #job.submit()
    #job.waitForCompletion()
    """
DelamResin=dict()
for ind,crack_len in enumerate([0.5*x for x in range(15)]):
    modelName=modelnameFormat%ind
    o=openOdb('Job-%s.odb'%modelName,readOnly=True)
    model1=mdb.models[modelName]
    sfo=o.steps.values()[-1].frames[-1].fieldOutputs['S']
    es=model1.rootAssembly.instances['Part-1-1'].elements.getByBoundingBox(yMax=ws[0]/2.0-crack_len+ACISEPS,yMin=ws[0]/2.0-crack_len-ws[0]/ws[1]-ACISEPS)
    elementCoord={e.label:lipeng.elementBound(e) for e in es}
    sortedElm=sorted(elementCoord.keys(),key=lambda label:elementCoord[label][0][2])
    botElmLabel=sortedElm[0]
    d={e.label:None for e in es}
    cnt=0
    RF1=o.steps.values()[-1].historyRegions['Node PART-2-1.1'].historyOutputs['RF1'].data[1][1]
    for value in sfo.values:
        if (value.elementLabel in d ) and d[value.elementLabel] is None:
            d[value.elementLabel]=(elementCoord[value.elementLabel],
                    value.data,value.localCoordSystem,
                    )
            cnt=cnt+1
            if(cnt==len(d)):
                break
    
    Delam[crack_len]=(RF1,d)
    print(modelName,crack_len,RF1,d[botElmLabel][1]#,elementCoord[botElmLabel],)
    """