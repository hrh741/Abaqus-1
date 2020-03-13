#! /bin/python2
# _*_ coding: UTF-8 _*_

"""
建立带有Cohesive 单元的层合板模型
"""

import os
				  
from abaqus import *
from abaqusConstants import *
import regionToolset
import mesh
import numpy as np
import visualization
import xyPlot
from odbAccess import *

if __name__ == '__main__' and __package__ is None:
    import inspect
    from os import sys, path
    #curpath=path.dirname(path.dirname(path.abspath(__file__)))
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    curpath = os.path.dirname(os.path.abspath(filename))
    pth=os.path.dirname(curpath)
    if pth not in sys.path:
        sys.path.append(pth)

from Abaqus.lipeng import edge2vector,labels2sequence
from Abaqus import Materials,UMATMaterial
			

def modelCohesiveLaminate(model1,ls,ws,plies,halfStructure,EquationOrDisplacement,ex=0.1):
    """
    在model1 中建立每层用Cohesive单元连接的层合板
    model1:  
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

    plysetkeys=[]
    for i in xrange(plynum):
        if plyProps[i]=="Cohesive":
            continue
        z=plyLowZs[i]+plyThicks[i]/2.0
        part1.Set(name="Ply-%1d"%(i),cells=part1.cells.findAt(((0,0,z),)) )
        r=regionToolset.Region(cells=part1.cells.findAt(((0,0,z),)))
        plyname='Ply-%1d-%3.1f'%(i,plyAngles[i])
        plysetkeys.append(plyname)
										
        compositeLayup.CompositePly(suppressed=False, plyName=plyname, region=r, 
        material=plyProps[i], thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
        orientationType=SPECIFY_ORIENT, orientationValue=plyAngles[i], 
        additionalRotationType=ROTATION_NONE, additionalRotationField='', 
        axis=AXIS_3, angle=0.0, numIntPoints=1) # must use 1 for integral point 

    ## Cohesive Section Assignment
    model1.CohesiveSection(name='Cohesive', material='Cohesive', response=TRACTION_SEPARATION, outOfPlaneThickness=None)
    coh_cells=labels2sequence([part1.cells.findAt((0,0,plyLowZs[i]+plyThicks[i]/2)).index for i in range(1,plynum) if plies[i][2]=="Cohesive"],part1.cells)
    coh_region=regionToolset.Region(cells=coh_cells)
    part1.SectionAssignment(region=coh_region, sectionName='Cohesive', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
    
    composite_cells=labels2sequence([part1.cells.findAt((0,0,plyLowZs[i]+plyThicks[i]/2)).index for i in range(0,plynum) if plies[i][2]!="Cohesive"],part1.cells)
    composite_region=regionToolset.Region(cells=composite_cells)

    ## Assembly
    rasm=model1.rootAssembly
    model1.rootAssembly.DatumCsysByDefault(CARTESIAN)
    myInstance1=model1.rootAssembly.Instance(dependent=ON, name='Part-1-1', 
      part=part1)
		   

    ## Step (Field Output & History Outpu)    
    model1.StaticStep(name='Step-1', noStop=OFF, previous='Initial', timeIncrementationMethod=AUTOMATIC, initialInc=0.01,maxNumInc=1000,minInc=1e-08, maxInc=0.01,nlgeom=ON)
    
    model1.fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'E', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'UR', 'RF', 'CF','CSTRESS', 'CDISP', 'SDV', 'STATUS','NFORC'))
    coh_set=rasm.Set(name="Cohesive",cells=myInstance1.cells.findAt([(0,0,plyLowZs[i]+plyThicks[i]/2) for i in range(1,plynum) if plies[i][2]=="Cohesive"]))
    model1.FieldOutputRequest(name='F-Output-Coh',createStepName='Step-1', variables=('SDEG', 'NE'), region=coh_set,sectionPoints=DEFAULT, rebar=EXCLUDE)
    compsite_set=rasm.Set(name="Composite",
        cells=myInstance1.cells.findAt(*[((0,0,plyLowZs[i]+plyThicks[i]/2),) for i in range(0,plynum) if plies[i][2]!="Cohesive"]))
    model1.FieldOutputRequest(name='F-Output-Composite',createStepName='Step-1', variables=('S', 'E','U','SDV', 'STH', 'STATUS','NFORC'), region=compsite_set,sectionPoints=DEFAULT, rebar=EXCLUDE)
     
    ## Mesh
    part1.setElementType(regions=composite_region,elemTypes=(mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD), mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD),
        mesh.ElemType(elemCode=C3D8R, elemLibrary=STANDARD, secondOrderAccuracy=OFF, kinematicSplit=AVERAGE_STRAIN, hourglassControl=ENHANCED, distortionControl=DEFAULT)))
    part1.setElementType(regions=coh_region,
        elemTypes=(mesh.ElemType(elemCode=COH3D8, elemLibrary=STANDARD,viscosity=0.0001,elemDeletion=ON),
            mesh.ElemType(elemCode=COH3D6, elemLibrary=STANDARD),
            mesh.ElemType(elemCode=UNKNOWN_TET, elemLibrary=STANDARD)))
    part1.assignStackDirection(cells=part1.cells, referenceRegion=part1.faces.getByBoundingBox(zMin=lam_hgt)[0])

    xedges=part1.edges.findAt(coordinates=map(lambda eg:eg.pointOn[0],
          filter(lambda eg:abs(abs(edge2vector(part1.vertices,eg)[0])-lam_len)<1e-6,part1.edges)))
    part1.seedEdgeByBias(biasMethod=DOUBLE, centerEdges=xedges, constraint=FINER, number=xdiv, ratio=xratio)

    yedges=part1.edges.findAt(coordinates=map(lambda eg:eg.pointOn[0],
          filter(lambda eg:abs(abs(edge2vector(part1.vertices,eg)[1])-lam_wid)<1e-6,part1.edges)))
    part1.seedEdgeByBias(biasMethod=DOUBLE, constraint=FINER, endEdges=yedges, number=ydiv, ratio=yratio)
    
    for i in range(plynum):
        zdiv=plies[i][3]
        zratio=plies[i][4]
        zedges=part1.edges.findAt(coordinates=map(lambda eg:eg.pointOn[0],
            filter(lambda eg:sum([abs(x) for x in edge2vector(part1.vertices,eg)][:2])<1e-6,
                part1.edges.getByBoundingBox(zMin=plyLowZs[i],zMax=plyLowZs[i-1] if i>0 else lam_hgt))))
        part1.seedEdgeByBias(biasMethod=DOUBLE, constraint=FINER, endEdges=zedges, number=zdiv,ratio=zratio)
    part1.generateMesh()

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
        #notEqNodeLabel.update([n.label for n in Z0.nodes])
        ## 不要不要将Z0放入不添加约束的集合里，U1,U2的约束没有的话计算结果会出现错误
        #print(len(notEqNodeLabel))
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
        #model1.Equation(name='Constraint-1', terms=((1.0, 'x-load',1), (-1.0, 'x-ref', 1)))
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
        ## 注意不要将边界条件都约束住 U2,U3 否则 会发现中间区域的应力与经典层合板理论不符合
        #model1.EncastreBC(name='X=0-ENCAST',createStepName='Initial',region=X0, localCsys=None)
        #model1.DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        #    distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='X=1-TENSION',
        #    region=X1, u1=ex*lam_len, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

        model1.DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='X=0-TENSION',
            region=X0, u1=0.0, u2=0, u3=0, ur1=UNSET, ur2=UNSET, ur3=UNSET)
        model1.DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='X=1-TENSION',
            region=X1, u1=ex*lam_len, u2=0, u3=0, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    
    ## Job    
    descinfo='%s\r\n %s\r\n %s\r\n'%(str((lam_len,lam_wid,lam_hgt)),",".join(map(str,(xdiv,ydiv))),str(plyAngles))
    suffix="".join(("XYZ"[i] if halfStructure[i]==1 else "" for i in xrange(3)))
    suffix=suffix+("" if len(suffix)==0 else "SYMM")
    #jobname="Job-%s-%s-%s"%(plyprefix,suffix,'Eq' if EquationOrDisplacement==0 else 'Disp')
    jobname="Job-%s"%model1.name
    job=mdb.Job(name=jobname,
        description=descinfo,model=model1,numCpus=1,numGPUs=0)
    
    return model1,job,notEqNodeLabel


if __name__=="__main__":
    ## 输入参数
    ls=(0.05,1,1.0) #  层合板长度(x方向) , 单元数目，以及double seed ratio (end/center)
    ws=(8.0,400,4.0)  #  层合板宽度(y方向) , 单元数目，以及double seed ratio (center/end)
    plies=[ (45,    0.125,    "T800-M21", 8,	4.0),
            (0,     0.00625, "Cohesive",    1,	1.0),
            (-45,   0.125,    "T800-M21", 8,	4.0),
            (0,     0.00625, "Cohesive",    1,	1.0),
            (-45,   0.125,    "T800-M21", 8,	4.0),
            (0,     0.00625, "Cohesive",    1,	1.0),
            (45,    0.125,    "T800-M21", 8,	4.0)] # 从顶部到底部的每层的角度、厚度、材料以及厚度
    halfStructure=[0,0,1] # 每个tuple分别包含是否在x,y,z方向上是否使用半结构 0否 1是
    # 注意宽度方向上的对称条件不知道为什么一直错误 所以宽度方向尽量不要使用半结构
    EquationOrDisplacement=0 #在拉伸边界上使用 参考点Equation约束(0) 或者 位移边界(1)
    ex=0.06
    PBC=(EquationOrDisplacement==0)
    modelName='ResinLayer45-Cohesive-3'

    model1=mdb.Model(name=modelName)
    model1.setValues(description='Longitude: %s\n Width: %s \n Plies: angle\t thickness \t material \t seedsize In height \n %s'%(repr(ls),repr(ws),repr(plies)))
    
    Materials(model1)
    UMATMaterial(model1)

    ## Uguen A, Zubillaga L, Turon A, et al. Comparison of cohesive zone models used to predict delamination initiated from free-edges: validation against experimental results[C]//ECCM-16TH European Conference on Composite Materials. 2014.
    El,Et,Ez=130.0e3,8.0e3,8.0e3
    nult,nulz,nutz=0.31,0.31,0.45
    Glt,Glz,Gtz=4.0e3,4.0e3,4.0e3

    GIc,GIIc,GIIIc=0.24,0.74,0.74
    sigma_zz0,tau_xz0,tau_yz0=46,75,75
    Hillerborg=Et*GIIc/(tau_xz0**2)
    M=0.21
    le=M*Hillerborg/10

    material3=model1.Material(name='Cohesive')     
    material3.Elastic(type=TRACTION, table=((1e9, 1e9, 1e9), ))
    material3.MaxsDamageInitiation(table=((sigma_zz0,tau_xz0,tau_yz0), ))
    material3.maxsDamageInitiation.DamageEvolution(type=ENERGY, mixedModeBehavior=BK, power=2.0,table=((GIc, GIIc, GIIIc), ))
    
    m,j,notEqNodeLabel=modelCohesiveLaminate(model1,ls,ws,plies,halfStructure,EquationOrDisplacement,ex=ex)
    model1.rootAssembly.regenerate()

    if PBC:
        from Abaqus import lipeng
        import time
        oldtime=time.time()
        #lipeng.Pipes_Pagano(m,ls[0])
        SetINP,EquationINP=lipeng.Pipes_Pagano_INP(model1,notInclude=notEqNodeLabel,byNode=False)
        print(len(SetINP),len(EquationINP))
        j.writeInput(consistencyChecking=OFF)
        
        jobname=j.name
        src=open("./%s.inp"%jobname,"r")
        fname="./%s_Modified.inp"%jobname
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
        job=mdb.JobFromInputFile(name='%s_PBC'%jobname, inputFileName=fname)
        print('Elapsed:',time.time()-oldtime,'s' )
