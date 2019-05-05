#!/usr/bin/python2
# -*- coding: UTF-8 -*-

"""
建立等厚度层合板模型

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

from Abaqus.lipeng import edge2vector

ACISEPS=1E-6

def modelLaminate(model1,ls,ws,plies,halfStructure,EquationOrDisplacement,ex=0.1):
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
          axis=AXIS_3, angle=0.0, numIntPoints=3) # must use 1 for integral point 
        print("ply ",i," end")

    ## Step
    model1.StaticStep(name='Step-1',noStop=OFF, previous='Initial', timeIncrementationMethod=AUTOMATIC, initialInc=0.1,maxNumInc=1000,minInc=1e-05)
    model1.fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'UR','RF', 'CF', 'CSTRESS', 'CDISP', 
    'SDV'))

    ## Mesh
    part1.setElementType(elemTypes=(mesh.ElemType(
    elemCode=C3D8R, elemLibrary=STANDARD, secondOrderAccuracy=OFF, 
    kinematicSplit=AVERAGE_STRAIN, hourglassControl=ENHANCED, 
    distortionControl=DEFAULT), mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD), 
    mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)), regions=(part1.cells,))
    part1.assignStackDirection(cells=part1.cells, referenceRegion=part1.faces.getByBoundingBox(zMin=lam_hgt)[0])

    xedges=part1.edges.findAt(coordinates=map(lambda eg:eg.pointOn[0],
          filter(lambda eg:abs(abs(edge2vector(part1.vertices,eg)[0])-lam_len)<ACISEPS,part1.edges)))
    part1.seedEdgeByBias(biasMethod=DOUBLE, centerEdges=xedges, constraint=FINER, number=xdiv, ratio=xratio)

    yedges=part1.edges.findAt(coordinates=map(lambda eg:eg.pointOn[0],
          filter(lambda eg:abs(abs(edge2vector(part1.vertices,eg)[1])-lam_wid)<ACISEPS,part1.edges)))
    part1.seedEdgeByBias(biasMethod=DOUBLE, constraint=FINER, endEdges=yedges, number=ydiv, ratio=yratio)
    
    for i in range(plynum):
        zdiv=plies[i][3]
        zratio=plies[i][4]
        zedges=part1.edges.findAt(coordinates=map(lambda eg:eg.pointOn[0],
            filter(lambda eg:sum([abs(x) for x in edge2vector(part1.vertices,eg)][:2])<ACISEPS,
                part1.edges.getByBoundingBox(zMin=plyLowZs[i]-ACISEPS,zMax=(plyLowZs[i-1] if i>0 else lam_hgt)+ACISEPS))))
        part1.seedEdgeByBias(biasMethod=DOUBLE, constraint=FINER, endEdges=zedges, number=zdiv,ratio=zratio)
    part1.generateMesh()

    ## Assembly
    rasm=model1.rootAssembly
    model1.rootAssembly.DatumCsysByDefault(CARTESIAN)
    myInstance1=model1.rootAssembly.Instance(dependent=ON, name='Part-1-1', 
      part=part1)

    ## BC
    region = rasm.Set(vertices=myInstance1.vertices.findAt(((0,0,0),)), name='Origin')
    model1.PinnedBC(name='X0Y0Z0-Pin', createStepName='Initial', region=region, localCsys=None)

    X0=rasm.Set(faces=myInstance1.faces.getByBoundingBox(xMax=0.0), name='X0')
    X1=rasm.Set(faces=myInstance1.faces.getByBoundingBox(xMin=lam_len), name='X1')
    Y0=rasm.Set(faces=myInstance1.faces.getByBoundingBox(yMax=0.0), name='Y0')
    Y1=rasm.Set(faces=myInstance1.faces.getByBoundingBox(yMin=lam_wid), name='Y1')
    Z0=rasm.Set(faces=myInstance1.faces.getByBoundingBox(zMax=0.0), name='Z0')
    Z1=rasm.Set(faces=myInstance1.faces.getByBoundingBox(zMin=lam_hgt), name='Z1')

    if halfStructure[0]==1:
        model1.XsymmBC(name='X=0-XSYMM', createStepName='Initial', region=X0, localCsys=None)    
    else:
        model1.EncastreBC(name='X=0-ENCAST',createStepName='Initial',region=X0, localCsys=None)
    if halfStructure[1]==1:
        #model1.YsymmBC(name='Y=0-YSYMM', createStepName='Initial', region=Y0, localCsys=None)
        model1.DisplacementBC(name='Y=0-LAMINATE_SYMM', createStepName='Initial', 
            region=Y0, u1=UNSET, u2=SET, u3=UNSET, ur1=SET, ur2=SET, ur3=UNSET, 
            amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
    if halfStructure[2]==1:
        model1.ZsymmBC(name='Z=0-ZSYMM', createStepName='Initial', region=Z0, localCsys=None)    

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

    ## Job    
    descinfo='%s\r\n %s\r\n %s\r\n'%(str((lam_len,lam_wid,lam_hgt)),",".join([str(x) for x in (xdiv,ydiv)]),str(plyAngles))
    suffix="".join(("XYZ"[i] if halfStructure[i]==1 else "" for i in xrange(3)))
    suffix=suffix+("" if len(suffix)==0 else "SYMM")
    #jobname="Job-%s-%s-%sX2"%(plyprefix,suffix,'Eq' if EquationOrDisplacement==0 else 'Disp')
    jobname="Job-%s"%model1.name
    job=mdb.Job(name=jobname,
        description=descinfo,model=model1,numCpus=1,numGPUs=0)
    
    return model1,job

if __name__=="__main__":
    
    model1=mdb.Model(name="VCCT")

    ## Material
    material1=model1.Material(name='Pipes-Pagano',description="1970 Pipes & Pagano\n (Psi)")
    material1.Elastic(type=ENGINEERING_CONSTANTS, table=((20000000.0, 2100000.0, 2100000.0, 0.21, 0.21, 0.21, 850000.0, 850000.0, 850000.0), ))
    
    material1=model1.Material(name='T300-7901',description="""FIBER:230000.00000000 15000.00000000 15000.00000000     0.03571429     0.20000000     0.20000000  7241.37931034 15000.00000000 15000.00000000  2500.00000000  2000.00000000     0.00000000
                                                              MATRIX:  3170.00000000  3170.00000000  3170.00000000     0.35500000     0.35500000     0.35500000  1169.74169742  1169.74169742  1169.74169742    85.10000000   107.00000000    52.60000000
                                                              VF=0.67""")
    material1.Elastic(type=ENGINEERING_CONSTANTS, table=((155.14610000000,9.54845283869,9.54845283869,0.25115000,0.25115000,0.27783858,5.43890821519,5.43890821519,3.73617332913), ))
    
    material2=model1.Material(name='Matrix',description=" 7901 (MPa)")
    material2.Elastic(type=ISOTROPIC, table=((3.17, 0.355),))
       
    ## Uguen A, Zubillaga L, Turon A, et al. Comparison of cohesive zone models used to predict delamination initiated from free-edges: validation against experimental results[C]//ECCM-16TH European Conference on Composite Materials. 2014.
    El,Et,Ez=130.0e3,8.0e3,8.0e3
    nult,nulz,nutz=0.31,0.31,0.45
    Glt,Glz,Gtz=4.0e3,4.0e3,4.0e3
    material3=model1.Material(name='T800-M21',description="2014 Uguen A & Zubillaga L & Turon A ")
    material3.Elastic(type=ENGINEERING_CONSTANTS, table=((El,Et,Ez,nult,nulz,nutz,Glt,Glz,Gtz), ))
    
    ## Lorriot T, Marion G, Harry R, et al. Onset of free-edge delamination in composite laminates under tensile loading[J]. Composites Part B: Engineering, 2003, 34(5): 459-471.
    El,Et,Ez=139.0e3,8.4e3,8.4e3
    nult,nulz,nutz=0.33,0.33,0.5
    Glt,Glz,Gtz=4.1e3,4.1e3,4.1e3 
    material3=model1.Material(name='T800-914',description="2003 Lorriot T, Marion G, Harry R, et al. Onset of free-edge delamination in composite laminates under tensile loading")
    material3.Elastic(type=ENGINEERING_CONSTANTS, table=((El,Et,Ez,nult,nulz,nutz,Glt,Glz,Gtz), ))
    
    ## 输入参数
    ls=(0.2,1,1.0) #  层合板长度(x方向) , 单元数目，以及double seed ratio (end/center)
    ws=(4.0,100,1.0)  #  层合板长度(y方向) , 单元数目，以及double seed ratio (center/end)
    plies=[ (45,    1,    "T300-7901", 4,	4.0),
            (0,     0.08, "Matrix",    1,	4.0),
            (-45,   1,    "T300-7901", 4,	4.0),
            (0,     0.08, "Matrix",    1,	4.0),
            (-45,   1,    "T300-7901", 4,	4.0),
            (0,     0.08, "Matrix",    1,	4.0),
            (45,    1,    "T300-7901", 4,	4.0)] # 从顶部到底部的每层的角度、厚度、材料以及厚度
    
    plies=[ (40,    1,    "T300-7901", 4,	4.0),
            (-45,   1,    "T300-7901", 4,	4.0),
            (-45,   1,    "T300-7901", 4,	4.0),
            (45,    1,    "T300-7901", 4,	4.0)] # 从顶部到底部的每层的角度、厚度、材料以及厚度

    plies=[ (40,    1,    "T300-7901", 4,	4.0),
            (-45,   1,    "T300-7901", 4,	4.0),
            (-45,   1,    "T300-7901", 4,	4.0),
            (45,    1,    "T300-7901", 4,	4.0)] # 从顶部到底部的每层的角度、厚度、材料以及厚度
    
    plies=[ (25,    0.125,    "T800-914", 8,	4.0),
            (25,    0.125,    "T800-914", 8,	4.0),
            (-25,   0.125,    "T800-914", 8,	4.0),
            (-25,   0.125,    "T800-914", 8,	4.0),
            (90,    0.125,    "T800-914", 8,	4.0),
            (90,    0.125,    "T800-914", 8,	4.0),
            (-25,   0.125,    "T800-914", 8,	4.0),
            (-25,   0.125,    "T800-914", 8,	4.0),
            (25,    0.125,    "T800-914", 8,	4.0),
            (25,    0.125,    "T800-914", 8,	4.0),
            ] # 从顶部到底部的每层的角度、厚度、材料以及厚度

    ws=(2.0,2,1.0)  #  层合板长度(y方向) , 单元数目，以及double seed ratio (center/end)
    plies=[ (45,    0.125,    "Pipes-Pagano", 1,    1.0),
            (-45,    0.125,   "Pipes-Pagano", 1,    1.0),
            (-45,   0.125,    "Pipes-Pagano", 1,    1.0),
            (45,    0.125,    "Pipes-Pagano", 1,    1.0),
            ] # 从顶部到底部的每层的角度、厚度、材料以及厚度

    halfStructure=[0,0,1] # 每个tuple分别包含是否在x,y,z方向上是否使用半结构 0否 1是
    # 注意宽度方向上的对称条件不知道为什么一直错误 所以宽度方向尽量不要使用半结构
    EquationOrDisplacement=0 #在拉伸边界上使用 参考点Equation约束(0) 或者 位移边界(1)

    ## 默认的model设置
    suffix="".join(("XYZ"[i] if halfStructure[i]==1 else "" for i in xrange(3)))
    suffix=suffix+("" if len(suffix)==0 else "SYMM")
    modelname='VCCT_%s_%s'%(suffix,'Eq' if EquationOrDisplacement==0 else 'Disp')

    model1.setValues(description='Longitude: %s\n Width: %s \n Plies: angle\t thickness \t material \t seedsize In height \n %s'%(repr(ls),repr(ws),repr(plies)))

    m,j=modelLaminate(model1,ls,ws,plies,halfStructure,EquationOrDisplacement)
    
    from Abaqus import lipeng
    import time
    print(time.time())
    lipeng.Pipes_Pagano(m,ls[0])
    print(time.time())
