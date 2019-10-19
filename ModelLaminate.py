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
    model1.StaticStep(name='Step-1',noStop=OFF, previous='Initial', timeIncrementationMethod=AUTOMATIC, initialInc=1.0,maxNumInc=1000,minInc=1e-05)
    model1.fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'UR','RF', 'CF', 'CSTRESS', 'CDISP', 
    'NFORC','SDV'))

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
    Origin = rasm.Set(vertices=myInstance1.vertices.findAt(((0,0,0),)), name='Origin')
    X0=rasm.Set(faces=myInstance1.faces.getByBoundingBox(xMax=0.0), name='X0')
    X1=rasm.Set(faces=myInstance1.faces.getByBoundingBox(xMin=lam_len), name='X1')
    Y0=rasm.Set(faces=myInstance1.faces.getByBoundingBox(yMax=0.0), name='Y0')
    Y1=rasm.Set(faces=myInstance1.faces.getByBoundingBox(yMin=lam_wid), name='Y1')
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
            frequency=1)
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
    descinfo='%s\r\n %s\r\n %s\r\n'%(str((lam_len,lam_wid,lam_hgt)),",".join([str(x) for x in (xdiv,ydiv)]),str(plyAngles))
    suffix="".join(("XYZ"[i] if halfStructure[i]==1 else "" for i in xrange(3)))
    suffix=suffix+("" if len(suffix)==0 else "SYMM")
    #jobname="Job-%s-%s-%sX2"%(plyprefix,suffix,'Eq' if EquationOrDisplacement==0 else 'Disp')
    jobname="Job-%s"%model1.name
    job=mdb.Job(name=jobname,
        description=descinfo,model=model1,numCpus=1,numGPUs=0)
    
    return model1,job,notEqNodeLabel

materials=[
    ## 吴庆欣
    ('T300-7901',(137.78e3,8.91e3,8.91e3,0.3,0.3,0.48,4.41e3,4.41e3,3.01e3),None,None), #
    ('Pipes-Pagano',(20.0, 2.1, 2.1, 0.21, 0.21, 0.21, 0.85, 0.85, 0.85),None,None), #MPsi
    ## Uguen A, Zubillaga L, Turon A, et al. Comparison of cohesive zone models used to predict delamination initiated from free-edges: validation against experimental results[C]//ECCM-16TH European Conference on Composite Materials. 2014.
    ### GIc,GIIc,GIIIc,sigma_zz0,tau_xz0,tau_yz0=0.24,0.74,0.74,46,75,75
    ('T800-M21',(130.0e3,8.0e3,8.0e3,0.31,0.31,0.45,4.0e3,4.0e3,4.0e3),None,None), #MPa
    ## 
    ('S2-SP250 GlasdEpoxy',None,None,None),
    ## Initiation of free-edge delamination in composite laminates
    ('G947-M18',(97.6e3, 8.0e3 ,8.0e3,0.37,0.37,0.5, 3.1e3, 3.1e3, 2.7e3),None,None),
    ## Lorriot T, Marion G, Harry R, et al. Onset of free-edge delamination in composite laminates under tensile loading[J]. Composites Part B: Engineering, 2003, 34(5): 459-471.
    ('T800-914',(159.0e3,8.4e3,8.4e3,0.31,0.31,0.45,4.0e3,4.0e3,4.0e3),None,None),
    ## 顾嘉杰 毕业论文
    ('IM7-914C',(153.6e3,10.2e3,10.2e3,0.27,0.27,0.46,5.7e3,5.7e3,3.5e3 ),None,None),
    ## 
    ('Expoxy 7901',(3.17e3, 0.355),None,None)
]


if __name__=="__main__":
    PBC=False
    # 几何参数

    """
    ws=(2.0,2,1.0)  #  层合板宽度(y方向) , 单元数目，以及double seed ratio (center/end)
    plies=[ (45,    0.125,    "Pipes-Pagano", 1,    1.0),
            (-45,    0.125,   "Pipes-Pagano", 1,    1.0),
            (-45,   0.125,    "Pipes-Pagano", 1,    1.0),
            (45,    0.125,    "Pipes-Pagano", 1,    1.0),
            ] # 从顶部到底部的每层的角度、厚度、材料以及厚度

    ls=(0.05,1,1.0) #  层合板长度(x方向) , 单元数目，以及double seed ratio (end/center)
    ws=(8.0,100,4.0)  #  层合板宽度(y方向) , 单元数目，以及double seed ratio (center/end)    
    # 从顶部到底部的每层的角度、厚度、材料以及厚度
    plies=[ (25,    0.125,    "IM7-914C", 8,	4.0),
            (25,    0.125,    "IM7-914C", 8,	4.0),
            (-25,   0.125,    "IM7-914C", 8,	4.0),
            (-25,   0.125,    "IM7-914C", 8,	4.0),
            (90,    0.0625,    "IM7-914C", 8,	4.0),
            (90,    0.0625,    "IM7-914C", 8,	4.0),
            (-25,   0.125,    "IM7-914C", 8,	4.0),
            (-25,   0.125,    "IM7-914C", 8,	4.0),
            (25,    0.125,    "IM7-914C", 8,	4.0),
            (25,    0.125,    "IM7-914C", 8,	4.0),
            ] 
    halfStructure=[0,0,1] # 每个tuple分别包含是否在x,y,z方向上是否使用半结构 0否 1是
    # 注意宽度方向上的对称条件不知道为什么一直错误 所以宽度方向尽量不要使用半结构
    EquationOrDisplacement=0 #在拉伸边界上使用 参考点Equation约束(0) 或者 位移边界(1)
    

    ls=(0.01,1,1.0) #  层合板长度(x方向) , 单元数目，以及double seed ratio (end/center)
    ws=(4.0,200,10.0)  #  层合板长度(x方向) , 单元数目，以及double seed ratio (center/end)
    plies=[ (20,    0.175,    "T800-M21", 	10,	8.0),
            #(0,     0.005,   "Cohesive",   	1 ,	8.0),
            (-20,   0.175,    "T800-M21", 	10,	8.0),
            (-20,   0.175,    "T800-M21", 	8 ,	8.0),
            #(0,     0.001,   "Cohesive",    1 ,	8.0),
            (20,    0.175,    "T800-M21", 	8 ,	8.0)] # 从顶部到底部的每层的角度、厚度、材料以及厚度
    halfStructure=[0,1,1] # 每个tuple分别包含是否在x,y,z方向上是否使用半结构 0否 1是
    # 注意宽度方向上的对称条件不知道为什么一直错误 所以宽度方向尽量不要使用半结构
    EquationOrDisplacement=0 #在拉伸边界上使用 参考点Equation约束(0) 或者 位移边界(1)
    
    ls=(0.01,1,1.0) #  层合板长度(x方向) , 单元数目，以及double seed ratio (end/center)
    ws=(20,50,20.0)  #  层合板长度(x方向) , 单元数目，以及double seed ratio (center/end)
    plies=[ (10,    0.125,    "Pipes-Pagano", 	5,	1.0),
            (45,   0.125,    "Pipes-Pagano", 	5,	1.0),
            (45,   0.125,    "Pipes-Pagano", 	5,	1.0),
            (10,    0.125,    "Pipes-Pagano", 	5,	1.0)] # 从顶部到底部的每层的角度、厚度、材料以及厚度
    halfStructure=[0,1,1] # 每个tuple分别包含是否在x,y,z方向上是否使用半结构 0否 1是
    # 注意宽度方向上的对称条件不知道为什么一直错误 所以宽度方向尽量不要使用半结构
    EquationOrDisplacement=0 #在拉伸边界上使用 参考点Equation约束(0) 或者 位移边界(1)
    ex=0.01
    PBC=(EquationOrDisplacement==0)
    modelName='Laminate-10-45-YZSYMM-5-2'

    ls=(200,100,10) #  层合板长度(x方向) , 单元数目，以及double seed ratio (end/center)
    ws=(20,50,20.0)  #  层合板长度(x方向) , 单元数目，以及double seed ratio (center/end)
    plies=[ (10,    0.125,    "Pipes-Pagano", 	5,	1.0),
            (45,   0.125,    "Pipes-Pagano", 	5,	1.0),
            (45,   0.125,    "Pipes-Pagano", 	5,	1.0),
            (10,    0.125,    "Pipes-Pagano", 	5,	1.0)] # 从顶部到底部的每层的角度、厚度、材料以及厚度
    halfStructure=[0,0,1] # 每个tuple分别包含是否在x,y,z方向上是否使用半结构 0否 1是
    # 注意宽度方向上的对称条件不知道为什么一直错误 所以宽度方向尽量不要使用半结构
    EquationOrDisplacement=1 #在拉伸边界上使用 参考点Equation约束(0) 或者 位移边界(1)
    ex=0.01
    PBC=(EquationOrDisplacement==0)
    modelName='Laminate-10-45-Wanzheng-5'

    ls=(200,50,10) #  层合板长度(x方向) , 单元数目，以及double seed ratio (end/center)
    ws=(20,50,20.0)  #  层合板长度(x方向) , 单元数目，以及double seed ratio (center/end)
    plies=[ (25,    0.125,    "Pipes-Pagano", 	6,	1.0),
            (-25,   0.125,    "Pipes-Pagano", 	6,	1.0),
            (90,   0.0625,    "Pipes-Pagano", 	3,	1.0),
            (90,   0.0625,    "Pipes-Pagano", 	3,	1.0),
            (-25,   0.125,    "Pipes-Pagano", 	6,	1.0),
            (25,    0.125,    "Pipes-Pagano", 	6,	1.0)] # 从顶部到底部的每层的角度、厚度、材料以及厚度
    halfStructure=[0,0,1] # 每个tuple分别包含是否在x,y,z方向上是否使用半结构 0否 1是
    # 注意宽度方向上的对称条件不知道为什么一直错误 所以宽度方向尽量不要使用半结构
    EquationOrDisplacement=1 #在拉伸边界上使用 参考点Equation约束(0) 或者 位移边界(1)
    ex=0.01
    PBC=(EquationOrDisplacement==0)
    modelName='Laminate-+-25-90-Wanzheng-1'

    ls=(0.05,1,1.0) #  层合板长度(x方向) , 单元数目，以及double seed ratio (end/center)
    ws=(8.0,400,5.0)  #  层合板宽度(y方向) , 单元数目，以及double seed ratio (center/end)
    plies=[ (45,    0.125,    "UMAT-Composite", 8,	4.0),
            (0,     0.00625, "UMAT-Matrix",    1,	1.0),
            (-45,   0.125,    "UMAT-Composite", 8,	4.0),
            (0,     0.00625, "UMAT-Matrix",    1,	1.0),
            (-45,   0.125,    "UMAT-Composite", 8,	4.0),
            (0,     0.00625, "UMAT-Matrix",    1,	1.0),
            (45,    0.125,    "UMAT-Composite", 8,	4.0)] # 从顶部到底部的每层的角度、厚度、材料以及厚度
    halfStructure=[0,0,1] # 每个tuple分别包含是否在x,y,z方向上是否使用半结构 0否 1是
    # 注意宽度方向上的对称条件不知道为什么一直错误 所以宽度方向尽量不要使用半结构
    EquationOrDisplacement=0 #在拉伸边界上使用 参考点Equation约束(0) 或者 位移边界(1)
    ex=0.06
    PBC=(EquationOrDisplacement==0)
    modelName='ResinLayer45-Quan-RF-3'
    """
    ls=(200,50,10) #  层合板长度(x方向) , 单元数目，以及double seed ratio (end/center)
    ws=(20,50,20.0)  #  层合板长度(x方向) , 单元数目，以及double seed ratio (center/end)
    plies=[ (0,    1.5,    "Pipes-Pagano", 	6,	1.0),
            (0,   1.5,    "Pipes-Pagano", 	6,	1.0),] # 从顶部到底部的每层的角度、厚度、材料以及厚度
    halfStructure=[0,0,1] # 每个tuple分别包含是否在x,y,z方向上是否使用半结构 0否 1是
    # 注意宽度方向上的对称条件不知道为什么一直错误 所以宽度方向尽量不要使用半结构
    EquationOrDisplacement=1 #在拉伸边界上使用 参考点Equation约束(0) 或者 位移边界(1)
    ex=0.01
    PBC=(EquationOrDisplacement==0)
    modelName='DCB-All'

    ## 生成模型
    model1=mdb.Model(name=modelName)
    model1.setValues(description='Longitude: %s\n Width: %s \n Plies: angle\t thickness \t material \t seedsize In height \n %s'%(repr(ls),repr(ws),repr(plies)))
    
    for mat in materials:
        if mat[1]:
            material1=model1.Material(name=mat[0])
            if len(mat[1])==9:
                material1.Elastic(type=ENGINEERING_CONSTANTS, table=(mat[1], ))
            elif len(mat[1])==2:
                material1.Elastic(type=ISOTROPIC, table=(mat[1], ))

    material1=model1.Material(name='UMAT-Composite',description='T300-7901')
    material1.Depvar(n=15)
    material1.UserMaterial(mechanicalConstants=
        (4100.0, 4100.0, 4100.0, 0.46, 0.46, 0.46, 1404.0, 1404.0, 1404.0,  # Matrix Elastic
        121.0, 210.0, 76.0, # Matrix : tensile compress shear 
        276000.0, 19000.0, 19000.0, 0.2, 0.2, 0.36, 27000.0, 27000.0, 6985.0, # Fiber Elastic
        4850.0, 3000.0, # fiber: tensile compress
        0.575, # Vf
        0.3, 0.3, # alpha,beta
        6.0, # MSEG 基体折线段数目 
        42.7,   53.7,   76.5,   101.5,  111.3,  125,# ETM(1,:)
        4.1E3,  2.5E3,  2E3,    1.4E3,  0.8E3,  0.41E3, # ETM(2,:)
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
        2.32225977 ,5.00655592, 1.65594384 , 1.44730762 , 1.86990463 , 1.68355870 ,1.70323652, 54.4, 0.01 ,1 #KT22,KTT22,KC22,K12,K23,KT22BI,KC22BI, 界面的临界脱粘强度, 衰减系数
    ))

    material1=model1.Material(name='UMAT-Matrix',description='7901')
    material1.Depvar(n=15)
    material1.UserMaterial(mechanicalConstants=
        (4100.0, 0.46, # Matrix Elastic
        121.0, 210.0, 76.0, # Matrix : tensile compress shear 
        6.0, # MSEG 基体折线段数目 
        42.7,   53.7,   76.5,   101.5,  111.3,  125,# ETM(1,:)
        4.1E3,  2.5E3,  2E3,    1.4E3,  0.8E3,  0.41E3, # ETM(2,:)
        0.001, # 衰减系数
        )
    )
    
    m,j,notEqNodeLabel=modelLaminate(model1,ls,ws,plies,halfStructure,EquationOrDisplacement,ex=ex)
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
        print 'Elapsed:',time.time()-oldtime,'s' 

    #lipeng.extractDataFromPath('Job-VCCT_PBC.odb','Path-z',['S11','S22','S33','S12','S13','S23'],prefix='+-25-90-Z-')
    #lipeng.extractDataFromPath('Job-VCCT_PBC.odb','Path-y_+25_-25',['S11','S22','S33','S12','S13','S23'],prefix='+25_-25-Y-')
    #lipeng.extractDataFromPath('Job-VCCT_PBC.odb','Path-y_-25-90',['S11','S22','S33','S12','S13','S23'],prefix='-25_90-Y-')
    #lipeng.extractDataFromPath('Job-VCCT_PBC.odb','Path-y',['S11','S22','S33','S12','S13','S23'],prefix='mid-Y-')

    #lipeng.extractDataFromPath('Job-+25_2_-25_2_90s_PBC.odb','Path-z',['S11','S22','S33','S12','S13','S23'],prefix='gjj-Z-')
    """
for i in range(1,14,1):
	session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=i)
	lipeng.extractDataFromPath('Job-ResinLayer45-Quan-1_PBC.odb','Path-3',['S11','S22','S33','S12','S13','S23','NFORC1','NFORC2','NFORC3'],prefix='Resin-Quan-2-',suffix='-F%1d'%i,csyname='CSYC-1')
    """