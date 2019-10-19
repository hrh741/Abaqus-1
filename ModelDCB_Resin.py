#!/usr/bin/python2
# -*- coding: UTF-8 -*-

"""
建立带基体层的DCB模型
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
from Abaqus import lipeng

ACISEPS=1E-4

def modelDCB_ZSYMM(model1,ls,a,ws,plies,resinName,hresin=0.02,delta=1):
    """
    ls: tuple (层合板长度,长度方向的单元数目,double seed 的长度比例)
    ws: tuple (层合板宽度,宽度方向上的单元数目,double seed 的长度比例)
    """
    l,xdiv,xratio=ls
    w,ydiv,yratio=ws
    lam_len,lam_wid=l,w

    plynum=len(plies)
    plyAngles=[plies[i][0] for i in xrange(plynum)]
    plyThicks=[plies[i][1] for i in xrange(plynum)]
    plyProps=[plies[i][2] for i in xrange(plynum)]

    ## Part
    hr=hresin/2.0
    lam_hgt=sum(plyThicks)+hr
    plyprefix=",".join(["%02.0f"%x for x in plyAngles])
    plyLowZs=np.cumsum(plyThicks,dtype='f')
    plyLowZs=(hr+plyLowZs[-1])-plyLowZs
    print(plyLowZs)

    part1=model1.Part(name="Part-1",dimensionality=THREE_D,type=DEFORMABLE_BODY)
    s=model1.ConstrainedSketch(name="Sketch-1",sheetSize=200.0)
    H=plyThicks[0]+plyLowZs[0]
    s.Line((a,0),(l,0))
    s.Line((l,0),(l,hr))
    s.Line((l,hr),(l,H))
    s.Line((l,H),(0,H))
    s.Line((0,H),(0,hr))
    s.Line((0,hr),(a,hr))
    s.Line((a,hr),(a,0))
    part1.BaseSolidExtrude(sketch=s,depth=w)
    
    #return 
    part1.PartitionCellByExtendFace(extendFace=part1.faces.getByBoundingBox(xMax=a+ACISEPS,yMin=hr/2,yMax=hr+ACISEPS)[0], cells=part1.cells) 
    for i in xrange(plynum-1):
        dt=part1.DatumPlaneByPrincipalPlane(offset=plyLowZs[i], principalPlane=XZPLANE)
        part1.PartitionCellByDatumPlane(cells=part1.cells,datumPlane=part1.datums[dt.id])
    part1.regenerate()				  
    
    ## Composite Laminate 
    Dt=part1.DatumCsysByThreePoints(name='XZ-Datum', coordSysType=CARTESIAN, origin=(
      0.0, 0.0, 0.0), line1=(1.0, 0.0, 0.0), line2=(0.0, 0.0, -1.0))
    layupOrientation = part1.datums[Dt.id]
    compositeLayup = part1.CompositeLayup(
      name="Laminate", description='/'.join(( str(x) for x in plyAngles)), elementType=SOLID, 
      symmetric=False, thicknessAssignment=FROM_SECTION)
    compositeLayup.ReferenceOrientation(orientationType=SYSTEM, 
    localCsys=layupOrientation, fieldName='', 
    additionalRotationType=ROTATION_NONE, angle=0.0, 
    additionalRotationField='', axis=AXIS_3, stackDirection=STACK_3)

    #print("part1.cells:",len(part1.cells))
    plysetkeys=[]
    for i in xrange(plynum):
        print("ply ",i," start")
        z=plyLowZs[i]+plyThicks[i]/2.0
        r=part1.Set(name="Ply-%1d"%(i),cells=part1.cells.findAt(((0,z,0),)) )
        plyname='Ply-%1d-%3.1f'%(i,plyAngles[i])
        plysetkeys.append(plyname)
        print("ply ",i," 3","region:",r)
        compositeLayup.CompositePly(suppressed=False, plyName=plyname, region=r, 
          material=plyProps[i], thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
          orientationType=SPECIFY_ORIENT, orientationValue=plyAngles[i], 
          additionalRotationType=ROTATION_NONE, additionalRotationField='', 
          axis=AXIS_3, angle=0.0, numIntPoints=3) # must use 1 for integral point 
        print("ply ",i," end")

    r=part1.Set(name="Ply-resin",cells=part1.cells.findAt((((a+l)/2.0,hr/2.0,0),)) )
    compositeLayup.CompositePly(plyName="Ply-resin", region=r, 
          material=resinName, thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
          orientationType=SPECIFY_ORIENT, orientationValue=0, 
          additionalRotationType=ROTATION_NONE, additionalRotationField='', 
          axis=AXIS_3, angle=0.0, numIntPoints=3)     
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
    part1.assignStackDirection(cells=part1.cells, referenceRegion=part1.faces.getByBoundingBox(yMin=lam_hgt-ACISEPS)[0])

    xedges=part1.edges.findAt(coordinates=map(lambda eg:eg.pointOn[0],
          filter(lambda eg:abs(abs(lipeng.edge2vector(part1.vertices,eg)[0])-lam_len)<ACISEPS,part1.edges)))
    #part1.seedEdgeByBias(biasMethod=DOUBLE, centerEdges=xedges, constraint=FINER, number=xdiv, ratio=xratio)
    part1.seedEdgeBySize(edges=xedges, size=l/xdiv, deviationFactor=0.1, constraint=FINER)

    xedges1=part1.edges.findAt(coordinates=map(lambda eg:eg.pointOn[0],
          filter(lambda eg:abs(abs(lipeng.edge2vector(part1.vertices,eg)[0])-lam_len)>ACISEPS,part1.edges)))
    #part1.seedEdgeByBias(biasMethod=DOUBLE, endEdges=xedges1, minSize=l/xdiv,maxSize=l/xdiv, constraint=FINER)
    part1.seedEdgeBySize(edges=xedges1, size=l/xdiv, deviationFactor=0.1, constraint=FINER)

    yedges=part1.edges.findAt(coordinates=map(lambda eg:eg.pointOn[0],
          filter(lambda eg:abs(abs(lipeng.edge2vector(part1.vertices,eg)[2])-lam_wid)<ACISEPS,part1.edges)))
    part1.seedEdgeByBias(biasMethod=DOUBLE, constraint=FINER, endEdges=yedges, number=ydiv, ratio=yratio)
    
    for i in range(plynum):
        zdiv=plies[i][3]
        zratio=plies[i][4]
        zedges=part1.edges.findAt(coordinates=map(lambda eg:eg.pointOn[0],
            filter(lambda eg:sum([abs(x) for x in lipeng.edge2vector(part1.vertices,eg)][0:3:2])<ACISEPS,
                part1.edges.getByBoundingBox(yMin=plyLowZs[i]-ACISEPS,yMax=(plyLowZs[i-1] if i>0 else lam_hgt)+ACISEPS))))
        part1.seedEdgeByBias(biasMethod=DOUBLE, constraint=FINER, endEdges=zedges, number=zdiv,ratio=zratio)
    part1.generateMesh()

    ## Assembly
    rasm=model1.rootAssembly
    model1.rootAssembly.DatumCsysByDefault(CARTESIAN)
    myInstance1=model1.rootAssembly.Instance(dependent=ON, name='Part-1-1', part=part1)
    
    ## BC
    Y0=rasm.Set(faces=myInstance1.faces.getByBoundingBox(yMax=0.0), name='Y0')
    X0Y1=rasm.Set(edges=myInstance1.edges.getByBoundingBox(xMax=0.0+ACISEPS,yMin=H-ACISEPS), name='X0Y1')

    model1.rootAssembly.regenerate()
    
    part2 = model1.Part(name='Part-2', dimensionality=THREE_D,
            type=DEFORMABLE_BODY)
    part2.ReferencePoint(point=(lam_len+1, 0.0, 0.0))
    myInstance2=model1.rootAssembly.Instance(dependent=ON, name='Part-2-1', 
        part=part2)
    r1 = myInstance2.referencePoints.findAt((lam_len+1.0,0,0))
    refPoints1=(r1,)
    xref=rasm.Set(referencePoints = refPoints1, name='RP1')
    model1.Equation(name='X0Y1-1', terms=((1.0, 'X0Y1', 1), (-1.0, 'RP1', 1)))
    model1.Equation(name='X0Y1-2', terms=((1.0, 'X0Y1', 2), (-1.0, 'RP1', 2)))
    model1.Equation(name='X0Y1-3', terms=((1.0, 'X0Y1', 3), (-1.0, 'RP1', 3)))
    
    model1.DisplacementBC(name='DCB-Tension', createStepName='Step-1',
        region=xref, u1=0.0, u2=delta, u3=0.0, ur1=UNSET, ur2=UNSET,ur3=UNSET, fixed=OFF, distributionType=UNIFORM,localCsys=None)
    model1.historyOutputRequests['H-Output-1'].setValues(frequency=1)
    model1.HistoryOutputRequest(name='H-Output-Ref',createStepName='Step-1', variables=('U', 'RF'), region=xref,sectionPoints=DEFAULT, rebar=EXCLUDE)
    
    model1.YsymmBC(name='BC-YSYMM', createStepName='Initial', region=Y0, localCsys=None)
    jobname="Job-%s"%model1.name
    job=mdb.Job(name=jobname,model=model1)
    return job


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
def Suo(Exx,Eyy,Gxy,nuxy,a,h,B):
    nuyx=nuxy*Eyy/Exx
    rho=np.sqrt(Exx*Eyy)/(2.0*Gxy)-np.sqrt(nuxy*nuyx)
    beta=(0.677+0.146*(rho-1)-0.0178*(rho-1)**2+0.00242*(rho-1)**3)/np.power(Eyy/Exx,1/4)
    a_h=a/h
    C=1/(24/(B*Exx)*(a_h**3/3+beta*a_h**2+beta**2*a_h))
    return rho,beta,C

if __name__=="__main__":
    ls=(150.0,300,1.0) #  层合板长度(x方向) , 单元数目，以及double seed ratio (end/center)
    ws=(25.0,10,1.0)  #  层合板宽度(y方向) , 单元数目，以及double seed ratio (center/end)
    plies=[ (0,    1.5,    "T300-7901", 	20,	1.0),]
    resinName='Expoxy 7901'
    hresin=0.02

    crack_list=range(50,110,5)

    modelnameFormat='DCB-05mm-%d-2'
    for crack_a in crack_list:
        print '%d/%d : '%(crack_a-50,60),crack_a
        modelName=modelnameFormat%crack_a
        model1=mdb.Model(name=modelName)
        model1.setValues(description='Longitude: %s\n Width: %s \n Plies: angle\t thickness \t material \t seedsize In height \n %s'%(repr(ls),repr(ws),repr(plies)))
        for mat in materials:
            if mat[1]:
                material1=model1.Material(name=mat[0])
                if len(mat[1])==9:
                    material1.Elastic(type=ENGINEERING_CONSTANTS, table=(mat[1], ))
                elif len(mat[1])==2:
                    material1.Elastic(type=ISOTROPIC, table=(mat[1], ))

        j=modelDCB_ZSYMM(model1,ls,crack_a,ws,plies,resinName,hresin=hresin)
        import time
        oldTime=time.time()
        j.submit()
        j.waitForCompletion()
        newTime=time.time()
        print 'Elapsed ',newTime-oldTime,' s'

if 'DCB' not in globals():
    DCB=dict()

for a in crack_list:
    modelName=modelnameFormat%a
    print modelName
    o=openOdb('Job-%s.odb'%modelName,readOnly=True)
    model1=mdb.models[modelName]
    sfo=o.steps.values()[-1].frames[-1].fieldOutputs['S']
    es=model1.parts['Part-1'].elements.getByBoundingBox(xMax=a+ls[0]/ls[1]+ACISEPS,xMin=a-ACISEPS,yMax=hresin/2.0+ACISEPS)
    
    elementCoord={e.label:lipeng.elementBound(e) for e in es}
    sortedElm=sorted(elementCoord.keys(),key=lambda label:elementCoord[label][0][2])
    midElmLabel=sortedElm[int(len(sortedElm)/2)]
    d={e.label:None for e in es}
    cnt=0
    for value in sfo.values:
        if (value.elementLabel in d ) and d[value.elementLabel] is None:
            d[value.elementLabel]=(elementCoord[value.elementLabel],
                    value.data,value.localCoordSystem,
                    )
            cnt=cnt+1
            if(cnt==len(d)):
                break
    RF2=o.steps.values()[-1].historyRegions['Node PART-2-1.1'].historyOutputs['RF2'].data[1][1]
    DCB[a]=(d,RF2,midElmLabel)
    o.close()

with open('E:/DCB-05mm.txt','w') as fp:
    fp.write('\t'.join(['ModelName','CrackSize','RF2','ElementLabel','width-Z','S11','S22','S33','S12','S13','S23','RF2']))
    fp.write('\n')
    for a in crack_list:
        modelName=modelnameFormat%a
        d,RF2,midElmLabel=DCB[a]
        fp.write('\t'.join([modelName,'%d'%a,'%f'%RF2]))
        fp.write('\n')
        for k,v in d.items():
            coordBound,stress,localCoord=v
            widthY=(coordBound[0][2]+coordBound[1][2])/2.0
            fp.write('\t'.join(['','','','%d'%k,'%f'%widthY,'\t'.join(('%f'%x for x in stress))]))
            fp.write('\n')
        print modelName,a,'\t'.join(['%f'%x for x in d[midElmLabel][1]]),RF2

