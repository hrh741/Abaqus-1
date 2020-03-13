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
from Abaqus.Materials import Materials,UMATMaterial

ACISEPS=1E-4
def modelDCBResin(model1,ls,a0,a,ws,plies,resinName,hresin=0.02,Nresin=2,delta=1,rPlastic=0.0):
    """
    ls: tuple (层合板长度,长度方向的单元数目,double seed 的长度比例)
    ws: tuple (层合板宽度,宽度方向上的单元数目,double seed 的长度比例)
    rPlastic: 裂尖partition的长度,用于画局部加密网格
    a0: 树脂富集区尖端长度
    a: 裂纹长度
    """
    #l,xdiv,xratio=ls
    l,XseedDict=ls
    w,ydiv,yratio=ws
    lam_len,lam_wid=l,w

    halfplynum=len(plies)/2
    halfplyAngles=[plies[i][0] for i in xrange(halfplynum)]
    halfplyThicks=[plies[i][1] for i in xrange(halfplynum)]

    ## Part
    hr=hresin/2
    halfplyLowZs=np.cumsum(halfplyThicks,dtype='f')
    halfplyLowZs=(hr+halfplyLowZs[-1])-halfplyLowZs
    lam_hgt=sum(halfplyThicks)+hr
    
    plyAngles=np.concatenate((halfplyAngles,[0,],np.flipud(halfplyAngles)))
    plyprefix=",".join(["%02.0f"%x for x in plyAngles])
    plyThicks=np.concatenate((halfplyThicks,[hresin,],np.flipud(halfplyThicks)))
    plyLowZs=np.concatenate((halfplyLowZs,[-hr,],-np.flipud(halfplyLowZs+halfplyThicks))) # abaqus numpy version limitation
    plies.insert(halfplynum,[0,hresin,resinName,Nresin,1.0])
    plynum=len(plies)
    plyProps=[plies[i][2] for i in  xrange(plynum)]
    print(plyLowZs)

    part1=model1.Part(name="Part-1",dimensionality=THREE_D,type=DEFORMABLE_BODY)
    s=model1.ConstrainedSketch(name="Sketch-1",sheetSize=200.0)
    H=plyThicks[0]+plyLowZs[0]
    s.Line((l,hr),(l,H))
    s.Line((l,H),(0,H))
    s.Line((0,H),(0,hr))
    s.Line((0,hr),(a0,hr))
    s.Line((a0,hr),(a0,-hr))
    s.Line((a0,-hr),(0,-hr))
    s.Line((0,-hr),(0,-H),)
    s.Line((0,-H),(l,-H))
    s.Line((l,-H),(l,-hr))
    s.Line((l,-hr),(l,hr))
    part1.BaseSolidExtrude(sketch=s,depth=w)
    
    #return 
    part1.PartitionCellByExtendFace(extendFace=part1.faces.getByBoundingBox(xMax=a0+ACISEPS,yMin=hr-ACISEPS,yMax=hr+ACISEPS)[0], cells=part1.cells) 
    part1.PartitionCellByExtendFace(extendFace=part1.faces.getByBoundingBox(xMax=a0+ACISEPS,yMin=-hr-ACISEPS,yMax=-hr+ACISEPS)[0], cells=part1.cells) 
    for i in xrange(plynum-1):
        dt=part1.DatumPlaneByPrincipalPlane(offset=plyLowZs[i], principalPlane=XZPLANE)
        try:
            part1.PartitionCellByDatumPlane(cells=part1.cells,datumPlane=part1.datums[dt.id])
        except:
            print('Partition y=',plyLowZs[i],' failed !')
    
    dt=part1.DatumPlaneByPrincipalPlane(offset=a0, principalPlane=YZPLANE)
    part1.PartitionCellByDatumPlane(cells=part1.cells,datumPlane=part1.datums[dt.id])    

    for x in ((a,a+rPlastic) if rPlastic>0 else (a,)):
        if x==a0:
            continue
        dt=part1.DatumPlaneByPrincipalPlane(offset=x, principalPlane=YZPLANE)
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
        #z=plyLowZs[i]+plyThicks[i]/2.0
        #r=part1.Set(name="Ply-%1d"%(i),cells=part1.cells.findAt(((0,z,0),)) )
        z1,z2=plyLowZs[i],plyLowZs[i]+plyThicks[i]
        z=(z1+z2)/2.0
        r=part1.Set(name="Ply-%1d"%(i),cells=part1.cells.getByBoundingBox(yMin=z1-ACISEPS,yMax=z2+ACISEPS))
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
    'NFORC','SDV'))
    ## Mesh
    part1.setElementType(elemTypes=(mesh.ElemType(
    elemCode=C3D8R, elemLibrary=STANDARD, secondOrderAccuracy=OFF, 
    kinematicSplit=AVERAGE_STRAIN, hourglassControl=ENHANCED, 
    distortionControl=DEFAULT), mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD), 
    mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)), regions=(part1.cells,))
    part1.assignStackDirection(cells=part1.cells, referenceRegion=part1.faces.getByBoundingBox(yMin=lam_hgt-ACISEPS)[0])

    #xedges=part1.edges.findAt(coordinates=map(lambda eg:eg.pointOn[0],
    #      filter(lambda eg:abs(abs(lipeng.edge2vector(part1.vertices,eg)[0])-lam_len)<ACISEPS,part1.edges)))
    ##part1.seedEdgeByBias(biasMethod=DOUBLE, centerEdges=xedges, constraint=FINER, number=xdiv, ratio=xratio)
    #part1.seedEdgeBySize(edges=xedges, size=l/xdiv, deviationFactor=0.1, constraint=FINER)
#
    #xedges1=part1.edges.findAt(coordinates=map(lambda eg:eg.pointOn[0],
    #      filter(lambda eg:abs(abs(lipeng.edge2vector(part1.vertices,eg)[0])-lam_len)>ACISEPS,part1.edges)))
    ##part1.seedEdgeByBias(biasMethod=DOUBLE, endEdges=xedges1, minSize=l/xdiv,maxSize=l/xdiv, constraint=FINER)
    #part1.seedEdgeBySize(edges=xedges1, size=l/xdiv, deviationFactor=0.1, constraint=FINER)
    for xl in XseedDict:
        if(xl==0):
            continue
        xdiv,xratio=XseedDict[xl]
        coords=map(lambda eg:eg.pointOn[0],
                    filter(lambda eg:abs(abs(lipeng.edge2vector(part1.vertices,eg)[0])-xl)<ACISEPS,
                            part1.edges))
        print(xl,xdiv,xratio,coords)
        if len(coords)>0:
            xedges=part1.edges.findAt(coordinates=coords)
            part1.seedEdgeByBias(biasMethod=DOUBLE, endEdges=xedges, constraint=FINER, number=int(xdiv), ratio=xratio)
            #rasm.seedEdgeBySize(edges=xedges, size=l/xdiv, deviationFactor=0.1, constraint=FINER)

    yedges=part1.edges.findAt(coordinates=map(lambda eg:eg.pointOn[0],
          filter(lambda eg:abs(abs(lipeng.edge2vector(part1.vertices,eg)[2])-lam_wid)<ACISEPS,part1.edges)))
    part1.seedEdgeByBias(biasMethod=DOUBLE, constraint=FINER, endEdges=yedges, number=ydiv, ratio=yratio)
    
    for i in range(plynum):
        zdiv=plies[i][3]
        zratio=plies[i][4]
        print(plyLowZs[i]-ACISEPS,(plyLowZs[i-1] if i>0 else lam_hgt)+ACISEPS)
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
    Y0=rasm.Set(faces=myInstance1.faces.getByBoundingBox(xMin=a,yMax=0.0), name='Y0')
    X0Y1=rasm.Set(edges=myInstance1.edges.getByBoundingBox(xMax=0.0+ACISEPS,yMin=H-ACISEPS), name='X0Y1')
    X0Y_1=rasm.Set(edges=myInstance1.edges.getByBoundingBox(xMax=0.0+ACISEPS,yMax=-H+ACISEPS), name='X0Y_1')
    model1.rootAssembly.regenerate()
    
    for (ind,refPoints),setName,u2 in zip(enumerate(((0, H+1, 0.0),(0, 0-1, 0.0))),
                                        ('X0Y1','X0Y_1'),(delta,-delta)):
        partName='Part-%d'%(ind+2)
        rpName='RP-%d'%(ind+1)
        part = model1.Part(name=partName, dimensionality=THREE_D,
                type=DEFORMABLE_BODY)
        part.ReferencePoint(point=refPoints)
        myInstance=model1.rootAssembly.Instance(dependent=ON, name='%s-1'%partName, 
            part=part)
        refPoints1=(myInstance.referencePoints.findAt(refPoints),)
        xref=rasm.Set(referencePoints = refPoints1, name=rpName)
        model1.Equation(name='%s-1'%setName, terms=((1.0, setName, 1), (-1.0, rpName, 1)))
        model1.Equation(name='%s-2'%setName, terms=((1.0, setName, 2), (-1.0, rpName, 2)))
        model1.Equation(name='%s-3'%setName, terms=((1.0, setName, 3), (-1.0, rpName, 3)))

        model1.DisplacementBC(name='DCB-Tension-%d'%(ind+1), createStepName='Step-1',
            region=xref, u1=0.0, u2=u2, u3=0.0, ur1=UNSET, ur2=UNSET,ur3=UNSET, fixed=OFF, distributionType=UNIFORM,localCsys=None)
        model1.HistoryOutputRequest(name='H-Output-Ref-%d'%(ind+1),createStepName='Step-1',
            variables=('U', 'RF'), region=xref,sectionPoints=DEFAULT, rebar=EXCLUDE)
    
    #model1.YsymmBC(name='BC-YSYMM', createStepName='Initial', region=Y0, localCsys=None)
    jobname="Job-%s"%model1.name
    job=mdb.Job(name=jobname,model=model1)
    return job

def modelDCB_ZSYMM(model1,ls,a0,a,ws,plies,resinName,hresin=0.02,Nresin=1,delta=1,rPlastic=0.0):
    """
    ls: tuple (层合板长度,长度方向的单元数目,double seed 的长度比例)
    ws: tuple (层合板宽度,宽度方向上的单元数目,double seed 的长度比例)
    rPlastic: 裂尖partition的长度,用于画局部加密网格
    a0: 树脂富集区尖端长度
    a: 裂纹长度
    """
    #l,xdiv,xratio=ls
    l,XseedDict=ls
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
    s.Line((a0,0),(l,0))
    s.Line((l,0),(l,hr))
    s.Line((l,hr),(l,H))
    s.Line((l,H),(0,H))
    s.Line((0,H),(0,hr))
    s.Line((0,hr),(a0,hr))
    s.Line((a0,hr),(a0,0))
    part1.BaseSolidExtrude(sketch=s,depth=w)
    
    #return 
    part1.PartitionCellByExtendFace(extendFace=part1.faces.getByBoundingBox(xMax=a0+ACISEPS,yMin=hr/2,yMax=hr+ACISEPS)[0], cells=part1.cells) 
    for i in xrange(plynum-1):
        dt=part1.DatumPlaneByPrincipalPlane(offset=plyLowZs[i], principalPlane=XZPLANE)
        part1.PartitionCellByDatumPlane(cells=part1.cells,datumPlane=part1.datums[dt.id])
    
    dt=part1.DatumPlaneByPrincipalPlane(offset=a0, principalPlane=YZPLANE)
    part1.PartitionCellByDatumPlane(cells=part1.cells,datumPlane=part1.datums[dt.id])    

    for x in ((a,a+rPlastic) if rPlastic>0 else (a,)):
        if x==a0:
            continue
        dt=part1.DatumPlaneByPrincipalPlane(offset=x, principalPlane=YZPLANE)
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
        #z=plyLowZs[i]+plyThicks[i]/2.0
        #r=part1.Set(name="Ply-%1d"%(i),cells=part1.cells.findAt(((0,z,0),)) )
        z1,z2=plyLowZs[i],plyLowZs[i]+plyThicks[i]
        z=(z1+z2)/2.0
        r=part1.Set(name="Ply-%1d"%(i),cells=part1.cells.getByBoundingBox(yMin=z1-ACISEPS,yMax=z2+ACISEPS))
        plyname='Ply-%1d-%3.1f'%(i,plyAngles[i])
        plysetkeys.append(plyname)
        print("ply ",i," 3","region:",r)
        compositeLayup.CompositePly(suppressed=False, plyName=plyname, region=r, 
          material=plyProps[i], thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
          orientationType=SPECIFY_ORIENT, orientationValue=plyAngles[i], 
          additionalRotationType=ROTATION_NONE, additionalRotationField='', 
          axis=AXIS_3, angle=0.0, numIntPoints=1) # must use 1 for integral point 
        print("ply ",i," end")

    r=part1.Set(name="Ply-resin",cells=part1.cells.getByBoundingBox(yMax=hr+ACISEPS))
    compositeLayup.CompositePly(plyName="Ply-resin", region=r, 
          material=resinName, thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
          orientationType=SPECIFY_ORIENT, orientationValue=0, 
          additionalRotationType=ROTATION_NONE, additionalRotationField='', 
          axis=AXIS_3, angle=0.0, numIntPoints=1)     
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

    #xedges=part1.edges.findAt(coordinates=map(lambda eg:eg.pointOn[0],
    #      filter(lambda eg:abs(abs(lipeng.edge2vector(part1.vertices,eg)[0])-lam_len)<ACISEPS,part1.edges)))
    ##part1.seedEdgeByBias(biasMethod=DOUBLE, centerEdges=xedges, constraint=FINER, number=xdiv, ratio=xratio)
    #part1.seedEdgeBySize(edges=xedges, size=l/xdiv, deviationFactor=0.1, constraint=FINER)
#
    #xedges1=part1.edges.findAt(coordinates=map(lambda eg:eg.pointOn[0],
    #      filter(lambda eg:abs(abs(lipeng.edge2vector(part1.vertices,eg)[0])-lam_len)>ACISEPS,part1.edges)))
    ##part1.seedEdgeByBias(biasMethod=DOUBLE, endEdges=xedges1, minSize=l/xdiv,maxSize=l/xdiv, constraint=FINER)
    #part1.seedEdgeBySize(edges=xedges1, size=l/xdiv, deviationFactor=0.1, constraint=FINER)
    for xl in XseedDict:
        if(xl==0):
            continue
        xdiv,xratio=XseedDict[xl]
        coords=map(lambda eg:eg.pointOn[0],
                    filter(lambda eg:abs(abs(lipeng.edge2vector(part1.vertices,eg)[0])-xl)<ACISEPS,
                            part1.edges))
        print xl,xdiv,xratio,coords
        if len(coords)>0:
            xedges=part1.edges.findAt(coordinates=coords)
            part1.seedEdgeByBias(biasMethod=DOUBLE, endEdges=xedges, constraint=FINER, number=int(xdiv), ratio=xratio)
            #rasm.seedEdgeBySize(edges=xedges, size=l/xdiv, deviationFactor=0.1, constraint=FINER)

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
    Y0=rasm.Set(faces=myInstance1.faces.getByBoundingBox(xMin=a,yMax=0.0), name='Y0')
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
    
    model1.DisplacementBC(name='DCB-Tension-1', createStepName='Step-1',
        region=xref, u1=0.0, u2=delta, u3=0.0, ur1=UNSET, ur2=UNSET,ur3=UNSET, fixed=OFF, distributionType=UNIFORM,localCsys=None)
    model1.historyOutputRequests['H-Output-1'].setValues(frequency=1)
    model1.HistoryOutputRequest(name='H-Output-Ref',createStepName='Step-1', variables=('U', 'RF'), region=xref,sectionPoints=DEFAULT, rebar=EXCLUDE)
    
    model1.YsymmBC(name='BC-YSYMM', createStepName='Initial', region=Y0, localCsys=None)
    jobname="Job-%s"%model1.name
    job=mdb.Job(name=jobname,model=model1)
    return job

if __name__=="__main__":
    """ 
    for n  in [0.115,0.25,0.5,1,1.5,2,3]:
        crack_a=50
        rPlastic=6
        ls=(150.0,{rPlastic:[int(rPlastic/n),1.0],crack_a:[int(crack_a/n),1.0],(150-rPlastic-crack_a):[10,1.0]}) #  层合板长度(x方向) , 单元数目，以及double seed ratio (end/center)
        ws=(25.0,10,1.0)  #  层合板宽度(y方向) , 单元数目，以及double seed ratio (center/end)
        plies=[ (0,    1.5,    "T300-7901", 	20,	1.0),]
        resinName='Epoxy7901'
        hresin=0.01
        criticalDelta={50 :2.77 ,52 :2.98 ,54 :3.20 ,56 :3.43 ,58 :3.67 ,
                    60 :3.91 ,62 :4.17 ,64 :4.43 ,66 :4.70 ,68 :4.97 ,
                    70 :5.26 ,72 :5.55 ,74 :5.85 ,76 :6.16 ,78 :6.47 ,
                    80 :6.80 ,82 :7.13 ,84 :7.47 ,86 :7.81 ,88 :8.17 ,
                    90 :8.53 ,92 :8.90 ,94 :9.28 ,96 :9.67 ,98 :10.07 ,
                    100:10.47 ,102:10.88 ,104:11.30 ,106:11.73 ,108:12.16 }

        modelName="DCB-Resin001-%d"%int(1000*n)
        model1=mdb.Model(name=modelName)
        model1.setValues(description='Longitude: %s\n Width: %s \n Plies: angle\t thickness \t material \t seedsize In height \n %s'%(repr(ls),repr(ws),repr(plies)))
        # 材料
        Materials(model1)

        modelDCB_ZSYMM(model1,ls,crack_a,ws,plies,resinName,
                        hresin=hresin,delta=criticalDelta[crack_a],rPlastic=rPlastic)

    """
    a0=50.0
    crack_a=52.0
    rPlastic=8.0
    ls=(150.0,{150:[150,1.0]}) #  层合板长度(x方向) , 单元数目，以及double seed ratio (end/center)
    # 对DCB样件宽度方向建议>=10
    ws=(25.0,10,1.0)  #  层合板宽度(y方向) , 单元数目，以及double seed ratio (center/end)
    plies=[ 
        (0,    1.5,    "T300-7901", 	20,	1.0),
        ]

    resinName='Epoxy7901'
    hresin=0.02
    modelnameFormat='DCB-New-%d-7'

    crack_list=range(50,100,5)
    Ntotal,Npermm=len(crack_list),4

    for ind,crack_a in enumerate(crack_list):
        print '%d/%d : '%(ind,Ntotal),crack_a
        ls=(150.0,{150:[150,1.0]})
        ls[1][a0]=[int(2*a0),1.0]
        ls[1][crack_a-a0]=[int(Npermm*(crack_a-a0)),1.0]
        ls[1][rPlastic]=[int(Npermm*rPlastic),1.0]
        ls[1][150-crack_a-rPlastic]=[10,1.0]

        modelName=modelnameFormat%crack_a
        model1=mdb.Model(name=modelName)
        model1.setValues(description='Longitude: %s\n Width: %s \n Plies: angle\t thickness \t material \t seedsize In height \n %s'%(repr(ls),repr(ws),repr(plies)))
        # 材料
        Materials(model1)

        j=modelDCB_ZSYMM(model1,ls,a0,crack_a,ws,plies,resinName,hresin=hresin,rPlastic=rPlastic)
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
        es=model1.parts['Part-1'].elements.getByBoundingBox(xMax=a+rPlastic/ls[1][rPlastic][0]+ACISEPS,xMin=a-ACISEPS,yMax=hresin/2.0+ACISEPS)
        
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

    with open('E:/%s.txt'%(modelnameFormat%(Npermm)),'w') as fp:
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
