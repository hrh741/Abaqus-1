#!/usr/bin/python2
# -*- coding: UTF-8 -*-

"""
建立DCB模型
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
from Abaqus import Materials,UMATMaterial

ACISEPS=1E-4
def modelDCB(model1,ls,a,ws,plies,delta=1,zsymm=True):
    """
    ls: tuple (层合板长度, 长度方向不同边长对应的seed数目和ratio的字典XseedDict)
    ws: tuple (层合板宽度,宽度方向上的单元数目,double seed 的ration)
    a : 裂纹长度
    plies : list of  各个单层的信息(铺排角,厚度,材料,seed数目,double ratio)
    zsymm : 是否使用Z 对称建模
    """
    l,XseedDict=ls
    w,ydiv,yratio=ws
    
    ## 根据半结构设置生成自定义参数 不要修改
    lam_len,lam_wid=l,w
    plynum=len(plies)//(2 if zsymm else 1)
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
    #if zsymm:
    #    part1.PartitionFaceByDatumPlane(faces=part1.faces.findAt(((lam_len/2.0,lam_wid/2.0,0))),datumPlane=part1.datums[dt.id])
    #else:
    #    part1.PartitionCellByDatumPlane(cells=part1.cells,datumPlane=part1.datums[dt.id])
    for x in (a,a+2):
        dt=part1.DatumPlaneByPrincipalPlane(offset=x, principalPlane=YZPLANE)
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
        print("ply ",i," start")
        z1,z2=plyLowZs[i],plyLowZs[i]+plyThicks[i]
        z=(z1+z2)/2.0
        r=part1.Set(name="Ply-%1d"%(i),cells=part1.cells.getByBoundingBox(zMin=z1-ACISEPS,zMax=z2+ACISEPS))
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
    'NFORC','SDV'))
    

    ## Assembly
    rasm=model1.rootAssembly
    model1.rootAssembly.DatumCsysByDefault(CARTESIAN)
    part11=model1.rootAssembly.Instance(dependent=OFF, name='Part-1-1', part=part1)
    
    if not zsymm:
        crack= rasm.Set(faces=part11.faces.getByBoundingBox(xMax=a+ACISEPS,zMax=lam_hgt/2+ACISEPS,zMin=lam_hgt/2.0-ACISEPS), name='Crack')
        rasm.engineeringFeatures.assignSeam(regions=crack)

    ## Mesh
    rasm.setElementType(elemTypes=(mesh.ElemType(
    elemCode=C3D8R, elemLibrary=STANDARD, secondOrderAccuracy=OFF, 
    kinematicSplit=AVERAGE_STRAIN, hourglassControl=ENHANCED, 
    distortionControl=DEFAULT), mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD), 
    mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)), regions=(part11.cells,))
    rasm.assignStackDirection(cells=part11.cells, referenceRegion=part11.faces.getByBoundingBox(zMin=lam_hgt-ACISEPS)[0])

    for xl in XseedDict:
        xdiv,xratio=XseedDict[xl]
        coords=map(lambda eg:eg.pointOn[0],
                    filter(lambda eg:abs(abs(lipeng.edge2vector(part11.vertices,eg)[0])-xl)<ACISEPS,
                            part11.edges))
        print(xl,xdiv,xratio,coords)
        if len(coords)>0:
            xedges=part11.edges.findAt(coordinates=coords)
            rasm.seedEdgeByBias(biasMethod=DOUBLE, endEdges=xedges, constraint=FINER, number=xdiv, ratio=xratio)
            #rasm.seedEdgeBySize(edges=xedges, size=l/xdiv, deviationFactor=0.1, constraint=FINER)

    yedges=part11.edges.findAt(coordinates=map(lambda eg:eg.pointOn[0],
          filter(lambda eg:abs(abs(lipeng.edge2vector(part11.vertices,eg)[1])-lam_wid)<ACISEPS,part11.edges)))
    rasm.seedEdgeByBias(biasMethod=DOUBLE, constraint=FINER, endEdges=yedges, number=ydiv, ratio=yratio)
    
    for i in range(plynum):
        zdiv=plies[i][3]
        zratio=plies[i][4]
        zedges=part11.edges.findAt(coordinates=map(lambda eg:eg.pointOn[0],
            filter(lambda eg:sum([abs(x) for x in lipeng.edge2vector(part11.vertices,eg)][0:2])<ACISEPS,
                part11.edges.getByBoundingBox(zMin=plyLowZs[i]-ACISEPS,zMax=(plyLowZs[i-1] if i>0 else lam_hgt)+ACISEPS))))
        rasm.seedEdgeByBias(biasMethod=DOUBLE, constraint=FINER, endEdges=zedges, number=zdiv,ratio=zratio)
    rasm.generateMesh(regions=(part11,))

    ## BC
    Z0=rasm.Set(faces=part11.faces.getByBoundingBox(xMin=a-ACISEPS,zMax=0.0), name='Z0')
    X0Z1=rasm.Set(edges=part11.edges.getByBoundingBox(xMax=0.0+ACISEPS,zMin=lam_hgt-ACISEPS), name='X0Z1')
    X0Z0=rasm.Set(edges=part11.edges.getByBoundingBox(xMax=0.0+ACISEPS,zMax=ACISEPS), name='X0Z0')
    model1.rootAssembly.regenerate()
    
    part2 = model1.Part(name='Part-2', dimensionality=THREE_D,
            type=DEFORMABLE_BODY)
    part2.ReferencePoint(point=(0,lam_wid/2,lam_hgt))
    myInstance2=model1.rootAssembly.Instance(dependent=ON, name='Part-2-1', 
        part=part2)
    r1 = myInstance2.referencePoints.findAt((0,lam_wid/2,lam_hgt))
    xref=rasm.Set(referencePoints = (r1,), name='RP1')
    model1.DisplacementBC(name='DCB-Tension-1', createStepName='Step-1',
        region=xref, u1=0.0, u2=0.0, u3=delta, ur1=UNSET, ur2=UNSET,ur3=UNSET, fixed=OFF, distributionType=UNIFORM,localCsys=None)
    model1.historyOutputRequests['H-Output-1'].setValues(frequency=1)
    model1.HistoryOutputRequest(name='H-Output-RF-1',createStepName='Step-1', variables=('U', 'RF'), region=xref,sectionPoints=DEFAULT, rebar=EXCLUDE)
    
    model1.Equation(name='X0Z1-1', terms=((1.0, 'X0Z1', 1), (-1.0, 'RP1', 1)))
    model1.Equation(name='X0Z1-2', terms=((1.0, 'X0Z1', 2), (-1.0, 'RP1', 2)))
    model1.Equation(name='X0Z1-3', terms=((1.0, 'X0Z1', 3), (-1.0, 'RP1', 3)))

    if zsymm:
        model1.ZsymmBC(name='BC-ZSYMM', createStepName='Initial', region=Z0, localCsys=None)
    else:
        part3 = model1.Part(name='Part-3', dimensionality=THREE_D,
                type=DEFORMABLE_BODY)
        part3.ReferencePoint(point=(0,lam_wid/2,0.0))
        myInstance3=model1.rootAssembly.Instance(dependent=ON, name='Part-3-1', 
            part=part3)
        r2 = myInstance3.referencePoints.findAt((0,lam_wid/2.0,0.0))
        xref1=rasm.Set(referencePoints = (r2,), name='RP2')
        model1.DisplacementBC(name='DCB-Tension-2', createStepName='Step-1',
            region=xref1, u1=0.0, u2=0.0, u3=-delta, ur1=UNSET, ur2=UNSET,ur3=UNSET, fixed=OFF, distributionType=UNIFORM,localCsys=None)
        model1.HistoryOutputRequest(name='H-Output-RF-2',createStepName='Step-1', variables=('U', 'RF'), region=xref1,sectionPoints=DEFAULT, rebar=EXCLUDE)

        model1.Equation(name='X0Z0-1', terms=((1.0, 'X0Z0', 1), (-1.0, 'RP2', 1)))
        model1.Equation(name='X0Z0-2', terms=((1.0, 'X0Z0', 2), (-1.0, 'RP2', 2)))
        model1.Equation(name='X0Z0-3', terms=((1.0, 'X0Z0', 3), (-1.0, 'RP2', 3)))

    jobname="Job-%s"%model1.name
    job=mdb.Job(name=jobname,model=model1)
    return job


if __name__=="__main__":
    crack_a=50
    ls=(150.0,{2:[5,1.0],50:[100,1.0],98:[10,1.0]}) #  层合板长度(x方向) , 单元数目，以及double seed ratio (end/center)
    ws=[25.0,5,1.0]  #  层合板宽度(y方向) , 单元数目，以及double seed ratio (center/end)
    plies=[ [0,    1.5,    "T300-7901", 	40,	1.0],
            (0,    1.5,    "T300-7901", 	40,	1.0),]
    zsymm=True
    modelName='DCB'
    
    seeds=[16,400,20,5,50]
    modelName='DCB-%d-%d-%d-%d-%d'%(seeds[0],seeds[1],seeds[2],seeds[3],seeds[4])
    print(modelName)
    ls[1][2][0]=seeds[0]
    ls[1][50][0]=seeds[1]
    ls[1][98][0]=seeds[2]
    ws[1]=seeds[3]
    plies[0][3]=seeds[4]

    model1=mdb.Model(name=modelName)
    model1.setValues(description='Longitude: %s\n Width: %s \n Plies: angle\t thickness \t material \t seedsize In height \n %s'%(repr(ls),repr(ws),repr(plies)))
    Materials(model1)

    j=modelDCB(model1,ls,crack_a,ws,plies,zsymm=zsymm)
    import time
    oldTime=time.time()
    j.submit()
    j.waitForCompletion()
    newTime=time.time()
    o=openOdb('Job-%s.odb'%modelName,readOnly=True)
    RF3=o.steps.values()[-1].historyRegions['Node PART-2-1.1'].historyOutputs['RF3'].data[1][1]
    print(modelName,RF3)
    o.close()
    print('Elapsed ',newTime-oldTime,' s')

    """
    for i in range(5):
        seeds=[10,20,30,10,20]
        for n in [2,80]:
            seeds[i]=n
            
            modelName='DCB-%d-%d-%d-%d-%d'%(seeds[0],seeds[1],seeds[2],seeds[3],seeds[4])
            print(modelName)
            ls[1][2][0]=seeds[0]
            ls[1][50][0]=seeds[1]
            ls[1][98][0]=seeds[2]
            ws[1]=seeds[3]
            plies[0][3]=seeds[4]

            model1=mdb.Model(name=modelName)
            model1.setValues(description='Longitude: %s\n Width: %s \n Plies: angle\t thickness \t material \t seedsize In height \n %s'%(repr(ls),repr(ws),repr(plies)))
            Materials(model1)
            j=modelDCB(model1,ls,crack_a,ws,plies,zsymm=zsymm)
            import time
            oldTime=time.time()
            j.submit()
            j.waitForCompletion()
            newTime=time.time()
            print('Elapsed ',newTime-oldTime,' s')

for i in range(5):
    seeds=[10,20,30,10,20]
    for n in [2,5,10,20,40,80]:
        seeds[i]=n
        
        modelName='DCB-%d-%d-%d-%d-%d'%(seeds[0],seeds[1],seeds[2],seeds[3],seeds[4])
        model1=mdb.models[modelName]
        es=model1.rootAssembly.instances['Part-1-1'].elements.getByBoundingBox(xMax=crack_a+2.0/seeds[0]+ACISEPS,xMin=crack_a-ACISEPS,yMax=1.5/seeds[-1]+ACISEPS)
        
        o=openOdb('Job-%s.odb'%modelName,readOnly=True)
        #sfo=o.steps.values()[-1].frames[-1].fieldOutputs['S']
        RF3=o.steps.values()[-1].historyRegions['Node PART-2-1.1'].historyOutputs['RF3'].data[1][1]
        print(modelName,RF3)
        o.close()

    model1=mdb.Model(name=modelName)
    model1.setValues(description='Longitude: %s\n Width: %s \n Plies: angle\t thickness \t material \t seedsize In height \n %s'%(repr(ls),repr(ws),repr(plies)))
    Materials(model1)
    modelDCB(model1,ls,crack_a,ws,plies,zsymm=zsymm)

    crack_list=range(50,51)
    modelnameFormat='DCB-1mm-%d'
    for crack_a in crack_list:
        print('%d/%d : '%(crack_a-50,60),crack_a)
        modelName=modelnameFormat%crack_a

        j=modelDCB(model1,ls,crack_a,ws,plies)
        import time
        oldTime=time.time()
        j.submit()
        j.waitForCompletion()
        newTime=time.time()
        print('Elapsed ',newTime-oldTime,' s')

    if 'DCB' not in globals():
        DCB=dict()

    for a in crack_list:
        modelName=modelnameFormat%a
        print(modelName)
        o=openOdb('Job-%s.odb'%modelName,readOnly=True)
        model1=mdb.models[modelName]
        sfo=o.steps.values()[-1].frames[-1].fieldOutputs['S']
        es=model1.parts['Part-1'].elements.getByBoundingBox(xMax=a+ls[0]/ls[1]+ACISEPS,xMin=a-ACISEPS,yMax=hresin/2.0+ACISEPS)
        
        elementCoord={e.label:lipeng.elementBound(e) for e in es}
        sortedElm=sorted(elementCoord.keys(),key=lambda label:elementCoord[label][0][0][2])
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

    with open('E:/DCB-4mm.txt','w') as fp:
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
            print(modelName,a,'\t'.join(['%f'%x for x in d[midElmLabel][1]]),RF2)

    """