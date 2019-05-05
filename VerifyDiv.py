#!/usr/bin/python2
# -*- coding: UTF-8 -*-

"""
验证基于对称结构的Abaqus计算的有效性
"""

import os
os.chdir(r"E:/temp")
from abaqus import *
from abaqusConstants import *
import regionToolset
import mesh
import numpy as np
import visualization
import xyPlot

l,w=40.0,4.0
plyHs=[1.0,1.0]
plyAngles=[45.0,-45.0]#[90.0,0.0]#
xyzdivs=[(40,40,8)]#,(40,40,2),(40,40,4),(40,40,8)
# x方向上single的ration
xydoubleratio=[8.0,8.0,8.0,8.0,1.0]

plyprefix=",".join(["%02.0f"%x for x in plyAngles])

lam_len,lam_wid,lam_hgt=l,w,sum(plyHs)
plyLowZs=np.cumsum(plyHs,dtype='f')
plyLowZs=plyLowZs[-1]-plyLowZs
plynum=len(plyAngles)
plythick=lam_hgt/plynum

UMAT_PATH="E:/1730895/Work/BridgeMatrix/src/UMAT_lipeng_AVER.for"

#Mdb()
for div_ind in xrange(len(xyzdivs)):
  modelname='Model-Symm-Div-%d'%xyzdivs[div_ind][2]
  model1=mdb.Model(name=modelname,description='')

  # Part
  part1=model1.Part(name="Part-1",dimensionality=THREE_D,type=DEFORMABLE_BODY)
  s=model1.ConstrainedSketch(name="Sketch-1",sheetSize=200.0)
  s.rectangle(point1=(0,0),point2=(lam_len,lam_wid))
  part1.BaseSolidExtrude(sketch=s,depth=lam_hgt)

  for i in xrange(plynum-1):
    dt=part1.DatumPlaneByPrincipalPlane(offset=plyLowZs[i], principalPlane=XYPLANE)
    part1.PartitionCellByDatumPlane(cells=part1.cells,datumPlane=part1.datums[dt.id])

  # Property
  material1=model1.Material(name='UMAT')
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

  material2=model1.Material(name='Pipe-Pagano')
  material2.Elastic(type=ENGINEERING_CONSTANTS, table=((20000000.0, 2100000.0, 2100000.0, 0.21, 0.21, 0.21, 850000.0, 850000.0, 850000.0), ))

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
    z=plyLowZs[i]+plyHs[i]/2.0
    part1.Set(name="Ply-%1d"%(i),cells=part1.cells.findAt(((0,0,z),)) )
    r=regionToolset.Region(cells=part1.cells.findAt(((0,0,z),)))
    plyname='Ply-%1d-%3.1f'%(i,plyAngles[i])
    plysetkeys.append(plyname)
    compositeLayup.CompositePly(suppressed=False, plyName=plyname, region=r, 
      material='Pipe-Pagano', thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
      orientationType=SPECIFY_ORIENT, orientationValue=plyAngles[i], 
      additionalRotationType=ROTATION_NONE, additionalRotationField='', 
      axis=AXIS_3, angle=0.0, numIntPoints=1) # must use 1 for integral point 

  # Step
  model1.StaticStep(initialInc=1, name='Step-1', noStop=OFF, 
    previous='Initial', timeIncrementationMethod=FIXED)
  model1.fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 'CSTRESS', 'CDISP', 
    'SDV'))

  # Mesh

  part1.setElementType(elemTypes=(mesh.ElemType(
    elemCode=C3D8R, elemLibrary=STANDARD, secondOrderAccuracy=OFF, 
    kinematicSplit=AVERAGE_STRAIN, hourglassControl=ENHANCED, 
    distortionControl=DEFAULT), mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD), 
    mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)), regions=(part1.cells,))
  part1.assignStackDirection(cells=part1.cells, referenceRegion=part1.faces.getByBoundingBox(zMin=lam_hgt)[0])

  def e2vec(ns,e):
    n1=ns[e.getVertices()[0]].pointOn[0]
    n2=ns[e.getVertices()[1]].pointOn[0]
    return [n1[i]-n2[i] for i in xrange(3)]

  xedges=part1.edges.findAt(coordinates=map(lambda eg:eg.pointOn[0],
          filter(lambda eg:abs(abs(e2vec(part1.vertices,eg)[0])-lam_len)<1e-6,part1.edges)))
  yedges=part1.edges.findAt(coordinates=map(lambda eg:eg.pointOn[0],
          filter(lambda eg:abs(abs(e2vec(part1.vertices,eg)[1])-lam_wid)<1e-6,part1.edges)))
  zedges=part1.edges.findAt(coordinates=map(lambda eg:eg.pointOn[0],
          filter(lambda eg:sum([abs(x) for x in e2vec(part1.vertices,eg)][:2])<1e-6,part1.edges)))

  xdiv,ydiv,zdiv=xyzdivs[div_ind]
  ratio=xydoubleratio[div_ind]
  part1.seedEdgeByBias(biasMethod=DOUBLE, centerEdges=xedges, constraint=FINER, number=xdiv, ratio=ratio)
  part1.seedEdgeByBias(biasMethod=DOUBLE, constraint=FINER, endEdges=yedges, number=ydiv, ratio=ratio)
  part1.seedEdgeByNumber(constraint=FINER, edges=zedges, number=zdiv)
  part1.generateMesh()
  
  # Assembly
  rasm=model1.rootAssembly
  model1.rootAssembly.DatumCsysByDefault(CARTESIAN)
  model1.rootAssembly.Instance(dependent=ON, name='Part-1-1', 
      part=part1)

  # BC
  X0=rasm.Set(faces=
      rasm.instances['Part-1-1'].faces.getByBoundingBox(xMax=0.0), name='X0')
  model1.EncastreBC(createStepName='Initial', localCsys=None, name=
      'X=0-ENCAST', region=X0)
  X1=rasm.Set(faces=
      rasm.instances['Part-1-1'].faces.getByBoundingBox(xMin=lam_len,xMax=lam_len), name='X1')
  model1.DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
      distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
      'X=1-TENSION', region=X1, u1=0.1*lam_len, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
  Z0=rasm.Set(faces=
      rasm.instances['Part-1-1'].faces.getByBoundingBox(zMax=0.0), name='Z0')
  model1.ZsymmBC(name='Z=0-ZSYMM', createStepName='Initial', region=Z0, localCsys=None)    
  #Y0=rasm.Set(faces=
  #    rasm.instances['Part-1-1'].faces.getByBoundingBox(yMax=0.0), name='Y0')
  #model1.DisplacementBC(name='Y=0-LAMINATE_SYMM', createStepName='Initial', 
  #    region=Y0, u1=UNSET, u2=SET, u3=UNSET, ur1=SET, ur2=SET, ur3=UNSET, 
  #    amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)


  # Job    
  descinfo='%s\r\n %s\r\n %s\r\n'%(str((lam_len,lam_wid,lam_hgt)),str(xyzdivs[div_ind]),str(plyAngles))
  job=mdb.Job(name="%s-ZDIV-%d-Symm"%(plyprefix,zdiv),description=descinfo,model=model1,userSubroutine=UMAT_PATH,numCpus=1,numGPUs=0)
  #job.submit()
  #job.waitForCompletion()

#mdb.saveAs(pathName='./%s'%plyprefix)
