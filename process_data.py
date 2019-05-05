#!/usr/bin/python2
# -*- coding: UTF-8 -*-

import os
os.chdir(r"E:/temp")
from abaqus import *
from abaqusConstants import *
import regionToolset
import mesh
import numpy as np
import visualization
import xyPlot

m1=mdb.models['Model-ZDIV-1']
asm=m1.rootAssembly
p1=asm.instances[asm.instances.keys()[0]].part
plyprefix=",".join(["%02.0f"%x.orientationValue for x in (p1.compositeLayups['Laminate'].plies[i] for i in xrange(len(p1.compositeLayups['Laminate'].plies)))])
datas={}
for div_ind in xrange(len(xyzdivs)):
  zdiv=xyzdivs[div_ind][2]
  datas[zdiv]={}
  m1=mdb.models['Model-ZDIV-%d'%zdiv]
  asm=m1.rootAssembly
  p1=asm.instances[asm.instances.keys()[0]]
  ns=p1.nodes.getByBoundingBox(xMin=lam_len/2.0-0.01,xMax=lam_len/2.0+0.01,zMin=plyLowZs[0]-0.01,zMax=plyLowZs[0]+0.01)
  sort_ns=sorted(ns,cmp=lambda x,y:1 if x.coordinates[1]>y.coordinates[1] else -1)
  pthy=session.Path(name="%s-ZDIV-%d-Y"%(plyprefix,zdiv),type=NODE_LIST,expression=[(p1.name.upper(),[n.label,]) for n in sort_ns])
  ns=p1.nodes.getByBoundingBox(xMin=lam_len/2.0-0.01,xMax=lam_len/2.0+0.01,yMax=0.0)
  sort_ns=sorted(ns,cmp=lambda x,y:1 if x.coordinates[2]>y.coordinates[2] else -1)
  pthz=session.Path(name="%s-ZDIV-%d-Z"%(plyprefix,zdiv),type=NODE_LIST,expression=[(p1.name.upper(),[n.label,]) for n in sort_ns])
  o = session.openOdb(name='./%s-ZDIV-%d.odb'%(plyprefix,zdiv),readOnly=True)
  session.viewports['Viewport: 1'].setValues(displayedObject=o)
  scratchOdb = session.ScratchOdb(o)
  scratchOdb.rootAssembly.DatumCsysByThreePoints(name='CSYS-1', 
      coordSysType=CARTESIAN, origin=(0.0, 0.0, 0.0), point1=(1.0, 0.0, 0.0), 
      point2=(0.0, 1.0, 0.0))
  dtm = scratchOdb.rootAssembly.datumCsyses['CSYS-1']
  session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(
      transformationType=visualization.USER_SPECIFIED, datumCsys=dtm)
  d={"y":{},"z":{}}
  for component in ('S11','S22','S33','S23','S13','S12','E11','E22','E33','E23','E13','E12','U1','U2','U3'):
    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
      outputPosition=NODAL if component[0]=='U' else INTEGRATION_POINT,
      variableLabel=component[0],  refinement=(COMPONENT, component))
    xy1=session.XYDataFromPath(name='%s-Y-%s-%d'%(plyprefix,component,zdiv,),path=pthy, includeIntersections=False, 
        projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
        projectionTolerance=0, shape=UNDEFORMED, labelType=TRUE_DISTANCE)
    d["y"][component]=xy1.data
    xy1=session.XYDataFromPath(name='%s-Z-%s-%d'%(plyprefix,component,zdiv, ),path=pthz, includeIntersections=False, 
        projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
        projectionTolerance=0, shape=UNDEFORMED, labelType=TRUE_DISTANCE)          
    d["z"][component]=xy1.data
  datas[zdiv]=d

fp=open('E:/XY-U-%s.txt'%(plyprefix),'w')
for zdiv in datas.keys():
  data=datas[zdiv]
  fp.write('ZDIV %d\n'%(zdiv))
  for vara in ('S','E'):
    for ind in ('11','22','33','12','13','23'):
      component=vara+ind
      fp.write(' Y \t '+component+'\t')
  fp.write('\n')
  for i in xrange(len(data["y"]['S11'])):
    for vara in ('S','E'):
      for ind in ('11','22','33','12','13','23'):
        component=vara+ind
        fp.write('%16f\t%16f\t'%(data["y"][component][min(i,len(data["y"][component])-1)]))
    fp.write('\n')
  for vara in ('S','E'):
    for ind in ('11','22','33','12','13','23'):
      component=vara+ind
      fp.write('Z \t '+component+'\t')
  fp.write('\n')
  for i in xrange(len(data["z"]['S11'])):
    for vara in ('S','E'):
      for ind in ('11','22','33','12','13','23'):
        component=vara+ind
        fp.write('%16f\t%16f\t'%(data["z"][component][min(i,len(data["z"][component])-1)]))
    fp.write('\n')
fp.close()

def elementBound(elm):
  ns=elm.getNodes()
  xmin=min((n.coordinates[0] for n in ns))
  xmax=max((n.coordinates[0] for n in ns))
  ymin=min((n.coordinates[1] for n in ns))
  ymax=max((n.coordinates[1] for n in ns))
  zmin=min((n.coordinates[2] for n in ns))
  zmax=max((n.coordinates[2] for n in ns))
  return ((xmin,ymin,zmin),(xmax,ymax,zmax))

def outputPlyAver(part1,plysetkeys):
  elnum=len(part1.elements)
  averarr=[None for i in range(elnum)]
  arr=[ 0 for i in range(elnum)]
  for setname in plysetkeys:
    if (setname in part1.allInternalSets):
      plyelements=part1.allInternalSets[setname].elements
    elif (setname in part1.allSets):
      plyelements=part1.allSets[setname].elements
    else:
      print setname,"not found in sets,Please Check"
      continue
    print setname
    for elm in plyelements:
      ind=elm.label-1
      averlist=averarr[ind]
      if averlist is None:
        bb=elementBound(elm)
        averes=plyelements.getByBoundingBox(xMin=bb[0][0],xMax=bb[1][0],yMin=bb[0][1],yMax=bb[1][1])
        averlist=list(e.label for e in averes)
        for elabel in averlist:
          averarr[elabel-1]=averlist
        for i in range(len(averlist)-1):
          arr[averlist[i]-1]=averlist[i+1]
        arr[averlist[-1]-1]=averlist[0]
  return averarr
def stressTrans(ls):
  """
    从新系应力分量转换到旧系应力分量
    
    ls[0],ls[1],ls[2]分别为新系的三个基在旧系的坐标
    inds 为6个应力分量的下标
  """
  inds=((1,1),(2,2),(3,3),(1,2),(1,3),(2,3))
  st=np.zeros((6,6),dtype='f')
  for i in range(6):
    i1,j1=inds[i]
    for j in range(6):
      i2,j2=inds[j]
      st[i,j]=(ls[i2-1][i1-1]*ls[j2-1][j1-1]+ls[j2-1][i1-1]*ls[i2-1][j1-1])
      if i2==j2:
        st[i,j]=st[i,j]/2
  return st
def strainTrans(ls):
  """
    从新系应变分量转换到旧系应变分量
    注意应变分量中剪切应变为2
    
    ls[0],ls[1],ls[2]分别为新系的三个基在旧系的坐标
    inds 为6个应变分量的下标
  """
  inds=((1,1),(2,2),(3,3),(1,2),(1,3),(2,3))
  st=np.zeros((6,6),dtype='f')
  for i in range(6):
    i1,j1=inds[i]
    for j in range(6):
      i2,j2=inds[j]
      st[i,j]=(ls[i2-1][i1-1]*ls[j2-1][j1-1]+ls[j2-1][i1-1]*ls[i2-1][j1-1])
      if i1==j1:
        st[i,j]=st[i,j]/2
  return st

aver_datas={}
for div_ind in xrange(len(xyzdivs)):
  zdiv=xyzdivs[div_ind][2]
  m1=mdb.models['Model-ZDIV-%d'%zdiv]
  asm=m1.rootAssembly
  p1=asm.instances[asm.instances.keys()[0]]
  o = session.openOdb(name='./%s-ZDIV-%d.odb'%(plyprefix,zdiv),readOnly=True)
  sfo=o.steps.values()[-1].frames[-1].fieldOutputs['S']
  es=p1.elements.getByBoundingBox(xMin=lam_len/2.0,xMax=lam_len/2.0+0.3,zMin=plyLowZs[0],zMax=plyLowZs[0]+plyHs[0]/zdiv)
  print(zdiv,len(es))
  plyelements=p1.sets['Ply-0'].elements
  sigtrans={ls:stressTrans(ls) for ls in set(map(lambda x:x.localCoordSystem,sfo.values))}
  arr=[]
  for elm in es:
    bb=elementBound(elm)    
    averes=plyelements.getByBoundingBox(xMin=bb[0][0],xMax=bb[1][0],yMin=bb[0][1],yMax=bb[1][1])
    arr.append((bb,np.array([np.dot(sigtrans[sfo.values[e.label-1].localCoordSystem],sfo.values[e.label-1].data) for e in averes]).mean(0)))
  aver_datas[zdiv]=arr

with open('E:/aver%s.txt'%(plyprefix),'w') as fp:
  for div_ind in xrange(len(xyzdivs)):
    zdiv=xyzdivs[div_ind][2]
    data=aver_datas[zdiv]
    fp.write('ZDIV %d\n'%zdiv)
    fp.write('x0\ty0\tz0\tx1\ty1\tz1\t')
    fp.write('\t'.join(['S11','S22','S33','S12','S13','S23']))
    fp.write('\n')
    for (bb1,bb2),s in data:
      fp.write('\t'.join(["%f"%x for x in bb1]))
      fp.write('\t')
      fp.write('\t'.join(["%f"%x for x in bb2]))
      fp.write('\t')
      fp.write('\t'.join(["%f"%x for x in s]))
      fp.write('\n')
