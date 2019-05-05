#! /usr/bin/python2
# -*- coding: UTF-8 -*-

import visualization
from abaqus import *
from abaqusConstants import *
import numpy as np

ACISEPS=1e-6

def importAbqusModule():
  import inspect
  import os,sys
  #curpath=path.dirname(path.dirname(path.abspath(__file__)))
  filename = inspect.getframeinfo(inspect.currentframe()).filename
  curpath = os.path.dirname(os.path.abspath(filename))
  pth=os.path.dirname(curpath)
  if pth not in sys.path:
      sys.path.append(pth)

def labels2sequence(labels,geomArr):
  """
  Abaqus中经常有需要用到Sequence的地方如在建立Set的时候cells,edges等参数
  我个人经常使用tuple来代替sequence,这时候就需要用这个来转换
  
  """
  seq=geomArr[0:0]
  for l in labels:
    seq+=geomArr[l:l+1]
  return seq

def feature2datum(f,part):
  """
  feature to datum
  part1.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=10) 返回feature对象
  但是part1.PartitionCellByDatumPlane(datumPlane=dt, cells=part1.cells)中的dt需要DatumPlane 对象
  """
  return part.datums[f.id]

def edge2vector(ns,e):
	"""
	返回Edge对象e 的 向量
	"""
	n1=ns[e.getVertices()[0]].pointOn[0]
	n2=ns[e.getVertices()[1]].pointOn[0]
	return [n1[i]-n2[i] for i in xrange(3)]


def prettyDir(obj):
  """
  members=["vertices","edges","faces","cells","nodes","elements"]
  d={}
  for mb in members:
    objArray=getattr(part1,mb)
    d[mb]=prettyDir(objArray)
    obj=objArray[0]
    d[mb[:-1]]=prettyDir(obj)
  for k,v in d.items():
    print k,":",v
  """
  return {"members":",".join(getattr(obj,'__members__',[])),"methods":"(),".join(getattr(obj,'__methods__',[]))+"()"}

def generateNodePath(model1,coord=None,pthname=None):
  """
  返回coord中定义的直线的节点和直线的方向"x","y","z"或None
  pth=lipeng.generate_path(m1,coord={"x":20,"z":1})
  """
  m1=model1
  eps=1e-6
  asm=m1.rootAssembly
  p1=asm.instances[asm.instances.keys()[0]]

  bounding={}
  for s in coord.keys():
    bounding["%sMin"%s]=coord[s]-eps
    bounding["%sMax"%s]=coord[s]+eps
  for k,s in enumerate("xyz"):
    if s not in coord:
      v=(s,k)
  ns=p1.nodes.getByBoundingBox(**bounding)
  
  f=lambda x:getattr(x,"coordinates")[v[1]]
  sort_ns=sorted(ns,cmp=lambda x,y: 1 if f(x)>f(y) else -1)

  name= ",".join(["%s=%.2f"%(s,x) for s,x in coord.items()]) if pthname is None else pthname
  session.Path(name=name,type=NODE_LIST,expression=[(p1.name.upper(),[n.label,]) for n in sort_ns])

  return session.paths[name]

def generateEdgePath(model1,coord=None,pthname=None):
  """
  返回coord中定义的Edge组成的path和直线的方向"x","y","z"或None
  pth=lipeng.generate_path(m1,coord={"x":20,"z":1})
  If type=EDGE_LIST, expression must be a sequence of sequences. Each inner sequence contains two items, the first item is a String specifying the name of the part instance, and the second item is a sequence of tuples of four Ints that uniquely identify an element edge. The four Ints are:

  The element label.

  The element face index (one-based).

  The face edge index (one-based).

  The edge direction. A positive number specifies that the edge direction runs from the edge start node to the edge end node. A negative number specifies the opposite.


  """
  m1=model1
  eps=1e-6
  asm=m1.rootAssembly
  p1=asm.instances[asm.instances.keys()[0]]

  bounding={}
  for s in coord.keys():
    bounding["%sMin"%s]=coord[s]-eps
    bounding["%sMax"%s]=coord[s]+eps
  for k,s in enumerate("xyz"):
    if s not in coord:
      v=(s,k)
  ns=p1.edges.getByBoundingBox(**bounding)
  
  f=lambda x:getattr(x,"coordinates")[v[1]]
  sort_ns=sorted(ns,cmp=lambda x,y: 1 if f(x)>f(y) else -1)
  for i in range(len(sort_ns)-1):
    n1=sort_ns[i]
    n2=sort_ns[i+1]
    e1=filter(lambda e: n2.label in set([n.label for n in e.getNodes()]),n1.getElemEdges(),)[0]
    [np.mean(np.array(map(lambda n:n.coordinates,f0.getNodes())),axis=0) for f0 in elm1.getElemFaces()]
  name= ",".join(["%s=%.2f"%(s,x) for s,x in coord.items()]) if pthname is None else pthname
  session.Path(name=name,type=NODE_LIST,expression=[(p1.name.upper(),[n.label,]) for n in sort_ns])

  return session.paths[name]

def extractDataFromPath(odb,pth,variables,suffix="",csyname=None):
  """
  在当前session中提取pth所对应的所有应力数据，并保存在XYData中名字后缀suffix
  lipeng.extractDataFromPath("Job-X1.odb",lipeng.generate_path(m1,coord={"x":20,"z":0.5}),suffix="X1")  
  """
  labelTypeSet={
    "x":TRUE_DISTANCE_X,
    "y":TRUE_DISTANCE_Y,
    "z":TRUE_DISTANCE_Z,
  }
  if isinstance(pth,str):
    pth=session.paths[pth]
  if isinstance(odb,str):
    o = session.openOdb(name=odb,readOnly=True)
  else:
    o = odb
  session.viewports['Viewport: 1'].setValues(displayedObject=o)
  if not (csyname is None):
    scratchOdb = session.ScratchOdb(o)
    scratchOdb.rootAssembly.DatumCsysByThreePoints(name='CSYS-1', 
        coordSysType=CARTESIAN, origin=(0.0, 0.0, 0.0), point1=(1.0, 0.0, 0.0), 
        point2=(0.0, 1.0, 0.0))
    dtm = scratchOdb.rootAssembly.datumCsyses['CSYS-1']
    session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(
        transformationType=visualization.USER_SPECIFIED, datumCsys=dtm)
  for var in variables:
    var=var.upper()
    try:
      displayDict={}
      if var[0] in ("U","S","E"):
        displayDict["variableLabel"]=var[0]
        displayDict["refinement"]=(COMPONENT,var)
      else:
        displayDict["variableLabel"]=var
      displayDict["outputPosition"]=INTEGRATION_POINT if var[0] in ('E','S','H') else NODAL
      session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(**displayDict)
      session.XYDataFromPath(name=var+suffix, path=pth, 
        includeIntersections=False, projectOntoMesh=False, pathStyle=PATH_POINTS, 
        numIntervals=10, projectionTolerance=0, shape=UNDEFORMED, 
        labelType=TRUE_DISTANCE)
    except Exception as e:
      print("extract %s error: %s"%(var,e))  
    
    
  #for sigma in ('S11','S22','S33','S23','S13','S12'):
  #  try:
  #    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
  #      variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 
  #      sigma))
  #    session.XYDataFromPath(name=sigma+suffix, path=pth, 
  #      includeIntersections=False, projectOntoMesh=False, pathStyle=PATH_POINTS, 
  #      numIntervals=10, projectionTolerance=0, shape=UNDEFORMED, 
  #      labelType=TRUE_DISTANCE)
  #  except Exception as e:
  #    print("extract %s error: %s"%(sigma,e))  

def elementBound(elm):
  ns=elm.getNodes()
  xmin=min((n.coordinates[0] for n in ns))
  xmax=max((n.coordinates[0] for n in ns))
  ymin=min((n.coordinates[1] for n in ns))
  ymax=max((n.coordinates[1] for n in ns))
  zmin=min((n.coordinates[2] for n in ns))
  zmax=max((n.coordinates[2] for n in ns))
  return ((xmin,ymin,zmin),(xmax,ymax,zmax))

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

def Pipes_Pagano(m1,l):
  """
  m1=mdb.models["Model-XZSYMM--Eq"]
  # 注意 要做之前把所有之前的Set和Equation 删除
  # 另外 这个函数的CAE操作很费时 建议直接使用下面的Pipe_Pagano_INP函数的返回结果修改
  """
  rootasm=m1.rootAssembly
  rootasm.regenerate()
  p1=rootasm.instances['Part-1-1']
  ns=p1.nodes
  x1ns=ns.getByBoundingBox(xMin=l-1e-6)
  for n0 in ns.getByBoundingBox(xMax=1e-6):
    xyz=n0.coordinates
    n1=x1ns.getByBoundingBox(yMin=xyz[1]-1e-6,yMax=xyz[1]+1e-6,zMin=xyz[2]-1e-6,zMax=xyz[2]+1e-6,)[0]
    s1name,s2name="Node-%d"%n0.label,"Node-%d"%n1.label
    s1=rootasm.Set(name=s1name,nodes=p1.nodes[(n0.label-1):n0.label])
    s2=rootasm.Set(name=s2name,nodes=p1.nodes[(n1.label-1):n1.label])
    m1.Equation(name='pp-%d-%d-1'%(n0.label,n1.label), terms=((1.0, s1name, 1), (-1.0, s2name, 1),(1,"x-ref",1)))
    m1.Equation(name='pp-%d-%d-2'%(n0.label,n1.label), terms=((1.0, s1name, 2), (-1.0, s2name, 2)))
    m1.Equation(name='pp-%d-%d-3'%(n0.label,n1.label), terms=((1.0, s1name, 3), (-1.0, s2name, 3)))  

ACISEPS=1e-6
def Pipes_Pagano_INP(m1):
  """
  返回两个字符串变量: SetINP,和EquationINP
  SetINP: 是每个节点都建立一个Set的inp语句 在inp文件中建议插入在Assembly节中End Instance后
  EquationINP： 是每个节点建立Equation的inp语句 在inp文件中建议插入在End Assembly之前
  """
  rasm=m1.rootAssembly
  part11=rasm.instances["Part-1-1"]
  lam_len,lam_width,lam_height=part11.cells.getBoundingBox()["high"][0]
  
  plyZs=set()
  for i in xrange(len(part11.cells)):
    bb=part11.cells[i:i+1].getBoundingBox()
    plyZs.add(bb["high"][2])
    plyZs.add(bb["low"][2])
  plyZs=sorted(list(plyZs),reverse=True)
  plynum=len(plyZs)-1
  SetINP=""
  EquationINP=""
  EquationNodeSet=set({})
  for i in xrange(plynum):
    zmin,zmax=plyZs[i+1],plyZs[i]
    ply=rasm.Set(name="PP_Ply%d"%i,cells=part11.cells.getByBoundingBox(zMin=zmin-ACISEPS,zMax=zmax+ACISEPS))
    ns=ply.nodes
    x1ns=ns.getByBoundingBox(xMin=lam_len-ACISEPS)
    for n0 in ns.getByBoundingBox(xMax=ACISEPS):
      if n0.label in EquationNodeSet:
        continue
      xyz=n0.coordinates
      n1=x1ns.getByBoundingBox(yMin=xyz[1]-ACISEPS,yMax=xyz[1]+ACISEPS,zMin=xyz[2]-ACISEPS,zMax=xyz[2]+ACISEPS,)[0]
      SetINP=SetINP+"*Nset, nset=Node-{0:d}, instance=Part-1-1\n{0:d},\n*Nset, nset=Node-{1:d}, instance=Part-1-1\n{1:d},\n".format(n0.label,n1.label)
      EquationINP=EquationINP+"*Equation\n3\nNode-{0:d}, 1, 1.\nNode-{1:d}, 1, -1.\nx-ref, 1, 1.\n*Equation\n2\nNode-{0:d}, 2, 1.\nNode-{1:d}, 2, -1.\n*Equation\n2\nNode-{0:d}, 3, 1.\nNode-{1:d}, 3, -1.\n".format(n0.label,n1.label)
      #EquationINP=EquationINP+"*Equation\n3\n{0:d}, 1, 1.\n{1:d}, 1, -1.\nx-ref, 1, 1.\n*Equation\n2\n{0:d}, 2, 1.\n{1:d}, 2, -1.\n*Equation\n2\n{0:d}, 3, 1.\n{1:d}, 3, -1.\n".format(n0.label,n1.label)
      EquationNodeSet.add(n0.label)  
  return SetINP,EquationINP