#! /usr/bin/python2
# -*- coding: UTF-8 -*-

import visualization
from abaqus import *
from abaqusConstants import *
import numpy as np
import textwrap
from itertools import chain

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

#######################################################################
##                              建模                                 ##
#######################################################################

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

def elementBound(elm):
  """
  返回单元几何界限
  ()
  """
  ns=elm.getNodes()
  xmin=min((n.coordinates[0] for n in ns))
  xmax=max((n.coordinates[0] for n in ns))
  ymin=min((n.coordinates[1] for n in ns))
  ymax=max((n.coordinates[1] for n in ns))
  zmin=min((n.coordinates[2] for n in ns))
  zmax=max((n.coordinates[2] for n in ns))
  return ((xmin,ymin,zmin),(xmax,ymax,zmax))

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

#######################################################################
##                              后处理                               ##
#######################################################################
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

def extractDataFromPath(odb,pth,variables,prefix="",suffix="",csyname=None):
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
  scratchOdb = session.ScratchOdb(o)
  if csyname is None:
    csyname='CSYS-1'
    if csyname not in scratchOdb.rootAssembly.datumCsyses:
      scratchOdb.rootAssembly.DatumCsysByThreePoints(name=csyname, 
          coordSysType=CARTESIAN, origin=(0.0, 0.0, 0.0), point1=(1.0, 0.0, 0.0), 
          point2=(0.0, 1.0, 0.0))
  dtm = scratchOdb.rootAssembly.datumCsyses[csyname]
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
      if var.startswith('NFORC'):
        displayDict["outputPosition"]=ELEMENT_NODAL
      session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(**displayDict)
      session.XYDataFromPath(name=prefix+var+suffix, path=pth, 
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

def AbaqusStressTrans(ls):
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

def AbaqusStrainTrans(ls):
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

def calMean(jobname,setname=None):
  """
  提取odb中setname最后一帧的应力的平均值
  如果不设置setname 或者
  本程序将只含有一个积分点的单元在积分点积分点的应力作为单元的平均应力
  然后按初始体积计算集合所有单元的应力平均值
  
  calMean("Job-1","FIBER") 
  calMean("Job-1","MATRIX")
  """
  getdata=lambda v:v.data
  odb=openOdb(jobname+".odb")
  laststep=odb.steps[odb.steps.keys()[-1]]
  lastframe=laststep.frames[-1]
  rasm=odb.rootAssembly
  inst=rasm.instances[rasm.instances.keys()[0]]
  if setname:
    if setname not in inst.elementSets:
      raise Exception,'no set %s in %s, please Check !'%(setname,inst.name)
    sigma=np.array(map(getdata,lastframe.fieldOutputs['S'].getSubset(region=inst.elementSets[setname]).values))
    evol=np.array(map(getdata,lastframe.fieldOutputs['EVOL'].getSubset(region=inst.elementSets[setname]).values))
  else:
    setname="ALL"
    sigma=np.array(map(getdata,lastframe.fieldOutputs['S'].values))
    evol=np.array(map(getdata,lastframe.fieldOutputs['EVOL'].values))
  closeOdb(odb)
  print '%s Mean:\n'%(setname),np.dot(evol,sigma)/np.sum(evol)


#######################################################################
##                              inp修改                              ##
#######################################################################
def Pipes_Pagano(m1,l):
  """
  只适用于角铺层
  m1=mdb.models["Model-XZSYMM--Eq"]
  # 注意 要做之前把所有之前的Set和Equation 删除 不然会更加费时
  # 这个函数执行CAE操作会很费时 建议直接使用下面的Pipe_Pagano_INP函数的返回结果修改inp文件
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
def Pipes_Pagano_INP(m1,instanceName='Part-1-1',notInclude=set(),byNode=False,refName='x-ref'):
  """
  只适用于角铺层
  对层合板的每一层建立Pipes Pagano模型的位移约束
  必要条件:
    1. 层厚度为z方向
  输入
    m1: model 对象
    instanceName: 层合板所在的instance 的名称
    notInclude: 不建立Equation的节点label,因为Equation约束的第一个节点位移不能再有其它约束(边界条件也不能有)
    byNode : Bool 对象，如果为True，每个Node 建立一个Set, 如果为False,将Equation的所有节点建立在Set内
    refName : 参考点集合名称
  返回
    SetINP: 字符串，每个节点都建立一个Set的inp语句 建议插入在inp文件中Assembly节中End Instance后
    EquationINP： 字符串，是每个节点建立Equation的inp语句 建议插入在inp文件中End Assembly之前
  """
  rasm=m1.rootAssembly
  part11=rasm.instances[instanceName]

  lam_len,lam_width,lam_height=part11.cells.getBoundingBox()["high"]
  
  tol=max(lam_len,lam_width,lam_height)*ACISEPS*10
  plyZs=set()
  for i in xrange(len(part11.cells)):
    bb=part11.cells[i:i+1].getBoundingBox()
    plyZs.add(bb["high"][2])
    plyZs.add(bb["low"][2])
  plyZs=sorted(list(plyZs),reverse=True)
  plynum=len(plyZs)-1
  EquationNodeSet=set({})
  EquationPair=[]
  formatStr=textwrap.dedent("""\
  3
  Node-{0:d}, 1, 1.,Node-{1:d}, 1, -1.,x-ref, 1, 1.
  2
  Node-{0:d}, 2, 1.,Node-{1:d}, 2, -1.
  2
  Node-{0:d}, 3, 1.,Node-{1:d}, 3, -1.
  """)
  for i in xrange(plynum):
    zmin,zmax=plyZs[i+1],plyZs[i]
    rasm.Set(name="PP_Ply%d"%i,cells=part11.cells.getByBoundingBox(zMin=zmin-tol,zMax=zmax+tol))
    ply=rasm.sets["PP_Ply%d"%i]
    ns=ply.nodes
    x1ns=ns.getByBoundingBox(xMin=lam_len-tol)
    for n0 in ns.getByBoundingBox(xMax=tol):
      if n0.label in EquationNodeSet or n0.label in notInclude:
        continue
      EquationNodeSet.add(n0.label)
      xyz=n0.coordinates
      n1=x1ns.getByBoundingBox(yMin=xyz[1]-tol,yMax=xyz[1]+tol,zMin=xyz[2]-tol,zMax=xyz[2]+tol,)[0]
      EquationNodeSet.add(n1.label)
      EquationPair.append((n0.label,n1.label))

  if byNode:
    SetINP=""
    EquationINP=""
    for label in  EquationNodeSet:
      SetINP=SetINP+"*Nset, nset=Node-{0:d}, instance={1:s}\n{0:d},\n".format(label,instanceName)  

    EquationINP+="*Equation\n"
    for pair in EquationPair:
      EquationINP+="3\nNode-%d, 1, 1.,Node-%d, 1, -1.,%s, 1, 1.\n"%(pair[0],pair[1],refName)
  else:
    N=15 # Nset的行数据最多为16
    SetINP="*Nset, nset=Pipes_Pagano_Set-0, instance=%s,UNSORTED\n"%instanceName
    for i in range(len(EquationPair)/N+1):
      SetINP=SetINP+",".join((str(x[0]) for x in EquationPair[N*i:min(len(EquationPair),N*i+N)]))+'\n'
    
    SetINP=SetINP+"*Nset, nset=Pipes_Pagano_Set-1, instance=%s,UNSORTED\n"%instanceName
    for i in range(len(EquationPair)/N+1):
      SetINP=SetINP+",".join((str(x[1]) for x in EquationPair[N*i:min(len(EquationPair),N*i+N)]))+'\n'
    
    EquationINP=textwrap.dedent("""\
    *Equation
    2
    Pipes_Pagano_Set-0, 2, 1.,Pipes_Pagano_Set-1, 2, -1.
    2
    Pipes_Pagano_Set-0, 3, 1.,Pipes_Pagano_Set-1, 3, -1.
    3
    Pipes_Pagano_Set-0, 1, 1.,Pipes_Pagano_Set-1, 1, -1., %s, 1, 1.
    """%refName)
  return SetINP,EquationINP

def PBC(m1,byNode=False):
  """
  输入:
    m1: abaqus Model 对象
    byNode: 是否按节点建立Equation (默认为按Set建立), 注意 按节点建立Equation需要很多个只包含一个节点的集合 所以会导致文件很大
  返回:
    SetInp: 用于生成Set 的inp文件字符串
    EquationInp: 用于生成 Equation 的inp文件字符串
        目前未对8个顶点进行处理,只对内部节点和边界的内部节点建立了Equation 
        互相平行的四条边 每个方向根据周期性边界条件建立了3个Equation,如果建立4个Abaqus会报错

    建议在inp文件的*End Assembly 行前插入 SetInp 和EquationINp
  描述:
    创建周期性边界条件
    想更加深入了解RVE的PBC怎么设置 可以参考
    Development of an ABAQUS plugin tool for periodic RVE homogenisation
    Omairey S L, Dunning P D, Sriramula S. Development of an ABAQUS plugin tool for periodic RVE homogenisation[J]. Engineering with Computers, 2019, 35(2): 567-577.
    将节点分为三类:
      1. 角点 innerNodes
      2. 面内部节点 edges
      2. 边内部点 corners: 是文章中的顺序 (0,0,0) (1,0,0) (1,0,1) (0,0,1)  (0,1,0) (1,1,0) (1,1,1) (0,1,1)  

  问题：
     这个程序建立Equation没有问题，但是求解出来的结果会在位移边界上出现抖动，不连续，目前还不清楚原因
     不建议使用
  ```
  m1=mdb.models["Model-1"]
  #SetINP,EquationINP,EquationPair,EquationPairDict=PeriodicBC(m1)
  SetINP,EquationINP=PBC(m1,byNode=False)
  with open("E:/PBC_EquationInp.txt","w") as fp:
    fp.write(SetINP)
    fp.write(EquationINP)
  SetINP,EquationINP=PBC(m1,byNode=True)
  with open("E:/PBC_Node_EquationInp.txt","w") as fp:
    fp.write(SetINP)
    fp.write(EquationINP)  
  ```
  """
  rasm=m1.rootAssembly
  instanceName=rasm.instances.keys()[0]
  part11=rasm.instances[instanceName]

  x_high,y_high,z_high=part11.cells.getBoundingBox()["high"]
  x_low,y_low,z_low=part11.cells.getBoundingBox()["low"]

  tol=ACISEPS*max(*(x1-x0 for x1,x0 in zip(part11.cells.getBoundingBox()["high"],part11.cells.getBoundingBox()["low"])))

  ns=part11.nodes
  x0ns=ns.getByBoundingBox(xMax=x_low+tol)
  y0ns=ns.getByBoundingBox(yMax=y_low+tol)
  z0ns=ns.getByBoundingBox(zMax=z_low+tol)
  x1ns=ns.getByBoundingBox(xMin=x_high-tol)
  y1ns=ns.getByBoundingBox(yMin=y_high-tol)
  z1ns=ns.getByBoundingBox(zMin=z_high-tol)

  faces=[x0ns,y0ns,z0ns,x1ns,y1ns,z1ns]
  #abbr=['Back','Right','Bottom','Front','Left','Top']
  abbr=["X0","Y0","Z0","X1","Y1","Z1"]
  edges={}
  for i in range(6):
    for j in range(6):
      if j!=i and abs(j-i)!=3:
        edges['%s-%s'%(abbr[i],abbr[j])]=set(n.label for n in faces[i]) & set(n.label for n in faces[j])

  corners=[0 for i in range(8)]
  corners[0]=(edges["X0-Y0"] & edges["X0-Z0"]).pop()
  corners[1]=(edges["X1-Y0"] & edges["X1-Z0"]).pop()
  corners[2]=(edges["X1-Y0"] & edges["X1-Z1"]).pop()
  corners[3]=(edges["X0-Y0"] & edges["X0-Z1"]).pop()
  corners[4]=(edges["X0-Y1"] & edges["X0-Z0"]).pop()
  corners[5]=(edges["X1-Y1"] & edges["X1-Z0"]).pop()
  corners[6]=(edges["X1-Y1"] & edges["X1-Z1"]).pop()
  corners[7]=(edges["X0-Y1"] & edges["X0-Z1"]).pop()

  innerNodes={}
  for i in range(6):
    innerNodes["Inner-%s"%abbr[i]]=set(n.label for n in faces[i]).difference(*edges.values())
  for k in edges.keys():
    edges[k].difference_update(corners)

  EquationPair={direction:{} for direction in 'XYZ'} # 按方向
  for ind,direction in enumerate('XYZ'):
    ns0,ns1=faces[ind],faces[ind+3]
    try:
      for n0 in ns0:
        bounding={}
        xyz=n0.coordinates
        for i in range(3) :
          if i==ind:
            continue
          bounding["%sMin"%'xyz'[i]]=xyz[i]-tol
          bounding["%sMax"%'xyz'[i]]=xyz[i]+tol
        n1=ns1.getByBoundingBox(**bounding)[0]
        EquationPair[direction][n0.label]=n1.label
        EquationPair[direction][n1.label]=n0.label
    except Exception as e:
      print 'Error in find:',bounding
      raise Exception
  
  with open('E:/PBC_Node_Verify.txt','w') as fp:
    for ind,direction in enumerate('XYZ'):
      fp.write(direction)
      for k,v in EquationPair[direction].items():
        c1,c2=ns[k-1].coordinates,ns[v-1].coordinates
        fp.write('%d ->  %d : %f %s -> %s\n'%(k,v,sum([abs(c1[i]-c2[i]) for i in range(3) if i!=ind]),repr(c1),repr(c2)))
  
  for ind,direction in enumerate('XYZ'):
    innerNodes["Inner-%s0"%direction]=list(innerNodes["Inner-%s0"%direction])
    innerNodes["Inner-%s1"%direction]=[EquationPair[direction][n] for n in innerNodes["Inner-%s0"%direction]]
    d1,d2='XYZ'[(ind+1)%3],'XYZ'[(ind+2)%3]
    edges['%s0-%s0'%(d1,d2)]=list(edges['%s0-%s0'%(d1,d2)])
    edges['%s1-%s0'%(d1,d2)]=[EquationPair[d1][n] for n in edges['%s0-%s0'%(d1,d2)]]
    edges['%s0-%s1'%(d1,d2)]=[EquationPair[d2][n] for n in edges['%s0-%s0'%(d1,d2)]]
    edges['%s1-%s1'%(d1,d2)]=[EquationPair[d1][n] for n in edges['%s0-%s1'%(d1,d2)]]
    
    edges['%s0-%s0'%(d2,d1)]=edges['%s0-%s0'%(d1,d2)]
    edges['%s1-%s0'%(d2,d1)]=edges['%s1-%s0'%(d1,d2)]
    edges['%s0-%s1'%(d2,d1)]=edges['%s0-%s1'%(d1,d2)]
    edges['%s1-%s1'%(d2,d1)]=edges['%s1-%s1'%(d1,d2)]
  
  SetINP=""
  EquationINP="*Equation\n"
  dataline=textwrap.dedent("""\
  3
  {0:s},1, {1:g},{2:s},1, {3:g}.,{4:s},1, {5:g}
  3
  {0:s},2, {1:g},{2:s},2, {3:g}.,{4:s},2, {5:g}
  3
  {0:s},3, {1:g},{2:s},3, {3:g}.,{4:s},3, {5:g}
  """)
  RPDict={'X':'RP4','Y':'RP5','Z':'RP6'}
  if byNode:
    for n in set([n.label for f in faces for n in f ]):
      SetINP+="*Nset,nset=Node-%d,instance=%s\n%d\n"%(n,instanceName,n)

    for ind,direction in enumerate('XYZ'):
      f0,f1=innerNodes["Inner-%s0"%direction],innerNodes["Inner-%s1"%direction]
      for i in xrange(len(f0)):
        EquationINP+=dataline.format("Node-%d"%f0[i],1.0,"Node-%d"%f1[i],-1.0,RPDict[direction],1.0)
      d1,d2='XYZ'[(ind+1)%3],'XYZ'[(ind+2)%3]
      # Y0Z0 Y1Z0 Y0Z1 Y1Z1
      e1,e2,e3,e4=[edges["%s%d-%s%d"%(d1,i,d2,j)] for j in range(2) for i in range(2) ]
      for i in xrange(len(e1)):
        EquationINP+=dataline.format("Node-%d"%e1[i],1.0,"Node-%d"%e2[i],-1.0,RPDict[d1],1.0)
        EquationINP+=dataline.format("Node-%d"%e2[i],1.0,"Node-%d"%e4[i],-1.0,RPDict[d2],1.0)
        EquationINP+=dataline.format("Node-%d"%e4[i],-1.0,"Node-%d"%e3[i],1.0,RPDict[d1],1.0)
        # 方程重复了 不删除会报错,nodes are missing degree of freedoms. The MPC/Equation/kinematic coupling constraints can not be formed.
        #EquationINP+=dataline.format("Node-%d"%e3[i],-1.0,"Node-%d"%e1[i],1.0,RPDict[d1],1.0)  
  else:
    # create Equation By Set  
    N=15 # Nset的行数据最多为16
    for k,v in innerNodes.items():
      SetINP=SetINP+"*Nset, nset=%s, instance=%s,UNSORTED\n"%(k,instanceName)
      for i in range(len(v)/N+1):
        SetINP=SetINP+",".join((str(x) for x in v[N*i:min(len(v),N*i+N)]))+'\n'

    for k,v in edges.items():
      SetINP=SetINP+"*Nset, nset=%s, instance=%s,UNSORTED\n"%(k,instanceName)
      for i in range(len(v)/N+1):
        SetINP=SetINP+",".join((str(x) for x in v[N*i:min(len(v),N*i+N)]))+'\n'

    for k in range(8):
      SetINP=SetINP+"*Nset, nset=C%d, instance=%s,UNSORTED\n%d\n"%(k+1,instanceName,corners[k])

    for ind,direction in enumerate('XYZ'):
      EquationINP+=dataline.format("Inner-%s0"%direction,1.0,"Inner-%s1"%direction,-1.0,RPDict[direction],1.0)
      d1,d2='XYZ'[(ind+1)%3],'XYZ'[(ind+2)%3]
      e1,e2,e3,e4=["%s%d-%s%d"%(d1,i,d2,j) for j in range(2) for i in range(2)]  # (0,0) (1,0) (0,1) (1,1)
      EquationINP+=dataline.format(e1,1.0,e2,-1.0,RPDict[d1],1.0)
      EquationINP+=dataline.format(e2,1.0,e4,-1.0,RPDict[d2],1.0)
      EquationINP+=dataline.format(e4,-1.0,e3,1.0,RPDict[d1],1.0)
      # 方程重复了 不删除会报错,nodes are missing degree of freedoms. The MPC/Equation/kinematic coupling constraints can not be formed.
      #EquationINP+=dataline.format(e3,1.0,e1,-1.0,RPDict[d2],-1.0)
  return SetINP,EquationINP

from .Composite import *

def CompactTension():
  """
  Compact Tension 紧凑拉伸实验用于确定材料I型断裂韧性
  
  """
  # ASTM  D5045 − 14
  f1=lambda x: (2.0+x)*(0.886+4.64*x-13.32*x**2+14.72*x**3-5.6*x**4)/(1.0-x)**1.5
  # 沈成康断裂力学 P107
  f2=lambda x: 29.6*x**0.5-185.5*x**1.5+655.7*x**2.5-1017.0*x**3.5+638.9*x**4.5
  # Stress Analysis of the Compact Specimen Including the Effects of Pin Loading
  f3=lambda x: 4.55-40.32*x+414.7*x**2.0-1698*x**3+3781*x**4-4287.0*x**5+2017*x**6 
  for x in xrange(3,7):
    print 0.1*x,f1(0.1*x),f2(0.1*x),f3(0.1*x),f4(0.1*x)

