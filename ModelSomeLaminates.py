#!/usr/bin/python2
# -*- coding: UTF-8 -*-

"""
建立等厚度层合板模型

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

from Abaqus.ModelLaminate import modelLaminate

## 输入参数
ls=(40.0,40,8.0) #  层合板长度(x方向) , 单元数目，以及double seed ratio (end/center)
ws=(8.0,40,8.0)  #  层合板长度(y方向) , 单元数目，以及double seed ratio (center/end)
plies=[ (45,    0.5,    "Pipes-Pagano", 4),
		(-45,   0.5,    "Pipes-Pagano", 4),
		(-45,   0.5,    "Pipes-Pagano", 4),
		(45,    0.5,    "Pipes-Pagano", 4)] # 从顶部到底部的每层的角度、厚度、材料以及厚度
halfStructure=[0,0,1] # 每个tuple分别包含是否在x,y,z方向上是否使用半结构 0否 1是
# 注意宽度方向上的对称条件不知道为什么一直错误 所以宽度方向尽量不要使用半结构
EquationOrDisplacement=0 #在拉伸边界上使用 参考点Equation约束(0) 或者 位移边界(1)

# x,y,z方向上使用对称性条件建模半结构简化计算 
halfs=[(x,y,z) for x in (0,1) for y in (0,1) for z in (0,1)] 	
for halfStructure in halfs:	
	## 默认的model设置
	suffix="".join(("XYZ"[i] if halfStructure[i]==1 else "" for i in xrange(3)))
	suffix=suffix+("" if len(suffix)==0 else "SYMM")

	modelname='Model-%s-Div-%d-%s'%(suffix,xyzdiv[2],'Eq' if EquationOrDisplacement==0 else 'Disp')
	model1=mdb.Model(name=modelname,description='')

	## Material
	material2=model1.Material(name='Pipes-Pagano',description="1970 Pipes & Pagano\n Psi")
	material2.Elastic(type=ENGINEERING_CONSTANTS, table=((20000000.0, 2100000.0, 2100000.0, 0.21, 0.21, 0.21, 850000.0, 850000.0, 850000.0), ))

	m,j=modelLaminate(model1,ls,ws,plies,halfStructure,EquationOrDisplacement)
	
	#job.submit()
	#print jobname,"Submited"
	#job.waitForCompletion()
	#print jobname,"Complete"

## 奇怪的bug 无法在job completed 后进行数据提取
for halfStructure in halfs:	
	## 生成自定义参数 不要修改
	lam_len,lam_wid=l/(2.0 if halfStructure[0]==1 else 1),w/(2.0 if halfStructure[1]==1 else 1)
	plynum=len(plies)//(2 if halfStructure[2]==1 else 1)
	plyAngles=[plies[i][0] for i in xrange(plynum)]
	plyThicks=[plies[i][1] for i in xrange(plynum)]
	
	lam_hgt=sum(plyThicks)
	plyprefix=",".join(["%02.0f"%x for x in plyAngles])
	plyLowZs=np.cumsum(plyThicks,dtype='f')
	plyLowZs=plyLowZs[-1]-plyLowZs	
	
	## 默认的model设置
	suffix="".join(("XYZ"[i] if halfStructure[i]==1 else "" for i in xrange(3)))
	suffix=suffix+("" if len(suffix)==0 else "SYMM")
	modelname='Model-%s-Div-%d-%s'%(suffix,xyzdiv[2],'Eq' if EquationOrDisplacement==0 else 'Disp')
	model1=mdb.models[modelname]
	p1=model1.rootAssembly.instances['Part-1-1']
	jobname="Job-%s-%s-ZDIV-%d-%s"%(plyprefix,suffix,xyzdiv[2],'Eq' if EquationOrDisplacement==0 else 'Disp')
	
	try:
		ns=p1.nodes.getByBoundingBox(xMin=lam_len/2-0.1,xMax=lam_len/2+0.1,zMin=plyLowZs[0]-0.01,zMax=plyLowZs[0]+0.01)
		sort_ns=sorted(ns,cmp=lambda x,y:1 if x.coordinates[1]>y.coordinates[1] else -1)
		session.Path(name=jobname,type=NODE_LIST,expression=[(p1.name.upper(),[n.label,]) for n in sort_ns])
			
		o=session.openOdb(name="%s.odb"%jobname,readOnly=True)
		session.viewports['Viewport: 1'].setValues(displayedObject=o) # 此处如果在complete后执行可能会报错
		scratchOdb = session.ScratchOdb(o)
		scratchOdb.rootAssembly.DatumCsysByThreePoints(name='CSYS-1', 
		  coordSysType=CARTESIAN, origin=(0.0, 0.0, 0.0), point1=(1.0, 0.0, 0.0), 
		  point2=(0.0, 1.0, 0.0))
		dtm = scratchOdb.rootAssembly.datumCsyses['CSYS-1']
		session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(
		  transformationType=visualization.USER_SPECIFIED, datumCsys=dtm)

		for sigma in ('S11','S22','S33','S23','S13','S12'):
		  session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
			variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, sigma))
		  pth = session.paths[jobname]
		  session.XYDataFromPath(name=sigma+suffix, path=pth, 
			includeIntersections=False, projectOntoMesh=False, pathStyle=PATH_POINTS, 
			numIntervals=10, projectionTolerance=0, shape=UNDEFORMED, 
			labelType=TRUE_DISTANCE_Y)
	except e:
		print e
