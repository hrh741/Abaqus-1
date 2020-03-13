# -*- coding: utf-8 -*-
import sys
sys.path.append('E:/1730895/Work/')

from  Abaqus.ModelLaminate import modelLaminate
from Abaqus.Materials import Materials
import numpy as np

ls=(0.1,1,1.0) #  层合板长度(x方向) , 单元数目，以及double seed ratio (end/center)
ws=(25.0,500,1.0)  #  层合板宽度(y方向) , 单元数目，以及double seed ratio (center/end)

halfStructure=[0,0,0] # 每个tuple分别包含是否在x,y,z方向上是否使用半结构 0否 1是
EquationOrDisplacement=0 #在拉伸边界上使用 参考点Equation约束(0) 或者 位移边界(1)

for plyAngles in [
    [0,30,-30,60,-60,90],
    [30,-30,60,-60,0,90],
    [0,90,0,90],]:
    plyThicks=0.135*np.ones_like(plyAngles)
    plies=[(plyAngles[i],plyThicks[i],'T300-7901',5,1.0) for i in range(len(plyAngles))]

    modelName='YeLin-'+'_'.join(['%d'%x for x in plyAngles])+'_Pipes'
    model1=mdb.Model(name=modelName)
    Materials(model1)
    modelLaminate(model1,(0.1,1,1.0),(25,500,1.0),
                        plies,halfStructure,0)

    modelName='YeLin-'+'_'.join(['%d'%x for x in plyAngles])+'_3D'
    model1=mdb.Model(name=modelName)
    Materials(model1)
    modelLaminate(model1,(180,20,5.0),(25,25,2.0),
                    plies,halfStructure,1)