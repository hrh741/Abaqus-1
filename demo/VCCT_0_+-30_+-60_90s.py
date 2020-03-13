# -*- coding: utf-8 -*-
import sys
sys.path.append('E:/1730895/Work/')

from  Abaqus.ModelVCCTUEL import modelVCCTUEL
from Abaqus.Materials import Materials
import numpy as np

for i in range(1,10,1):
    ## 输入参数
    ls=(0.05,1,1.0) #  层合板长度(x方向) , 单元数目，以及double seed ratio (end/center)
    ws=(25.0,500,1.0)  #  层合板长度(y方向) , 单元数目，以及double seed ratio (center/end)
    plies=[ (0,     0.135,    "T300-7901",  8,    1.0),
            (30,    0.135,    "T300-7901",  8,    1.0),
            (0,    0.001,    "Epoxy914", 	1,	1.0),
            (0,    0.001,    "Epoxy914", 	1,	1.0),
            (-30,   0.135,    "T300-7901",  8,    1.0),
            (60,    0.135,    "T300-7901",  8,    1.0),
            (-60,   0.135,    "T300-7901",  8,    1.0),                
            (90,    0.135,    "T300-7901",  8,    1.0),
            (90,    0.135,    "T300-7901",  8,    1.0),                
            (-60,   0.135,    "T300-7901",  8,    1.0),                
            (60,    0.135,    "T300-7901",  8,    1.0),
            (-30,   0.135,    "T300-7901",  8,    1.0),
            (0,    0.001,    "Epoxy914", 	1,	1.0),
            (0,    0.001,    "Epoxy914", 	1,	1.0),
            (30,    0.135,    "T300-7901",  8,    1.0),
            (0,     0.135,    "T300-7901",  8,    1.0),
            ] # 从顶部到底部的每层的角度、厚度、材料以及厚度
    
    halfStructure=[0,1,0] # 每个tuple分别包含是否在x,y,z方向上是否使用半结构 0否 1是
    # 注意宽度方向上的对称条件不知道为什么一直错误 所以宽度方向尽量不要使用半结构
    EquationOrDisplacement=0 #在拉伸边界上使用 参考点Equation约束(0) 或者 位移边界(1)
    crack_len=0.05*i
    cracks=[(2,crack_len,1),]
    
    model1=mdb.Model(name="VCCT-11-20-1000a=%d"%(int(1000*crack_len)))
    Materials(model1)     
    job=modelVCCTUEL(model1,ls,ws,plies,halfStructure,EquationOrDisplacement,ex=0.01,cracks=cracks)
    
    job.submit()
    job.waitForCompletion()

import os,glob
Gs=[]
dst=open("E:/result.txt","w")
for fn in glob.glob("./Job-VCCT-11-20-1000*.log"):
    with open(fn,"r") as fp:
        for line in fp.readlines():
            if line.find("GI,GII,GIII") !=-1:
                G=line
        dst.write("%s : %s"%(fn,G))
        Gs.append((fn,[float(x) for x in G.split()[1:]]))

dst.close()