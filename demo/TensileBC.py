# -*- coding: utf-8 -*-
import sys
sys.path.append('E:/1730895/Work/')

from  Abaqus.ModelLaminate import modelLaminate
from Abaqus.Materials import Materials
import numpy as np

for modelSuf in ['X1_U1']:
    ls=(180,20,5.0)
    ws=(25,40,2.0)

    modelName='45_45-%s'%modelSuf
    plies=[(45,0.125,'Pipes-Pagano',2,1.0),
        (-45,0.125,'Pipes-Pagano',  2,1.0),
        (-45,0.125,'Pipes-Pagano',  2,1.0),
        (45,0.125,'Pipes-Pagano',   2,1.0),]
    model1=mdb.Model(name=modelName)
    Materials(model1)
    modelLaminate(model1,ls,ws,plies,[0,0,0],1,ex=0.01)