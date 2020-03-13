#!/usr/bin/python2
# -*- coding: UTF-8 -*-

"""
建立等厚度层合板模型

"""
import sys
sys.path.append('E:/1730895/Work/')

from Abaqus import modelLaminate,Materials,UMATMaterial,modelVCCTUEL

ls=(0.5,1,1.0)
ws=(25.0,100,1.0)
plies=[ (0,0.125,'T300-7901',  4,1.0),
        (0,0.001,'Epoxy7901',1,1.0),
        (30,0.125,'T300-7901',	4,1.0),
        (0,0.001,'Epoxy7901',1,1.0),
        (-30,0.125,'T300-7901',  4,1.0),
        (0,0.001,'Epoxy7901',1,1.0),
		(0,0.001,'Epoxy7901',1,1.0),
        (60,0.125,'T300-7901',	4,1.0),
        (0,0.001,'Epoxy7901',1,1.0),
        (-60,0.125,'T300-7901',  4,1.0),
        (0,0.001,'Epoxy7901',1,1.0),
        (90,0.125,'T300-7901',  4,1.0),
        (0,0.001,'Epoxy7901',1,1.0),
        (0,0.001,'Epoxy7901',1,1.0),
        (90,0.125,'T300-7901',  4,1.0),
        (0,0.001,'Epoxy7901',1,1.0),
        (-60,0.125,'T300-7901',  4,1.0),
        (0,0.001,'Epoxy7901',1,1.0),
        (60,0.125,'T300-7901',  4,1.0),
        (0,0.001,'Epoxy7901',1,1.0),
		(0,0.001,'Epoxy7901',1,1.0),
        (-30,0.125,'T300-7901',  4,1.0),
        (0,0.001,'Epoxy7901',1,1.0),
        (30,0.125,'T300-7901',  4,1.0),
        (0,0.001,'Epoxy7901',1,1.0),
        (0,0.125,'T300-7901',  4,1.0),]

"""
halfStructure=[0,0,0]
modelName='Laminate0_+-30_+-60_90s'
model1=mdb.Model(name=modelName)
Materials(model1)
UMATMaterial(model1)
m,job,_=modelLaminate(model1,ls,ws,plies,halfStructure,0,ex=0.02)
"""

halfStructure=[0,1,0]
modelName='EdgeDelaminated0_+-30_+-60_90s_90_90'
model1=mdb.Model(name=modelName)
Materials(model1)
UMATMaterial(model1)
#cracks=[(6,2.0,1)]# 
cracks=[(12,2.0,1)]# 
modelVCCTUEL(model1,ls,ws,plies,halfStructure,0,ex=0.02,cracks=cracks)
