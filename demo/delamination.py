# -*- coding: UTF-8 -*-

from __future__ import print_function

import sys
PackagePath='E:/1730895/Work/'
if PackagePath not in sys.path:
    sys.path.append('E:/1730895/Work/')

from abaqus import *
from abaqusConstants import *
from Abaqus import modelLaminate,modelDCB,UMATMaterial,Materials,Laminate,Compliance

"""
ls=(180,20,4.0)
ws=(25,25,1.0)

plies=[ (0,0.125,'UMAT-Composite',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (30,0.125,'UMAT-Composite',2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (-30,0.125,'UMAT-Composite',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (60,0.125,'UMAT-Composite',2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (-60,0.125,'UMAT-Composite',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (90,0.125,'UMAT-Composite',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (90,0.125,'UMAT-Composite',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (-60,0.125,'UMAT-Composite',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (60,0.125,'UMAT-Composite',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (-30,0.125,'UMAT-Composite',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (30,0.125,'UMAT-Composite',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (0,0.125,'UMAT-Composite',  2,1.0),]

halfStructure=[0,0,1]
modelName='Delamination_I'

model1=mdb.Model(name=modelName)
UMATMaterial(model1)
m,job,_=modelLaminate(model1,ls,ws,plies,[0,0,0],1,ex=0.02)

ls=(1,1,1.0)
ws=(1,1,1.0)

plies=[ (30,    0.125,  'T300-7901',   2,1.0),
        (0,     0.001,  'UMAT-Matrix',      1,1.0),
        (-30,   0.125,  'T300-7901',   2,1.0),
        (0,     0.001,  'UMAT-Matrix',      1,1.0),
        (-30,   0.125,  'T300-7901',   2,1.0),
        (0,     0.001,  'UMAT-Matrix',      1,1.0),
        (30,    0.125,  'T300-7901',   2,1.0),]

plies=[ (0,0.125,'T300-7901',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (30,0.125,'T300-7901',2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (-30,0.125,'T300-7901',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (60,0.125,'T300-7901',2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (-60,0.125,'T300-7901',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (90,0.125,'T300-7901',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (90,0.125,'T300-7901',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (-60,0.125,'T300-7901',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (60,0.125,'T300-7901',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (-30,0.125,'T300-7901',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (30,0.125,'T300-7901',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (0,0.125,'T300-7901',  2,1.0),]

halfStructure=[0,0,0]
modelName='testProgressive0_+-30_+-60_90s'
model1=mdb.Model(name=modelName)
Materials(model1)
UMATMaterial(model1)
m,job,_=modelLaminate(model1,ls,ws,plies,halfStructure,1,ex=0.02)

ls=(1,1,1.0)
ws=(1,1,1.0)
plies=[(0,1,'UMAT-Matrix',1,1.0),]
halfStructure=[0,0,0]
modelName='OneElement'
model1=mdb.Model(name=modelName)
Materials(model1)
UMATMaterial(model1)
m,job,_=modelLaminate(model1,ls,ws,plies,halfStructure,1,ex=0.1)

plies=[ (0,0.125,'T300-7901',  10,1.0),
        (30,0.125,'T300-7901',10,1.0),
        (-30,0.125,'T300-7901',  10,1.0),
        (60,0.125,'T300-7901',10,1.0),
        (-60,0.125,'T300-7901',  10,1.0),
        (90,0.125,'T300-7901',  10,1.0),
        (90,0.125,'T300-7901',  10,1.0),
        (-60,0.125,'T300-7901',  10,1.0),
        (60,0.125,'T300-7901',  10,1.0),
        (-30,0.125,'T300-7901',  10,1.0),
        (30,0.125,'T300-7901',  10,1.0),
        (0,0.125,'T300-7901',  10,1.0),]
plies=[ (30,    0.125,  'T300-7901',   2,1.0),
        (0,     0.001,  'UMAT-Matrix',      1,1.0),
        (-30,   0.125,  'T300-7901',   2,1.0),
        (0,     0.001,  'UMAT-Matrix',      1,1.0),
        (-30,   0.125,  'T300-7901',   2,1.0),
        (0,     0.001,  'UMAT-Matrix',      1,1.0),
        (30,    0.125,  'T300-7901',   2,1.0),]

##################################################################

ls=(180,20,1.0)
ws=(25,25,1.0)


plies=[ (0,0.125,'UMAT-Composite',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (30,0.125,'UMAT-Composite', 2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (-30,0.125,'UMAT-Composite',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (60,0.125,'UMAT-Composite', 2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (-60,0.125,'UMAT-Composite',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (90,0.125,'UMAT-Composite',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (90,0.125,'UMAT-Composite',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (-60,0.125,'UMAT-Composite',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (60,0.125,'UMAT-Composite',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (-30,0.125,'UMAT-Composite',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (30,0.125,'UMAT-Composite',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (0,0.125,'UMAT-Composite',  2,1.0),]

halfStructure=[0,0,0]
modelName='Delamiantion3D_I-1mm'
model1=mdb.Model(name=modelName)
Materials(model1)
UMATMaterial(model1)
m,job,_=modelLaminate(model1,ls,ws,plies,halfStructure,1,ex=0.02)

##################################################################
ls=(0.5,1,1.0)
ws=(25,25,1.0)


plies=[ (0,0.125,'UMAT-Composite',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (30,0.125,'UMAT-Composite', 2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (-30,0.125,'UMAT-Composite',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (60,0.125,'UMAT-Composite', 2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (-60,0.125,'UMAT-Composite',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (90,0.125,'UMAT-Composite',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (90,0.125,'UMAT-Composite',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (-60,0.125,'UMAT-Composite',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (60,0.125,'UMAT-Composite',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (-30,0.125,'UMAT-Composite',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (30,0.125,'UMAT-Composite',  2,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (0,0.125,'UMAT-Composite',  2,1.0),]

halfStructure=[0,0,0]
modelName='Delamiantion_I-1mm'
model1=mdb.Model(name=modelName)
Materials(model1)
UMATMaterial(model1)
m,job,_=modelLaminate(model1,ls,ws,plies,halfStructure,0,ex=0.02)
##################################################################
ls=(0.1,1,1.0)
ws=(25,200,5.0)


plies=[ (0,0.125,'UMAT-Composite',  10,1.0),
        (30,0.125,'UMAT-Composite', 10,1.0),
        (-30,0.125,'UMAT-Composite',  10,1.0),
        (60,0.125,'UMAT-Composite', 10,1.0),
        (-60,0.125,'UMAT-Composite',  10,1.0),
        (90,0.125,'UMAT-Composite',  10,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (0,0.001,'UMAT-Matrix',1,1.0),
        (90,0.125,'UMAT-Composite',  10,1.0),
        (-60,0.125,'UMAT-Composite',  10,1.0),
        (60,0.125,'UMAT-Composite',  10,1.0),
        (-30,0.125,'UMAT-Composite',  10,1.0),
        (30,0.125,'UMAT-Composite',  10,1.0),
        (0,0.125,'UMAT-Composite',  10,1.0),]

halfStructure=[0,0,0]
modelName='Delamination_I-Only90-0_125mm'
model1=mdb.Model(name=modelName)
Materials(model1)
UMATMaterial(model1)
m,job,_=modelLaminate(model1,ls,ws,plies,halfStructure,0,ex=0.02)
"""
ls=(10,10,1.0)
ws=(20,40,1.0)


plies=[ (30,    0.125,  'T300-7901',    2,  1.0),
        (-30,   0.125,  'T300-7901',    2,  1.0),
        (-30,   0.125,  'T300-7901',    2,  1.0),
        (30,    0.125,  'T300-7901',    2,  1.0),]

halfStructure=[0,0,1]
modelName='Angle30'
model1=mdb.Model(name=modelName)
Materials(model1)
UMATMaterial(model1)
mat1=model1.Material(name='Cohesive-T300-7901')
mat1.Elastic(type=TRACTION, table=((1e4, 1e4, 1e4), ))
mat1.QuadsDamageInitiation(table=((10, 50.0, 50.0), ))
#m1.quadsDamageInitiation.DamageEvolution(type=DISPLACEMENT, table=((0.01, ), ))
mat1.quadsDamageInitiation.DamageEvolution(type=ENERGY, mixedModeBehavior=BK, power=2.0,table=((0.3, 1.0, 1.0), ))

m,job,_=modelLaminate(model1,ls,ws,plies,halfStructure,1,ex=0.01)
