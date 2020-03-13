# -*- coding: utf-8 -*-
if __name__ == '__main__' and __package__ is None:
    """
    添加当前模块的路径
    """
    import inspect
    from os import sys, path
    #curpath=path.dirname(path.dirname(path.abspath(__file__)))
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    curpath = os.path.dirname(os.path.abspath(filename))
    pth=os.path.dirname(curpath)
    if pth not in sys.path:
        sys.path.append(pth)

from  Abaqus.ModelCrackLaminate import modelCrackLaminate
from Abaqus.Materials import Materials
import numpy as np

ls=(0.1,1,1.0) #  层合板长度(x方向) , 单元数目，以及double seed ratio (end/center)
ws=(25.0,500,1.0)  #  层合板宽度(y方向) , 单元数目，以及double seed ratio (center/end)

halfStructure=[0,1,1] # 每个tuple分别包含是否在x,y,z方向上是否使用半结构 0否 1是
EquationOrDisplacement=0 #在拉伸边界上使用 参考点Equation约束(0) 或者 位移边界(1)

n=1
ind=1
import time
modelNameFormat='Lorrit25-%d'

for n in [1,2,3,4,5]:
    plyAngles=np.array([25,-25,90,0,0,90,-25,25])
    plyThicks=np.array([0.125,0.125,0.0625,0.001,0.001,0.0625,0.125,0.125])
    plies=[((plyAngles[i],plyThicks[i]*n,'T800-914',5*n,1.0) if (plyAngles[i]!=0) else (plyAngles[i],plyThicks[i],'Epoxy914',1,1.0))
            for i in range(len(plyAngles))]
    for crack_len in [0,]:
        modelName=modelNameFormat%ind
        ind=ind+1

        cracks=[(3,crack_len,1),]

        model1=mdb.Model(name=modelName)
        ## Material
        Materials(model1)
        modelCrackLaminate(model1,ls,ws,plies,halfStructure,EquationOrDisplacement,ex=0.01,cracks=cracks)
for j in mdb.jobs.values():
    if j.status is None:
        j.submit()
        j.waitForCompletion()

from odbAccess import openOdb

ind=1
fp=open('E:/lorrit25-1.txt','w')
for n in [1,2,3,4,5]:
    plyAngles=np.array([25,-25,90,0,0,90,-25,25])
    plyThicks=np.array([0.125,0.125,0.0625,0.001,0.001,0.0625,0.125,0.125])
    plies=[((plyAngles[i],plyThicks[i]*n,'T800-914',5*n,1.0) if (plyAngles[i]!=0) else (plyAngles[i],plyThicks[i],'Epoxy914',1,1.0))
            for i in range(len(plyAngles))]
    for crack_len in [0,]:
        modelName=modelNameFormat%ind
        jobname='Job-%s'%modelName
        print modelName
        ind=ind+1
        model1=mdb.models[modelName]
        part11=model1.rootAssembly.instances['Part-1-1']
        ymax=12.5-crack_len
        elms=part11.sets['Ply-3'].elements.getByBoundingBox(yMin=ymax-0.5,yMax=ymax)
        if(len(elms)!=10):
            print 'there is some error in ',modelName
        o=openOdb('Job-%s.odb'%modelName,readOnly=True)
        sfo=o.steps.values()[-1].frames[-1].fieldOutputs['S']
        d={e.label: None for e in elms}
        cnt=0
        for value in sfo.values:
            if value.elementLabel in d and d[value.elementLabel] is None:
                d[value.elementLabel]=value.data
                cnt=cnt+1
                if(cnt==len(d)):
                    break
        RF1=o.steps.values()[-1].historyRegions['Node PART-2-1.1'].historyOutputs['RF1'].data[1][1]
        fp.write('\t'.join([
            modelName,'[±%d_%d]s'%(theta,n),
            '%.2f'%crack_len,
            '%.6f'%RF1,
            #'\t'.join(['%d\t%s'%(k,v) for k,v in d.items()]).replace('\n',''),
            repr(np.mean(np.array(d.values()),axis=0))[7:-7].replace(',','\t')
        ]))
        fp.write('\n')

fp.close()