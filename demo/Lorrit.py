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

from  Abaqus import modelCrackLaminate,Materials
import numpy as np

ls=(0.1,1,1.0) #  层合板长度(x方向) , 单元数目，以及double seed ratio (end/center)
ws=(25.0,500,1.0)  #  层合板宽度(y方向) , 单元数目，以及double seed ratio (center/end)

halfStructure=[0,1,1] # 每个tuple分别包含是否在x,y,z方向上是否使用半结构 0否 1是
EquationOrDisplacement=0 #在拉伸边界上使用 参考点Equation约束(0) 或者 位移边界(1)

n=1
ind=120
import time
modelNameFormat='Lorrit-%d'
for H in [0.001,0.0005,0.0001]:
    for theta in [20,]:
        plyAngles=np.array([theta,0,0,-theta,0,0,-theta,0,0,theta])
        plyThicks=0.125*n*np.ones_like(plyAngles)
        plies=[((plyAngles[i],plyThicks[i],'T800-914',5,1.0) if (plyAngles[i]!=0) else (plyAngles[i],H,'Epoxy914',1,1.0))
                for i in range(len(plyAngles))]
        for crack_len in [0,]:
            modelName=modelNameFormat%ind
            ind=ind+1

            cracks=[(1,crack_len,1  ),]

            model1=mdb.Model(name=modelName)
            ## Material
            Materials(model1)
            modelCrackLaminate(model1,ls,ws,plies,halfStructure,EquationOrDisplacement,ex=0.01,cracks=cracks)

for j in mdb.jobs.values():
    if j.status is None:
        j.submit()
        j.waitForCompletion()

from odbAccess import openOdb
ind=120
fp=open('E:/lorrit2.txt','w')
for H in [0.01,0.005,0.001]:
    for theta in [20,]:
        plyAngles=np.array([theta,0,0,-theta,0,0,-theta,0,0,theta])
        plyThicks=0.125*n*np.ones_like(plyAngles)
        plies=[((plyAngles[i],plyThicks[i],'T800-914',5,1.0) if (plyAngles[i]!=0) else (plyAngles[i],H,'Epoxy914',1,1.0))
                for i in range(len(plyAngles))]
        for crack_len in [0,]:
            modelName=modelNameFormat%ind
            jobname='Job-%s'%modelName
            print modelName
            ind=ind+1
            model1=mdb.models[modelName]
            part11=model1.rootAssembly.instances['Part-1-1']
            ymax=12.5-crack_len
            elms=part11.sets['Ply-1'].elements.getByBoundingBox(yMin=ymax-0.25,yMax=ymax)
            if(len(elms)!=5):
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