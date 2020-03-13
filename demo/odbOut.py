from regionToolset import Region
from odbAccess import *
import numpy as np
from abaqusConstants import *

a,b,h=10,10,100
r,h1,h2=2,25,75
rd=0.5
seedsize=1
SETMATRIX="SET-1"
SETFIBER="SET-2"

m1=mdb.models['Model-1']
p1=m1.parts['Part-1']
alles=p1.elements
x,y=1.85,1.85
p1.Set(name="WITHFIBER",elements=alles.getByBoundingCylinder(center1=(x,y,0),center2=(x,y,2.8),radius=1.0))
p1.Set(name="WITHOUTFIBER",elements=alles.getByBoundingCylinder(center1=(x,y,0),center2=(x,y,2),radius=1.0))
#p1.Set(name="Set-FES",elements=fes)
matelements=p1.sets['Set-1'].elements
fes=p1.sets["Set-2"].elements
mes=p1.cells.findAt((a/2.0+r+rd/2.0,b/2.0,(h1+h2)/2.0)).getElements()
mes+=p1.cells.findAt((a/2.0+r+rd/2.0,b/2.0,0)).getElements()
mes+=p1.cells.findAt((a/2.0+r+rd/2.0,b/2.0,h)).getElements()
p1.Set(name="Set-MES",elements=mes)
p1.Set(name="Set-Line",nodes=ns.getByBoundingBox(xMin=2.75,xMax=2.75,yMin=1.75,yMax=1.75))




jobname="TEST6"
job=mdb.Job(name=jobname,model="Model-1")
job.submit()
job.waitForCompletion()


meanSRs=[]
if job.status==COMPLETED:
    getdata=lambda v:v.data
    odb=openOdb(jobname+".odb")
    inst=odb.rootAssembly.instances['PART-1-1']
    # 纤维和基体的初始体积
    firststep=odb.steps[odb.steps.keys()[0]]
    firstframe=firststep.frames[0]
    evolfo=firstframe.fieldOutputs['EVOL']
    evolF=evolfo.getSubset(region=inst.elementSets[SETFIBER])
    volumeF=reduce(lambda x,y:x+y,(v.data for v in evolF.values))
    evolM=evolfo.getSubset(region=inst.elementSets[SETMATRIX])
    volumeM=reduce(lambda x,y:x+y,(v.data for v in evolM.values))
    volumeR=volumeM+volumeF
    # 纤维和基体的体平均应力
    laststep=odb.steps[odb.steps.keys()[-1]]
    lastframe=laststep.frames[-1]
    sfo=lastframe.fieldOutputs['S']
    evolfo=lastframe.fieldOutputs['EVOL']
    sigmaF=sfo.getSubset(region=inst.elementSets[SETFIBER])
    evolF=evolfo.getSubset(region=inst.elementSets[SETFIBER])
    sumSF=reduce(lambda x,y:x+y,(s.data*v.data for s,v in zip(sigmaF.values,evolF.values)))
    meanSF=sumSF/volumeF
    sigmaM=sfo.getSubset(region=inst.elementSets[SETMATRIX])
    evolM=evolfo.getSubset(region=inst.elementSets[SETMATRIX])
    sumSM=reduce(lambda x,y:x+y,(s.data*v.data for s,v in zip(sigmaM.values,evolM.values)))
    meanSM=sumSM/volumeM
    meanSRs.append((sumSF+sumSM)/(volumeM+volumeF))
    print jobname
    print 'fiber volume:',volumeF,'matrix volume:',volumeM,'total volume:',volumeR
    print 'last volume:',reduce(lambda x,y:x+y,(v.data for v  in evolfo.values))
    print 'fiber Mean:\n',str(meanSF).strip('\n')
    print 'matrix Mean:\n',str(meanSM).strip('\n')
    print 'RVE Mean:\n',str((sumSF+sumSM)/(volumeM+volumeF)).strip('\n')
    sigma=np.array(map(getdata,lastframe.fieldOutputs['S'].getSubset(region=inst.elementSets["SET-MES"]).values))
    evol=np.array(map(getdata,lastframe.fieldOutputs['EVOL'].getSubset(region=inst.elementSets["SET-MES"]).values))
    print 'wall Mean:\n',np.dot(evol,sigma)/np.sum(evol)
    sigma=np.array(map(getdata,lastframe.fieldOutputs['S'].getSubset(region=inst.elementSets["WITHFIBER"]).values))
    evol=np.array(map(getdata,lastframe.fieldOutputs['EVOL'].getSubset(region=inst.elementSets["WITHFIBER"]).values))
    print 'WITHFIBER Mean:\n',np.dot(evol,sigma)/np.sum(evol)
    sigma=np.array(map(getdata,lastframe.fieldOutputs['S'].getSubset(region=inst.elementSets["WITHOUTFIBER"]).values))
    evol=np.array(map(getdata,lastframe.fieldOutputs['EVOL'].getSubset(region=inst.elementSets["WITHOUTFIBER"]).values))
    print 'WITHOUTFIBER Mean:\n',np.dot(evol,sigma)/np.sum(evol)

from regionToolset import Region
from odbAccess import *
import numpy as np
from abaqusConstants import *
def calMean(jobname,setname=None):
    getdata=lambda v:v.data
    odb=openOdb(jobname+".odb")
    laststep=odb.steps[odb.steps.keys()[-1]]
    lastframe=laststep.frames[-1]    
    inst=odb.rootAssembly.instances['PART-1-1']
    if setname:
        sigma=np.array(map(getdata,lastframe.fieldOutputs['S'].getSubset(region=inst.elementSets[setname]).values))
        evol=np.array(map(getdata,lastframe.fieldOutputs['EVOL'].getSubset(region=inst.elementSets[setname]).values))
    else:
        setname="ALL"
        sigma=np.array(map(getdata,lastframe.fieldOutputs['S'].values))
        evol=np.array(map(getdata,lastframe.fieldOutputs['EVOL'].values))
    closeOdb(odb)
    print '%s Mean:\n'%(setname),np.dot(evol,sigma)/np.sum(evol)


m1=mdb.models['Model-1']
p1=m1.parts['Part-1']
alles=p1.elements
p1.Set(name="WITHFIBER",elements=alles.getByBoundingCylinder(center1=(1.75,1.75,0),center2=(1.75,1.75,2.6),radius=1.04))
print "WITHFIBER",len(p1.sets["WITHFIBER"].elements)
p1.Set(name="WITHOUTFIBER",elements=alles.getByBoundingCylinder(center1=(1.75,1.75,0),center2=(1.75,1.75,2),radius=1.04))
print "WITHOUTFIBER",len(p1.sets["WITHOUTFIBER"].elements)

jobname="Job-1"
job1=mdb.Job(name=jobname,model="Model-1")
job1.submit()
job1.waitForCompletion()

calMean(jobname,"WITHFIBER")
calMean(jobname,"WITHOUTFIBER")
calMean(jobname,"SET-1")
calMean(jobname,"SET-2")
calMean(jobname)


