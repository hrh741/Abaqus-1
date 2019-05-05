#!/usr/bin/python2

from odbAccess import *
import sys
if __name__=="__main__":
    odb=openOdb('Job-1.odb');
    steps=odb.steps;
    laststep=steps[steps.keys()[-1]];
    lastframe=laststep.frames[-1];
    fieldOutputs=lastframe.fieldOutputs
    if sys.argv[1]:
        argv=sys.argv[1]
    else:
        argv='EVOL'
    out=fieldOutputs[argv]
    with open('resulet.txt','w') as fp:
        for f in out.values:
            fp.write("%d\t%f\n"%(f.elementLabel,f.data))

from odbAccess import *
import sys
odb=openOdb('E:\Work\Abaqus/temp/a.odb');
steps=odb.steps;
laststep=steps[steps.keys()[-1]];
for frm in laststep.frames:
    sfo=frm.fieldOutputs['S']
    setnames=[]
    for name in setnames:
        sbfo=sfo.getSubset(region=odb.rootAssembly.elementsets[name])
        A,B
        for v in sbfo.values:
            fprint(v.elemenLabel,v.data) 
