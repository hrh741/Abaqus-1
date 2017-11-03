#! /usr/bin/python2


from odbAccess import *

odb=openOdb('Job-1.odb')
lastStep=odb.steps[odb.steps.keys()[-1]]
lastframe=lastStep.frames[-1]
out=lastframe.fieldOutputs['EVOL']
#for output in out.values:
with open('result.txt','w') as fp:
    fp.write('Index\tMaterial\tEVOL\n')
    for fop in out.values:
        sn=fop.instance.sectionAssignments[0].sectionName
        mt=odb.sections[sn].material
        fp.write('%d\t%s\t%f\n'%(fop.elementLabel,mt,fop.data))

    