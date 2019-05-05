# -*- coding: utf-8 -*-
from abaqus import *
from abaqusConstants import *
from regio  nToolset import Region

rvc_params=(1.6,1.6,1.538*48)
rvc_len,rvc_wid,rvc_hgt=rvc_params

fibcenter=(rvc_len/2,rvc_wid/2)
fibr=0.5
fibz1,fibz2=0.269*48,1.269*48

directs=['x','y','z']
subds=[(1,1),(2,2),(3,3),(2,3),(1,3),(1,2)]
for i in range(3):
di,dj=subds[i]
modelname='s%d%d'%(di,dj) 
model1=mdb.Model(name=modelname)
sketch1=model1.ConstrainedSketch(name='Sketch-1', sheetSize=200.0)
sketch1.rectangle(point1=(0.0, 0.0),point2=(rvc_len, rvc_wid))
part1=model1.Part(dimensionality=THREE_D, name='Part-1', type=DEFORMABLE_BODY)
part1.BaseSolidExtrude(depth=rvc_hgt, sketch=sketch1)

topface=part1.faces.findAt((rvc_len/2,rvc_wid,rvc_hgt/2))
zedge=part1.edges.findAt((rvc_len,rvc_wid,rvc_hgt/2))
tf1=part1.MakeSketchTransform(origin=(0.0,rvc_wid,rvc_hgt),
    sketchPlane=topface, sketchPlaneSide=SIDE1,
    sketchUpEdge=zedge, sketchOrientation=RIGHT)
sketch2=model1.ConstrainedSketch(gridSpacing=1.0, 
    name='Sketch-2',sheetSize=20.0, transform=tf1)
part1.projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=sketch2)
sketch2.rectangle(point1=(0,fibz1), point2=(rvc_len,fibz2))
part1.PartitionCellBySketch(cells=part1.cells, 
    sketch=sketch2, sketchPlane=topface, sketchUpEdge=zedge)


zedge=part1.edges.findAt((rvc_len,rvc_wid,rvc_hgt/2))
fiberEdge1=part1.edges.findAt((rvc_len/2,rvc_wid,fibz1))
fiberEdge2=part1.edges.findAt((rvc_len/2,rvc_wid,fibz2))
part1.PartitionCellByPlanePointNormal(cells=part1.cells,
    normal=zedge, point=part1.vertices[zedge.getVertices()[0]])

zedge=part1.edges.findAt((rvc_len,rvc_wid,rvc_hgt/2))
part1.PartitionCellByPlanePointNormal(cells=part1.cells,
    normal=zedge, point=part1.vertices[zedge.getVertices()[1]])

## fiber partition 
fiberBtmFace=part1.faces.findAt((rvc_len/2,rvc_wid/2,fibz1))
fiberBtmRight=part1.edges.findAt((0,rvc_wid/2,fibz1))
tf2 = part1.MakeSketchTransform(sketchPlane=fiberBtmFace,
    sketchOrientation=RIGHT,sketchUpEdge=fiberBtmRight,sketchPlaneSide=SIDE1,
    origin=(0.0, 0.0, 2.15))
sketch3=model1.ConstrainedSketch(gridSpacing=1, name='Sketch-3',sheetSize=20.0, transform=tf2)
part1.projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=sketch3)
sketch3.CircleByCenterPerimeter(center=(-rvc_len/2,rvc_wid/2), point1=(-rvc_len/2+fibr,rvc_wid/2))

fiberCells=(part1.cells.findAt((rvc_len/2,rvc_wid/2,(fibz1+fibz2)/2)),)
part1.PartitionCellBySketch(cells=fiberCells, sketch=sketch3,
    sketchPlane=fiberBtmFace,sketchUpEdge=fiberBtmRight)
zedge=part1.edges.findAt((rvc_len,rvc_wid,rvc_hgt/2))
fiberCircle=part1.edges.findAt((rvc_len/2+fibr,rvc_wid/2,fibz1))
part1.PartitionCellByExtrudeEdge(cells=fiberCells, edges=fiberCircle, line=zedge, sense=FORWARD)

# Assembly
fiberCell=part1.cells.findAt((rvc_len/2,rvc_wid/2,(fibz1+fibz2)/2))
fiberLabels=(fiberCell.index,)
matrixLabels=tuple(filter(lambda label:label not in fiberLabels,
    map(lambda cell:cell.index,part1.cells)))

fiberCells=part1.cells[fiberLabels[0]:fiberLabels[0]+1]
matrixCells=part1.cells[matrixLabels[0]:matrixLabels[0]+1]
for i in range(1,len(matrixLabels)):
    l=matrixLabels[i]
    matrixCells+=part1.cells[l:l+1]

model1.Material(name='Material-1').Elastic(table=((1e9,0.38), ))
model1.HomogeneousSolidSection(material='Material-1', name='Section-1')
part1.Set(cells=matrixCells, name='SET-MATRIXCELL')
part1.SectionAssignment(region=part1.sets['SET-MATRIXCELL'], sectionName='Section-1')

model1.Material(name='Material-2').Elastic(table=((3e10, 0.2), ))
model1.HomogeneousSolidSection(material='Material-2', name='Section-2')
part1.Set(cells=fiberCells, name='SET-FIBERCELL')
part1.SectionAssignment(region=part1.sets['SET-FIBERCELL'], sectionName='Section-2',)

# Step
model1.StaticStep(name='Step-1', previous='Initial')
model1.fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 'NFORC', 'CSTRESS', 
    'CDISP','EVOL'))

# Assembly
rootasm=model1.rootAssembly
rootasm.DatumCsysByDefault(CARTESIAN)
inst1=rootasm.Instance(dependent=ON, name='Part-1-1', part=part1)
inst1faces=inst1.faces
# BC
## z=0
z0faces=inst1faces.getByBoundingBox(zMin=0,zMax=0)
## z=rvc_hgt
zhfaces=inst1faces.getByBoundingBox(zMin=rvc_hgt,zMax=rvc_hgt)
## y=0 
y0faces=inst1faces.getByBoundingBox(yMin=0,yMax=0)
## y=rvc_wid
ywfaces=inst1faces.getByBoundingBox(yMin=rvc_wid,yMax=rvc_wid)
# x=0 
x0faces=inst1faces.getByBoundingBox(xMin=0,xMax=0)
# x=rvc_len
xlfaces=inst1faces.getByBoundingBox(xMin=rvc_len,xMax=rvc_len)

sixfaces=[(x0faces,xlfaces),(y0faces,ywfaces),(z0faces,zhfaces)]
symBCs={1:model1.XsymmBC,2:model1.YsymmBC,3:model1.ZsymmBC}

if di==dj:
    model1.DisplacementBC(name='BC-1',createStepName='Step-1',region=Region(faces=sixfaces[di-1][0]),**{'u%d'%(di):0.0})
    model1.DisplacementBC(name='BC-2',createStepName='Step-1',region=Region(faces=sixfaces[di-1][1]),**{'u%d'%(di):rvc_params[di-1]})
    for i in range(1,4):
        if i!=di:
            symbc=symBCs[i]
            symbc(name='BC-%d'%(i+2),createStepName='Step-1',
                region=Region(faces=sixfaces[i-1][0]+sixfaces[i-1][1]))
else:
    if di+dj==5:
        model1.DisplacementBC(name='BC-1',createStepName='Step-1',region=Region(faces=y0faces),u2=0.0)
        model1.DisplacementBC(name='BC-2',createStepName='Step-1',region=Region(faces=x0faces),u2=0.0)
        model1.ZsymmBC(name="BC-5",createStepName="Step-1",region=Region(faces=z0faces+zhfaces))
    elif di+dj==4:
        pass
    else:
        pass            
# Mesh
part1.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=0.1)
part1.generateMesh()
rootasm.regenerate()

# Job
mdb.Job(name=modelname,model=model1)

for jobname in mdb.jobs.keys():
job1=mdb.jobs[jobname]
job1.submit()
job1.waitForCompletion()

# 读取输出数据
from odbAccess import *
import numpy as np

SETMATRIX="SET-1"
SETFIBER="SET-2"
meanSRs=[]
for jobname in mdb.jobs.keys():
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
    print(jobname)
    print('fiber volume:',volumeF,'matrix volume:',volumeM,'total volume:',volumeR)
    print('last volume:',reduce(lambda x,y:x+y,(v.data for v  in evolfo.values)))
    print('fiber Mean:\n',str(meanSF))
    print('matrix Mean:\n',str(meanSM))
    print('RVE Mean:\n',str((sumSF+sumSM)/(volumeM+volumeF)))
meanSRs=np.array(meanSRs)
res=np.linalg.inv(meanSRs[0:3,0:3])
print('E11:',1/res[2,2])
