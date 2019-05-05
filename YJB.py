# -*- coding: utf-8 -*-
from abaqus import *
from abaqusConstants import *
from caeModules import *

rec_len,rec_wid,crack_len=64,60,42;


model1=mdb.Model(name="Model-1");
part1 = model1.Part(name='Part-1', dimensionality=THREE_D, type=DEFORMABLE_BODY)
sketch1 = model1.ConstrainedSketch(name='Sketch-1', sheetSize=200.0)
sketch1.setPrimaryObject(option=STANDALONE)
sketch1.rectangle(point1=(0.0, rec_wid), point2=(rec_len,0))
shell1=part1.BaseShell(sketch=sketch1)
sketch1.unsetPrimaryObject()

t = part1.MakeSketchTransform(sketchPlane=part1.faces[0], sketchUpEdge=part1.edges.findAt((rec_len,rec_wid/2,0)),
    sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0, 0, 0.0))
sketch2=model1.ConstrainedSketch(name='Sketch-2', sheetSize=175.45, gridSpacing=4.38, transform=t)
#sketch2.setPrimaryObject(option=SUPERIMPOSE)
part1.projectReferencesOntoSketch(sketch=sketch2,filter=COPLANAR_EDGES)
sketch2.Line(point1=(0.0, rec_wid/2), point2=(crack_len, rec_wid/2))
sketch2.Line(point1=(crack_len, rec_wid/2), point2=(crack_len, rec_wid))
sketch2.Line(point1=(crack_len, rec_wid/2), point2=(crack_len, 0))

pickedFaces = part1.faces.getSequenceFromMask(mask=('[#1 ]', ), )
part1.PartitionFaceBySketch(faces=part1.faces.findAt((0,0,0)),sketch=sketch2)
sketch2.unsetPrimaryObject()

mdb.models['Model-1'].Material(name='Material-1').Elastic(
    type=ENGINEERING_CONSTANTS, table=((76599.5, 76599.5, 9701.0, 0.0302, 
    0.3222, 0.3222, 4820.0, 4100.0, 4100.0), ))
mdb.models['Model-1'].HomogeneousShellSection(name='Section-1', 
    preIntegrate=OFF, material='Material-1', thicknessType=UNIFORM, 
    thickness=1.0, thicknessField='', idealization=NO_IDEALIZATION, 
    poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, 
    useDensity=OFF, integrationRule=SIMPSON, numIntPts=5)

faces = part1.faces.getSequenceFromMask(mask=('[#7 ]', ), )
region = part1.Set(faces=faces, name='SA_SET-1')
part1.SectionAssignment(region=region, sectionName='Section-1', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)

a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
a.Instance(name='Part-1-1', part=part1, dependent=OFF)
mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial')
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 'NFORC', 'CSTRESS', 
    'CDISP'))

e1 = a.instances['Part-1-1'].edges
edges1 = e1.getSequenceFromMask(mask=('[#8 ]', ), )
pickedRegions = a.Set(edges=edges1, name='SEAM_Set-1')
mdb.models['Model-1'].rootAssembly.engineeringFeatures.assignSeam(
    regions=pickedRegions)

partInstances =(a.instances['Part-1-1'], )
a.seedPartInstance(regions=partInstances, size=1.0, deviationFactor=0.1, 
    minSizeFactor=0.1)
elemType1 = mesh.ElemType(elemCode=S4R, elemLibrary=STANDARD, 
    secondOrderAccuracy=OFF, hourglassControl=DEFAULT)
elemType2 = mesh.ElemType(elemCode=S3, elemLibrary=STANDARD)

f1 = a.instances['Part-1-1'].faces
faces1 = f1.getSequenceFromMask(mask=('[#7 ]', ), )
pickedRegions =(faces1, )
a.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))

a.generateMesh(regions=partInstances)
a = mdb.models['Model-1'].rootAssembly
n1 = a.instances['Part-1-1'].nodes
nodes1 = n1.getSequenceFromMask(mask=('[#0:36 #200000 ]', ), )
a.Set(nodes=nodes1, name='Set-2')
#: The set 'Set-2' has been created (1 node).
a = mdb.models['Model-1'].rootAssembly
n1 = a.instances['Part-1-1'].nodes
nodes1 = n1.getSequenceFromMask(mask=('[#0:73 #10000000 ]', ), )
a.Set(nodes=nodes1, name='Set-3')
#: The set 'Set-3' has been created (1 node).
a = mdb.models['Model-1'].rootAssembly
region = a.sets['Set-2']
mdb.models['Model-1'].ConcentratedForce(name='Load-1', createStepName='Step-1', 
    region=region, cf2=1.0, distributionType=UNIFORM, field='', localCsys=None)
a = mdb.models['Model-1'].rootAssembly
region = a.sets['Set-3']
mdb.models['Model-1'].DisplacementBC(name='BC-1', createStepName='Initial', 
    region=region, u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET, 
    amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
mdb.Job(name='Job-1', model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, 
    numGPUs=0)
mdb.jobs['Job-1'].submit(consistencyChecking=OFF)
#: The job input file "Job-1.inp" has been submitted for analysis.
#: Job Job-1: Analysis Input File Processor completed successfully.
#: Job Job-1: Abaqus/Standard completed successfully.
#: Job Job-1 completed successfully. 
o3 = session.openOdb(name='Job-1.odb')
#: Model: D:/SIMULIA/Temp/Job-1.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       3
#: Number of Node Sets:          5
#: Number of Steps:              1
mdb.saveAs(pathName='E:/29.cae')
#: The model database has been saved to "D:\SIMULIA\Temp\29.cae".
