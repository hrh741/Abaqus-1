from abaqus import *
from abaqusConstants import *
import math
import part
import sketch
  
if __name__=='__main__':
  mymodel=mdb.Model(name="Model-1")
  
  # create part_1
  part1=mymodel.Part(name='Part-1',dimensionality=THREE_D,type=DEFORMABLE_BODY)
  sketch1=mymodel.ConstrainedSketch(name='Sketch-1',sheetSize=200.0)
  
  Height=40
  linesxy=((6,Height-8),(6,Height),(0,Height),(0,0),(40,0),(40,8),(12,8))
  for i in range(len(linesxy)-1):
    sketch1.Line(point1=linesxy[i],point2=linesxy[i+1])
  sketch1.ArcByStartEndTangent(point1=linesxy[0],point2=(8,22),vector=(0,-1))
  sketch1.Line(point1=(8,22),point2=(8,12))
  
  # doc says the ArcByCenterEnds create arc from point1 to point2 clockwisely
  # but my result shows that the arc is created  counterclockwisely
  # what a fuck
  sketch1.ArcByCenterEnds(center=(8+4,8+4),point2=linesxy[-1],point1=(8,12))

  part1.BaseSolidExtrude(sketch=sketch1,depth=15)

  e1=part1.edges.findAt(coordinates=(3,40,0))
  part1.Round(radius=5,edgeList=(e1,))

  f1=part1.faces.findAt((20,8,7.5))
  e1=part1.edges.findAt(coordinates=(40,8,7.5))
  tf=part1.MakeSketchTransform(sketchPlane=f1,
    origin=(linesxy[-2]+(15,)),sketchOrientation=LEFT,
    sketchUpEdge=e1)
  sketch2=mymodel.ConstrainedSketch(name='Sketch-2',sheetSize=100,gridSpacing=1,transform=tf)
  sketch2.CircleByCenterPerimeter(center=(10,0), point1=(5,0))
  # 
  #part1.projectReferencesOntoSketch(sketch=sketch2,filter=COPLANAR_EDGES)
  part1.CutExtrude(sketchPlane=f1, sketchUpEdge=e1, sketchPlaneSide=SIDE1, 
        sketchOrientation=LEFT, sketch=sketch2, flipExtrudeDirection=OFF)

  # create Material
  steel=mymodel.Material(name='Steel')
  steel.Elastic(table=((21e4,0.3),))

  mymodel.HomogeneousSolidSection(name='Section-1',material='Steel',thickness=None)


  # region is a set of sequence object 
  # each sequence is a collection of edge , face, cell
  part1.SectionAssignment(region=(part1.cells[0:1],),sectionName='Section-1')

  rootAssembly=mymodel.rootAssembly
  rootAssembly.DatumCsysByDefault(CARTESIAN)
  rootAssembly.Instance(name='Part-1-1',part=part1,dependent=ON)

  # 
  e1=part1.edges.findAt(coordinates=(18,8,15))
  part1.PartitionCellByPlanePointNormal(point=part1.InterestingPoint(edge=e1,rule=MIDDLE),normal=e1,cells=part1.cells)

  e1=part1.edges.findAt(coordinates=(6,Height-3,15))
  v1=part1.vertices.findAt(coordinates=linesxy[0]+(15,))
  part1.PartitionCellByPlanePointNormal(point=v1,normal=e1,cells=part1.cells)

  part1.seedPart(size=2.0)
  sqrt2=math.sqrt(2)
  part1.seedEdgeByNumber(edges=part1.edges.findAt(((12-2*sqrt2,12-2*sqrt2,15),)),
  number=12,constraint=FINER)

  import mesh
  part1.setMeshControls(regions=part1.cells,algorithm=MEDIAL_AXIS)
  elmtype1=mesh.ElemType(elemCode=C3D20,elemLibrary=STANDARD)
  elmtype2=mesh.ElemType(elemCode=C3D15,elemLibrary=STANDARD)
  elmtype3=mesh.ElemType(elemCode=C3D10,elemLibrary=STANDARD)
  part1.setElementType(regions=(part1.cells,),elemTypes=(elmtype1,elmtype2,elmtype3))
  part1.generateMesh() 

  # set Load
  stepLoad1=mymodel.StaticStep(name="Step-Load-1",previous="Initial",initialInc=0.2)
  rootAssembly.ReferencePoint(point=(30,20,15))
  rp=rootAssembly.referencePoints.findAt((30,20,15))
  setPoint=rootAssembly.Set(referencePoints=(rp,),name='Set-Point')
  rootFaces=rootAssembly.instances['Part-1-1'].faces
  surfHole=rootAssembly.Surface(side1Faces=rootFaces.findAt(((30,8,10),)),name="Surf-Hole")

  mymodel.Coupling(name='Constraint-Hole', controlPoint=setPoint, 
    surface=surfHole, influenceRadius=WHOLE_SURFACE, couplingType=DISTRIBUTING, 
    weightingMethod=UNIFORM, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, 
    ur2=ON, ur3=ON)
  mymodel.TabularAmplitude(name='Amp-1', timeSpan=STEP, 
    smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (0.2, 1.5), (0.4, 2.0), (1.0, 1.0)))
  mymodel.ConcentratedForce(name='Load-1', 
    createStepName='Step-Load-1', region=setPoint, cf2=-1000.0, 
    amplitude='Amp-1', distributionType=UNIFORM, field='', localCsys=None)

  stepLoad2=mymodel.StaticStep(name='Step-Load-2',previous='Step-Load-1')
  surfLoad=rootAssembly.Surface(side1Faces=rootFaces.findAt(((40,4,5),)),name='Surf-Load')
  mymodel.SurfaceTraction(name='Load-2', 
    createStepName='Step-Load-2', region=surfLoad, magnitude=36.0, 
    directionVector=((0.0, 0.0, 0.0), (0.0, -10.0, 0.0)), 
    distributionType=UNIFORM, field='', localCsys=None)

  # set BC
  setFix=rootAssembly.Set(faces=rootFaces.findAt(((0,10,10),),((0,Height-2,10),)),name='Set-Fix')
  mymodel.EncastreBC(name='BC-Fix', createStepName='Initial', 
    region=setFix, localCsys=None)

  wantCentoid=tuple(map(lambda f:f.getCentroid()[0],filter(lambda f:f.getCentroid()[0][2]==15,rootFaces)))
  setSym=rootAssembly.Set(faces=rootFaces.findAt(coordinates=wantCentoid),name='Set-Sym') 
  mymodel.ZsymmBC(name='BC-Sym', createStepName='Initial', 
    region=setSym, localCsys=None)

  mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 'CSTRESS', 'CDISP', 
    'EVOL'), region=MODEL, exteriorOnly=OFF, sectionPoints=DEFAULT, 
    rebar=EXCLUDE)

  mdb.Job(name='Job-1', model='Model-1',
      description='', type=ANALYSIS, 
      atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
      memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
      explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
      modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
      scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, 
      numGPUs=0)
  mdb.jobs['Job-1'].submit(consistencyChecking=OFF)
    