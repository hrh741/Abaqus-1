from abaqus import *
from abaqusConstants import *

if __name__=="__main__":
    model1=mdb.Model(name="Model-1")
    part1=model1.Part(name="Part-1",type=DEFORMABLE_BODY,dimensionality=THREE_D)

    sketch1=model1.ConstrainedSketch(name="Sketch-1",sheetSize=200.0)
    sketch1.rectangle(point1=(0,0),point2=(40,20))

    part1.BaseShell(sketch=sketch1)

    import material

    t800_914=model1.Material(name="T800-914")
    GPa=1e9
    t800_914.Elastic(type=ENGINEERING_CONSTANTS,table=((139*GPa,8.4*GPa,8.4*GPa,0.33,0.5,0.5, 4.1*GPa,4.1*GPa,4.1*GPa),))

    


