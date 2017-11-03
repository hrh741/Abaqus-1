from abaqus import *
from abaqusConstants import *


if __name__=='__main__':
    model1=mdb.Model(name="Model-1")
    part1=model1.Part(name="Part-1",dimensionality=THREE_D,type=DEFOMABLE_BODY)
    sketch1=model1.ConstrainedSketch(name='Shell-Sketch',sheetSize=200.0)
    sketch1.rectangle(point1=(0,0),point2=(40,20))
    part1.BaseShell(sketch=sketch1)

    
    t800_917=model1.Material(name="T800/917")
    t800_917.Elastic(type=ENGINEERING_CONSTANTS,table=((139e9,8.4e9,8.4e9,
    0.33,0.5,0.5,4.1,4.1,4.1)))

