#! /usr/bin/python2
#-*-coding:utf-8 -*-
from abaqus import *
from abaqusConstants import *
import math
import  regionToolset

"""
建立减层结构模型

# 层合板长度和宽度
laminate_l,laminate_w=140,40
# 每个元组包含为角度和厚度
Plies=[(90,1),(45,1),(0,1),(0,1),(-45,1),
    (0,1),(45,1),(0,1),(-45,1),(0,1),
    (0,1),(-45,1),(0,1),(45,1),(0,1),
    (-45,1),(0,1),(0,1),(45,1),(90,1),]

# 每个元组分别表示在减层开始的x坐标、第几层(在所有层的下标)，减连续的多少层，以及过渡区长度
# 注意 按照第一个元素从小到大排序 且 两个减层x区间不存在交叉
Drops=[
    (50,3,1,20),
    (70,9,1,20),
    (90,5,1,20),
    (110,6,1,20)]

"""
## 建模开始
mymodel=mdb.Model(name="Model-2")

## 输入参数

#### From : Delamination in asymmetrically tapered composites loaded in tension
#laminate_l,laminate_w=240,12.5 # 240,25
#Plies=[
#    (-45,0.125),(-45,0.125),(45,0.125),(45,0.125),(0,0.125),(0,0.125),(0,0.125),(0,0.125),
#    (0,0.125),(0,0.125),(0,0.125),(0,0.125),(45,0.125),(45,0.125),(-45,0.125),(-45,0.125),
#    (-45,0.125),(-45,0.125),(45,0.125),(45,0.125),(0,0.125),(0,0.125),(0,0.125),(0,0.125),
#    (0,0.125),(0,0.125),(0,0.125),(0,0.125),(45,0.125),(45,0.125),(-45,0.125),(-45,0.125),]
#Drops=[(120,9,8,2.5)] # 

#laminate_l,laminate_w=40,8 # 240,25
#Plies=[(45,1),(-45,1)]
#Drops=[] # []代表无递减层 也就是普通等厚度层合板

### From : Delamination in asymmetrically tapered composites loaded in tension
laminate_l,laminate_w=240,12.5 # 240,25
Plies=[
    (-45,0.125,"Composite"),(-45,0.125,"Composite"),(45,0.125,"Composite"),
    (45,0.125,"Composite"),(0,0.125,"Composite"),(0,0.125,"Composite"),
    (0,0.125,"Composite"),(0,0.125,"Composite"),(0,0.125,"Composite"),
    (0,0.125,"Composite"),(0,0.125,"Composite"),(0,0.125,"Composite"),
    (45,0.125,"Composite"),(45,0.125,"Composite"),(-45,0.125,"Composite"),
    (-45,0.125,"Composite"),(-45,0.125,"Composite"),(-45,0.125,"Composite"),
    (45,0.125,"Composite"),(45,0.125,"Composite"),(0,0.125,"Composite"),
    (0,0.125,"Composite"),(0,0.125,"Composite"),(0,0.125,"Composite"),
    (0,0.125,"Composite"),(0,0.125,"Composite"),(0,0.125,"Composite"),
    (0,0.125,"Composite"),(45,0.125,"Composite"),(45,0.125,"Composite"),
    (-45,0.125,"Composite"),(-45,0.125,"Composite"),]
Drops=[(120,9,8,2.5)] # 

## From : Strain-Energy-Release Rate Analysis of Delamination in a Tapered Laminate Subjected to Tension Load
h=0.125
laminate_l,laminate_w=180*h,20*h # 240,25
Plies=[
    (0,0.875,"Composite"),
    (45,0.125,"Composite"),(-45,0.125,"Composite"),
    (0,0.01,"Resin"),
    (45,0.125,"Composite"),(-45,0.125,"Composite"),
    (45,0.125,"Composite"),(-45,0.125,"Composite"),
    (45,0.125,"Composite"),(-45,0.125,"Composite"),
    (0,0.01,"Resin"),
    (0,0.125,"Composite"),
    (45,0.125,"Composite"),(-45,0.125,"Composite"),
    (0,0.125,"Composite"),]
Drops=[(60*h,5,2,20*h),
        (80*h,7,2,20*h),
        (100*h,9,2,20*h),] # X起始坐标,第i层开始,中断层数,过渡区长度


mymodel.Material(name='Composite',description="S2/SP250 GlasdEpoxy").Elastic(
    type=ENGINEERING_CONSTANTS, table=((7.3, 2.1, 2.1, 0.275, 
    0.275, 0.275, 0.88, 0.88, 0.88), ))
mymodel.Material(name='Resin').Elastic(table=((0.59, 0.33), ))
mymodel.HomogeneousSolidSection(name='Resin_Pocket', material='RESIN', thickness=None)


###  create part_1
part1=mymodel.Part(name='Part-1',dimensionality=THREE_D,type=DEFORMABLE_BODY)

sketch1=mymodel.ConstrainedSketch(name='Sketch-1',sheetSize=200.0)
point1=(0,0)
thick_h=sum([h for _,h,_ in Plies])
point2=(0,thick_h)
sketch1.Line(point1=point1,point2=point2)
for drop in Drops:
    point1=point2
    point2=(drop[0],point1[1])
    sketch1.Line(point1=point1,point2=point2)
    point1=point2
    x,y=point1
    y=y-sum([Plies[i][1] for i in range(drop[1]-1,drop[1]+drop[2]-1)])
    point2=(x+drop[3],y)
    sketch1.Line(point1=point1,point2=point2)
sketch1.Line(point1=point2,point2=(laminate_l,point2[1]))
sketch1.Line(point1=(laminate_l,point2[1]),point2=(laminate_l,0))
sketch1.Line(point1=(0,0),point2=(laminate_l,0))
part1.BaseSolidExtrude(sketch=sketch1,depth=laminate_w) 

## 建立单层
f=part1.faces.findAt((Drops[0][1]/2.0 if len(Drops)!=0 else laminate_l/2.0,thick_h/2.0,0))
t = part1.MakeSketchTransform(sketchPlane=f, sketchUpEdge=part1.edges.findAt((0,thick_h/2.0,0)), 
    sketchPlaneSide=SIDE2, origin=(0,0,0))
s= mymodel.ConstrainedSketch(name='cut',sheetSize=200,transform=t)
part1.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)

points=[] # 每个元素为 四边形底部横线中点坐标和所在层 构成的tuple
drop_points=[] # ply drop竖线 与  树脂聚集区域底部横线 的 中点
pocket_corner=[] # 树脂聚集区角点坐标
n_plys,n_drops=len(Plies),len(Drops)
for ply_ind in range(1,n_plys):
    # 在s上绘制第ply_ind层底部线段以及中间的竖线
    x=0
    y=thick_h-sum([Plies[i][1] for i in range(ply_ind)])
    ply_h=Plies[ply_ind-1][1]
    dropped=False # 本层是否被dropped
    for drop_i,drop in enumerate(Drops):
        # Drop 从第st个层开始到第ed 层结束（包括第ed层）
        st,ed=drop[1],drop[1]+drop[2]-1
        drop_h=sum([Plies[i][1] for i in range(st-1,ed)])
        drop_l=drop[3]
        if(ply_ind<=ed):
            # 注意判断顺序： 一定要先判断ed
            if(ply_ind==ed):
                # 注意这段是被dropped的Ply的最底层
                point2=(drop[0]+drop_l,y)
                drop_points.append((drop[0]+drop_l/2.0,y))
                pocket_corner.append((drop[0],y))
                dropped=True
            elif(ply_ind>=st and ply_ind<ed): 
                # 被dropped的ply 中间层
                point2=(drop[0],y)
                dropped=True
            else:
                point2=(drop[0],y)
            
            s.Line(point1=(x,y),point2=point2)
            if abs(x-drop[0])>1e-6:
                points.append(((x+drop[0])/2.0,y,ply_ind-1,drop_i)) 
                # 沿着厚度方向的竖线
                s.Line(point1=(drop[0],y),point2=(drop[0],y+ply_h))
                drop_points.append((drop[0],y+ply_h/2.0))            

            x,y=point2
            if(dropped):
                break
            point2=(drop[0]+drop_l,y-drop_h)
            s.Line(point1=(x,y),point2=point2)
            points.append(((x+point2[0])/2.0,(y+point2[1])/2.0,ply_ind-1))
            x,y=point2
            # 沿着厚度方向的竖线
            s.Line(point1=(x,y),point2=(x,y+ply_h))
            drop_points.append((x,y+ply_h/2.0))            

    if(not dropped):
        # 一般来说 不会有 树脂聚集区端部在层合板端部
        s.Line(point1=(x,y),point2=(laminate_l,y))
        points.append(((x+laminate_l)/2.0,y,ply_ind-1))
part1.PartitionFaceBySketch(faces=f,sketch=s)

for point in points:
    direction=part1.edges.findAt((0,0,laminate_w/2.0))
    es=part1.edges.findAt(coordinates=(point[0],point[1],0))
    #part1.PartitionCellByExtrudeEdge(edges=es, cells=part1.cells, line=direction,sense=FORWARD)
    part1.PartitionCellBySweepEdge(sweepPath=direction, cells=part1.cells,edges=(es,))

for point in drop_points:
    direction=part1.edges.findAt((0,0,laminate_w/2.0))
    es=part1.edges.findAt(coordinates=(point[0],point[1],0))
    #part1.PartitionCellByExtrudeEdge(edges=es, cells=part1.cells, line=direction,sense=FORWARD)
    try:
        part1.PartitionCellBySweepEdge(sweepPath=direction, cells=part1.cells,edges=(es,))
    except:
        pass

datum={}
points.append([laminate_l/2.0,0,len(Plies)-1])
for i,point in enumerate(points):
    angle,ply_h,prop=Plies[point[2]]
    angle=int(angle)
    es=part1.edges.findAt(coordinates=(point[0],point[1],0))
    n1,n2=es.getVertices()
    v1,v2=part1.vertices[n1],part1.vertices[n2]
    o,x=(v1,v2) if v1.pointOn[0][0]<v2.pointOn[0][0] else (v2,v1)
    #slope=(x.pointOn[0][1]-o.pointOn[0][1])/(x.pointOn[0][0]-o.pointOn[0][0])
    #if(abs(slope)<1e-6):
    #    dt=None
    #else:
    #    key='x=%.6f,tan=%.6f'%(point[0],slope)
    #    if not datum.has_key(key):
    #        f=part1.DatumCsysByThreePoints(name='PLY_%d_%d'%(point[2],i), coordSysType=CARTESIAN, 
    #        origin=o, point1=x, point2=(o.pointOn[0][0], o.pointOn[0][1], -1.0))
    #        datum[key]=part1.datums[f.id]
    f=part1.DatumCsysByThreePoints(name='PLY_%d_%d'%(point[2],i), coordSysType=CARTESIAN, 
        origin=o, point1=x, point2=(o.pointOn[0][0], o.pointOn[0][1], -1.0))    
    dt=part1.datums[f.id]
    name="PLY-%2d-%2d"%(point[2]+1,i)
    r=regionToolset.Region(cells=part1.cells.findAt(((point[0],point[1]+ply_h/2.0,0.0),)))
    layup=part1.CompositeLayup(
        name=name, description='PLY %d'%(point[2]+1), elementType=SOLID, 
        symmetric=False, thicknessAssignment=FROM_SECTION)
    layup.ReferenceOrientation(orientationType=SYSTEM, 
        localCsys=dt, fieldName='', 
        additionalRotationType=ROTATION_NONE, angle=0.0, 
        additionalRotationField='', axis=AXIS_3, stackDirection=STACK_3)
    layup.CompositePly(suppressed=False, plyName='Ply-%d'%i, region=r, 
        material=prop, thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
        orientationType=SPECIFY_ORIENT, orientationValue=angle, 
        additionalRotationType=ROTATION_NONE, additionalRotationField='', 
        axis=AXIS_3, angle=0.0, numIntPoints=1)


for i,corner in enumerate(pocket_corner):
    corner_x,corner_y=corner
    st=drop[1]
    ed=drop[1]+drop[2]-1
    drop_h=sum([Plies[k][1] for k in range(st-1,ed)])
    drop_l=drop[3]
    
    region = part1.Set(cells=part1.cells.findAt(((corner_x+drop_l/4.0,corner_y+drop_h/4.0,0),)), name='Pocket-%d'%i)
    part1.SectionAssignment(region=region, sectionName='Resin_Pocket', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)

for drop in Drops:
    dt=part1.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=drop[0])
    try:
        part1.PartitionCellByDatumPlane(datumPlane=part1.datums[dt.id], cells=part1.cells)
    except Exception as e:
        print e
    dt=part1.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=drop[0]+drop[3])
    try:
        part1.PartitionCellByDatumPlane(datumPlane=part1.datums[dt.id], cells=part1.cells)
    except Exception as e:
        print e

# Assembly
rootasm=mymodel.rootAssembly
rootasm.DatumCsysByDefault(CARTESIAN)
inst1=rootasm.Instance(dependent=ON, name='Part-1-1', part=part1)

# Step
mymodel.StaticStep(name='Step-1', previous='Initial',timeIncrementationMethod=AUTOMATIC, initialInc=1.0)
mymodel.fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 'NFORC', 'CSTRESS', 
    'CDISP','EVOL','SDV'))
