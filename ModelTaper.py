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
## 输入参数
#### From : Delamination in asymmetrically tapered composites loaded in tension
#laminate_l,laminate_w=240,12.5 # 240,25
#Plies=[
#    (-45,0.125),(-45,0.125),(45,0.125),(45,0.125),(0,0.125),(0,0.125),(0,0.125),(0,0.125),
#    (0,0.125),(0,0.125),(0,0.125),(0,0.125),(45,0.125),(45,0.125),(-45,0.125),(-45,0.125),
#    (-45,0.125),(-45,0.125),(45,0.125),(45,0.125),(0,0.125),(0,0.125),(0,0.125),(0,0.125),
#    (0,0.125),(0,0.125),(0,0.125),(0,0.125),(45,0.125),(45,0.125),(-45,0.125),(-45,0.125),]
#Drops=[(120,9,8,2.5)] # 

laminate_l,laminate_w=40,8 # 240,25
Plies=[(45,1),(-45,1)]
Drops=[] # []代表无递减层 也就是普通等厚度层合板

## 建模开始
mymodel=mdb.Model(name="Model-1")

###  create part_1
part1=mymodel.Part(name='Part-1',dimensionality=THREE_D,type=DEFORMABLE_BODY)

sketch1=mymodel.ConstrainedSketch(name='Sketch-1',sheetSize=200.0)
point1=(0,0)
thick_h=sum([h for _,h in Plies])
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


# UGENS_PROPS=(0.,1.,1.,6.,0.,0.,6.,7.,-1.,1.,0.,0.,0.,0.,0.,0.,0.3,0.3,276000.,0.2,19000.,0.36,27000.,4850.,3000.,0.,0.,0.,0.46,121.,210.,0.,0.,0.,42.7,53.7,76.5,101.5,111.3,125.,4105.8,2500.,2000.,1396.6,803.3,410.2,45.6,88.,128.4,153.2,173.7,193.4,4108.1,2395.5,1603.2,1097.3,699.7,310.2,0.575,-0.5,0.5,1.,0.,15.,100.,100.,100.,121.,210.,121.,210.,76.,12.,12.,12.,12.,12.,0.01,64.,5.083,1.447,)
# UMAT_PROPS=(4100.,  4100.,  4100.,   0.46,   0.46,   0.46,  1400.,  1400.,  1400.,   121.,   210.,    76.,276000., 19000., 19000.,    0.2,    0.2,   0.36, 27000., 27000.,  6980.,  4850.,  3000.,  0.575,    0.3,    0.3,     6.,   42.7,   53.7,   76.5,  101.5,  111.3,   125.,  4100.,  2500.,  2000.,  1400.,   800.,   410.,     1.,     2.,     3.,     4.,     5.,     6.,     7.,     8.,     9.,    10.,    11.,    12.,     0.,  2.322,  5.083,  1.656,  1.447,   1.87,   54.4,   0.01）
UMAT_NAME='IM7-8552'
UMAT_PROPS = (#基体的九个弹性常数加上三个拉压剪模量 E11,E22,E33,V12,V13,V23,G12,G13,G23
            4100., 4100., 4100., 0.46, 0.46,   0.46,  1400.,  1400., 1400.,121.,210.,76.,
            #纤维的九个弹性常数加上两个拉压模量 E11,E22,E33,V12,V13,V23,G12,G13,G23
            276000.,19000., 19000., 0.2, 0.2, 0.36, 27000.,27000., 6980., 4850., 3000.,
            # VF ALPHA BETA
            0.575,  0.3,    0.3,
            #MSEG
            6., 
            #ETM 基体塑性阶段折线强度和模量
            42.7,   53.7,   76.5,  101.5,  111.3,   125.,  
            4100.,  2500.,  2000.,  1400.,   800.,   410.,   
            #STATEV(13): SSF(6) SSM(6) 衰减次数
            0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,        
            #KT22  KTT22 KC22 K12 K23 KT22BI KC22BI
            2.34,5.083,1.656,1.447,1.87,1.72,1.69,            
            #2.322,  5.083,  1.656,  1.447,1.87,   
            #LMISE RF INDFAIL
            54.4,   0.01,  1
            )
UMAT_NAME='T300-914C'
UMAT_PROPS = (#基体的九个弹性常数加上三个拉压剪模量 E11,E22,E33,V12,V13,V23,G12,G13,G23
            4000.0, 4000.0, 4000.0, 0.35, 0.35, 0.35, 1481.0, 1481.0, 1481.0, 75.0,150.0, 70.0,
            #纤维的九个弹性常数加上两个拉压模量 E11,E22,E33,V12,V13,V23,G12,G13,G23
            230000.0, 15000.0, 15000.0, 0.2, 0.2, 0.07, 15000.0, 15000.0,7447.0, 2500.0, 2000.0,
            # VF ALPHA BETA
            0.6, 0.3, 0.3, 
            #MSEG
            6.0, 
            #ETM 基体塑性阶段折线强度和模量
            110.0, 120.0, 130.0, 140.0,150.0, 160.0, 
            4000.0, 4000.0, 4000.0, 4000.0, 4000.0, 4000.0,
            #STATEV(13): SSF(6) SSM(6) 衰减次数
            0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,        
            #KT22  KTT22 KC22 K12 K23 KT22BI KC22BI
            2.14, 5.07, 1.57,1.43, 1.73,1.75943,1.60124,
            #LMISE RF INDFAIL
            60,   0.01,  1
            )
usrmat=mymodel.Material(name=UMAT_NAME, description='UMAT')
usrmat.UserMaterial(mechanicalConstants=UMAT_PROPS)
# 设置nstav值
usrmat.Depvar(n=14)

#mymodel.Material(name='IM7-8552').Elastic(
#    type=ENGINEERING_CONSTANTS, table=((76599.5, 76599.5, 9701.0, 0.0302, 
#    0.3222, 0.3222, 4820.0, 4100.0, 4100.0), ))
mymodel.Material(name='RESIN').Elastic(table=((4000.0, 0.35), ))
mymodel.HomogeneousSolidSection(name='Resin_Pocket', material='RESIN', thickness=None)

datum={}
points.append([laminate_l/2.0,0,len(Plies)-1])
for i,point in enumerate(points):
    angle,ply_h=Plies[point[2]]
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
        material=UMAT_NAME, thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
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
    part1.PartitionCellByDatumPlane(datumPlane=part1.datums[dt.id], cells=part1.cells)
    dt=part1.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=drop[0]+drop[3])
    part1.PartitionCellByDatumPlane(datumPlane=part1.datums[dt.id], cells=part1.cells)

# Assembly
rootasm=mymodel.rootAssembly
rootasm.DatumCsysByDefault(CARTESIAN)
inst1=rootasm.Instance(dependent=ON, name='Part-1-1', part=part1)

# Step
mymodel.StaticStep(name='Step-1', previous='Initial',timeIncrementationMethod=FIXED, initialInc=0.01)
mymodel.fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 'NFORC', 'CSTRESS', 
    'CDISP','EVOL','SDV'))
