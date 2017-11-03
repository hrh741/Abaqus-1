# Abaqus

import os
os.chdir('E:/CodeStudy/Abaqus')
execfile('E:/CodeStudy/Abaqus/3-1.py')

## 名词解释

1. Part

2. feature

2. Section
3. Assemble
4. Mesh
5. Load
6. BC
7. JOB
8. Visualization
9. cell
10. faces
11. edges
12. region
13. 

### Abaqus 模型数据库的基本结构

每个Abaqus CAE 都能创建一个模型数据库，
但是每个模型数据库可以包含多个模型(model)
每个模型只能包含一个装配件(assembly),
每个装配件(assembly)由一个或多个实体(instance)组成,
每个实体(instance)是部件(part)在装配件(assembly)中的映射,
每个部件(part)可以对应一个或多个实体(instance)，
材料和截面属性(section)定义在部件(part)上，
相互作用(interaction)、边界条件(BC boundary condition)、载荷(Load)等定义在实体上，
网格(mesh)可以定义在部件(part)或者实体(instance)上，
对求解过程的控制定义在整个模型(model)上

## Abaqus 基础几何概念

### cell

Cells are volumetric regions of geometry.
cell 是三维区域

### face

### edge
### vertex
### feature
### node
### element

### 对象
#### Cell Object
方法:

1. getSize()
返回cell体积
2. getFaces
3. getEdges
4. getVertices
5. getNodes
6. getElements
7. getAdjacentCells

属性:


### abaqus 中颜色的含义
1. 灰色
> Part(部件)中
1. 蓝色
> Assembly(装配)
1. 橙色
> 无法使用默认的网格划分技术来划分网格
2. 黄色
> 可以使用扫掠划分网格
3. 绿色
> 可以使用结构化网格划分技术


## 算例学习
### 1.4.7 Debonding behavior of a double cantilever beam
