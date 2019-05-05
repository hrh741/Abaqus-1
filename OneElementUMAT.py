#! /usr/python2
# _*_ coding: utf-8 _*_

"""
创建UMAT用于测试
"""
if __name__ == '__main__' and __package__ is None:
    """

    """
    import inspect
    from os import sys, path
    #curpath=path.dirname(path.dirname(path.abspath(__file__)))
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    curpath = os.path.dirname(os.path.abspath(filename))
    pth=os.path.dirname(curpath)
    if pth not in sys.path:
        sys.path.append(pth)

    from Abaqus.ModelLaminate import modelLaminate

    ## 输入参数
    ls=(1,10,2.0) #  层合板长度(x方向) , 单元数目，以及double seed ratio (end/center)
    ws=(10,50,4.0)  #  层合板长度(x方向) , 单元数目，以及double seed ratio (center/end)
    plies=[	(0,     0.125,   "UMAT",   1,	1),
    		(90,     0.125,   "UMAT",   1,	1),
    		(90,     0.125,   "UMAT",   1,	1),
    		(0,     0.125,   "UMAT",   1,	1),] # 从顶部到底部的每层的角度、厚度、材料以及厚度
    halfStructure=[0,0,1] # 每个tuple分别包含是否在x,y,z方向上是否使用半结构 0否 1是
    # 注意宽度方向上的对称条件不知道为什么一直错误 所以宽度方向尽量不要使用半结构
    EquationOrDisplacement=0 #在拉伸边界上使用 参考点Equation约束(0) 或者 位移边界(1)

    ## 默认的model设置
    modelname='CrossPly'
    model1=mdb.Model(name=modelname,description='')

    material1=model1.Material(name='UMAT')
    material1.Depvar(n=15)
    material1.UserMaterial(mechanicalConstants=
        (4100.0, 4100.0, 4100.0, 0.46, 0.46, 0.46, 1404.0, 1404.0, 1404.0,  # Matrix Elastic
        121.0, 210.0, 76.0, # Matrix : tensile compress shear 
        276000.0, 19000.0, 19000.0, 0.2, 0.2, 0.36, 27000.0, 27000.0, 6985.0, # Fiber Elastic
        4850.0, 3000.0, # fiber: tensile compress
        0.575, # Vf
        0.3, 0.3, # alpha,beta
        6.0, # MSEG 基体折线段数目 
        42.7,   53.7,   76.5,   101.5,  111.3,  125,# ETM(1,:)
        4.1E3,  2.5E3,  2E3,    1.4E3,  0.8E3,  0.41E3, # ETM(2,:)
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
        2.32225977 ,5.00655592, 1.65594384 , 1.44730762 , 1.86990463 , 1.68355870 ,1.70323652, 54.4, 0.01 ,1 #KT22,KTT22,KC22,K12,K23,KT22BI,KC22BI, 界面的临界脱粘强度, 衰减系数
        ))

    m,j=modelLaminate(model1,ls,ws,plies,halfStructure,EquationOrDisplacement,ex=0.02)

    m.setValues(description='Longitude: %s\n Width: %s \n Plies: angle\t thickness \t material \t seedsize In height \n %s'%(repr(ls),repr(ws),repr(plies)))
    
    j.setValues(userSubroutine='E:/1730895/Work/BridgeMatrix/src/UMAT_lipeng.for')