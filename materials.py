#! /usr/bin/python3
# -*- coding: UTF-8 -*-


"""
材料属性为Abaqus 顺序
断裂韧性单位为KJ/m^2
厚度单位为mm
"""
materials={
    ## 吴庆欣
    'T300-7901':((137.78e3,8.91e3,8.91e3,0.3,0.3,0.48,4.41e3,4.41e3,3.01e3),None,None), #
    'Pipes-Pagano':((20.0, 2.1, 2.1, 0.21, 0.21, 0.21, 0.85, 0.85, 0.85),None,None), #MPsi
    ## Uguen A, Zubillaga L, Turon A, et al. Comparison of cohesive zone models used to predict delamination initiated from free-edges: validation against experimental results[C]//ECCM-16TH European Conference on Composite Materials. 2014.
    ### GIc,GIIc,GIIIc,sigma_zz0,tau_xz0,tau_yz0=0.24,0.74,0.74,46,75,75
    'T800-M21':((130.0e3,8.0e3,8.0e3,0.31,0.31,0.45,4.0e3,4.0e3,4.0e3),None,None), #MPa
    ## Initiation of free-edge delamination in composite laminates
    'G947-M18':((97.6e3, 8.0e3 ,8.0e3,0.37,0.37,0.5, 3.1e3, 3.1e3, 2.7e3),None,None),
    ## Lorriot T, Marion G, Harry R, et al. Onset of free-edge delamination in composite laminates under tensile loading[J]. Composites Part B: Engineering, 2003, 34(5): 459-471.
    'T800-914':((159.0e3,8.4e3,8.4e3,0.31,0.31,0.45,4.0e3,4.0e3,4.0e3),None,None),
    ## 顾嘉杰 毕业论文
    'IM7-914C':((153.6e3,10.2e3,10.2e3,0.27,0.27,0.46,5.7e3,5.7e3,3.5e3 ),None,None),
    ## Diaz A D, Caron J. Prediction of the onset of mode III delamination in carbon-epoxy laminates[J]. Composite Structures, 2006, 72(4): 438-445.
    'T700/CTE1-15':((153.82e3,10.61e3,10.61e3,0.315,0.315,0.315,5.58e3,5.58e3,5.58e3),None,None),
    ## 教研室数据
    'Epoxy7901':((3.17e3, 0.355),None,None),
    ## http://www.matweb.com/search/datasheet.aspx?MatGUID=c3ace011772449e28fa7644ec6c701c8
    'EpoxyM18':((3.5e3,0.38),None,None),
    ## http://www.matweb.com/search/DataSheet.aspx?MatGUID=10199622b9ba408c83347d0dc63bf686
    'Epoxy914':((3.9e3,0.41),None,None),
    ## A Fast Numerical Methodology for Delamination Growth Initiation Simulation
    'AS4/3501-6':((147.0e3,9.0e3,9.0e3,0.33,0.33,0.42,5.0e3,5.0e3,3.0e3),None,{'GI':0.175,'GII':0.532,'GIII':0.532,'t':0.188}),
    ## Role of Matrix Resin in Delamination Onset and Growth in Composite Laminates
    'T300-648':((137e3,9.1e3,9.1e3,0.31,0.31,0.31,5.3e3,5.3e3,5.3e3),None,None),
    'T300-634':((133e3,7.7e3,7.7e3,0.33,0.33,0.33,4.2e3,4.2e3,4.2e3),None,None),
}

def Materials(model1):
    """
    材料库
    """
    from abaqusConstants import ENGINEERING_CONSTANTS,ISOTROPIC
    for matName,props in materials.items():
        if props[0]:
            material1=model1.Material(name=matName)
            if len(props[0])==9:
                material1.Elastic(type=ENGINEERING_CONSTANTS, table=(props[0], ))
            elif len(props[0])==2:
                material1.Elastic(type=ISOTROPIC, table=(props[0], ))

def UMATMaterial(model1):
    material1=model1.Material(name='UMAT-Composite',description='T300-7901')
    material1.Depvar(n=15)
    material1.UserMaterial(mechanicalConstants=
        (3170.0, 3170.0, 3170.0, 0.355, 0.355, 0.355, 1169.74, 1169.74, 1169.74,  # Matrix Elastic
        85.1, 107, 52.6, # Matrix : tensile compress shear 
        230000.0, 15000.0, 15000.0, 0.2, 0.2, 0.07142857, 15000.0, 15000.0, 7000.0, # Fiber Elastic
        3500.0, 2000.0, # fiber: tensile compress
        0.59, # Vf
        0.3, 0.3, # alpha,beta
        1.0, # MSEG 基体折线段数目 
        42.7, # ETM(1,:)
        3170, # ETM(2,:)
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
        #KT22,KTT22,KC22,K12,K23,KT22BI,KC22BI, 界面的临界脱粘强度, 衰减系数
        1.0 ,1.0, 1.0, 1.0 , 1.0, 1.0,1.0, 54.4, 0.01 ,1 ,
    ))
    material1=model1.Material(name='UMAT-Matrix',description='7901')
    material1.Depvar(n=16,deleteVar=16)
    material1.UserMaterial(mechanicalConstants=
        (3170.0, 0.355, # Matrix Elastic
        85.1, 107, 52.6,  # Matrix : tensile compress shear 
        1.0, # MSEG 基体折线段数目 
        42.7, # ETM(1,:)
        3170,  # ETM(2,:)
        0.01, # 衰减系数
        )
    )
    #material1=model1.Material(name='UMAT-Composite',description='IM7-8552')
    #material1.Depvar(n=15)
    #material1.UserMaterial(mechanicalConstants=
    #    (4100.0, 4100.0, 4100.0, 0.46, 0.46, 0.46, 1404.0, 1404.0, 1404.0,  # Matrix Elastic
    #    121.0, 210.0, 76.0, # Matrix : tensile compress shear 
    #    276000.0, 19000.0, 19000.0, 0.2, 0.2, 0.36, 27000.0, 27000.0, 6985.0, # Fiber Elastic
    #    4850.0, 3000.0, # fiber: tensile compress
    #    0.575, # Vf
    #    0.3, 0.3, # alpha,beta
    #    6.0, # MSEG 基体折线段数目 
    #    42.7,   53.7,   76.5,   101.5,  111.3,  125,# ETM(1,:)
    #    4.1E3,  2.5E3,  2E3,    1.4E3,  0.8E3,  0.41E3, # ETM(2,:)
    #    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
    #    2.32225977 ,5.00655592, 1.65594384 , 1.44730762 , 1.86990463 , 1.68355870 ,1.70323652, 54.4, 0.01 ,1 #KT22,KTT22,KC22,K12,K23,KT22BI,KC22BI, 界面的临界脱粘强度, 衰减系数
    #))

    #material1=model1.Material(name='UMAT-Matrix',description='8552')
    #material1.Depvar(n=15)
    #material1.UserMaterial(mechanicalConstants=
    #    (4100.0, 0.46, # Matrix Elastic
    #    121.0, 210.0, 76.0, # Matrix : tensile compress shear 
    #    6.0, # MSEG 基体折线段数目 
    #    42.7,   53.7,   76.5,   101.5,  111.3,  125,# ETM(1,:)
    #    4.1E3,  2.5E3,  2E3,    1.4E3,  0.8E3,  0.41E3, # ETM(2,:)
    #    0.001, # 衰减系数
    #    )
    #)