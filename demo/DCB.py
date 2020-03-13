#! -*- coding: UTF-8 -*-
from __future__ import print_function

from abaqus import mdb

import sys
sys.path.append('E:/1730895/Work/')
from abaqusConstants import *
from Abaqus import Materials,UMATMaterial,modelDCB_ZSYMM,modelDCBResin

if __name__=="__main__":
    """
    crack_a=50
    ls=(150.0,{2:[5,1.0],50:[100,1.0],98:[10,1.0]}) #  层合板长度(x方向) , 单元数目，以及double seed ratio (end/center)
    ws=[25.0,5,1.0]  #  层合板宽度(y方向) , 单元数目，以及double seed ratio (center/end)
    plies=[ [0,    1.5,    "T300-7901", 	40,	1.0],
            (0,    1.5,    "T300-7901", 	40,	1.0),]
    zsymm=True
    modelName='DCB'
    
    seeds=[16,400,20,5,50]
    modelName='DCB-%d-%d-%d-%d-%d'%(seeds[0],seeds[1],seeds[2],seeds[3],seeds[4])
    print(modelName)
    ls[1][2][0]=seeds[0]
    ls[1][50][0]=seeds[1]
    ls[1][98][0]=seeds[2]
    ws[1]=seeds[3]
    plies[0][3]=seeds[4]

    model1=mdb.Model(name=modelName)
    model1.setValues(description='Longitude: %s\n Width: %s \n Plies: angle\t thickness \t material \t seedsize In height \n %s'%(repr(ls),repr(ws),repr(plies)))
    Materials(model1)

    j=modelDCB(model1,ls,crack_a,ws,plies,zsymm=zsymm)
    """
    a0=50.0
    crack_a=50.0
    rPlastic=60.0
    Npermm=2 # 每mm网格数

    ls=(150.0,{150:[150,1.0]})#  层合板长度(x方向) , 单元数目，以及double seed ratio (end/center)
    ls[1][a0]=[int(Npermm*a0),1.0]
    ls[1][crack_a-a0]=[int(Npermm*(crack_a-a0)),1.0]
    ls[1][rPlastic]=[int(Npermm*rPlastic),1.0]
    ls[1][150-crack_a-rPlastic]=[int(150-crack_a-rPlastic),1.0]
   
    # 对DCB样件宽度方向建议>=10
    ws=(25.0,10,1.0)  #  层合板宽度(y方向) , 单元数目，以及double seed ratio (center/end)
    plies=[ 
        (0,   0.75,    "T300-7901", 	10,	1.0),
        #(0,   0.001,    "UMAT-Matrix", 	10,	1.0),
        (90,   0.75,    "T300-7901", 	10,	1.0),
        #(0,   1.5,    "T300-7901", 	5,	1.0),
        ]
    plies=[ 
        (0,   0.125,    "T300-7901", 	1,	1.0),
        (90,   0.125,    "T300-7901", 	1,	1.0),
        (0,   0.125,    "T300-7901", 	1,	1.0),
        (90,   0.125,    "T300-7901", 	1,	1.0),
        (0,   0.125,    "T300-7901", 	1,	1.0),
        (90,   0.125,    "T300-7901", 	1,	1.0),
        (0,   0.125,    "T300-7901", 	1,	1.0),
        (90,   0.125,    "T300-7901", 	1,	1.0),
        (0,   0.125,    "T300-7901", 	1,	1.0),
        (90,   0.125,    "T300-7901", 	1,	1.0),
        (0,   0.125,    "T300-7901", 	1,	1.0),
        (90,   0.125,    "T300-7901", 	1,	1.0),
        ]

    resinName='UMAT-Matrix'#'Epoxy7901'
    hresin=0.02
    
    delta=5 #张开位移
    
    modelName='DCB90A-Z-N%d-Kc183-U10'%(Npermm)
    model1=mdb.Model(name=modelName)
    model1.setValues(description='Longitude: %s\n Width: %s \n Plies: angle\t thickness \t material \t seedsize In height \n %s'%(repr(ls),repr(ws),repr(plies)))
    # 材料
    Materials(model1)
    UMATMaterial(model1)
    mat1=model1.Material(name='Cohesive-T300-7901')
    mat1.Elastic(type=TRACTION, table=((1e5, 1e5, 1e5), ))
    mat1.QuadsDamageInitiation(table=((25, 50.0, 50.0), ))
    #m1.quadsDamageInitiation.DamageEvolution(type=DISPLACEMENT, table=((0.01, ), ))
    mat1.quadsDamageInitiation.DamageEvolution(type=ENERGY, mixedModeBehavior=BK, power=2.0,table=((0.34, 1.0, 1.0), ))

    #j=modelDCBResin(model1,ls,a0,crack_a,ws,plies,resinName,hresin=hresin,delta=delta,rPlastic=rPlastic)
    j=modelDCB_ZSYMM(model1,ls,a0,crack_a,ws,plies,resinName,hresin=hresin,delta=delta,rPlastic=rPlastic)
    
    if True:
        """
        model1.CohesiveSection(name='Section-Cohesive', 
            material='Cohesive-T300-7901', response=TRACTION_SEPARATION, 
            outOfPlaneThickness=None)
        p = model1.parts['Part-1']
        region=p.sets['Ply-resin'] if 'Ply-resin' in p.sets else p.sets['Ply-1']
        p.SectionAssignment(region=region, sectionName='Section-Cohesive', offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', 
            thicknessAssignment=FROM_SECTION)
        compositeLayup = p.compositeLayups['Laminate']
        compositeLayup.plies[1].suppress()
        
        import mesh
        elemType1 = mesh.ElemType(elemCode=COH3D8, elemLibrary=STANDARD, 
            elemDeletion=ON, viscosity=1E-6)
        p.setElementType(regions=region, elemTypes=(elemType1,))
        model1.rootAssembly.regenerate()
        """

        model1.steps['Step-1'].setValues(initialInc=0.1, maxInc=0.1)
        model1.StaticStep(name='Step-2', previous='Step-1', maxNumInc=3000, 
                            initialInc=0.001, maxInc=0.002)

        model1.boundaryConditions['DCB-Tension-1'].setValuesInStep(stepName='Step-1', u2=2.0)    
        model1.boundaryConditions['DCB-Tension-1'].setValuesInStep(stepName='Step-2', u2=10.0)    
        try:
            model1.boundaryConditions['DCB-Tension-2'].setValuesInStep(stepName='Step-1', u2=-2.0)    
            model1.boundaryConditions['DCB-Tension-2'].setValuesInStep(stepName='Step-2', u2=-10.0)    
        except:
            pass

        model1.fieldOutputRequests['F-Output-1'].setValues(
            variables=('S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'UR', 'RF', 'CF', 
            'NFORC', 'CSTRESS', 'CDISP', 'SDEG', 'SDV','STATUS'),frequency=10)
        umat=model1.materials['UMAT-Matrix']
        umat.depvar.setValues(deleteVar=16, n=16)
        Ks={4:1.03,2:1.83,1:2.98}
        umat.userMaterial.setValues(
            mechanicalConstants=(3170.0, 0.355, 
                85.1, 107.0, 52.6, 
                1.0, 42.7, 3170.0, 
                0.01, Ks[Npermm], 1.0))    
        j.setValues(userSubroutine='E:/1730895/Work/Abaqus/UMAT/UMAT_lipeng.for', numCpus=4, numDomains=4)
