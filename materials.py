from abaqus import *
from abaqusConstants import *

materials=[
    ## 吴庆欣
    ('T300-7901',(137.78e3,8.91e3,8.91e3,0.3,0.3,0.48,4.41e3,4.41e3,3.01e3),None,None), #
    ('Pipes-Pagano',(20.0, 2.1, 2.1, 0.21, 0.21, 0.21, 0.85, 0.85, 0.85),None,None), #MPsi
    ## Uguen A, Zubillaga L, Turon A, et al. Comparison of cohesive zone models used to predict delamination initiated from free-edges: validation against experimental results[C]//ECCM-16TH European Conference on Composite Materials. 2014.
    ### GIc,GIIc,GIIIc,sigma_zz0,tau_xz0,tau_yz0=0.24,0.74,0.74,46,75,75
    ('T800-M21',(130.0e3,8.0e3,8.0e3,0.31,0.31,0.45,4.0e3,4.0e3,4.0e3),None,None), #MPa
    ('S2/SP250 GlasdEpoxy',None,None,None),
    ('G947/M18',(97.6e3, 8.0e3 ,8.0e3,0.37,0.37,0.5, 3.1e3, 3.1e3, 2.7))
]

m1=mdb.models['DCB']

for mat in materials:
    if mat[1]:
        material1=m1.Material(name=mat[0])
        material1.Elastic(type=ENGINEERING_CONSTANTS, table=(mat[1], ))
