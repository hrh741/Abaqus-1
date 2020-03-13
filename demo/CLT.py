#! -*- Coding: UTF-8 -*-
import sys
sys.path.append('E:/1730895/Work/')

from Abaqus import Laminate,Compliance
from Abaqus.Materials import materials
import numpy as np

T300_7901=Compliance(materials['T300-7901'][0],abaqus=True)
plyAngles=[0,30,-30,60,-60,90,90,-60,60,-30,30,0]
lam=Laminate(T300_7901,
            [x/180.0*np.pi for x in plyAngles],
            t=0.125)

print(lam.DelaminateStiffness(),lam.DelaminateStiffness(cracks=[5,]))

CLTstf,Cmps=lam.CLTStiffness(isPlaneStress=False)
for cmp in Cmps:
    print(np.linalg.inv(cmp)[:,0])

np.set_printoptions(linewidth=np.inf)
formula=0
if sys.argv[1:]:
    formula=int(sys.argv[1])
Pipes=[20.0, 2.1, 2.1, 0.21, 0.21, 0.21, 0.85, 0.85, 0.85]

T300_7901=[ 136999.70000000,8570.38221423,8570.38221423,
            0.26355000,0.26355000,0.33647467,
            4494.49736748,4494.49736748,3206.33918429]

T300_7901=[137.78e3,8.91e3,8.91e3,
            0.3,0.3,0.48,
            4.41e3,4.41e3,3.01e3]
T300_7901_Vf593=[137680.19000000 ,8605.27093889, 8605.27093889 , 
                0.26308500  , 0.26308500  ,0.33505172 , 
                4525.95447127 , 4525.95447127,3222.82306907]

T300_7901=[138e3,8.9e3,8.9e3,
            0.27,0.27,0.48,
            4.8e3,4.8e3,3.01e3]

T300_7901=[137.78e3,8.04e3,8.04e3,
            0.3,0.3,0.48,
            4.41e3,4.41e3,3.01e3]

mat=T300_7901
lamS=Compliance(mat,abaqus=True)
plyAngles=np.pi/180*np.array([38,-38,-38,38])
plyAngles=np.pi/180*np.array([10,-10,-10,10])
plyAngles=np.pi/180*np.array([45,-45,-45,45])
plyThicks=0.135*np.ones_like(plyAngles)
lam0=Laminate(lamS,plyAngles,plyThicks)
print(lam0.UniformAxiaExtentionSolution(ex=1.0,verbose=True),lam0.DelaminateStiffness(),lam0.DelaminateStiffness(cracks=[0,2]))

"""
plyAngles=np.pi/180*np.array([0,60,-60,-60,60,0])
plyAngles=np.pi/180*np.array([0,60,-30,30,-60,90,90,-60,30,-30,60,0])
plyAngles=np.pi/180*np.array([0,30,-30,60,-60,90,90,-60,60,-30,30,0])
plyThicks=0.135*np.ones_like(plyAngles)
lam0=Laminate(lamS,plyAngles,plyThicks)
E0=lam0.DelaminateStiffness()
E1=lam0.DelaminateStiffness(cracks=[5,])
E2=lam0.DelaminateStiffness(cracks=[2,8])
E3=lam0.DelaminateStiffness(cracks=[1,9])
print(E0,E1,E2,E3,E1/E0,E2/E0,E3/E0)


lamSs=[Compliance(T300_7901,abaqus=True),Compliance([3.17e3,0.355])]
plyAngles=np.pi/180*np.array([0,0,30,0,-30,0,
                             60,0,-60,0,90,0,
                             0,90,0,-60,0,60,
                             0,-30,0,30,0,0])
plyThicks=[ 0.125,0.001,0.125,0.001,0.125,0.001,
            0.125,0.001,0.125,0.001,0.125,0.001,
            0.001,0.125,0.001,0.125,0.001,0.125,
            0.001,0.125,0.001,0.125,0.001,0.125,]
lam0=Laminate([lamSs[0] if t==0.125 else lamSs[1] for t in plyThicks],plyAngles,plyThicks=plyThicks)
print(lam0.UniformAxiaExtentionSolution(ex=0.02))
"""