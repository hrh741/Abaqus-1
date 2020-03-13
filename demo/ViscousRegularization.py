# -*- coding: UTF-8 -*- 

"""
展示Viscous Regularization 的原理

"""
import sys
sys.path.append("E:/1730895/Work/")

import math
from Abaqus.RungaKutta import rk4
from functools import partial
from matplotlib import pyplot as plt

def showViscousRegularization(d,eta,n=1000):
    """
    $d_v'=(d-d_v)/\eta$
    """
    f=partial(lambda t,y,d,eta: (d(t)-y)/eta,d=d,eta=eta)
    vx,vy=rk4(f,0,0,1,n)
    analy=list(map(d,vx))
    print(analy[-1],vy[-1],)
    plt.plot(vx,vy)
    plt.plot(vx,analy)
    plt.legend(['Numerical','Analytical'])
    plt.show()


def UMAT(u1,du,sgm0,k,u_c,u_f,eta=0.0):
    """
    u1:   u at the end of the step
    du:   u increment of the step
    sgm0: sgm at the start of step
    """
    if u1<=u_c:
        sgm=k*u1
        ddsdde=k
    elif u1<=u_f:            
        sgm=k*u_c/(u_f-u_c)*(u_f-u1)
        
        ddsdde=-u_c*k/(u_f-u_c)  #real
        #ddsdde=k*u_c*u_f/(u_f-u_c)*(u_f/u1-1)
        #print(u,sgm)
    elif u1>u_f:
        sgm=0.0
        ddsdde=0.0

    return ddsdde,sgm

def UMATofDamage(u1,du,sgm0,k,u_c,u_f,statev,eta=0.0):
    """

    statev: [dmgvold,]
    """
    dmgvold=statev[0]
    if u1<=u_c:
        dmg=0.0
        ddsdde=k
    elif u1<=u_f:
        dmg=1.0-u_c/(u_f-u_c)*(u_f-u1)/u1
        
        dmg=dmg*du/(eta+du)+dmgvold*eta/(eta+du)
        #ddsdde=-u_c*k/(u_f-u_c)  #real tan
        ddsdde=1.0  #viscous 
        #ddsdde=k*u_c*u_f/(u_f-u_c)*(u_f/u1-1) # secant
        #print(u,sgm)
    elif u1>u_f:
        dmg=1.0
        ddsdde=0.0
    
    sgm=(1.0-dmg)*k*u1

    statev[0]=dmg
    return ddsdde,sgm

"""
eta=0.01
d=lambda t: math.sqrt(t) # t-eta+eta*exp(-t/eta)
n=1000
showViscousRegularization(d=d,eta=eta,n=n)
"""

# Bilinear 
k,dlt_f,dlt_c,k0=2,2,1,1.0

## viscous regularization
eta=0.001

## 加载长度
u_c=10
du=0.001

u,dlt,sgm=0.0,0.0,0.0
ddsdde=1.0
sdeg=0.0

sol=[(u,dlt,sgm)]
statev=[0,]
while u<u_c:
    # increment
    dp=k0*du
    iterCnt=0

    while abs(dp)>1e-6:
        if iterCnt>100:
            print('convergence error in ',u,dlt,ddsdde,dp,statev)
            sys.exit(1)    

        # iteration 
        presgm,predlt,preddsdde=sgm,dlt,ddsdde

        ddlt=dp/(preddsdde+k0)
        """
        if abs(ddlt)<1e-8:
            u=u+du-dp/k0
            break
        """
        dlt=dlt+ddlt
        #ddsdde,sgm=UMAT(dlt,ddlt,presgm,k,dlt_c,dlt_f)
        ddsdde,sgm=UMATofDamage(dlt,ddlt,presgm,k,dlt_c,dlt_f,statev)
        dp=dp-(sgm-presgm)-k0*(dlt-predlt)
        iterCnt+=1
    u=u+du
    sol.append((u,dlt,sgm))

import matplotlib.pyplot as plt
import numpy as np

sol=np.array(sol)
plt.subplot('121')
plt.plot(sol[:,0],sol[:,2])
plt.xlabel('$u$')
plt.ylabel('$\sigma$')

plt.subplot('122')
plt.plot(sol[:,1],sol[:,2])
plt.xlabel('$\delta$')
plt.ylabel('$\sigma$')

plt.show()
