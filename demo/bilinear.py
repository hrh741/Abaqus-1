# -*-  Coding: UTF-8 -*-
#  刚度曲线| bilinear |         k0 |
# |<-^^^^^^^^^^^->|<-^^^^^^^^^^^->

from collections import deque
import time

k,dlt_f,dlt_c,k0=2,2,1,1

## viscous regularization
eta=0.01

def UMAT(u1,du,sgm0,k,u_c,u_f):
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
## 加载长度
u_c=10

du=0.001

"""
t0=time.time()
u,dlt,sgm=0.0,0.0,0.0
sol=[]
sol.append((u,dlt,sgm))
ddsdde=1.0
while u<u_c:
    # increment
    dp=k0*du
    iterCnt=0

    while abs(dp)>1e-8:
        # iteration 
        presgm,predlt,preddsdde=sgm,dlt,ddsdde

        ddlt=dp/(preddsdde+k0)
        if abs(ddlt)<1e-8:
            u=u+du-dp/k0

            break
        
        dlt=dlt+ddlt
        if dlt<=dlt_c:
            sgm=k*dlt
            ddsdde=k
        elif dlt<=dlt_f:            
            sgm=k*dlt_c/(dlt_f-dlt_c)*(dlt_f-dlt)
            
            #ddsdde=-dlt_c*k/(dlt_f-dlt_c)  #real
            #ddsdde=k*dlt_c*dlt_f/(dlt_f-dlt_c)*(dlt_f/dlt-1)
            ddsdde= # viscous regularization
            #print(dlt,sgm)
        elif dlt>dlt_f:
            sgm=0.0
            ddsdde=1.0
        dp=dp-(sgm-presgm)-k0*(dlt-predlt)
        iterCnt+=1
    
    u=u+du
    sol.append((u,dlt,sgm))
"""
t0=time.time()
u,dlt,sgm=0.0,0.0,0.0
sol=[]
sol.append((u,dlt,sgm))
ddsdde=1.0
sdeg=0.0
while u<u_c:
    # increment
    dp=k0*du
    iterCnt=0

    while abs(dp)>1e-8:
        # iteration 
        presgm,predlt,preddsdde=sgm,dlt,ddsdde

        ddlt=dp/(preddsdde+k0)
        if abs(ddlt)<1e-8:
            u=u+du-dp/k0
            break
        dlt=dlt+ddlt
        ddsdde,sgm=UMAT(dlt,ddlt,presgm,k,dlt_c,dlt_f)
        dp=dp-(sgm-presgm)-k0*(dlt-predlt)
        iterCnt+=1
    
    u=u+du
    sol.append((u,dlt,sgm))

print(time.time()-t0)
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
