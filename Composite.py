# -*- coding: utf-8 -*-
import numpy as np
from itertools import chain

def Ortho_Compliance(El,Et,Ez,nutz,nulz,nult,Gtz,Glz,Glt):
  """
  各项异性柔度矩阵,传统顺序 11 22 33 23 13 12
  """
  lamS=np.array([1/El,-nult/El,-nulz/El,0,0,0,
                -nult/El,1/Et,-nutz/Et,0,0,0,
                -nulz/El,-nutz/Et,1/Ez,0,0,0,
                0,0,0,1/Gtz,0,0,
                0,0,0,0,1/Glz,0,
                0,0,0,0,0,1/Glt]).reshape((6,6))
  return lamS

def Iso_Compliance(E,nu):
  G=E/(2.0+2.0*nu)
  return Ortho_Compliance(E,E,E,nu,nu,nu,G,G,G)

def Compliance(props,abaqus=False):
  """
  如果abaqus 为True,代表props中的顺序为abaqus的顺序, 11 22 33 12 13 23 12 13 23
  """
  if len(props)==2:
    return Iso_Compliance(props[0],props[1])
  elif len(props)==9:
    if abaqus:
      return Ortho_Compliance(props[0],props[1],props[2],
                              props[5],props[4],props[3],
                              props[8],props[7],props[6])
    else:
      return Ortho_Compliance(*props)

def transformers(theta=None,axis=None,rotate=None,abaqus=False):
    """
    根据绕着axis轴转角为theta 或者 rotate 给定 柔度矩阵的 变换
    
    输入:
      theta : 旋转角 弧度制
      axis :  [None,'x','y','z']
      rotate : 的每一行都是 新坐标系下轴在旧坐标系下的向量
    返回:
      Tepl  : 应变变换矩阵
      Tsgm  : 
    transformers(theta=np.pi/2,axis='x')
    transformers(rotate=np.array([[1,0,0],[0,1,0],[0,0,1]]))
    
    abaqus Odb中element的应力先从局部坐标系转换为全局坐标系
    transformers(rotate=np.array(sfo.values[0].localCoordSystem).T)
    """
    axis_d={'x':0,'y':1,'z':2}
    if not (theta is None or axis is None):
        c,s=np.cos(theta),np.sin(theta)
        rotate=np.zeros((3,3))
        if axis not in axis_d:
            raise ValueError("axis must be in ['x','y','z']")
        k=axis_d[axis]
        rotate[k,k]=1
        i,j=(k+1)%3,(k+2)%3
        rotate[i,i]=c
        rotate[i,j]=s
        rotate[j,j]=c
        rotate[j,i]=-s
    elif rotate is None :
        raise ValueError('pleast give at least theta and axis or rotate')
    
    Tsgm,Tepl=np.zeros((6,6)),np.zeros((6,6))
    if abaqus:
      vogit_seq=((1,1),(2,2),(3,3),(1,2),(1,3),(2,3))
    else:
      vogit_seq=((1,1),(2,2),(3,3),(2,3),(1,3),(1,2))
    for n1 in range(6):
        i,j=vogit_seq[n1][0]-1,vogit_seq[n1][1]-1
        for n2 in range(6):
            k,l=vogit_seq[n2][0]-1,vogit_seq[n2][1]-1
            Tsgm[n1,n2]=(rotate[i,k]*rotate[j,l]+rotate[i,l]*rotate[j,k])/2 #if k==l else rotate[i,k]*rotate[j,l]+rotate[i,l]*rotate[j,k]

    route=np.diag([1,1,1,2,2,2,])
    route_inv=np.diag([1,1,1,0.5,0.5,0.5,])
    Tepl=route.dot(Tsgm)
    Tsgm=Tsgm.dot(route)
    return Tsgm,Tepl

def transform_Compliance(S,theta,axis):
    """
    计算在绕axis转动theta后的局部坐标系中柔度矩阵

    Parameters
    ----------
    S : 
        柔度矩阵
    theta: 
        弧度值
    axis : 
      0,'x',1,'y',2,'z' 


    Returns
    ----------
      S1
      
    Examples
    ----------


    """
    _,Tepl=transformers(theta=theta,axis=axis)
    return Tepl.dot(S).dot(Tepl.T)

def Suo(Exx,Eyy,Gxy,nuxy,a,h,B,planeStress=False,nuyz=0.0):
  """
  来自: Xu W , Guo Z Z . A simple method for determining the mode I interlaminar fracture toughness of composite without measuring the growing crack length[J]. Engineering Fracture Mechanics, 2018, 191:476-485.
  
  Exx,Eyy,Gxy,nuxy: 单层板材料属性
  a : 裂纹长度
  h : 层合板厚度的一半
  B : 层合板宽度
  返回:
  C : DCB的刚度
  G : 单位张开位移时的能量释放率，注意 G与张开位移平方成正比
  
  注意DCB处于平面应变状态，这里的默认计算结果为平面应变状态的
  """
  if not planeStress:
    if nuyz==0:
      raise Exception('不正确的nuyz %f'%nuyz)
    Exx,Eyy,nuxy=Exx/(1-nuxy**2*(Eyy/Exx)),Eyy/(1-nuyz**2),(nuxy+nuxy*nuyz)/Exx
    nuxy=Exx*nuxy
  nuyx=nuxy*Eyy/Exx
  rho=np.sqrt(Exx*Eyy)/(2.0*Gxy)-np.sqrt(nuxy*nuyx)
  beta=(0.677+0.146*(rho-1)-0.0178*(rho-1)**2+0.00242*(rho-1)**3)/np.power(Eyy/Exx,0.25)
  a_h=a/h
  C=(B*Exx)/(24.0*(a_h**3/3+beta*a_h**2+beta**2*a_h))
  G=(12.0*C**2/(B**2*Exx*h))*(a_h**2+2*beta*a_h+beta**2)
  #print Exx,Eyy,Gxy,nuxy,rho,beta
  return C,G

def linear_fit(x0,y0):
  """
  ax+b拟合x0,y0
  """
  res0=np.linalg.lstsq(np.vstack([x0,np.ones(len(x0))]).T,y0)
  a,b=res0[0][0],res0[0][1]
  return a,b

def TsaiWu(stressOfAbaqus,strengths):
  t,c,s=strengths
  F1,F11,F66=1./t-1./c,1./(t*c),1./s**2
  sg1,sg2,sg3,sg12,sg13,sg23=stressOfAbaqus
  tw=F1*(sg1+sg2+sg3)+F11*(sg1**2+sg2**2+sg3**2-sg1*sg2-sg1*sg3-sg2*sg3)+F66*(sg12**2+sg13**2+sg23**2)
  return tw