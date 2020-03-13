#! /usr/bin/python3
# -*- coding: UTF-8 -*-
try:
    from .Composite import transform_Compliance,Compliance
except:
    import sys
    sys.path.append('E:/1730895/Work/')
    from Abaqus import transform_Compliance,Compliance

from itertools import chain
from collections import namedtuple
import numpy as np


Lamina=namedtuple('Lamina',['Compliance','angle','thickness','z0','z1',])
class Laminate:
    """
    复合材料层合板

        MatCompiances: compliance array or list of compliance array
            层合板材料体系的柔度矩阵
        plyAngles:
            弧度值，层合板铺层角度列表
        plyThicks:
            层合板厚度列表
        zmid:
            中面坐标，默认为0
        t: 
            单层厚度，如果plyThicks为0时使用

    """
    def __init__(self,MatCompliances,plyAngles,plyThicks=None,zmid=0.0,t=0.125):
        """
        """
        if(plyThicks is None):
            plyThicks=t*np.ones_like(plyAngles)
        else:
            assert len(plyAngles)==len(plyThicks),"please make len(plyAngles)==len(plyThicks)!" 
                
        # the list of compliance matrix of materials in each ply
        if isinstance(MatCompliances,np.ndarray):
            MatCompliances=[MatCompliances for _ in range(len(plyAngles))]

        s=np.cumsum(plyThicks)
        plyLowZs=(zmid+s[-1]/2.0)-s
        
        # the list of nametuple `Lamina` contain information of ply 
        # material,angle(in rad),thick,z0,z1
        self.plies=[Lamina(MatCompliances[i],
                            plyAngles[i],plyThicks[i],plyLowZs[i],plyLowZs[i]+plyThicks[i])
            for i in range(len(plyThicks))]

        self.H=sum((ply[2] for ply in self.plies))
        self.Z=zmid
        self.N=len(self.plies)
    
    def __getitem__(self,item):
        """
        注意，不变动子层坐标
        """
        if isinstance(item,slice):
            start,stop,step=item.start,item.stop,item.step
            
            return Laminate(MatCompliances=[self.plies[i][0] for i in range(start,stop)],
                            plyAngles=[self.plies[i][1] for i in range(start,stop)],
                            plyThicks=[self.plies[i][2] for i in range(start,stop)],
                            zmid=(self.plies[start][4]+self.plies[stop-1][3])/2.0)

    def __repr__(self):
        return '['+','.join(('(%3d,%.4f,%.4f)'%(int(ply[1]/np.pi*180),ply[3],ply[4]) for ply in self.plies))+']'

    def CLTStiffness(self,toMid=False,isPlaneStress=True):
        """
        计算经典层合板理论的刚度矩阵
                            [   A | B
                                B | D]
        注意经典层板理论处于平面应力状态
        注意 默认刚度矩阵的积分是相对于Z=0的

        toMid: True 如果向中面计算
        """
        Z=self.Z if toMid else 0.0
        A,B,D=np.zeros((3,3),dtype=np.float64),np.zeros((3,3),dtype=np.float64),np.zeros((3,3),dtype=np.float64)
        globalCompliance=[None for _ in range(self.N)]
        for i,ply in enumerate(self.plies):
            lamS,angle,t,z0,z1=ply

            compliance=transform_Compliance(lamS,-angle,axis='z')
            globalCompliance[i]=compliance
            ixgrid=np.ix_([0,1,5],[0,1,5])
            if isPlaneStress:
                planeStressCompliance=compliance[ixgrid]
                Q=np.linalg.inv(planeStressCompliance)
            else:
                Q=np.linalg.inv(compliance)[ixgrid]
            A=A+Q*(z1-z0)
            B=B+Q*((z1-Z)**2-(z0-Z)**2)/2
            D=D+Q*((z1-Z)**3-(z0-Z)**3)/3
        
        CLTstf=np.vstack((np.hstack((A,B)),
                        np.hstack((B,D))))
        return CLTstf,globalCompliance
    
    def CLTSolution(self,CLTLoad):
        """
        return the Stress and Strain in the middle of each ply
        """
        CLTstf,globalCompliance=self.CLTStiffness()
        
        CLTek=np.linalg.solve(CLTstf,CLTLoad)
        emid,kappa=CLTek[:3],CLTek[3:]

        CLTSgm,CLTEpl=np.zeros((self.N,6)),np.zeros((self.N,6))
        for i,ply in enumerate(self.plies):
            e=emid+kappa*((ply[3]+ply[4])/2.0)
            CLTSgm[i,[0,1,5]]=np.linalg.inv(globalCompliance[i][[0,1,5]][:,[0,1,5]]).dot(e)
            CLTEpl[i,:]=globalCompliance[i].dot(CLTSgm[i,:])

        return CLTSgm,CLTEpl
    
    def UniformAxiaExtentionSolution(self,ex,verbose=0):
        """
        利用经典层合板理论计算均匀轴向应变时（轴向应变为ex)时各个单层应力和应变

        ex: 
            轴向应变
        verbose: optional
            输出 中面应变、曲率、层间平均应力和 合力和合力矩
        """
        CLTstf,globalCompliance=self.CLTStiffness()
        
        """
        ## Ny=Nxy=Mx=My=Mxy=0
        CLTstf1=CLTstf
        A1,B1,D1=CLTstf1[:3,:3],CLTstf1[:3,3:],CLTstf1[3:,3:]
        invD=np.linalg.inv(D1)
        S=np.linalg.inv((A1-np.dot(B1,np.dot(invD,B1.T)))/H) 
        emid=S[:,0]*(ex/S[0,0])
        kappa=-np.dot(invD,np.dot(B1.T,emid))    

        ## Ny=Nxy=kx=My=Mxy=0
        CLTstf1=CLTstf-np.dot(CLTstf[:,[3,]]/CLTstf[3,3],CLTstf[[3,],:])
        CLTstf1[:,3]=CLTstf[3,:]/CLTstf[3,3]
        CLTstf1[3,:]=CLTstf[3,:]/CLTstf[3,3]
        CLTstf1[3,3]=-1./CLTstf[3,3]
        A1,B1,D1=CLTstf1[:3,:3],CLTstf1[:3,3:],CLTstf1[3:,3:]
        invD=np.linalg.inv(D1)
        S=np.linalg.inv((A1-np.dot(B1,np.dot(invD,B1.T)))/H) # Note : B1 is may not symmetry
        emid=S[:,0]*(ex/S[0,0])
        kappa=-np.dot(invD,np.dot(B1.T,emid))
        kappa[0]=0.0

        ## Ny=Nxy=kxy=My=Mxy=0
        ixgrid=np.ix_((0,1,2,4,5),(0,1,2,3,4))
        CLTstf1=CLTstf[ixgrid]
        S=np.linalg.inv(CLTstf1)
        emid=S[:3,0]*(ex/S[0,0])
        kappa=np.zeros((3,))
        kappa[:2]=S[3:5,0]*(ex/S[0,0])
        kappa[2]=0.0
        """

        ## Ny=Nxy=My=kx=kxy=0
        ixgrid=np.ix_((0,1,2,4),(0,1,2,4))
        CLTstf1=CLTstf[ixgrid]
        S=np.linalg.inv(CLTstf1)
        emid=S[:3,0]*(ex/S[0,0])
        kappa=np.zeros((3,))
        kappa[1]=S[3,0]*(ex/S[0,0])

        NM=CLTstf.dot(np.vstack((emid,kappa)).reshape(-1,1))
        if verbose>0:
            print('eps:', emid,
                '\nkappa:',kappa,
                '\n global Compliance\n',globalCompliance,
                '\nCLT stiffness:',CLTstf,
                '\nCLT Load:',NM)
        
        CLTSgm,CLTEpl=np.zeros((self.N,6)),np.zeros((self.N,6))
        for i,ply in enumerate(self.plies):
            e=emid+kappa*((ply[3]+ply[4])/2.0)
            CLTSgm[i,[0,1,5]]=np.linalg.inv(globalCompliance[i][[0,1,5]][:,[0,1,5]]).dot(e)
            CLTEpl[i,:]=globalCompliance[i].dot(CLTSgm[i,:])

        return CLTSgm,CLTEpl
    
    def DelaminateStiffness(self,cracks=[],formula=3):
        """
        参考CHARACTERIZATION OF DELAMINATION ONSET AND GROWTH IN A COMPOSITE LAMINATE
        根据O'Brien的混合法则计算完全分层后的层合板的轴向模量Ex
        """
        np.set_printoptions(linewidth=np.inf)
                
        st,Ex=0,0.0
        for i in chain(cracks,[self.N-1,]):
            ed=i+1
            
            lam0=self[st:ed]
            CLTstf=lam0.CLTStiffness()[0]
            A1,B1,D1=CLTstf[:3,:3],CLTstf[:3,3:],CLTstf[3:,3:]

            if(formula==0):
                ## O'Brien Model Method
                ## Note A,B,D in Obrien's Model Must be calculated in middle of sublaminate
                dH=lam0.Z-self.Z
                D1=D1+dH**2*A1-(2*dH)*B1
                B1=B1-dH*A1
                S=np.linalg.inv(A1-np.dot(B1,np.dot(np.linalg.inv(D1),B1)))
                Ei=1/S[0,0]
            elif(formula==1):
                # My own thought
                # Nyy=Nxy=0 κx=0 Myy=Mxy=0
                CLTstf1=CLTstf-np.dot(CLTstf[:,[3,]]/CLTstf[3,3],CLTstf[[3,],:])
                CLTstf1[:,3]=CLTstf[3,:]/CLTstf[3,3]
                CLTstf1[3,:]=CLTstf[3,:]/CLTstf[3,3]
                CLTstf1[3,3]=-1./CLTstf[3,3]
                A,B,D=CLTstf1[:3,:3],CLTstf1[:3,3:],CLTstf1[3:,3:],
                S=np.linalg.inv(A-np.dot(B,np.dot(np.linalg.inv(D),B.T)))
                Ei=1/S[0,0]
            elif(formula==2):
                ## Ye Lin: Ny=My=κx=γxy=0
                Ei=A1[0,0]-((A1[0,1]*(A1[0,1]*D1[1,1]-B1[0,1]*B1[1,1]))+(B1[0,1]*(A1[1,1]*B1[0,1]-A1[0,1]*B1[1,1])))/(A1[1,1]*D1[1,1]-B1[1,1]*B1[1,1])
            elif(formula==3):
                # Ny=Nxy=My=κx=κxy=0
                ixgrid=np.ix_((0,1,2,4),(0,1,2,4))
                CLTstf1=CLTstf[ixgrid]
                S=np.linalg.inv(CLTstf1)
                Ei=1/S[0,0]
            Ex=Ex+Ei
            st=ed
        return Ex/self.H

    def ObrienErr(self,cracks=[],ex=1.0,formula=3):
        """
        根据O'Brien的混合法则计算在cracks中的层间出现分层时的总能量释放率
        !!! 注意: 如果出现多个裂纹，计算结果是这几个裂纹一同扩展时的总能量释放率
        """
        E0=self.DelaminateStiffness(formula=formula)
        E1=self.DelaminateStiffness(cracks=cracks,formula=formula)
        return (E0-E1)*self.H/2.0*ex**2
    
    def KimErr(self,thcrack,ex,plyMats=None):
        """
        Generalized theoretical analysis method for free-edge delaminations in composite laminates
        只支持单个裂纹能量释放率计算
        """
        N=self.N
        wholeStf=np.zeros((3,3))
        regions=[]
        ixgrid1,ixgrid2,ixgrid3=np.ix_([0,3,5],[0,3,5]),np.ix_([0,3,5],[1,2,4]),np.ix_([1,2,4],[1,2,4])
        for regioni0,regioni1 in [(0,N),(0,thcrack+1),(thcrack+1,N)]:
            stf,_=self[regioni0:regioni1].CLTStiffness()
            # Epsilon_yy,Gamma_xy,Kyy
            y1=-np.linalg.solve(stf[ixgrid3],stf[0,[1,2,4]])
            regions.append([stf,y1])
            print(y1)
            #print(y1,stf.dot(np.array([1,y1[0],y1[1],0,y1[2],0])))
        
        y1=regions[0][1]
        ny,nx,myy=stf[[1,2,4],:].dot(np.array([1,y1[0],y1[1],0,y1[2],0]))
        myy=myy-ny*(self.H/2-t)
        print(myy,ny,nx,regions[1][1][2]-regions[2][1][2])
        GI=myy*(regions[1][1][2]-regions[2][1][2])/2
        GII=ny*(regions[2][1][0]-regions[1][1][0])/2
        GIII=nx*(regions[2][1][1]-regions[1][1][1])/2
        return GI,GII,GIII


import unittest
class TestLaminate(unittest.TestCase):
    def testUniformAxiaExtentionSolution_Pipes1970(self):
        """
        Pipes, R B. and N. J. Pagano. "Interlaminar Stresses in Composite Materials under Uniform Axial Extension," J Composite Materials, 4:538-548 (1970).
        Fig 3
        """
        Pipes=[20.0, 2.1, 2.1, 0.21, 0.21, 0.21, 0.85, 0.85, 0.85]
        lamS=Compliance(Pipes,abaqus=True)
        plyAngles=np.pi/180*np.array([45,-45,-45,45])
        plyThicks=0.125*np.ones_like(plyAngles)

        lam0=Laminate(lamS,plyAngles,plyThicks)
        sgm,epl=lam0.UniformAxiaExtentionSolution(ex=1.0)

        self.assertAlmostEqual(sgm[0][0],2.96,places=2)
        self.assertAlmostEqual(sgm[0][5],1.15,places=2)
        self.assertAlmostEqual(epl[0][0],1.0, places=2)

    def testgetitem(self):
        Pipes=[20.0, 2.1, 2.1, 0.21, 0.21, 0.21, 0.85, 0.85, 0.85]
        lamS=Compliance(Pipes,abaqus=True)
        plyAngles=np.pi/180*np.array([45,-45,-45,45])
        plyThicks=0.125*np.ones_like(plyAngles)

        lam0=Laminate(lamS,plyAngles,plyThicks)

        lam1=lam0[0:2]
        for i in range(2):
            self.assertEqual(lam1.plies[i].angle,plyAngles[i])
            self.assertEqual(lam1.plies[i].thickness,plyThicks[i])
        #print(lam1)

    def testOBrien(self):
        """
        O'Brien T K. Characterization of delamination onset and growth in a composite laminate,1982.
        """
        OBrienMat=[138e3,15e3,15e3,0.21,0.21,0.21,5.9e3,5.9e3,5.9e3]
        lamS=Compliance(OBrienMat,abaqus=True)
        plyAngles=np.pi/180*np.array([30,-30,30,-30,90,90,90,-30,30,-30,30])
        plyThicks=1.51/11*np.ones_like(plyAngles)
        lam0=Laminate(lamS,plyAngles,plyThicks)

        E0=lam0.DelaminateStiffness(cracks=[])
        E1=lam0.DelaminateStiffness(cracks=[3,6],formula=0)
        E1_1=lam0.DelaminateStiffness(cracks=[3,6],formula=3)
        
        self.assertAlmostEqual(E1/E0,0.69,delta=0.004)
        self.assertAlmostEqual(E1_1/E0,0.742,delta=0.003)
        self.assertAlmostEqual(lam0.ObrienErr(cracks=[3,6],ex=0.00347,formula=3)*1000,137,delta=2)
     
    def testSSwang(self):
        """
        Wang S S. Edge Delamination in Angle-Ply Composite Laminates[J]. AIAA journal, 1984, 22(2): 256-264.
        """
        Pipes=[20.0, 2.1, 2.1, 0.21, 0.21, 0.21, 0.85, 0.85, 0.85]
        lamS=Compliance(Pipes,abaqus=True)

        for theta,Exp in zip([15,30,45,60,75],[8.1076,4.0506,0.5740,0.0138,0.0036]):
            plyAngles=np.pi*theta/180*np.array([1,-1,-1,1])
            plyThicks=1.0*np.ones_like(plyAngles)
            lam0=Laminate(lamS,plyAngles,plyThicks)
            self.assertAlmostEqual(lam0.ObrienErr(cracks=[0,2])/Exp,2.0,delta=0.05)

    def testKweonKim(self):
        """
        (1999) Kim I K, Kang H G, Onohara K. Ultimate strength of composite laminates with free-edge delamination[J]. KSME International Journal, 1999, 13(7): 569-574.
        """
        #mats=[94.1e3,5690,5690,0.3,0.01995,0.33,2190,2110,2110]
        mats=[94.1e3,5690,5690,0.3,0.3,0.34834,2190,2110,2110]
        lamS=Compliance(mats)

        plyAngles=np.pi/180*np.array([30,30,-30,-30,90,90,-30,-30,30,30])
        plyThicks=np.ones_like(plyAngles)

        lam0=Laminate(lamS,plyAngles,plyThicks)

        #for crack,stiff in zip(([],[4,],[3,],[1,],[3,5],[1,7]),
        #                        (3880,3530,3480,3410,3370,2980)):
        #    print(lam0.DelaminateStiffness(cracks=crack)*1/3+lam0.DelaminateStiffness()*2/3,stiff)
        
        self.assertAlmostEqual(lam0.DelaminateStiffness(),3880)
        self.assertAlmostEqual(lam0.DelaminateStiffness([4,]),3530)
        self.assertAlmostEqual(lam0.DelaminateStiffness([3,]),3480)
        self.assertAlmostEqual(lam0.DelaminateStiffness([1,]),3410)
        self.assertAlmostEqual(lam0.DelaminateStiffness([3,5]),3370)
        self.assertAlmostEqual(lam0.DelaminateStiffness([1,7]),2980)
        
    def testQuasiIsotropic(self):
        """
        test isotropic laminate Stiffness
        """
        for mat in ([137.78e3,8.91e3,8.91e3,0.3,0.3,0.48,4.41e3,4.41e3,3.01e3],
                    [94.1e3,5690,5690,0.3,0.3,0.34834,2190,2110,2110]):
            E11,E22,nu12,G12=mat[0],mat[1],mat[3],mat[6]
            q=1./(1-nu12**2/E11*E22)
            Q11,Q12,Q22,Q66=q*E11,q*E22*nu12,q*E22,G12
            nu=(Q11+Q22-4*Q66+6*Q12)/(3*Q11+3*Q22+2*Q12+4*Q66)
            Ex=(Q11+Q22-2*Q12+4*Q66)*(1+nu)*2/8
            
            lamS=Compliance(mat,abaqus=True)
            for k in range(3,8):
                plyAngles=np.array([np.pi/k*i for i in range(k)])
                plyAngles=np.concatenate((plyAngles,np.flip(plyAngles,axis=0)))
                plyThicks=np.ones_like(plyAngles)
                lam0=Laminate(lamS,plyAngles,plyThicks=plyThicks)
                self.assertAlmostEqual(lam0.DelaminateStiffness(),Ex)

if __name__=='__main__':
    import doctest
    doctest.testmod()

    unittest.main(verbosity=4)
