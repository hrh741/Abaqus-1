import sys
sys.path.append('E:/1730895/Work/')

from Abaqus.Composite import  *

np.set_printoptions(linewidth=np.inf)
formula=0
if sys.argv[1:]:
    formula=int(sys.argv[1])

# O'Brien T K. Characterization of delamination onset and growth in a composite laminate,1982.
OBrienMat=[138e3,15e3,15e3,0.21,0.21,0.21,5.9e3,5.9e3,5.9e3]
lamS=Compliance(OBrienMat,abaqus=True)

plyAngles=np.pi/180*np.array([30,-30,30,-30,90,90,90,-30,30,-30,30])
plyThicks=1.51/11*np.ones_like(plyAngles)
E0=delaminatedEx(lamS,plyAngles,plyThicks,[])
## 137=(0.258*58211.1)*1.51/2.0*0.00347**2
## 直接利用Orien的公式是得不到正确结果137的
print("%.3f J/m^2"%(OBrienERR(lamS,plyAngles,1.51/11*np.ones_like(plyAngles),[3,6],0.00347,formula=0)*10**3)) # 单位 Mpa * mm = 10^3 J/m^2)

## Obrien在±30_2的层板计算结果进行了修正
plyAngles=np.array([ x*np.pi/180.0 for x in (30,-30,30,-30)])
plyThicks=1.51/11*np.ones_like(plyAngles)
E30_1=delaminatedEx(lamS,plyAngles,plyThicks,[],formula=0)
plyAngles=np.array([ x*np.pi/180.0 for x in (30,-30,-30,30)])
plyThicks=1.51/11*np.ones_like(plyAngles)
E30_2=delaminatedEx(lamS,plyAngles,plyThicks,[],formula=0)
E90=delaminatedEx(lamS,[x*np.pi/180.0 for x in (90,90,90)],1.51/11*np.ones((3,)),[],formula=0)
## O'brien 的论文中 (30/-30)_2 为0.69 不考虑(30/-30) 的拉弯耦合为0.743
## 0.258约等于1-0.743
print(E0,(E30_1*8+E90*3)/(E0*11),(E30_2*8+E90*3)/(E0*11))

## 实际上基于kxy=0的假设可以直接得到很好的结果
plyAngles=np.array([ x*np.pi/180.0 for x in (30,-30,30,-30)])
plyThicks=1.51/11*np.ones_like(plyAngles)
E30_1=delaminatedEx(lamS,plyAngles,plyThicks,[],formula=2)
plyAngles=np.array([ x*np.pi/180.0 for x in (30,-30,-30,30)])
plyThicks=1.51/11*np.ones_like(plyAngles)
E30_2=delaminatedEx(lamS,plyAngles,plyThicks,[],formula=2)
E90=delaminatedEx(lamS,[x*np.pi/180.0 for x in (90,90,90)],1.51/11*np.ones((3,)),[],formula=0)
print(E0,(E30_1*8+E90*3)/(E0*11),(E30_2*8+E90*3)/(E0*11))
