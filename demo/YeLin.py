import sys
sys.path.append('E:/1730895/Work/')

from Abaqus.Composite import  *

np.set_printoptions(linewidth=np.inf)
formula=0
if sys.argv[1:]:
    formula=int(sys.argv[1])

#(1987)Role of Matrix Resin in Delamination Onset and Growth in Composite Laminates.pdf
#  El,Et,Ez,nutz,nulz,nult,Gtz,Glz,Glt
T300_648=[137e3,9.1e3,9.1e3,0.31,0.31,0.31,5.3e3,5.3e3,5.3e3]
T300_634=[133e3,7.7e3,7.7e3,0.33,0.33,0.33,4.2e3,4.2e3,4.2e3]
for plyAngles in [[30,-30,30,-30,30,-30,90,90,90,90,-30,30,-30,30,-30,30],
                [45,-45,45,-45,0,0,90,90,90,90,0,0,-45,45,-45,45]]:
    plyAngles=np.pi/180*np.array(plyAngles)
    plyThicks=0.135*np.ones_like(plyAngles)
    for mat in [T300_648,T300_634]:
        lamS=Ortho_Compliance(*mat)
        print(  delaminatedEx(lamS,plyAngles,plyThicks,[],formula=formula),
                delaminatedEx(lamS,plyAngles,plyThicks,[5,],formula=formula),
                delaminatedEx(lamS,plyAngles,plyThicks,[5,9],formula=formula))
