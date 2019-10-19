inpfile="E:/Abaqus/Taper-2.inp"
## 从Seam的Set中查找ns和elms

ns=[] # Seam的节点集合
elms= [] # Seam的单元集合

Nodes={}
Elements={}
fp=open(inpfile,"r")
inNode=0 # 三种状态 0: 未开始 1: 开始 2: 结束
inElement=0 #  三种状态 0: 未开始 1: 开始 2: 结束
inSeamNset=0
inSeamElset=0
for line in fp:
    s=line.strip().upper()
    if inNode==1:
        if s.startswith("*"):
            inNode=2
        else:
            nodeLabel,x,y,z=s.split(",")
            nodeLabel,x,y,z=int(nodeLabel),float(x),float(y),float(z)
            Nodes[nodeLabel]=(x,y,z)
    if inSeamNset==1:
        if s.startswith("*"):
            inSeamNset=2
        else:
            ns.extend([int(n) for n in s.split(",")])
    if inSeamElset==1:
        if s.startswith("*"):
            inSeamElset=2
        else:
            elms.extend([int(n) for n in s.split(",")])
    if inElement==1:
        if s.startswith("*"):
            inElement=2
        else:
            elmdata=s.split(",")
            Elements[int(elmdata[0])]=[int(n) for n in elmdata[1:]]
    if s.startswith("*NODE") and inNode==0:
        inNode=1
    elif s.startswith("*ELEMENT,") and inElement==0:
        inElement=1
    elif s.startswith("*NSET, NSET=SEAM,") and inSeamNset==0:
        inSeamNset=1
    elif s.startswith("*ELSET, ELSET=SEAM,") and inSeamElset==0:
        inSeamElset=1        
fp.close()
#for n in ns:
#    print("\t",n,":",Nodes[n])
#for elm in elms:
#    print(elm)
#    for i in Elements[elm]:
#        print("\t",i,":",Nodes[i])

## 查找ns中坐标相同的点
sameCoord={}
sortedNs=sorted(sorted(ns,key=lambda x:Nodes[x][2]),key=lambda x:Nodes[x][0])
for n in sortedNs:
    print(n,':',Nodes[n],sep='\t')
## 创建裂纹前沿为z方向的裂纹

## seam的端点，是seam的节点集中不会出现两个相同坐标的点
endPoints=set()#set({19 ,462,20 ,67,1164,65 })
sameCoord={}

i,n=0,len(sortedNs)
while i<n-1:
    x1,y1,z1=Nodes[sortedNs[i]]
    x2,y2,z2=Nodes[sortedNs[i+1]]
    if x1==x2 and y1==y2 and z1==z2:
        sameCoord[sortedNs[i]]=sortedNs[i+1]
        sameCoord[sortedNs[i+1]]=sortedNs[i]
        i=i+2
    else :
        endPoints.add(sortedNs[i])
        i=i+1
endPoints.add(sortedNs[n-1])
print("坐标相同的节点映射",sameCoord)
print("seam端点:",endPoints)

#裂纹前沿上的点
cracktips=[] #[(428, 4850, ),(4850,  426 , ),(818 , 5356, ),(815 , 5356, )]
i,n=0,len(sortedNs)
while i<n-1:
    x1,y1,z1=Nodes[sortedNs[i]]
    x2,y2,z2=Nodes[sortedNs[i+1]]
    if x1==x2 and y1==y2 and z1==z2:
        if Nodes[sortedNs[i+2]][2]>z1:
            cracktips.append((sortedNs[i],sortedNs[i+2]))
        else:
            break
        i=i+2
    else :
        i=i+1
i=n-1
while i>=1:
    x1,y1,z1=Nodes[sortedNs[i]]
    x2,y2,z2=Nodes[sortedNs[i-1]]
    if x1==x2 and y1==y2 and z1==z2:
        if Nodes[sortedNs[i-2]][2]<z1:
            cracktips.append((sortedNs[i],sortedNs[i-2]))
        else:
            break
        i=i-2
    else :
        i=i-1

print("裂纹尖端前沿节点对数组:",cracktips)

print("""** UEL of Interface element for VCCT
*USER ELEMENT, TYPE=U1, NODES=8, COORD=3, PROPERTIES=1, VARIABLES=3
1,2,3
*ELEMENT, TYPE=U1,ELSET=SEAM_UEL""")
for ind,cracktip in enumerate(cracktips):
    # cracktip为裂纹前沿上的点
    uelNodes=[]
    #print(cracktip)
    cracktip=set({sameCoord[n] for n in cracktip}) | set(cracktip)
    for elm in elms:
        nset=set(Elements[elm])
        if len(nset & set(cracktip))>1 and len(nset & endPoints)==0: # 选择seam的中间单元
            #print("\t",elm,nset,nset & set(ns))
            uelNodes.extend(sorted(sorted(list(nset & set(ns)),key=lambda x:Nodes[x][0]),key=lambda x:Nodes[x][2])) # 分别用x,z坐标就地排序
    #print(uelNodes)
    if uelNodes[0] in cracktip:
        for i in range(4):
            uelNodes[2*i],uelNodes[2*i+1]=uelNodes[2*i+1],uelNodes[2*i]
    #for n in uelNodes:
    #    print('\t',n,':',Nodes[n],sep='\t')
    UELformat="%d,%s"%(1000000+ind,",".join([str(i) for i in uelNodes])) # 起始单元号为1000000 可以修改成其它不再elements中出现的单元号
    for n in uelNodes:
        print(Nodes[n],end=" ")
    print("")
    print(UELformat)
print("*UEL PROPERTY, ELSET=SEAM_UEL\n1.E6,")
