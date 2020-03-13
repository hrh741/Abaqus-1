#! /usr/python2
# -*- coding: UTF-8 -*- 

"""
渐进式删除单元实现分层破坏分析
给定inp文件，注意
1. inp文件要添加EL file 
2. 开启restart,Restart, write, frequency=1
3. 将step 调整好
"""
import sys,os,subprocess

def abaqusSetDataLine(s,maxEntries=16):
    """
    将编号转换为abaqus 的数据行,每行最多maxEntries个编号
    """
    n=len(s)
    dataLine=''
    elmlst=list(s)
    for start in range(0,n,maxEntries):
        dataLine+=','.join(('%d'%x for x in elmlst[start:min(start+maxEntries,n)]))+'\n'
    return dataLine

restartInpFormat="""*Heading
** 
*Preprint, echo=NO, model=NO, history=NO, contact=NO
*Restart,read,step={kstep},inc={kinc},end step
**
**DUMMY STEP delte previous failed element
**
*STEP,name=DUMMY-{stepName}
*Elset, elset=Assembly.{elementSetName},internal, instance=Part-1-1
{failedElements}
*STATIC
*MODEL CHANGE,REMOVE,TYPE=ELEMENT
Assembly.{elementSetName}
*END STEP
** 
** STEP: {stepName}
** 
*Step, name={stepName}, nlgeom=NO, inc=1000
*Static,direct
0.01, 1., 1e-05, 1.
** 
** BOUNDARY CONDITIONS
** 
** Name: X=1-Tension Type: Displacement/Rotation
*Boundary
x-ref, 1, 1, 0.01
** 
** OUTPUT REQUESTS
** 
*Restart, write, frequency=1
** 
** FIELD OUTPUT: F-Output-1
** 
*Output, field
*Node Output
CF, RF, U, UR
*Element Output, directions=YES
LE, NFORC, PE, PEEQ, PEMAG, S, SDV
*Contact Output
CDISP, CSTRESS
** 
** HISTORY OUTPUT: H-Output-Ref
** 
*Output, history
*Node Output, nset=x-ref
RF1, RF2, RF3, RM1, RM2, RM3, U1, U2
U3, UR1, UR2, UR3
** 
** HISTORY OUTPUT: H-Output-1
** 
*Output, history, variable=PRESELECT
*El file
SDV
*End Step
"""
DCBrestartInpFormat="""
*Heading
** 
*Preprint, echo=NO, model=NO, history=NO, contact=NO
*Restart,read,step={kstep},inc={kinc},end step
**
**DUMMY STEP delte previous failed element
**
*STEP,name=DUMMY-{stepName}
*Elset, elset=Assembly.{elementSetName},internal, instance=Part-1-1
{failedElements}
*STATIC
*MODEL CHANGE,REMOVE,TYPE=ELEMENT
Assembly.{elementSetName}
*END STEP
** 
** STEP: {stepName}
** 
*Step, name={stepName}, nlgeom=NO
*Static,direct
0.02, 1., 1e-05, 1.
** 
** BOUNDARY CONDITIONS
** 
** Name: DCB-Tension Type: Displacement/Rotation
*Boundary
RP1, 1, 1
RP1, 2, 2, 5
RP1, 3, 3
** 
** OUTPUT REQUESTS
** 
*Restart, write, frequency=1
** 
** FIELD OUTPUT: F-Output-1
** 
*Output, field
*Node Output
CF, RF, U, UR
*Element Output, directions=YES
LE, NFORC, PE, PEEQ, PEMAG, S, SDV
*Contact Output
CDISP, CSTRESS
** 
** HISTORY OUTPUT: H-Output-Ref
** 
*Output, history
*Node Output, nset=RP1
RF1, RF2, RF3, RM1, RM2, RM3, U1, U2
U3, UR1, UR2, UR3
** 
** HISTORY OUTPUT: H-Output-1
** 
*Output, history, variable=PRESELECT
*El file
SDV
*End Step
"""
def getRestartInp(**args):
    return restartInpFormat.format(**args)

inpfname='E:/Abaqus/RestartAnalysis/0_±30_±60_90s/0.25mm-K=1/Job-Delamiantion_I-quatermm_Modified.inp' #'E:/Abaqus/RestartAnalysis/0_±30_±60_90s/3D/1mm/Job-Delamiantion3D_I-1mm.inp'#
outdir=os.path.dirname(inpfname)
os.chdir(outdir)
jobname=os.path.basename(inpfname).split('.')[0]
restartfnameFormat=outdir+'/'+jobname+'-restart-%d.inp'

ind=int(sys.argv[1]) if sys.argv[1:] else 0
while True:
    ind+=1
    restartfname=restartfnameFormat%(ind)
    jobname=os.path.basename(restartfname).split('.')[0]
    if ind>30:
        break
    if ind==1:
        with open(restartfname,'w') as fp:
            fp.write(open(inpfname).read())
            # DUMMY-STEP-1 中删除单元1，如果不再第一个inp文件中添加modelChange的语句会报错
            fp.write('''\n**\n*STEP,name=DUMMY-STEP-1\n*ELset,elset=d1,internal,instance=Part-1-1\n1\n*STATIC\n*MODEL CHANGE,REMOVE,TYPE=ELEMENT\nd1\n*END STEP''')

        print('Starting ',jobname)
        p1=subprocess.Popen(['E:/Program/SIMULIA/Abaqus/Commands/abaqus.bat',
                            'job='+jobname,
                            'user=E:/1730895/Work/Abaqus/UMAT/UMAT_lipeng.for',
                            'interactive'],
                            stdin=subprocess.PIPE,stdout=open('./'+jobname+'.log','w'))
        out,err=p1.communicate('y'.encode())
        if p1.returncode==0:
            print(restartfname,' caculation completed!')
        else:
            print('caculation error,please check logfile!')
            exit(1)
    else:
        oldjob=os.path.basename(restartfnameFormat%(ind-1)).split('.')[0]
        fiberFailFile='E:/UMAT-OUTPUT/'+oldjob+'-fiber.txt'
        # verify whether fiber failure occur
        fiberFailure= os.stat(fiberFailFile).st_size>0 #True
        # get the kstep and kinc from previous caculation
        if fiberFailure:
            print ('Fiber fail in '+oldjob,' !,stop analysis')
            break

        kstep,kinc=1,1
        # extract failed interface element from result file 
        # in the step of kstep and increment of kinc
        failedInterfaceElements=set()
        with open('E:/UMAT-OUTPUT/'+oldjob+'-resinLayer.txt','r') as fp:
            for line in fp:
                if line.find('Resin_Layer_Broken') !=-1:
                    linesep=line.strip().split()
                    kstep,kinc=int(linesep[4]),int(linesep[6])
                    failedInterfaceElements.add(int(linesep[2]))
        if len(failedInterfaceElements)==0:
            print('no interface resin element fail in this analysis! stop analysis!')
            break
        print('starting ',jobname,'from ',oldjob,' kstep: ',kstep,' kinc:',kinc)
        
        elementSetName='delete%d'%ind
        dataLine=abaqusSetDataLine(failedInterfaceElements)
        with open(restartfname,'w') as fp:
            fp.write(getRestartInp(kstep=kstep,kinc=kinc,
                elementSetName=elementSetName,failedElements=dataLine,
                stepName='Step-%d'%(ind)))
        
        p1=subprocess.Popen(['E:/Program/SIMULIA/Abaqus/Commands/abaqus.bat',
                            'job='+jobname,
                            'oldjob='+oldjob,
                            'user=E:/1730895/Work/Abaqus/UMAT/UMAT_lipeng.for',
                            'interactive'],
                            stdin=subprocess.PIPE,stdout=open('./'+jobname+'.log','w'))
        out,err=p1.communicate('y'.encode())
        if p1.returncode==0:
            print(restartfname,' caculation completed!')
        else:
            print('caculation error,please check logfile!')
            exit(1)

"""
import glob
from abaqus import *
from abaqusConstants import *
from odbAccess import openOdb

for fname in glob.glob('*.odb'):
    o=openOdb(fname,readOnly=True)

realtime=0.0
for stpname in o.steps.keys():
    if stpname.startswith('DUMMY'):
        continue
        
    stp=o.steps[stpname]
    rfs=stp.historyRegions['Node PART-2-1.1'].historyOutputs['RF1'].data
    us=stp.historyRegions['Node PART-2-1.1'].historyOutputs['U1'].data
    for i in range(len(rfs)):
        print stpname,rfs[i][0],rfs[i][0]*(1-realtime)+realtime,us[i][1],rfs[i][1]
    realtime=rfs[i][0]*(1-realtime)+realtime

"""