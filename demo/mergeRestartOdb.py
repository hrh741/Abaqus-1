#! -*- Coding: UTF-8 -*-
import glob,os,subprocess,re

"""
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
def mergeRestartOdb(fnameFormat,directory):
    os.chdir(directory)
    flist=flist=sorted(glob.glob('%s/%s'%(directory,fnameFormat)), 
        key=lambda s:int(re.findall('\d+',s)[-1])) #key=os.path.getmtime

    if(len(flist)<2):
        print('必须提供2个以上满足条件的odb!')
        return 
    f1,f2=flist[:2]
    if not (os.path.isfile(f1) or os.path.isfile(f2)):
        print(f1,f2,' are not file')
        return 
    print('start join restart odb!')
    ret=subprocess.call(['E:/Program/SIMULIA/Abaqus/Commands/abaqus.bat', 
                        'restartjoin',
                        'originalodb='+os.path.basename(f1),
                        'restartodb='+os.path.basename(f2),
                        '-history','-copyoriginal'
                        ],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
    if ret==0:
        print('success join %s and %s'%(f1,f2))
    else:
        print('failed to join %s and %s'%(f1,f2))
    
    for fname in flist[2:]:
        if not os.path.isfile(fname):
            break
        print('joining ',fname)
        ret=subprocess.call(['E:/Program/SIMULIA/Abaqus/Commands/abaqus.bat', 
                            'restartjoin',
                            'originalodb=Restart_'+os.path.basename(f1),
                            'restartodb='+os.path.basename(fname),
                            '-history'#,'-copyoriginal'
                            ],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
        if ret==0:
            print('success join %s'%fname)
        else:
            print('failed to join %s'%fname)
            exit(1)

if __name__=="__main__":
    mergeRestartOdb('Job-Delamiantion_I-1mm_Modified-restart-*.odb','E:/Abaqus/RestartAnalysis/0_±30_±60_90s/1mm')
