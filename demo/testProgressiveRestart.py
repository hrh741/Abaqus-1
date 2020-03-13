#! /usr/python2
# -*- coding: UTF-8 -*- 

import sys,os,subprocess

inpfname='E:/Abaqus/Restart/Pipes-1970-First.inp'
inpstr=open(inpfname).read()

outdir=os.path.dirname(inpfname)
os.chdir(outdir)

p1=subprocess.Popen(['E:/Program/SIMULIA/Abaqus/Commands/abaqus.bat','job='+os.path.basename(inpfname),'interactive'],
                    stdin=subprocess.PIPE,stdout=subprocess.PIPE)
out,err=p1.communicate('y'.encode())
print(out.decode())