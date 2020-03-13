import sys
sys.path.append('E:/1730895/Work/')

import numpy as np
from Abaqus import  Compliance,Laminate
import os

fname=os.path.join(os.getcwd(),sys.argv[1])
print(os.getcwd(),fname,os.path.realpath(fname),os.path.isfile(fname))
