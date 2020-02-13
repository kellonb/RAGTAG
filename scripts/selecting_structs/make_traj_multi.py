import warnings
import pytraj as pt
import seaborn as sns
from pytraj import matrix
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import os 

#usage /mnt/raidc2/haichit/anaconda3/bin/python
#home/haichit/anaconda3/lib/python3.5/

#this show the functions in the module pytraj (pt)
#print (dir(pt))
#print (help(pt))

# used for random frame
traj_num=np.loadtxt("random_frame", dtype='int')

#python starts from 0 cpptraj start from 1 is offset by 1
traj_num1 = [x-1 for x in traj_num]

top = "../../ALY.ff14SB.MM.parm7" 
crd = "./prod.skip100.nc"
traj = pt.iterload(crd, top)

# write trajectory
# write for random frame
pt.write_traj('./output/500traj.nc', traj, frame_indices=traj_num1, overwrite=True)


# multiple
#pt.write_traj('./500_alpha/frame.pdb', traj, frame_indices=traj_num1, overwrite=True, options='multi')
