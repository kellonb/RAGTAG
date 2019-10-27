import numpy as np 
import genA as gA
import os
import time
import sys
import pytraj as pt
# any function with gA come from the genA module/library that I created 

# debug by using debug flag 
def debug(x):
  if DEBUG:
     print (x)

DEBUG=False # change to false if not debugging  
top = '../../../01.Make_structures_fromLeap/02.heat_pick_relax/alpha_noparam/CNX.ff14SB.MM.parm7'# topology path 

nconf = 500 # number of conformation
#convert xyz to pdb 
olist=np.arange(1, nconf+1, 1)
backbone = ['alpha', 'opt']
# for each backbone
for bb in backbone:
 # for each structure in a given backbone
  for om in olist:
    string = "struct%d" %(om)
    debug(string)
    # cd into the directory created
    os.chdir("./%s_struct%d" %(bb,om))

    # convert xyz to rst7 file 
    gA.XYZ_to_rst7("%s.xyz" %(string), "%s.rst7" %(string))
    debug("converted xyz to rst7")

    # cd out of the directory 
    os.chdir("../")


for bb in backbone:
  # make a directory to hold QM trajectory
  os.mkdir("%s_QMtraj" %(bb))
  os.chdir("%s_QMtraj" %(bb))

  # make a list of the QM structures for backbone
  alist = open("QM_list", "w")
  for om in olist: 
    alist.write("../%s_struct%d/struct%d\n" %(bb,om,om))
  alist.close()
 
  
  # Make a trajectory of only alpha bb or only beta bb conformation
  gA.create_traj("../%s"%(top), "QM500", nconf, 1, "./QM_list")
  # cd out of the directory 
  os.chdir("../")


