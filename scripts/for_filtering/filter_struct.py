import linecache
import re
import genA as gA
import sys
import os
import numpy as np
import time

# debug by using debug flag 
def debug(x):
  if DEBUG:
     print (x)

# make a directory alpha and opt
os.system("mkdir alpha")
os.system("mkdir opt")


DEBUG=True # change to false if not debugging  
top = 'CNX.ff14SB.MM.parm7'# topology path 
nconf = 500 # number of conformation
#nconf = 3 # number of conformation
backbone = ['opt', 'alpha']
olist = np.arange(1, nconf+1, 1)

# Convert all xyz to rst7 
# for each structure in each backbone
for bb in backbone:
  # open a run script 
  vdwout = open("%s/%s_filter.sh" %(bb,bb), "w")
  vdwout.write("#!/bin/bash\n")
  vdwout.write("../filter_vdwelec ../%s " %(top))

  # list of all structures
  all_struct = []
  # cd into the alpha or opt directory
  os.chdir("./%s" %(bb))
  # for each structure
  for om in olist:
    # copy over the .xyz file from the QM directory which is one directory up
    string = "%s_struct%s/struct%s" %(bb,om,om)
    debug(string)

    # copy over xyz from QM directory to this directory
    os.system("cp ../../%s.xyz ./" %(string))
    
    # convert xyz to rst7 file 
    gA.XYZ_to_rst7("struct%d.xyz" %(om), "struct%d.rst7" %(om))
    debug("converted xyz to rst7")

    # continue writing each rst7 to the file 
    vdwout.write("struct%d.rst7 " %(om))
    all_struct.append("struct%d.rst7" %(om))
  vdwout.write("\n")
  vdwout.close()
  os.system("bash %s_filter.sh" %(bb))
  # get the suspect list and write out files 
  time.sleep(2) # sleep for 2 seconds 
  suspect_list = gA.suspect_list('Suspect_structures.vdw.dat', 'Suspect_structures.elec.dat', '../%s.suspect_list' %(bb))
  print (suspect_list)
  # get the non suspect list: 
  non_suspect_list = []
  for struct in all_struct:
    if struct not in suspect_list:
      non_suspect_list.append(struct)
  # make trajectories of suspected list and nonsuspect to look at in vmd or to do futher analysis in cpptraj or pytraj
  gA.create_traj("../%s" %(top), "../%s.suspect_structures" %(bb), nconf, 2, suspect_list)
  gA.create_traj("../%s" %(top), "../%s.nonsuspect_structures" %(bb), nconf, 2, non_suspect_list)
  # print a file that just had the suspect and nonsuspect
  out = open("../%s.nonSuspectlist" %(bb), "w")
  for i in non_suspect_list:
    out.write("%s\n" %(i))
 
  # print number of suspects
  print("For backbone %s: You have %s Suspects and %s nonSuspects" %(bb, len(suspect_list), len(non_suspect_list)))
      
  # cd out of the bb directory 
  os.chdir("../")
 

