import os
import linecache
import re
import sys
import time
import numpy as np

# this script uses python2
# usage python2 makeorca.py molecule_name number_of_atoms
#       python2 makeorca_NC.py alphabbMet 29

def create_esp(mol_name):
  # make a directory based on molecule name for charges
  os.system("mkdir %s_charge" %(mol_name))
  # cd into the directory
  os.chdir("./%s_charge" %(mol_name))
  # copy mol_name.esp to the working directory
  os.system("cp ../%s_esp/*.esp ./" %(mol_name))
  # copy the pdb to the working directory
  os.system("cp ../%s/*.pdb ./" %(mol_name))
  # convert pdb to .ac file using the antechamber on our cluster
  os.system("/opt/amber/bin/antechamber -i %s.pdb -fi pdb -o %s.ac -fo ac" %(mol_name,mol_name) )
  # make all .esp files into one
  os.system("cat %s.esp >> ../charge.esp" %(mol_name))
  #go back to main directory
  os.chdir("../")

def run_resp(first_mol_name, num_mol):
  # only need the first .ac file for the atom name and index
  # create a charge directory
  os.system("mkdir charge_dir")
  # cd into the directory
  os.chdir("./charge_dir")
  # copy the .ac file from the first directory on your dirlist
  os.system("cp ../%s_charge/*.ac ./" %(first_mol_name))
  # create the input files for resp
  os.system("/opt/amber/bin/respgen -i %s.ac -o step1.respin1 -f resp1 -a ../qin.dat -n %s" %(first_mol_name, num_mol))
  os.system("/opt/amber/bin/respgen -i %s.ac -o step2.respin2 -f resp2 -n %s" %(first_mol_name, num_mol))
  # run resp
  os.system("/opt/amber/bin/resp -O -i step1.respin1 -o step1.respout1 -e ../charge.esp -q QIN -t qout_stage1")
 
  os.system("/opt/amber/bin/resp -O -i step2.respin2 -o step2.respout2 -e ../charge.esp -q qout_stage1 -t qout_stage2")


# Everything start here

#get the number of molecule you are fitting
num_mol=sys.argv[1]
#print num_mol

# read the dirlist file
dirlist=np.genfromtxt("dirlist", dtype=str)

# remove charge.esp incase it is there from before
os.system("rm ./charge.esp")
 
# Incase the list only have one structure 
if num_mol == "1":
  mol_name=dirlist
  #print mol_name
  create_esp(mol_name) 
else:
  for mol_name in dirlist:
    #print mol_name
    create_esp(mol_name)

if num_mol == "1":
  first_mol_name=dirlist 
  run_resp(first_mol_name, num_mol)
else:
  first_mol_name=dirlist[0] #only need the first .ac file on your list
  run_resp(first_mol_name, num_mol)

os.system("/opt/amber/bin/antechamber -i ../%s_charge/%s.ac -fi ac -o %s_resp.ac -fo ac -c rc -cf charge_dir/qout_stage2" %(first_mol_name,first_mol_name,first_mol_name))
