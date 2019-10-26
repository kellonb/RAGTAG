import linecache
import re
import sys

# this script uses python2
# usage python2 makeorca.py pdbfile molecule_name

# remember orca index start from 0 so minus 1 from pdb index

# If you don't need to do constraints the delete psi phi, chi1, chi2 and chi3 line

import sys
import linecache
import re
import os

# this script uses python2
# usage python makeorca.py pdbfile molecule_name charge multiplicity 
#       python2 makeorca.py alphabb.pdb alphabb 0 1
# remember orca index start from 0 so minus 1 from pdb index

def make_orca(pdbfile, mol_name, charge, multiplicity):
    f1 = open('%s.inp' %(mol_name),'w')
    l1="! RHF 6-31G(d) TightSCF TightOpt KeepDens"
    l2="%MaxCore 4000"
    l3='%base' + "\"%s\"" %(mol_name)
    l4="%geom Constraints"
    psi="      { D 6 8 30 32 C }"
    phi="      { D 4 6 8 30 C }"
    l5="      end"
    l6="    end"
    l7=" "
    l8="* xyz %s %s" %(charge, multiplicity)
    f1.write('{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n'.format(l1,l2,l3,l4,psi,phi,l5,l6,l7,l8))
    f2 = open(pdbfile,'r').readlines()
    for line in f2:
        if re.match('TER',line):
            break
        else:
            xyz=line.split()[5:8]
            xyz='   '.join(xyz)
            atom_name=line.split()[2][0]
            atom_name=''.join(atom_name)
            f1.write(atom_name+'   '+xyz+'\n')
    f1.write("*\n")


def create_orcarun(mol_name):
  orcarun_template='''#!/bin/bash

#SBATCH -n 1
#SBATCH --partition=980
#SBATCH --job-name $dir.{mol_name} 

export PATH=$PATH:/opt/orca 

#this command run orca
/opt/orca/orca {mol_name}.inp > {mol_name}.log
'''
  molecule={
    "mol_name":mol_name
  }
  with open("run_orca.sh", "w") as f:
    f.write(orcarun_template.format(**molecule))
    f.close()
#...............................................................#
# Everything start here

# load your pdb file
pdbfile=sys.argv[1]
# load your molecule name
mol_name=sys.argv[2]
# what is your molecule charge
charge=sys.argv[3]
# what is your molecule multiplicity
multiplicity=sys.argv[4]
# make a directory based on molecule name
os.system("mkdir %s" %(mol_name))
# cd into the directory
os.chdir("./%s" %(mol_name))
# function that generate your orca input file
make_orca("../%s" %(pdbfile), mol_name, charge, multiplicity)
# function that create the run script
create_orcarun(mol_name)
# execute the run script
os.system("sbatch run_orca.sh")
#go back to main directory
os.chdir("../")
print "Your Job has been submitted"

