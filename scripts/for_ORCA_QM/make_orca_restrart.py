import linecache
import re

import sys
import os
import numpy as np


def make_orca(pdbfile, mol_name, charge, multiplicity):
    f1 = open('%s.inp' %(mol_name),'w')
    l1="! HF 6-31G(d) TightSCF TightOpt KeepDens"
    l2="%MaxCore 4000"
    l3='%id' + "\"%s\"" %(mol_name)
    l4="%geom Constraints"
    f1.write('{}\n{}\n{}\n{}\n'.format(l1,l2,l3,l4))
    k1 = open("../orca_index", 'r').readlines()
    for line in k1:
       f1.write(line)
    l5="      end"
    l6="    end"
    l7=" "
    l8="* xyz %s %s" %(charge, multiplicity)
    f1.write('{}\n{}\n{}\n{}\n'.format(l5,l6,l7,l8))
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

orca=/cavern/kbelfon/orca_4.0.1.2/orca 
#orca=/u/sciteam/belfon/orca/orca

qsub << EOF
#$ -N {mol_name}
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -q cpu_short
#$ -P carlosprj


#this command run orca
$orca {mol_name}.inp > {mol_name}.log

exit 0
EOF
'''
  molecule={
    "mol_name":mol_name
  }
  with open("run_orca.sh", "w") as f:
    f.write(orcarun_template.format(**molecule))
    f.close()





# path to where your MM structures are stored
#path="/u/sciteam/belfon/nma/omega"
path="/cavern/kbelfon/GLY"
name="GLY"

phi = int(sys.argv[1])
psi = int(sys.argv[2])

string="%s.phi%s.psi%s" %(name,phi,psi)
print (string)
os.system("mkdir phi%s_psi%s" %(phi,psi))
os.chdir("./phi%s_psi%s" %(phi,psi))
os.system("cp ../01.structure/%s.pdb ./" %(string))
make_orca("%s.pdb" %(string), "%s"%(string), 0, 1)
create_orcarun("%s"%(string))
os.system("bash run_orca.sh")
#os.system("/cavern/kbelfon/orca_4.0.1.2/orca %s.inp > %s.log" %(struct,struct))
os.chdir("../")


