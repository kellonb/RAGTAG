import linecache
import re

import sys
import os
import numpy as np


def make_orca_redo(mol_name, charge, multiplicity):
    f1 = open('%s.redo.inp' %(mol_name),'w')
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
    l8="* xyzfile %s %s %s.xyz" %(charge, multiplicity,mol_name)
    f1.write('{}\n{}\n{}\n{}\n'.format(l5,l6,l7,l8))

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
$orca {mol_name}.redo.inp > {mol_name}.log

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

#begin = int(sys.argv[1])
#end = int(sys.argv[2])
#philist = np.arange(-180, 171, 10)
#philist = np.arange(-180, -10, 10)
#philist = np.arange(-10, 171, 10)
#philist = np.arange(120, 171, 10)
#psilist = np.arange(-180, 171, 10)


phi= int(sys.argv[1])
psi= int(sys.argv[2])
string="%s.phi%s.psi%s" %(name,phi,psi)
print (string)
os.chdir("./phi%s_psi%s" %(phi,psi))
make_orca_redo("%s"%(string), 0, 1)
create_orcarun("%s"%(string))
os.system("bash run_orca.sh")
#os.system("/cavern/kbelfon/orca_4.0.1.2/orca %s.inp > %s.log" %(struct,struct))
os.chdir("../")


