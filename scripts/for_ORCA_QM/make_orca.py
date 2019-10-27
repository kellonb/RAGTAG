import linecache
import re

import sys
import os
import numpy as np


def make_orca(pdbfile, mol_name, charge, multiplicity):
    f1 = open('%s.inp' %(mol_name),'w')
    l1="! UHF 6-31G(d) TightSCF TightOpt KeepDens"
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
path="/cavern/kbelfon/CNX"
name="CNX"

backbone = ['opt', 'alpha']
#bb = 'alpha'
#olist = np.arange(1, 501, 1)
olist = np.arange(1, 100, 1)
#tlist = np.arange(10, 171, 10)

for bb in backbone:
  for om in olist:
    string="frame.pdb.%d" %(om)
    print (string)
    os.system("mkdir %s_struct%d" %(bb,om))
    os.chdir("./%s_struct%d" %(bb,om))
    #os.system("cp ../500_%s1/%s ./" %(bb,string))
    os.system("cp ../%s/struct.relax.pdb.%s ./" %(bb,om))
    make_orca("struct.relax.pdb.%s" %(om), "struct%s"%(om), 0, 2)
    create_orcarun("struct%s"%(om))
    os.system("bash run_orca.sh")
    #os.system("/cavern/kbelfon/orca_4.0.1.2/orca %s.inp > %s.log" %(struct,struct))
    os.chdir("../")


#begin = int(sys.argv[1])
#end = int(sys.argv[2])
