import linecache
import re

import sys
import os
import numpy as np


def orca_spEnergy(mol_name, charge, multiplicity):
  spEnergy_template='''! MP2 6-31+G(d,p) TightSCF ENERGY
%id "{mol_name}"
%MaxCore 4000
%scf
HFTyp ROHF
end
* xyzfile {charge} {multiplicity} {mol_name}.xyz
'''
  molecule={
    "mol_name":mol_name,
    "charge":charge,
    "multiplicity":multiplicity
  }
  with open("Energy.inp", "w") as f1:
    f1.write(spEnergy_template.format(**molecule))
    f1.close()


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
$orca Energy.inp > Energy.log

exit 0
EOF
'''
  molecule={
    "mol_name":mol_name
  }
  with open("run_orcaEnergy.sh", "w") as f:
    f.write(orcarun_template.format(**molecule))
    f.close()



# path to where your MM structures are stored
#path="/u/sciteam/belfon/nma/omega"
path="/cavern/kbelfon/CNX"
name="CNX"

backbone = ['opt', 'alpha']
olist = np.arange(500, 501, 1)

for bb in backbone:
  for om in olist:
    string="struct%s" %(om)
    print (string)
    os.chdir("./%s_struct%d_rohf" %(bb,om))
    orca_spEnergy("%s"%(string), 0, 2)
    create_orcarun("%s"%(string))
    os.system("bash run_orcaEnergy.sh")
    os.chdir("../")


