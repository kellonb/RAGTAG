import linecache
import re

import sys
import os
import numpy as np


def orca_spEnergy(mol_name, charge, multiplicity):
  spEnergy_template='''! MP2 6-31+G(d,p) TightSCF ENERGY
%id "{mol_name}"
%MaxCore 4000
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


def create_orcarun(mol_name, bb):
  orcarun_template='''#!/bin/bash

path=/cavern/kbelfon/CNX
orca=/cavern/kbelfon/orca_4.0.1.2/orca 

qsub << EOF
#$ -N {bb}_{mol_name}
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -q cpu_short
#$ -P carlosprj

cd \$TMPDIR
mkdir {bb}_{mol_name}/
cd {bb}_{mol_name}/
cp $path/{bb}_{mol_name}/Energy.inp ./
cp $path/{bb}_{mol_name}/{mol_name}.xyz ./
#this command run orca
$orca Energy.inp > Energy.log
cp Energy.log $path/{bb}_{mol_name}/Energy.log
cp Energy.gbw $path/{bb}_{mol_name}/Energy.gbw
cp Energy.prop $path/{bb}_{mol_name}/Energy.prop
cp Energy_property.txt $path/{bb}_{mol_name}/Energy_property.txt
rm -rf ../{bb}_{mol_name}/
exit 0
EOF
'''
  molecule={
    "mol_name":mol_name,
    "bb":bb
  }
  with open("run_orca.sh", "w") as f:
    f.write(orcarun_template.format(**molecule))
    f.close()





# path to where your MM structures are stored
#path="/u/sciteam/belfon/nma/omega"
path="/cavern/kbelfon/CNX/"
name="CNX"

backbone = ['opt', 'alpha']
#olist = np.arange(1, 501, 1)
begin = int(sys.argv[1])
end = int(sys.argv[2])
olist = np.arange(begin,end, 1)
#tlist = np.arange(10, 171, 10)

for bb in backbone:
  for om in olist:
    string="frame.pdb.%d" %(om)
    print (string)
    #os.system("mkdir %s_struct%d" %(bb,om))
    os.chdir("./%s_struct%d" %(bb,om))
    #make_orca("%s" %(string), "struct%s"%(om), 0, 1)
    orca_spEnergy("struct%s"%(om), 0, 2)
    create_orcarun("struct%s"%(om), bb)
    os.system("bash run_orca.sh")
    #orca_spEnergy("struct%s"%(om), 0, 1)
    #create_orcarun("struct%s"%(om))
    #os.system("bash run_orcaEnergy.sh")
    #os.system("/cavern/kbelfon/orca_4.0.1.2/orca %s.inp > %s.log" %(struct,struct))
    os.chdir("../")

