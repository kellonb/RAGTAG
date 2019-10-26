import os
import linecache
import re
import sys
import time

# this script uses python2
# usage python2 make_orcavpot.py molecule_name number_of_atoms
#       python2 make_orcavpot.py alphabbMet 29

def create_orca_vpot_run(mol_name):
  orcarunvpot_template='''#!/bin/bash

sbatch << EOF
#!/bin/bash 

#SBATCH -n 1
#SBATCH --partition=980
#SBATCH --job-name orcavpot.{mol_name}

export PATH=$PATH:/opt/orca

#this command run orca
/opt/orca/orca_vpot {mol_name}.gbw {mol_name}.scfp ../{mol_name}_esppoints/esppoints_14 esppoints_14.vpot.out
/opt/orca/orca_vpot {mol_name}.gbw {mol_name}.scfp ../{mol_name}_esppoints/esppoints_16 esppoints_16.vpot.out
/opt/orca/orca_vpot {mol_name}.gbw {mol_name}.scfp ../{mol_name}_esppoints/esppoints_18 esppoints_18.vpot.out
/opt/orca/orca_vpot {mol_name}.gbw {mol_name}.scfp ../{mol_name}_esppoints/esppoints_20 esppoints_20.vpot.out
EOF
'''
  molecule={
    "mol_name":mol_name
  }
  with open("run_orcavpot.sh", "w") as f:
    f.write(orcarunvpot_template.format(**molecule))
    f.close()

def create_esp(mol_name, num_atoms):
  createEsp_template='''#!/bin/bash

#editing the .xyz file from ang to bohr
awk '{{printf "%16s%16.7E%16.7E%16.7E\\n","", $2/0.52917721092, $3/0.52917721092, $4/0.52917721092}}' {mol_name}.xyz > {mol_name}.tmp
#delete line 1-2 in .xyz file
sed -i '1,2d' {mol_name}.tmp
for i in `seq 7 10`
do
 I=$((2*i))
 # move column 4(vpot) to column 1 so respgen can read the file
 awk '{{printf "%16.7E%16.7E%16.7E%16.7E\\n",$4, $1, $2, $3}}' esppoints_$I.vpot.out > esppoints.$I.out 
 #editing the .vpot.out, delete line 1 in vpot.out
 sed -i '1d' esppoints.$I.out
done

#combine xyz and esp to one file
cat esppoints.14.out esppoints.16.out esppoints.18.out esppoints.20.out > esppoints.tmp
cat {mol_name}.tmp esppoints.tmp > {mol_name}.esppoints.tmp
t_lines=$(cat esppoints.tmp | wc -l)
echo {num_atom} $t_lines 0 > line1.txt
awk '{{printf "   %s %s    %s\\n",$1, $2, $3}}' line1.txt > line.txt
cat line.txt {mol_name}.esppoints.tmp > {mol_name}.esp
#remove all unncessary files
rm -f blank.txt
rm -f line1.txt
rm -f line.txt
rm -f {mol_name}.esppoints.tmp
rm -f *.tmp
rm -f esppoints.tmp
echo done with .esp
'''
  molecule={
    "mol_name":mol_name,
    "num_atom":num_atom
  }
  with open("create_esp.sh", "w") as f:
    f.write(createEsp_template.format(**molecule))
    f.close()
#..........................................................#
# load your molecule name
mol_name=sys.argv[1]
# load the number of atoms in your molecule
num_atom=sys.argv[2]
# make a directory based on molecule name
os.system("mkdir %s_esp" %(mol_name))
# cd into the directory
os.chdir("./%s_esp" %(mol_name))
# copy your .gbw and .scfp file to the working directory to calculate the vpot
os.system("cp ../%s/*.gbw ./" %(mol_name))
os.system("cp ../%s/*.scfp ./" %(mol_name))
os.system("cp ../%s/*.xyz ./" %(mol_name))
# function that generate your orca input file to calculate the electrostatic potential (vpot)
create_orca_vpot_run(mol_name)
# execute the run script
os.system("bash run_orcavpot.sh")
# wait 30 secs for the command above to be done
# if you get "awk: fatal: cannot open file" error extend time.sleep to 60s
time.sleep(100)
# edit the files to get the write format for fortran
create_esp(mol_name, num_atom)
os.system("bash create_esp.sh")
#go back to main directory
os.chdir("../")
