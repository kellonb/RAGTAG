import os
import linecache
import re
import sys 

# this script uses python3 because James wrote the generate script
# usage python make_esppoints.py molecule_name density 
#       python make_esppoints.py alphabbMet 30
#...................................................................# 
# load your molecule name
mol_name=sys.argv[1]
# load the density value to control the number of grid points
# remember the total hardcoded limit of 99,999 in RESP
density=sys.argv[2]
# make a directory based on molecule name
os.system("mkdir %s_esppoints" %(mol_name))
# cd into the directory
os.chdir("./%s_esppoints" %(mol_name))
# run james python script to generate the gridpoints. It uses the QM log file 
os.system("python3 ../jenerate_resp_points.py ../%s/%s.log %s" %(mol_name,mol_name,density))
#go back to main directory
os.chdir("../")
