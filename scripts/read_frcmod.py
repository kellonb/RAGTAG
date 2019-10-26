import numpy as np
import sys
import pandas as pd

'''
usage: 
python read_frcmod.py ala.frcmod
'''

# Get dihedral parameters 
data = []
with open(sys.argv[1]) as file_handle:
    for line in file_handle:
        if line.strip() == 'DIHE': #Begin line where I want data
            break
    for line in file_handle:
        if line.strip() == 'IMPROPER': #end line where I want data
            break 
        data.append(line)
#print (data)

# Get the specific information, dependent on frcmod file, if code break it is here
dih_atomtype=[];divider=[];amplitude=[]
phase=[];periodicity=[]
#CT-CX-C -O     1     0.00000000    0.000   2.0    SCEE=1.2 SCNB=2.0
for line in data:
  a = line.strip()
  dih_at = a[0:15]
  dih_atomtype.append(dih_at)
  d = a[15:21]
  divider.append(d)
  V = a[21:33]
  amplitude.append(V)
  ph = a[33:42]
  phase.append(ph)
  p = a[42:48]
  periodicity.append(p)

#print (dih_atomtype)
#print (dih_atomtype[0])
#print (dih_atomtype[0][3:8])
#print (dih_atomtype[0][6:8])

# index array 
ind = np.arange(0, len(dih_atomtype), 1)

# mid is the two middle atom 
mid = []
for ptr in ind:
  mid1 = dih_atomtype[ptr][3:8]
  mid.append(mid1)
  #print (mid1)

# create a dataframe to work with
genA_df = pd.DataFrame({
          'index' : ind,
          'atomtype' : dih_atomtype,
          'divider' : divider, 
          'amplitude' : amplitude,
          'phase' : phase,
          'periodicity' : periodicity,
          'centralatom' : mid,
          })

#print (genA_df) 
# sort pandas dataframe by the central atom
genA_df = genA_df.sort(['centralatom'])
#print (genA_df)

# create arrays grouped based on central atoms
dihedral={}
for name, group in genA_df.groupby(['centralatom']):
  #print(name)
  #print(group)
  dihedral["{0}".format(name)] = group

# dihedral[centralatom] is the fitting groups
# if a fitting dihedral has a central atom then all dihedrals with that central atom should be fitted 

# examples
print (dihedral['N -CX'])
# get the atomtypes, or any onf the column name in genA_df e.g periodicity, phase etc
print (dihedral['N -CX'].atomtype)

#write amplitude parameters 

