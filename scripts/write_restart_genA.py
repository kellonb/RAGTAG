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
  ph = a[33:42]
  phase.append(ph)
  V = a[21:33]
  if int(ph) == 0:
    V = float(V)
    amplitude.append(V)
  elif int(ph) == 180:
    V = float(V)
    V = -V
    amplitude.append(V)
  p = a[42:48]
  periodicity.append(p)

print (V)
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

print (genA_df) 

#write amplitude parameters 
#create arrays grouped based on atomtype
#genA_df = genA_df.sort(['periodicity'])
dihedral={}
for name, group in genA_df.groupby(['atomtype']):
  #print(name)
  #print(group)
  dihedral["{0}".format(name)] = group

print (genA_df)
#-0.2823 0.23434 0.0786192 -0.658793 -0.077958  0.128863 0.0425491 0.0940626 0.0164876


'''
N- CX-2C-2C       1    0.282300     180   -1
N- CX-2C-2C       1    0.234340       0   -2
N- CX-2C-2C       1    0.078619       0    3
CX-2C-2C-SE       1    0.658793     180   -1
CX-2C-2C-SE       1    0.077958     180   -2
CX-2C-2C-SE       1    0.128863       0    3
2C-2C-SE-CT       1    0.042549       0   -1
2C-2C-SE-CT       1    0.094063       0   -2
2C-2C-SE-CT       1    0.016488       0    3
'''

nchrom = sys.argv[2]

print (dihedral['N -CX-2C-2C'].amplitude)

f = open("amber_prev_ff.restart", "w")
for i in np.arange(0, int(nchrom), 1):
  #for j in genA_df.amplitude:
  for j in dihedral[dih_atomtype].amplitude:
    f.write(str(j) + " ")
  f.write("\n")


'''
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

'''
