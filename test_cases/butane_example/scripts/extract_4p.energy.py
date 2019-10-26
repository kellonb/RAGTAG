import numpy as np
import matplotlib.pyplot as plt
import sys
import pytraj as pt
import os
from itertools import islice


def extract_QMenergy(filename, option, nconf):
  # nconf is the number of conf, so number of lines of energy we take
  # option 1: from orca scan
  # option 2: from individual orca files
  if option == 1:
    data=[]
    with open(filename) as file_handle:
      for line in file_handle:
        if line.rstrip() == "The Calculated Surface using the MP2 energy":
          qme = (list(islice(file_handle, nconf)))
          for n in qme:
            data.append(n)  #append lines with both energy and number
          break

    # get only the energy
    QM_eng=[]
    for conf in np.arange(0, nconf, 1):
      temp = data[conf].split()
      temp1 = float(temp[1]) * 627.509608031
      QM_eng.append(temp1)
    return QM_eng

def extract_MMenergy(mol_name, nconf):
  # their mol_name. eg mol_name.000 to name.xxx where xxx is nconf
  crd=[]
  EAMBER=[]
  for n in np.arange(1, nconf+1, 1):
    n = "%03d" %(n,) #make n 000 001 002 etc
    string="%s%s.energy" %(mol_name, n)
    crd.append(string)

  data=[]
  for struct in crd:
    with open("%s.info"%(struct)) as file_handle:
      for line in file_handle:
        data.append(line)

  for line in data:
    if "Etot" in line:
      # split lines
      info_4 = line.strip();info_4 = line.split()

      # populate the total energy
      EAMBER_1 = info_4[2]
      if EAMBER_1 == "=":
        EAMBER_1 = info_4[3];EAMBER.append(float(EAMBER_1))
      elif EAMBER_1 == "*************":
        EAMBER.append('NaN')
      else:
        EAMBER.append(float(EAMBER_1))

  return (EAMBER)

#for p in np.arange(1,5,1): # for four periodicity
  #
f = open("MM_energy4p.dat", "w")
p=4
mol_name="butane"
mol_name3="../05.mm_after_genA/%sp%s_mm_" %(p,mol_name)
nconf=12
MMeng = extract_MMenergy("%s"%(mol_name3), 12)
for i in MMeng:
  f.write(str(i))
  f.write("\n")


