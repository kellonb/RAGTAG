import numpy as  np
import pytraj as pt
import sys
import os
from itertools import islice

#This program used pytraj, 
#This program needs a trajectory, file with MM energies and a file with QM energies
#note your trajectory frames should match the order of your energies
#use the multiple m command if you are fitting different molecules, or conformation

#functions
#load MM structures or QM structures to calculate dihedral (structures must be loaded at trajectory)
def cal_dih(top, trajectory, dih_mask):
  #dih_mask is first and last atomname
  #maybe add more functionality
  dihout_array=[]
 
   # load trajectory
  traj=pt.Trajectory(trajectory, top) 
  
  #calculate dihedrals
  dih_out_array = pt.dihedral(traj,  dih_mask) #range360=True, if I need it to be 0 to 360deg

  dih = []
  for dih1 in dih_out_array:
    dih.append(dih1)
  return dih

def get_all_dihedral(top):
  # create a parmed input to get dihedral mask 
  with open("dih.in", "w") as f:
    f.write("parm %s" %top + "\n")
    f.write("printDihedrals" + "\n")
    f.write("quit")
    f.close()

  #get atoms index for all dihedral from parmed -n ensure the logo is not printed
  os.system("parmed -i dih.in -n > dih_out.log")
  os.system("rm dih.in") # remove unncessaryu file

  # load the file, skip the first six lines and last line (Done!)
  dihinfo = np.genfromtxt('dih_out.log', dtype=str, skip_header=6, skip_footer=1)
  length, width = dihinfo.shape

  # get the atom index for each dihedral including the improper ones
  atom1_index = [dihinfo[x][0] for x in np.arange(0,length,1)]
  atom2_index = [dihinfo[x][4] for x in np.arange(0,length,1)]
  atom3_index = [dihinfo[x][8] for x in np.arange(0,length,1)]
  atom4_index = [dihinfo[x][12] for x in np.arange(0,length,1)]
  dih_index1 = np.array([atom1_index, atom2_index, atom3_index, atom4_index])
  dih_index1 = dih_index1.transpose()
  #print dih_index1

  # Here I remove redundancy in the dihedral index
  dih_index=[]
  temp = [tuple(row) for row in dih_index1]
  for i in temp:
    if i not in dih_index:
      dih_index.append(i)

  #print dih_index
  # get the atom name for each dihedral including the improper ones
  atom1_name = [dihinfo[x][1] for x in np.arange(0,length,1)]
  atom2_name = [dihinfo[x][5] for x in np.arange(0,length,1)]
  atom3_name = [dihinfo[x][9] for x in np.arange(0,length,1)]
  atom4_name = [dihinfo[x][13] for x in np.arange(0,length,1)]
  dih_name = np.array([atom1_name, atom2_name, atom3_name, atom4_name])
  dih_name = dih_name.transpose()
  #print dih_name

  # get the atom type for each dihedral including the improper ones
  atom1_type = [dihinfo[x][3] for x in np.arange(0,length,1)]
  atom2_type = [dihinfo[x][7] for x in np.arange(0,length,1)]
  atom3_type = [dihinfo[x][11] for x in np.arange(0,length,1)]
  atom4_type = [dihinfo[x][15] for x in np.arange(0,length,1)]
  dih_type = np.array([atom1_type, atom2_type, atom3_type, atom4_type])
  dih_type = dih_type.transpose()
  #print dih_type

  #calculate all dihedral at a time using atom index
  '''
  dih_data=[]
  for dih in dih_index:
    data = pt.dihedral(traj,'@%s @%s @%s @%s' %(dih[0],dih[1],dih[2],dih[3]) ) #, range360=True)
    dih_data.append(data)
  dih_data=np.vstack(dih_data)
  '''
  return dih_type or dih_name 

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

  #if option == 2:   
  #  data=[]
  #  for struct in orca_struct_list:
  #    with open("orca.log") as file_handle:
  #      for line in file_handle:


'''
The Calculated Surface using the MP2 energy
   1.00000000 -157.96305502
   2.00000000 -157.96053810
   3.00000000 -157.95765966
   4.00000000 -157.96026532
   5.00000000 -157.96209966
   6.00000000 -157.95817518
   7.00000000 -157.95354356
   8.00000000 -157.95817361
   9.00000000 -157.96209994
  10.00000000 -157.96028663
  11.00000000 -157.95766300
  12.00000000 -157.96053381


Timings for individual modules:
'''    

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

'''
-chi1 N- CX-2C-2C
-chi2 CX-2C-2C-SE
-chi3 2C-2C-SE-CT
MSEalpha chi1 chi2 chi3
 -60.0000 -151.3518 -172.2707 -1865195.7795 -20.7497
 +60.0000 -133.4966 +162.7864 -1865192.2533 -17.6159
 -20.2001 -175.8127 -19.7117 -1865192.8184 -18.2698
/
MSEopt chi1 chi2 chi3
 -43.2101 +107.8337 -168.6532 -1865197.6927 -21.6032
 +165.1010 +159.6474 -95.6575 -1865199.4515 -23.8997
 +135.1019 +148.1114 +30.8011 -1865197.5170 -21.9696
/
'''

#create input file
def make_genA_input(num_dihedral, dihname, dihatom, dihedrals, QMeng, MMeng):
  #num_dihedral is the number of dihedral
  with open("genA_input", "w") as genA:
    if num_dihedral == 1:
      genA.write("-%s" %(dihname) + " " + "%s" %(dihatom) + "\n")
      genA.write("%s" %(mol_name) + " " + "%s" %(dihname) + "\n")
      for dihedral, QM, MM in zip(dihedral1, QMeng, MMeng):
        #print (str(dihedral) + " " + QM + " " + MM)
        genA.write('%10.4f %16.4f %10.4f\n' %(float(dihedral), QM, float(MM)) )
    elif num_dihedral > 1:
      for dih, atom in zip(dihname, dih_atom): 
         genA.write("-%s" %(dih) + " " + "%s" %(atom) + "\n")
         genA.write("%s" %(mol_name) + " " + "%s " %(dihname) + "\n")

#***********************************************************************************#
#Program start here

# get name of molecule
mol_name="butane_mm0_"


# get name of dihedral
dihname="dih1" #e.g phi, psi, etc
dihatom="c3-c3-c3-c3" # write a function to get this from frcmod file

#for more than one dihedral dihname and dihatom is an array of string
#dihname=['dih1', 'dih2', 'phi1']
#dihatom=['c3-c3-c3-c3', 'c2-c2-c2-c3', 'c1-c2-c3-c3']


# get QM and MM energy    
QMeng = extract_QMenergy("orca.log", 1, 12)
MMeng = extract_MMenergy("%s"%(mol_name), 12)

#print (QMeng)
#print (MMeng)


#calculate dihedrals based on mask that user give
# this function does it for all  frame and return an array
dihedral1 = cal_dih("butane.parm7", "%s.nc" %(mol_name), '@2 @1 @6 @9')
#dihedral2 = cal_dih("butane.parm7", "%s.nc" %(mol_name), '@3 @4 @4 @4')
#print (dihedral1)

make_genA_input(1, dihname, dihatom, dihedral1, QMeng, MMeng)  


