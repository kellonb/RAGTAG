import numpy as  np
import pytraj as pt
import sys
import os
from itertools import islice
import warnings
import math 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
def circular_avg(angle_array):
  ''' calculate the average of a dataset of dihedral angles which is circular'''
  a = b = 0 
  conv = 3.14000000 / 180.00000000 # pi/180
  i_conv = 180.00000000 / 3.14000000 # 180/pi 
  for angle in angle_array:
    a += np.sin(angle*conv) 
    b += np.cos(angle*conv)
  a /= len(angle_array)
  b /= len(angle_array)
  c = 90.0 - np.arctan2(b, a) * i_conv
  if (c < 0): 
     c += 360 
  return math.fmod(c, 360) 

def make_chirst(atom_num_array, r1_array, r23_array, r4_array, rk23_array, filename, flag):
  ''' assuming r2 and r3 equals and rk2 and rk3 equals
      flag == S means only a single restraint dihedral angle, flag == M means multiple dihedral restrain angles'''
  fin = open(filename, "w")
  if flag == 'S':   # only one dihedral for chi rst 
    fin.write("&rst iat=%s,r1=%s,r2=%s,r3=%s,r4=%s,rk2=%s,rk3=%s,\n" 
            %(atom_num_array,r1_array,r23_array,r23_array,r4_array,rk23_array,rk23_array))
    fin.write("&end\n")
  elif flag == 'M':
    for atom_num, r1, r23, r4, rk23, in zip (atom_num_array, r1_array, r23_array, r4_array, rk23_array):
      fin.write("&rst iat=%s,r1=%s,r2=%s,r3=%s,r4=%s,rk2=%s,rk3=%s,\n" 
               %(atom_num,r1,r23,r23,r4,rk23,rk23))
      fin.write("&end\n")

def XYZ_to_rst7(filename_xyz, filename_rst7):
  ''' filename_xyz is the name of your xyz file, this function 
      will create a .rst7 file with the name filename_rst7'''
  fin = open('%s' %(filename_xyz), "r")
  fout = open('%s' %(filename_rst7), "w")
  data = fin.readlines()
  tot_lines = len(data) # total lines in the file
  # write second line from xyz as first line for rst7
  fout.write("%s\n" %(data[1].strip()))
  # write first line from xyz as second line for rst7
  fout.write("%6s\n" %(data[0].strip()))
  # loop through every 2 lines and write them as one line in rst7
  for line_num in np.arange(2, tot_lines, 1): #12.7f
    # let us get the line and strip and split it
    # if even number then write on the left   
    if (float(line_num) % 2 == 0):
      l1 = data[line_num].strip().split()
      element1 = l1[0]; X_coord1 = float(l1[1]); Y_coord1 = float(l1[2]); Z_coord1 = float(l1[3])
      fout.write("%12.7f%12.7f%12.7f" %(X_coord1,Y_coord1,Z_coord1))
    # now get the odd num line and do the same
    else: 
      l2 = data[line_num].strip().split()
      element2 = l2[0]; X_coord2 = float(l2[1]); Y_coord2 = float(l2[2]); Z_coord2 = float(l2[3])
      fout.write("%12.7f%12.7f%12.7f\n" %(X_coord2,Y_coord2,Z_coord2))
  fin.close()
  fout.close()
#load MM structures or QM structures to calculate dihedral (structures must be loaded at trajectory)
def cal_dih(top, trajectory, dih_mask, dihnum):
  '''
  top is an amber topology file
  trajectory is an amber trajectorty, it can be both netcdf or ascii
  dih_mask is the dihedral atom number mask  e.g @2 @3 @5 @8
  dihnum is the number of dihedral being fitted
  cal_dih() returns dihedral values for a given mask
  '''
  # first load trajectory

  traj=pt.Trajectory(trajectory, top)
  # if only one dihedral then 
  if dihnum == 1:
    dih_data = pt.dihedral(traj, '%s' %(dih_mask))
    return dih_data
  # we have multiple dihedrals
  else:  
    dih_data=[]
    #calculate all dihedrals
    for mask in dih_mask:
      data = pt.dihedral(traj,  '%s' %(mask)) 
      dih_data.append(data)   
    return dih_data

def get_all_dihedral(top, flag): 
  ''' flag controls what this function will return. It can return 
      dihedral index using dih_index, atom types and atom names'''
  # create a parmed input to get dihedral mask
  with open("dih.in", "w") as f:
    f.write("parm %s" %top + "\n")
    f.write("printDihedrals" + "\n")
    f.write("quit")
    f.close()
  
  #get atoms index for all dihedral from parmed -n ensure the logo is not printed
  os.system("parmed -i dih.in -n > dih_out.log")
  os.system("rm dih.in") # remove unncessary file
  os.system('sed -i  "s/M/ /" dih_out.log') # remove M in the file  
  #os.system('sed -i  "s/I/ /" dih_out.log') # remove I in the file, uncomment if you want the index for improper dihedrals  
  os.system('sed -i  "/I/d" dih_out.log') # delete lines with I in the file  

  # load the file, skip the first five lines and last two line (Done!)
  dihinfo = np.genfromtxt('dih_out.log', dtype=str, skip_header=5, skip_footer=1)
  length, width = dihinfo.shape

  # get the atom index for each dihedral including the improper ones
  atom1_index = [dihinfo[x][0] for x in np.arange(0,length,1)]
  atom2_index = [dihinfo[x][4] for x in np.arange(0,length,1)]
  atom3_index = [dihinfo[x][8] for x in np.arange(0,length,1)]
  atom4_index = [dihinfo[x][12] for x in np.arange(0,length,1)]
  dih_index1 = np.array([atom1_index, atom2_index, atom3_index, atom4_index])
  dih_index1 = dih_index1.transpose()
  #print (dih_index1)

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

  # make dictionarys with dih_index1 as key 
  dih_dict=[]
  for key, value1, value2 in zip(dih_index1, dih_name, dih_type):
    dih_dict.append("%s:%s:%s" %(key, value1, value2))
  #print (dih_dict)  

  # Remove redundancy in the dihedral index
  #dih_index=[]
  #temp = [tuple(row) for row in dih_index1]
  #for i in temp:
  #  if i not in dih_index:
  #    dih_index.append(i)
  #print (dih_index)

  # Remove redundacy and get only dih_name and dih_type based on dih_index 
  dih_name1=[];dih_index=[];dih_type1=[];temp=[] 
  for comb in dih_dict:
    t = comb.split(":") 
    # here we are removing redundacy and stripping unncessary things 
    if t[0] not in temp:
      temp.append(t[0])
      t[0] = (t[0].strip('[]')).replace('\'','')
      t[1] = (t[1].strip('[]')).replace('\'','')
      t[2] = ((t[2].strip('[]')).replace('\'','')).replace(')','')
      dih_index.append(t[0].split())
      dih_name1.append(t[1].split())
      dih_type1.append(t[2].split())
 
 
  # return an array based on the flag
  if flag == 'dih_type': 
    return dih_type1 
  elif flag == 'dih_name':
    return dih_name1
  elif flag == 'dih_index':
    return dih_index


def extract_QMenergy(filename, nconf, Flag):
  '''
  nconf is the number of conf, so number of lines of energy we take
  filename is the name of the file with the QM energies, or with path to QM files 
  Flag options: 
    0: get energy from an orca relaxed scan
    1: get energy from a filename where each line has a path to an orca.log file with QM energy 
    2: get the energies from one file with format molecule_index Qmenergy, but QM energy is in hatrees
    3: get the energies from one file with format molecule_index Qmenergy, but QM energy is kcal/mol
  '''
  
  nconf=int(nconf)
  # get MP2 energy from a relaxed orca scan.
  if Flag == 0:
    data=[]
    with open(filename) as file_handle:
      for line in file_handle:
        # May need to change this line if you you are not using MP2 energy
        if line.rstrip() == "The Calculated Surface using the MP2 energy":
          qme = (list(islice(file_handle, nconf))) # extract nconf of lines after this
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

  elif (Flag == 1 and nconf == 1): # a file name that is an orca enegy .log file
    data = []
    QM_eng = []
    with open(filename) as file_handle:
      for line in file_handle:
        # orca example MP2 TOTAL ENERGY  -494.619011792790
        line = line.rstrip() # remove
        if "MP2 TOTAL ENERGY" in line:
          info = line.split()
          qme = info[3] 
          qme = float(qme) * 627.509608031 # convert to kcal/mol
          QM_eng.append(qme)

    return QM_eng

  elif (Flag == 1 and nconf > 1): # a file name was given with strings matching structure name
    data = []
    QM_eng = []
    file_list = np.genfromtxt(filename, dtype=str, usecols=0)
    for string in file_list:
      file_orca = "%s/orca.log" %(string)
      #file_orca = "%s" %(string)
      with open(file_orca) as file_handle:
        for line in file_handle:
          # orca example MP2 TOTAL ENERGY  -494.619011792790
          line = line.rstrip() # remove
          if "MP2 TOTAL ENERGY" in line:
            info = line.split()
            qme = info[3] 
            qme = float(qme) * 627.509608031 # convert to kcal/mol
            QM_eng.append(qme)
          #  =x86-Linux-G98RevA.7\HF=-968.4663305\MP2=-970.4209847\RMSD=6.597e-09\P
          # \\Version=x86-Linux-G98RevA.7\HF=-968.4697675\MP2=-970.4253704\RMSD=5.
          if "\MP2" in line:
            info =  line.split("MP2",1)[1]                    
            qme = info.split("\\")
            qme = qme[0][1:-1]
            #print (qme)
            qme = float(qme) * 627.509608031 # convert to kcal/mol
            QM_eng.append(qme)

    return QM_eng


  elif Flag == 2: # a file name was given with energies listed in hatree header #structure qm_eng
    QM_eng = np.genfromtxt(filename, dtype='f', usecols=1)
    QM_eng *= 627.509608031
    return QM_eng

  elif Flag == 3: # a file name was given with energies in kcal/mol with header #structure qm_energy
    QM_eng = np.genfromtxt(filename, dtype='f', usecols=1)
    return (QM_eng)


def extract_MMenergy(mol_name, nconf, Flag, filename=None, file_output='info'):
  '''
  mol_name is the molecule name e.g butanol 
  nconf is the number of conformations
  filename is a file with MM energies, default is none
  file_output is the mm output of the file used to get the MM energies. for example an .info or .out file
  Flag options: 
    0: get MM energy using the mol_name. The structures has to be of format
         mol_name.000.energy to mol_name.nconf.energy (e.g butane000.energy to 
                                       1butane012.energy when nconf = 12 
                                       mol_name.nconf (mol_name.012)

    1: get energy from a filename where each line has a path to an MM file with MM energy 
    2: get the energies from one file with format molecule_index MM energy (kcal/mol)
 
  '''
  crd=[]; data=[];  EAMBER=[]
  nconf=int(nconf)

  #---------------------------------------------------------------------------------------#
  if Flag == 0: #structure name was created as mol_name.000 to mol_name.nconf(mol_name.012)
    for n in np.arange(1, nconf+1, 1):
      n = "%03d" %(n,) #make n 000 001 002 etc
      string="%s%s.energy" %(mol_name, n)
      crd.append(string)
    if file_output == 'info':
      for struct in crd:
        with open("%s.info"%(struct)) as file_handle:
          for line in file_handle:
            data.append(line)
      # use the info file 
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
        elif "EAMBER" in line: 
          #split lines
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
    elif file_output == 'out':
      for struct in crd:
        with open("%s.out"%(struct)) as file_handle:
          for line in file_handle:
            line = line.rstrip()
            if "      A V E R A G E S   O V E R" in line:
              mmdata = (list(islice(file_handle, 8))) # extract 8 lines after this
              for n in mmdata:
                data.append(n)  #append lines with both energy and number
              break
      # use the .out file 
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
  #--------------------------------------------------------------------------------------------#

  elif (Flag == 1 and nconf == 1): # a file name was given that is the md out or md info file
    # info file 
    if file_output == 'info':
      with open("%s.info"%(filename)) as file_handle:
        for line in file_handle:
          data.append(line)

      # use the info file 
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
        elif "EAMBER" in line: 
          #split lines
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
    # output file 
    elif file_output == 'out':
      with open("%s.out"%(filename)) as file_handle:
        for line in file_handle:
          data.append(line)
      # use the .out file      
      for line in data:
        if "EAMBER" in line:
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
      
  #--------------------------------------------------------------------------------------------#

  elif (Flag == 1 and nconf > 1): # a file name was given with strings matching structure name
    file_list = np.genfromtxt(filename, dtype=str, usecols=0)
    for string in file_list:
      crd.append(string)
    if file_output == 'info':
      for struct in crd:
        with open("%s.info"%(struct)) as file_handle:
          for line in file_handle:
            data.append(line)
      # use the info file 
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
        elif "EAMBER" in line: 
          #split lines
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
    elif file_output == 'out':
      for struct in crd:
        with open("%s.out"%(struct)) as file_handle:
          for line in file_handle:
            line = line.rstrip()
            if "      A V E R A G E S   O V E R" in line:
              mmdata = (list(islice(file_handle, 8))) # extract 8 lines after this
              for n in mmdata:
                data.append(n)  #append lines with both energy and number
              break
      # use the .out file 
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

  #--------------------------------------------------------------------------------------------#
  elif Flag == 2: # a file name was given with energies with header #structure mm_energy
    MMeng = np.genfromtxt(filename, dtype='f', usecols=1)
    return (MMeng)

#e.g of genA_input, write a function to print this
def print_example():
  example_genA_input = '''
-dih1_name dih1_atom_types
-dih2_name dih2_atom_types
-dihN_name dihN_atom_types
molecule_name dih1_name dih2_name dihN_name
conf1_dihedral_value1 conf1_dihedral_value2 conf1_dihedral_valueN QM_eng(kcal/mol) MM_energy(kcal/mol)
conf2_dihedral_value1 conf2_dihedral_value2 conf2_dihedral_valueN QM_eng(kcal/mol) MM_energy(kcal/mol)
confN_dihedral_value1 confN_dihedral_value2 confN_dihedral_valueN QM_eng(kcal/mol) MM_energy(kcal/mol)
/ seperate different set of structures e.g different backbones
conf1_dihedral_value1 conf1_dihedral_value2 conf1_dihedral_valueN QM_eng(kcal/mol) MM_energy(kcal/mol)
conf2_dihedral_value1 conf2_dihedral_value2 conf2_dihedral_valueN QM_eng(kcal/mol) MM_energy(kcal/mol)
confN_dihedral_value1 confN_dihedral_value2 confN_dihedral_valueN QM_eng(kcal/mol) MM_energy(kcal/mol)
/

HERE IS AN EXAMPLE:

-chi1 N- CX-2C-2C
-chi2 CX-2C-2C-SE
-chi3 2C-2C-SE-CT
MSEalpha chi1 chi2 chi3
-60.0000 -151.3518 -172.2707 -1865195.7795 -20.7497
60.0000 -133.4966 162.7864 -1865192.2533 -17.6159
-20.2001 -175.8127 -19.7117 -1865192.8184 -18.2698
/
MSEbeta chi1 chi2 chi3
-43.2101  107.8337 -168.6532 -1865197.6927 -21.6032
165.1010  159.6474 -95.6575 -1865199.4515 -23.8997
135.1019 +148.1114 +30.8011 -1865197.5170 -21.9696
/
'''
  print (example_genA_input)   


# function that create the input file

#helper functions
def string_format_dihname(dihname):
  dns = " %" + "%d" %(len(dihname)) + "s" # e.g %6s
  return (dns %(dihname))
def string_format_molname(molname): 
  mns = "%" + "%d" %(len(molname)) + "s"
  return (mns %(molname))
def string_format_atype(dihname, atomtype):
  das = "-%s %-11s" %(dihname, atomtype)
  return (das)

# main function 
def make_genA_input(num_dihedral, mol_name, dihname, dihatom, dihedrals, QMeng, MMeng, nbreaks):
  '''
  num_dihedral is the number of dihedrals you are fitting
  mol_name includes directory , only need the molecule name
  dihname is the dihedral name, it is an array 
  dihatom is the atom mask, it is an array
  dihedrals is an array with the dihedral values, obtain from the cal_dih() function 
  QMeng is the QM energy
  MMeng is the MM energy 
  nbreaks is the number of breaks to separate your different datasets 
  '''
  # for cases when you have directory as part of mol name, normally when dealing with 1 dataset
  if nbreaks == 1:
    a = mol_name.split("/")
    mol_name = a[-1]
  
  # open a file name genA_input, this is the input file for RAGTAG 
  with open("genA_input", "w") as genA:

    # for only one dihedral and one dataset (nbreaks=1)
    if num_dihedral == 1 and nbreaks == 1:
      genA.write("%s\n" %(string_format_atype(dihname, dihatom)))
      genA.write("%s" %(string_format_molname(mol_name)))
      genA.write("%s\n" %(string_format_dihname(dihname)))
      for dihedral, QM, MM in zip(dihedrals, QMeng, MMeng):
        genA.write('%10.4f %16.4f %10.4f\n' %(float(dihedral), float(QM), float(MM)) )
      genA.write("/\n")
   
    # for only 1 dihedral, but multiple datasets (nbreaks)
    if num_dihedral == 1 and nbreaks > 1:
      for dih, atom in zip(dihname[0:num_dihedral], dihatom):
        # write the dihname and atom type  
        genA.write("%s\n" %(string_format_atype(dih, atom)))
      # get length of molname strings for the input to make the format perfect
      genA.write("%s" %(string_format_molname(mol_name[0])))
      # loop through and write the dihedrals name
      print (dihname)
      for dihn in dihname[0:num_dihedral]:
        genA.write("%s" %(string_format_dihname(dihn)))
      genA.write("\n")
      data = {}
      for j in np.arange(0, (len(dihedrals)), 1): #dihedrals is an arrray of arrays holding dihedrals
         data[j] = dihedrals[j]
      data[(len(dihedrals))] = QMeng
      data[(len(dihedrals))+1] = MMeng

      # loop through dihedrals
      all_values=[] # hold the dihedral values, QM and MM0
      # dihedral values are separated based on datasets so using for i in nbreaks
      for nb in np.arange(0, nbreaks, 1):
         all_values.append(data[nb]) # append the dihedral values
      all_values.append(data[nbreaks]) # append the QM which
      all_values.append(data[nbreaks+1]) # append the MM
    
      # for every dataset 
      for nb in np.arange(0, nbreaks, 1): 
        for d1 in np.arange(0, (len(all_values[nb])), 1): 
           genA.write('%10.4f ' %(all_values[nb][r1][d1])) #dihedrals
           genA.write('%16.4f ' %(all_values[nbreaks][nb][d1])) #QM
           genA.write('%10.4f' %(all_values[nbreaks+1][nb][d1])) #MM
        genA.write("\n")
        genA.write("/\n")
        if nb != (nbreaks -1): # dont write it for the last nbreaks
          genA.write("%s" %(string_format_molname(mol_name[nb+1])))
          for dihn in dihname[num_dihedral*(nbreaks-1):num_dihedral*nbreaks]:
            genA.write("%s" %(string_format_dihname(dihn)))
          genA.write("\n")

    # for multiple dihedrals, but only one dataset (nbreaks)
    if num_dihedral > 1 and nbreaks == 1:
      for dih, atom in zip(dihname, dihatom):
        # write the dihname and atom type  
        genA.write("%s\n" %(string_format_atype(dih, atom)))
      # get length of molname strings for the input to make the format perfect
      genA.write("%s" %(string_format_molname(mol_name)))
      for dihn in dihname:
        genA.write("%s" %(string_format_dihname(dihn)))
      genA.write("\n")
      data = {}
      for j in np.arange(0, (len(dihedrals)), 1): #dihedrals is an arrray of arrays holding dihedrals
        data[j] = dihedrals[j]
      data[(len(dihedrals))] = QMeng 
      data[(len(dihedrals))+1] = MMeng
      # loop through dihedrals and print 
      all_values=[] # hold the dihedral values, QM and MM0 
      for dval in np.arange(0, num_dihedral, 1):
         all_values.append(data[dval]) # append the dihedral values 
      all_values.append(data[num_dihedral]) # append the QM which
      all_values.append(data[num_dihedral+1]) # append the MM 
      #get the number of columns and number of rows  
      num_cols = (len(all_values))
      num_rows = (len(all_values[0]))
      # for each row in each columns
      for r1 in np.arange(0, num_rows, 1):
        for c1 in np.arange(0, num_dihedral, 1):
           genA.write('%10.4f ' %(all_values[c1][r1])) #dihedrals
        genA.write('%16.4f ' %(all_values[num_dihedral][r1])) #QM
        genA.write('%10.4f' %(all_values[num_dihedral+1][r1]))
        genA.write("\n")

    #----------------------------------------------------------------------------------------#
    # for multiple dihedrals, but with multiple dataset
    if num_dihedral > 1 and nbreaks > 1:
      for dih, atom in zip(dihname[0:num_dihedral], dihatom):
        # write the dihname and atom type  
        genA.write("%s\n" %(string_format_atype(dih, atom)))
      # get length of molname strings for the input to make the format perfect
      genA.write("%s" %(string_format_molname(mol_name[0])))
      # loop through and write the dihedrals name 
      print (dihname) 
      for dihn in dihname[0:num_dihedral]:
        genA.write("%s" %(string_format_dihname(dihn)))
      genA.write("\n")
      data = {}
      for j in np.arange(0, (len(dihedrals)), 1): #dihedrals is an arrray of arrays holding dihedrals
         data[j] = dihedrals[j]
      data[(len(dihedrals))] = QMeng 
      data[(len(dihedrals))+1] = MMeng
 
      # loop through dihedrals 
      all_values=[] # hold the dihedral values, QM and MM0 
      # dihedral values are separated based on datasets so using for i in nbreaks
      for nb in np.arange(0, nbreaks, 1):
         all_values.append(data[nb]) # append the dihedral values 
      all_values.append(data[nbreaks]) # append the QM which
      all_values.append(data[nbreaks+1]) # append the MM 

      #get the number of columns and number of rows  
      num_cols = (len(all_values)) # number of columns is dih for nbreaks, then QM then MM 
      num_rows1 = (len(all_values[0])) # number of dihderals in a dataset (nbreak)
      num_rows2 = (len(all_values[nbreaks])) # number of QM and MM values which is based on the number of dataset (nbreak)
      num_dih_rows = (len(all_values[0][0]))
      # for example 2 datasets (nbreaks = 2) with 3 set of dihedrals will have 4 columns (dataset1 dataset2 QM and MM) 
      # the first two cols will have 3 rows each(num_rows1=3), where each row have all of the dihedrals (row1 is dihedral 1, row2 is dihedrals)
      # the last two cols (QM and MM) will have 2 rows(num_rows2=2), where each row have the QM and MM valuses for each dataset

      for nb in np.arange(0, nbreaks, 1): 
        for d1 in np.arange(0, (len(all_values[nb][0])), 1):
          for r1 in np.arange(0, num_rows1, 1):
            genA.write('%10.4f ' %(all_values[nb][r1][d1])) #dihedrals 
          genA.write('%16.4f ' %(all_values[nbreaks][nb][d1])) #QM 
          genA.write('%10.4f' %(all_values[nbreaks+1][nb][d1])) #MM
          genA.write("\n")
        genA.write("/\n")
        if nb != (nbreaks -1): # dont write it for the last nbreaks 
          genA.write("%s" %(string_format_molname(mol_name[nb+1])))
          for dihn in dihname[num_dihedral*(nbreaks-1):num_dihedral*nbreaks]:
            genA.write("%s" %(string_format_dihname(dihn)))
          genA.write("\n")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


# If we do this we have to write the restart file for users based on
# their mol_name. eg mol_name.000 to name.xxx where xxx is nconf
def create_traj(top, mol_name, nconf, Flag, filename=None):
  ''' top is topology
      mol_name is molecule name, must be set of structures number mol_name001 to mol_namexn1
      nconf is conformation number
      Flag:
        0 if using index moleculename.000 to moleculename.nconf (moleculename.012)
        1 filename where each line has a file path to the rst7 file
  '''
  crd=[] # holds the filename for list of structures
  nconf=int(nconf)

  if Flag == 0:
    for n in np.arange(1, nconf+1, 1):
      n = "%03d" %(n,) #make n 000 001 002 etc
      string="%s%s.energy.rst7" %(mol_name, n)
      crd.append(string)
  # read from a file name that have list of names with no extension 
  elif Flag == 1:
    file_list = np.genfromtxt(filename, dtype=str, usecols=0)
    for string in file_list:
      string = "%s.rst7" %(string)
      crd.append(string)
  # read from a list of structures
  elif Flag == 2:
    for string in filename:
      crd.append(string)

  # write trajectory
  traj = pt.iterload(crd, top)
  pt.write_traj('%s.nc' %(mol_name), traj, overwrite=True)
  print ("trajectory name %s.nc is completed" %(mol_name))

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
def get_QM_2Dscan(nconf, orcalog_path):
  QMeng=[]
  int(nconf)
  for k in np.arange(1, nconf+1, 1):
    k  = "%03d" %(k,)
    QMeng_temp = extract_QMenergy("%s.%s/orca.log" %(k), 1, 12)
    QMeng.append(QMeng_temp)
  #put in 1D space
  QMeng = np.concatenate(QMeng, axis=0)

def get_MM_2Dscan(nconf, MMpath, mol_name):
  MMeng=[]
  int(nconf)
  for k in np.arange(1, nconf+1, 1):
    k  = "%03d" %(k,)
    MMeng_temp = extract_MMenergy("%s%s_"%(mol_name,k), 12)
    MMeng.append(MMeng_temp)
  MMeng = np.concatenate(MMeng, axis=0)

def run_MMenergy(top, maxcyc, chirst_name, infilename, outfilename, min_filename):
  #infilename is the name for your amber input eg. infilename.rst7, 
  #outfilename is the name for your amber outputfile eg. filename.rst7, filename.out
  min_template='''minimizing QM structure
&cntrl
imin = 1, maxcyc = {maxcyc}, ntb = 0, igb=0,
ntx = 1, ntwr = 100, ntpr = 100 
cut = 999.0
nmropt = 1,
/
&wt type='END'   /
DISANG=./{chirst_name}
END
'''
  info={
      "chirst_name":chirst_name,
      "maxcyc":maxcyc 
  }
  with open("%s.in" %(min_filename), "w") as f:
    f.write(min_template.format(**info))
    f.close()
  os.system("$AMBERHOME/bin/sander -O -i ./%s.in -p ./%s -c ./%s.rst7 -ref ./%s.rst7 -o ./%s.out -r ./%s.rst7 -inf ./%s.info" %(min_filename,top,infilename,infilename,outfilename, outfilename,outfilename))
 
# specific for proteins relaxing the back bone to an average chi.rst (chirst_name)
def relax_backbone(filename, chirst_name):
    relaxbb_template='''relax backbone to average min.in
&cntrl
imin = 1, ntx = 1, maxcyc = 100000,
ntwr=100, ntpr=100,
cut = 99.0,
ntb=0, igb = 0, saltcon = 0.0,
nmropt = 1,
/
&wt TYPE='REST',istep1=0,istep2=499,value1=100.,value2=100.,
/
 &wt type='REST',istep1=500,istep2=999,value1=500.,value2=500.,
/
 &wt TYPE='REST',istep1=1000,istep2=1499,value1=1000.,value2=1000.,
/
 &wt TYPE='REST',istep1=1500,istep2=1999,value1=5000.,value2=5000.,
/
 &wt TYPE='REST',istep1=2000,value1=10000.,
/
&wt type='END'  /
DISANG=./{chirst_name}
END
'''
    info={
      "chirst_name":chirst_name
    }
    with open("%s.in" %(filename), "w") as f:
      f.write(relaxbb_template.format(**info))
      f.close()

#specific for sorting fileterd structures from the filter_vdwelec.c program
def sort_name(name, struct_array, nconf, ext):
  # if you have a set of structures name1 to nameN 
  # for example struct1 to struct500
  # nconf is number of structures
  # extension such as rst7
  sorted_list = []
  for num in np.arange(1, nconf+1, 1):
    for struct in struct_array:
      # remove redundancy first
      if struct in sorted_list:
        continue
      elif "%s%d.%s"%(name,num,ext) == struct:
        sorted_list.append(struct)
  return (sorted_list)

def suspect_list(vdw_file, elec_file, output_filename):

  vdwin = np.genfromtxt(vdw_file, usecols=0, dtype=str)
  elecin = np.genfromtxt(elec_file, usecols=0, dtype=str)

  # output file with a list of suspect structures
  output = open(output_filename , "w")
  output.write("#These structures are suspects\n")
  output.write("#%18s %6s %6s\n" %("Structure", "vdw", "elec"))
  y = "yes"; n ="no"
  suspect_array = []
  
  if (elecin.size == 0 and vdwin.size == 0): 
     print("There is no suspects, Use all structures")

  # No suspects for elecin and more than one suspects for vdwin 
  elif (elecin.size == 0 and vdwin.size > 1):
    for struct1 in vdwin:
       output.write("%18s %6s %6s\n" %(struct1, y, n))
       suspect_array.append(struct1)
  # No susepect for elecin and vdwin only have one suspect
  elif (elecin.size == 0 and vdwin.size == 1):
     struct1 = vdwin 

  # No suspects for vdwin and more than one suspects for elec
  elif (vdwin.size == 0 and elecin.size > 1) :
    for struct2 in elecin:
       output.write("%18s %6s %6s\n" %(struct2, n, y))
       suspect_array.append(struct2)
  # No suspects for vdwin and only one suspects for elec
  elif (vdwin.size == 0 and elecin.size == 1) :
     struct2 = elecin

  # If there is 1 suspect structure in elec but multiple supect structures in vdw
  elif (elecin.size == 1 and vdwin.size > 1): # only 1 structure so cannot loop
    struct2 = elecin
    print (struct2)
    for struct1 in vdwin: 
       if struct1 == struct2: 
            output.write("%18s %6s %6s\n" %(struct1, y, y))
            suspect_array.append(struct1)
       else:
            if (struct1 not in suspect_array):
              # failed the vdw criteria
              output.write("%18s %6s %6s\n" %(struct1, y, n))
              suspect_array.append(struct1)
            elif (struct2 not in suspect_array):
              # failed the elec criteria
              output.write("%18s %6s %6s\n" %(struct2, n, y))
              suspect_array.append(struct2)

  # if there is 1 suspect structure in vdw and multiple in elec
  elif (vdwin.size == 1 and elecin.size > 1): # only 1 structure so cannot loop
    struct1 = vdwin
    print (struct1)
    for struct2 in elecin: 
       if struct2 == struct1: 
            output.write("%18s %6s %6s\n" %(struct1, y, y))
            suspect_array.append(struct1)
       else:
            if (struct1 not in suspect_array):
              # failed the vdw criteria
              output.write("%18s %6s %6s\n" %(struct1, y, n))
              suspect_array.append(struct1)
            elif (struct2 not in suspect_array):
              # failed the elec criteria
              output.write("%18s %6s %6s\n" %(struct2, n, y))
              suspect_array.append(struct2)


  # multiple structures in all      
  if ( elecin.size > 1 and vdwin.size > 1):
    for struct1 in vdwin: 
      # loop through the suspect list based on elec
      for struct2 in elecin:
         # if the same structures exist in both vdw and elec file then the reason 
         # it is a suspect structure is because it failed both vdw and elec criteria
         if struct1  == struct2: 
            output.write("%18s %6s %6s\n" %(struct1, y, y))
            suspect_array.append(struct1)
         else:
            if (struct1 not in suspect_array):
              # failed the vdw criteria
              output.write("%18s %6s %6s\n" %(struct1, y, n))
              suspect_array.append(struct1)
            elif (struct2 not in suspect_array):
              # failed the elec criteria
              output.write("%18s %6s %6s\n" %(struct2, n, y))
              suspect_array.append(struct2)

  output.close()
  # load the output file and sort and make a new file
  list_a = []; list_b = []  #list for dictionary 
  with open(output_filename) as f1:
    for line in f1:
      line = line.strip()
      # if line do not have a # sign in front
      if not line.startswith("#"):
        dict_key = line.split()[0]
        dict_elem = "%6s %6s" %(line.split()[1], line.split()[2])
        list_a.append(dict_key); list_b.append(dict_elem)  
  dict_struct = dict(zip(list_a, list_b))
  sortoutput = open("%s_sorted" %(output_filename), 'w')
  sortoutput.write("#These structures are sorted suspects\n")
  sortoutput.write("#%18s %6s %6s\n" %("Structure", "vdw", "elec"))
  new_key = sort_name('struct',list_a, 500, 'rst7')
  for key in new_key:
    value = dict_struct.get(key)
    sortoutput.write("%18s %s\n" %(key, value))
  sortoutput.close()
  # return a sorted list
  sorted_suspect_array = sort_name('struct', suspect_array, 500, 'rst7')
  return (sorted_suspect_array)


# create lib files 
def create_lib(mol_name, leaprc, frcmod, top_name): 
  '''
  mol_name is the molecule name
  leaprc is the leaprc for the forcefield 
  frcmod we will use to create the top
  top_name is the topology file name
  '''
  lib_template='''source {leaprc}
{RES} = loadmol2 ./{mol_name}.mol2
loadAmberParams ./{frcmod}
check {RES}
saveoff {RES} {mol_name}.lib
saveamberparm {RES} {top_name}.parm7 {top_name}.rst7
quit
'''
  info={
      "mol_name":mol_name, 
      "top_name":top_name, 
      "frcmod":frcmod, 
      "leaprc":leaprc,
      "RES":mol_name[0:3]
  }
  with open("make_off.in", "w") as f:
    f.write(lib_template.format(**info))
    f.close()
  os.system("$AMBERHOME/bin/tleap -f make_off.in")

# minimization step before heating 
def run_min_gas(top, maxcyc, inpcrd, chirst_name=None):
  #infilename is the name for your amber input eg. infilename.rst7, 
  #outfilename is the name for your amber outputfile eg. filename.rst7, filename.out
  if chirst_name != None: 
    min_template='''minimizing initial structure
&cntrl
imin = 1, maxcyc = {maxcyc}, igb=0, drms=0.001,
ntx = 1, ntwr = 100, ntpr = 100, ntwx=100,  
ntc = 1, ntf = 1, ntb = 0,  
cut = 999.0, irest=0, ig = -1, 
ntt = 3, gamma_ln = 1.0, temp0 = 500.0, 
ntr = 1, restraint_wt = 10, restraintmask="!@H"
nmropt = 1,
/
&ewald
  eedmeth=5 
/
&wt type='END'   /
DISANG=./{chirst_name}
END
'''
    info={
      "chirst_name":chirst_name,
      "maxcyc":maxcyc 
    }
    with open("1min.in", "w") as f:
      f.write(min_template.format(**info))
      f.close()
    os.system("$AMBERHOME/bin/sander -O -i ./1min.in -p ./%s -c ./%s -ref ./%s -o ./1min.out -r 1min.r" %(top,inpcrd,inpcrd))
  else: 
    min_template='''minimizing initial structure
&cntrl
imin = 1, maxcyc = {maxcyc}, igb=0, drms=0.001,
ntx = 1, ntwr = 100, ntpr = 100, ntwx=100,  
ntc = 1, ntf = 1, ntb = 0,  
cut = 999.0, irest=0, ig = -1, 
ntt = 3, gamma_ln = 1.0, temp0 = 500.0, 
ntr = 1, restraint_wt = 10, restraintmask="!@H",
/
&ewald
  eedmeth=5 
/
'''
    info={
      "maxcyc":maxcyc 
    }
    with open("1min.in", "w") as f:
      f.write(min_template.format(**info))
      f.close()
    os.system("$AMBERHOME/bin/sander -O -i ./1min.in -p ./%s -c ./%s -ref ./%s -o ./1min.out -r 1min.r" %(top,inpcrd,inpcrd))

def run_heat_gas(top, nstlim, inpcrd, chirst_name=None):
  #infilename is the name for your amber input eg. infilename.rst7, 
  #outfilename is the name for your amber outputfile eg. filename.rst7, filename.out
  if chirst_name != None: 
    heat_template='''heat to 500
&cntrl
imin = 0, nstlim = {nstlim}, dt = 0.001,
ntx = 1, ntwr = 500, ntpr = 500, ntwx=500,  
ntc = 2, ntf = 2, ntb = 0, igb = 0, 
cut = 999.0, irest=0, ig = -1, dielc = 4,
ntt = 3, gamma_ln = 1.0, tempi = 100.0, temp0 = 500.0, 
ntr = 1, restraint_wt = 10, restraintmask="!@H"
nmropt = 1,
/
&ewald
  eedmeth=5 
/
&wt type='TEMP0', istep1=0, istep2={half}
        value1=100.0, value2=500.0,             ! Varies target temperature
/
&wt type='TEMP0', istep1={half}, istep2={nstlim},
        value1=500.0, value2=500.0,             ! Varies target temperature
/
&wt type='END'   /
DISANG=./{chirst_name}
END
'''
    info={
      "chirst_name":chirst_name,
      "nstlim":nstlim,
      "half":nstlim/2
    }
    with open("2mdheat.in", "w") as f:
      f.write(heat_template.format(**info))
      f.close()
    os.system("$AMBERHOME/bin/sander -O -i ./2mdheat.in -p ./%s -c 1min.r -ref ./1min.r -o ./2md.out -r 2md.r" %(top))
  else: 
    heat_template='''heat to 500K
&cntrl
imin = 0, nstlim = {nstlim}, dt = 0.001,
ntx = 1, ntwr = 500, ntpr = 500, ntwx=500,  
ntc = 2, ntf = 2, ntb = 0, igb = 0,  
cut = 999.0, irest=0, ig = -1, dielc = 4,  
ntt = 3, gamma_ln = 1.0, tempi = 100.0, temp0 = 500.0, 
ntr = 1, restraint_wt = 10, restraintmask="!@H",
/
&ewald
  eedmeth=5 
/
&wt type='TEMP0', istep1=0, istep2={half}
        value1=100.0, value2=500.0,             ! Varies target temperature
/
&wt type='TEMP0', istep1={half}, istep2={nstlim},
        value1=500.0, value2=500.0,             ! Varies target temperature
/
&wt type='END'   /
'''
    info={
      "nstlim":nstlim,
      "half":nstlim/2
    }
    with open("2mdheat.in", "w") as f:
      f.write(heat_template.format(**info))
      f.close()
    os.system("$AMBERHOME/bin/sander -O -i ./2mdheat.in -p ./%s -c 1min.r -ref ./1min.r -o ./2md.out -r 2md.r" %(top))

def run_hightemp_gas(top, nstlim, inpcrd, chirst_name=None):
  #infilename is the name for your amber input eg. infilename.rst7, 
  #outfilename is the name for your amber outputfile eg. filename.rst7, filename.out
  if chirst_name != None: 
    prod_template='''heat to 500
&cntrl
imin = 0, nstlim = {nstlim}, dt = 0.001,
ntx = 5, irest = 1, ntwr = 10000, ntpr = 1000, ntwx=1000,  
ntc = 2, ntf = 2, ntb = 0, igb = 0, 
cut = 999.0, ig = -1, dielc = 4,
ntt = 3, gamma_ln = 1.0, temp0 = 500.0, 
ntr = 1, restraint_wt = 10, restraintmask="!@H"
nmropt = 1,
ioutfm = 1, ntwx = 2, nscm = 1000, 
/
&ewald
  eedmeth=5 
/
DISANG=./{chirst_name}
END
'''
    info={
      "chirst_name":chirst_name,
      "nstlim":nstlim,
    }
    with open("prod.in", "w") as f:
      f.write(prod_template.format(**info))
      f.close()
    os.system("$AMBERHOME/bin/sander -O -i ./prod.in -p ./%s -c 2md.r -o ./prod.out -x prod.nc -r prod.r" %(top))
  else: 
    prod_template='''heat to 500K
&cntrl
imin = 0, nstlim = {nstlim}, dt = 0.001,
ntx = 5, irest = 1, ntwr = 10000, ntpr = 1000, ntwx=1000,  
ntc = 2, ntf = 2, ntb = 0, igb = 0,  
cut = 999.0, ig = -1, dielc = 4,  
ntt = 3, gamma_ln = 1.0, temp0 = 500.0,
ioutfm = 1, ntwx = 2, nscm = 1000, 
/
&ewald
  eedmeth=5 
/
'''
    info={
      "nstlim":nstlim,
    }
    with open("prod.in", "w") as f:
      f.write(prod_template.format(**info))
      f.close()
    os.system("$AMBERHOME/bin/sander -O -i ./prod.in -p ./%s -c 2md.r -o ./prod.out -x prod.nc -r prod.r" %(top))

