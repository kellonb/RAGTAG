import numpy as np 
import genA as gA
import os
import time
import sys
import pytraj as pt
# any function with gA come from the genA module/library that I created 

# debug by using debug flag 
def debug(x):
  if DEBUG:
     print (x)

DEBUG=False # change to false if not debugging  
top = '../../ALY.ff14SB.MM.parm7'# topology path 
crd = '../output/500traj.nc'  # trajectory with the 500 structures 
nconf = 500 # number of conformation
backbone = bb = 'alpha'

olist = np.arange(1, nconf+1, 1)
# first convert trajectory to restart file 
traj = pt.iterload(crd, top)
pt.write_traj('./struct.picked.rst7', traj, overwrite=True, options='multi,keepext')

# get the dihedral name and index for each dihedral  
dih_name = gA.get_all_dihedral(top, 'dih_name')
dih_index = gA.get_all_dihedral(top, 'dih_index')
dih_type = gA.get_all_dihedral(top, 'dih_type')
#print (dih_name)
#print (dih_index)
#print (dih_type)

# File that has the atom mask of all dihedrals
fout = open("dihedrals_avg_info.dat", "w")
fout1 = open("dihedrals_all_info.dat", "w")
fout.write("#%-15s%-16s%-16s%-16s%-16s%-16s\n" %('Atom Name','Atom Number','Atom Type','Dihedral Avg','Dihedral min','Dihedral max'))
fout1.write("#%-15s%-16s%-16s " %('Atom Name','Atom Number','Atom Type'))
for frame in np.arange(1, nconf+1, 1):
  fout1.write("Frame%-4d" %(frame))
fout1.write("\n")
alpha_dihedrals=[] # arrays will hold all the dihedral for alpha and opt

for atomname, atomnum, atomtype in zip(dih_name,dih_index,dih_type): 
  # write atom name and atom number index for a given dihedral,
  atna = '%s %s %s %s' %(atomname[0],atomname[1],atomname[2],atomname[3])
  atnu = '%s %s %s %s' %(atomnum[0],atomnum[1],atomnum[2],atomnum[3])
  aty  = '%s %s %s %s' %(atomtype[0],atomtype[1],atomtype[2],atomtype[3])
  fout.write("%-16s%-16s%-16s" %(atna,atnu,aty))
  fout1.write("%-16s%-16s%-16s" %(atna,atnu,aty))
    
  # calculate the dihedral values to get the averages, 1 is the dihnum  
  alpha_data = gA.cal_dih(top, crd, '@%s @%s @%s @%s'%(atomnum[0],atomnum[1],atomnum[2],atomnum[3]), 1)
  # convert from -180 to 180 scale to 0  to 360 scale for doing averages 
  # np.where(condition, x, y) if condition is true give x, if false give y
  # so if alpha_data negative, add 360 else return alpha_data
  a_data = np.where(alpha_data < 0, alpha_data+360, alpha_data)
  min_val = min(a_data);max_val = max(a_data)
  # convert back to -180 to 180
  min_val = np.where(min_val > 180, min_val-360, min_val) 
  max_val = np.where(max_val > 180, max_val-360, max_val) 
  # print some information for debug purpose
  debug(atna);debug(atnu);debug(aty)
  for dihed in alpha_data: 
    fout1.write("%-16.4f" %(dihed))
  
  fout1.write("\n")
  alpha_dihedrals.append(alpha_data) # populate alpha dihedrals 

  # calculate circular averages 
  a_avg = gA.circular_avg(alpha_data)
  debug(alpha_data)
  ## convert averages back into -180 to 180 
  a_avg = np.where(a_avg > 180, a_avg-360, a_avg) 
  fout.write("%-16.4f%-16.4f%-16.4f\n" %(a_avg, max_val, min_val))
fout.close()
fout1.close()


# write chi.rst 

# first load dihedral info file
# since file has columns of width 16, then delimiter is 16
dih_data = np.genfromtxt("dihedrals_avg_info.dat", delimiter=16, dtype=str, autostrip=True)
debug(dih_data)
atom_name = dih_data[:,0]; atom_num = dih_data[:,1] # get the 1st and second column 
atom_type = dih_data[:,2]; avg_dih = dih_data[:,3]; 

debug(len(atom_name))
alpha_dihedrals = np.vstack(alpha_dihedrals)
debug(alpha_dihedrals)
debug(alpha_dihedrals[0])

# these files will have the dihedrals used and not used check them out
outAme = open('ACE_methylrotation.dat', "w")
outNme = open('NME_methylrotation.dat', "w")
outSc = open('SC_methylrotation.dat', "w")
outphi = open('phi.dat', "w")
outpsi = open('phi.dat', "w")
outom = open('omega.dat', "w")
outchi = open('chi.dat', "w")

# files that has the constraint index for orca and gaussian 
orca_index = open('orca_index', "w")
gaus_index = open('gaussian_index', "w")

# Now loop through all structures and write chi.rst
for om in olist:
   # array to hold the r2 and r3 values
   alpha_r23 = []
   # array to hold the atom number
   adih_index = []
   # Below is the actual dihedral values for that structure 
   actual_alpha = alpha_dihedrals[:,om-1] 
   # below is the avergage dihedral values for that structure 
   avg_dih 

   # choose whether to use average or actual dihedral values
   # if it is backbone dihedral then use average
   for atn, ati, avg, act in zip(atom_name, atom_num, avg_dih, actual_alpha):
      atn = atn.split() # split string of 4 atoms into individual atoms for if loops
      ati = ati.split() # split string of 4 atoms number
      ati_rst = '%s,%s,%s,%s' %(ati[0],ati[1],ati[2],ati[3])
      atn_rst = '%s,%s,%s,%s' %(atn[0],atn[1],atn[2],atn[3])

      # middle two atoms will have index 1 and 2 
      # let us start with the ACE methyl rotation 
      if (atn[1] == 'CH3' and atn[2] == 'C') or (atn[1] =='C' and atn[2] == 'CH3'): # methyl rotation 
        # do not include any methyl rotation with H2 or H3 or with an O only using H1
        if (atn[0] == 'H2' or atn[0] == 'H3' or atn[3] == "O"):
          outAme.write("Struct%s notused %-16s %-16s\n" %(om,atn_rst,ati_rst))
        else:
          alpha_r23.append(-5.00) # or use -5 H1 CH3 C N  1 2 5 7 
          adih_index.append(ati_rst)
          outAme.write("Struct%s used %-19s %-16s %-16s\n" %(om,atn_rst,ati_rst,"-5.00"))
          if (om == 1): # only write the index once 
             orca_index.write("{ D %d %d %d %d C }\n" %(int(ati[0])-1, int(ati[1])-1, int(ati[2])-1, int(ati[3])-1))
             gaus_index.write("%s %s %s %s F\n" %(ati[0],ati[1],ati[2],ati[3]))
   # repeating the for loop to keep things sorted, need to figure out a better faster way 
   # let us do omega1 and omega2, we will use avg for all and the omega in the sidechain  
   for atn, ati, avg, act in zip(atom_name, atom_num, avg_dih, actual_alpha):
      atn = atn.split() 
      ati = ati.split()
      ati_rst = '%s,%s,%s,%s' %(ati[0],ati[1],ati[2],ati[3])
      atn_rst = '%s,%s,%s,%s' %(atn[0],atn[1],atn[2],atn[3])
      #if (atn[1] == 'C' and atn[2] == 'N') or (atn[1] == 'N' and atn[2] == 'C'): 
      if (atn[1] == 'C' and atn[2] == 'N') or (atn[1] == 'N' and atn[2] == 'C') or (atn[1] == 'CH' and atn[2] == 'NZ') or (atn[1] == 'NZ' and atn[2] == 'CH'): 
        alpha_r23.append(avg)
        adih_index.append(ati_rst)
        outom.write("Struct%s used %-16s %-16s %-16s\n" %(om,atn_rst,ati_rst,avg))
        if (om == 1): # only write the index once 
           orca_index.write("{ D %d %d %d %d C }\n" %(int(ati[0])-1, int(ati[1])-1, int(ati[2])-1, int(ati[3])-1))
           gaus_index.write("%s %s %s %s F\n" %(ati[0],ati[1],ati[2],ati[3]))
 
   # let us do phi 
   for atn, ati, avg, act in zip(atom_name, atom_num, avg_dih, actual_alpha):
      atn = atn.split() 
      ati = ati.split()
      ati_rst = '%s,%s,%s,%s' %(ati[0],ati[1],ati[2],ati[3])
      atn_rst = '%s,%s,%s,%s' %(atn[0],atn[1],atn[2],atn[3])
      if (atn[1] == 'N' and atn[2] == 'CA') or (atn[1] == 'CA' and atn[2] == 'N'): # phi
        alpha_r23.append(avg)
        adih_index.append(ati_rst)
        outphi.write("Struct%s used %-16s %-16s %-16s\n" %(om,atn_rst,ati_rst,avg))
        if (om == 1): # only write the index once 
           orca_index.write("{ D %d %d %d %d C }\n" %(int(ati[0])-1, int(ati[1])-1, int(ati[2])-1, int(ati[3])-1))
           gaus_index.write("%s %s %s %s F\n" %(ati[0],ati[1],ati[2],ati[3]))
 
   # let us do psi
   for atn, ati, avg, act in zip(atom_name, atom_num, avg_dih, actual_alpha):
      atn = atn.split() 
      ati = ati.split()
      ati_rst = '%s,%s,%s,%s' %(ati[0],ati[1],ati[2],ati[3])
      atn_rst = '%s,%s,%s,%s' %(atn[0],atn[1],atn[2],atn[3])
      if (atn[1] == 'CA' and atn[2] == 'C' ) or (atn[1] == 'C' and atn[2] == 'CA'): # psi central two
        alpha_r23.append(avg)
        adih_index.append(ati_rst)
        outpsi.write("Struct%s used %-16s %-16s %-16s\n" %(om,atn_rst,ati_rst,avg))
        if (om == 1): # only write the index once 
           orca_index.write("{ D %d %d %d %d C }\n" %(int(ati[0])-1, int(ati[1])-1, int(ati[2])-1, int(ati[3])-1))
           gaus_index.write("%s %s %s %s F\n" %(ati[0],ati[1],ati[2],ati[3]))
 
   # let us do Nme methyl rotation 
   for atn, ati, avg, act in zip(atom_name, atom_num, avg_dih, actual_alpha):
      atn = atn.split() 
      ati = ati.split()
      ati_rst = '%s,%s,%s,%s' %(ati[0],ati[1],ati[2],ati[3])
      atn_rst = '%s,%s,%s,%s' %(atn[0],atn[1],atn[2],atn[3])
      if (atn[1] == 'N' and atn[2] == 'CH3') or (atn[1] == 'CH3' and atn[2] == 'N'): # NME methyl rotation 
        if (atn[3] == 'HH32' or atn[3] == 'HH33' or atn[0] == 'H'):
          outNme.write("Struct%s notused %-16s %-16s\n" %(om,atn_rst,ati_rst))
        else:
          alpha_r23.append(60.00) # or use 60 C N CH3 HH31 
          adih_index.append(ati_rst)
          outNme.write("Struct%s used %-19s %-16s %-16s\n" %(om,atn_rst,ati_rst,"60.00"))
          if (om == 1): # only write the index once 
             orca_index.write("{ D %d %d %d %d C }\n" %(int(ati[0])-1, int(ati[1])-1, int(ati[2])-1, int(ati[3])-1))
             gaus_index.write("%s %s %s %s F\n" %(ati[0],ati[1],ati[2],ati[3]))
  # let us do side chain methyl rotation 
   for atn, ati, avg, act in zip(atom_name, atom_num, avg_dih, actual_alpha):
      atn = atn.split()
      ati = ati.split()
      ati_rst = '%s,%s,%s,%s' %(ati[0],ati[1],ati[2],ati[3])
      atn_rst = '%s,%s,%s,%s' %(atn[0],atn[1],atn[2],atn[3])
      if (atn[1] == 'CH' and atn[2] == 'CH3') or (atn[1] == 'CH3' and atn[2] == 'CH'): # NME methyl rotation
        if (atn[3] == 'HHA2' or atn[3] == 'HHA3' or atn[0] == 'OH'):
          outSc.write("Struct%s notused %-16s %-16s\n" %(om,atn_rst,ati_rst))
        else:
          alpha_r23.append(60.00) # or use 60 C N CH3 HH31
          adih_index.append(ati_rst)
          outSc.write("Struct%s used %-19s %-16s %-16s\n" %(om,atn_rst,ati_rst,"60.00"))
          if (om == 1): # only write the index once
             orca_index.write("{ D %d %d %d %d C }\n" %(int(ati[0])-1, int(ati[1])-1, int(ati[2])-1, int(ati[3])-1))
             gaus_index.write("%s %s %s %s F\n" %(ati[0],ati[1],ati[2],ati[3]))

 
   # Now side chains chi1 and chi2. Only do the ones you are fitting here 
   chi_1st = ['N', 'CA', 'CB', 'CG', 'CD', 'OH'] # all the first atoms for each chi you are fitting
   chi_4th = ['CG', 'CD', 'CE', 'NZ', 'CH', 'HZ'] # all the fourth atoms for each chi you are fitting
   for chi_1,chi_4 in zip(chi_1st,chi_4th):

     for atn, ati, avg, act in zip(atom_name, atom_num, avg_dih, actual_alpha):
        atn = atn.split() 
        ati = ati.split()
        ati_rst = '%s,%s,%s,%s' %(ati[0],ati[1],ati[2],ati[3])
        atn_rst = '%s,%s,%s,%s' %(atn[0],atn[1],atn[2],atn[3])
        # chi 
        if (atn[0] == chi_1 and atn[3] == chi_4) or (atn[0] == chi_4 and atn[3] == chi_1): # only so example N and SG, then CA and SD, 
          alpha_r23.append("%.4f"%(act)) # using actual since avg is already formatted
          adih_index.append(ati_rst)
          outchi.write("Struct%s used %-16s %-16s %-16.4f chi1\n" %(om,atn_rst,ati_rst,act))
          if (om == 1): # only write the index once 
             orca_index.write("{ D %d %d %d %d C }\n" %(int(ati[0])-1, int(ati[1])-1, int(ati[2])-1, int(ati[3])-1))
             gaus_index.write("%s %s %s %s F\n" %(ati[0],ati[1],ati[2],ati[3]))
 
 
   # make r1 and r4 arrays, Here I am only using -360 and 360 but if you want to customize u can do it in the loop below
   array_length = len(atom_name)
   r1_array = np.full((array_length), -360)
   r4_array = np.full((array_length), 360)
   rk23_array = np.full((array_length), 1.)
   #rk23 is 1 in this case since we are using the force constant in the input file and not in the chi.rst file
   debug(r1_array)
   debug(r4_array)
   debug(rk23_array)
   # make chi.rst
   gA.make_chirst(adih_index, r1_array, alpha_r23, r4_array, rk23_array, 'struct%s.chi.rst' %(om))
   # make input file for relaxation, 1st argument is filename and second arguments is
   # chirst name 
   gA.relax_backbone("struct%s"%(om), 'struct%s.chi.rst' %(om))
   # run the relaxation 
   os.system("/opt/amber/bin/sander -O -i ./struct%s.in -p ./%s -c ./struct.picked.rst7.%s -ref ./struct.picked.rst7.%s -o ./struct.relax.%s.out -r ./struct.relax.%s.rst7" %(om,top,om,om,om,om))
# now take the relax rst7, make a new traj and then same the pdbs for qm optimization 
alist = open("relaxrst7_list", "w")
for om in olist: 
  alist.write("struct.relax.%s\n" %(om))
alist.close()
# Make a trajectory of only alpha bb or only beta bb conformation  
# Make a trajectory of only alpha bb or only beta bb conformation
gA.create_traj(top, "relax500", nconf, 1, "./relaxrst7_list")
relax_traj = pt.iterload("relax500.nc", top)
pt.write_traj('struct.relax.pdb', relax_traj, overwrite=True, options='multi')
