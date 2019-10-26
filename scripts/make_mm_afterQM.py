import numpy as np 
import genA as gA
import os
import time
import sys
import pytraj as pt 

# debug by using debug flag 
def debug(x):
  if DEBUG:
     print (x)

DEBUG=False # change to false if not debugging  
QM_direc = '../02.QM/CNX_cavern/CNX/' # QM directory with the xyz files
top = 'CNX.ff14SB.MM.parm7'# topology path 
top_zero = 'CNX.ff14SB.CYXc1c2.parm7'# topology path c1c2 14SB
nconf = 500 # number of conformation
#nconf = 3 # number of conformation
backbone = ['opt', 'alpha']
olist = np.arange(1, nconf+1, 1)
# maximum cycle for minimization of QM structures
maxcyc = 10000
#maxcyc = 1
# number of dihedral fitting 
num_dih = 4
chi3_at = '2C-S -S -CT'  ; chi3_am = '@11 @14 @15 @16'  ; chi3_an = 'CB SG SD CE'
chi4_at = 'S -S -CT-MC'  ; chi4_am = '@14 @15 @16 @19'  ; chi4_an = 'SG SD CE CZ'
chi5_at = 'S -CT-MC-MC'  ; chi5_am = '@15 @16 @19 @40'  ; chi5_an = 'SD CE CZ CK'
chi5p_at = 'S -CT-MC-YC'  ; chi5p_am = '@15 @16 @19 @20'  ; chi5p_an = 'SD CE CZ CH'
dihname = ['chi3', 'chi4', 'chi5', 'chi5p']
dihatom_fit = [chi3_at, chi4_at, chi5_at, chi5p_at]
dihatom_mask = [chi3_am, chi4_am, chi5_am, chi5p_am]
dihatom_name = [chi3_an, chi4_an, chi5_an, chi5p_an]

# make the dihedral_info file for genA
dih_info = open("dihedral_info", "w")
dih_info.write("%-10s%-18s%-18s\n" %("#dih_name", "dih_mask", "dih_atom_type"))
# for each dataset 
for bb in backbone: 
  for dih_name, mask, atype in zip(dihname, dihatom_mask, dihatom_fit):
    dih_info.write("%-10s%-18s%-18s\n" %("%s,"%(dih_name), "%s,"%(mask), "%s,"%(atype)))

# for each backbone
for bb in backbone:
 # for each structure in a given backbone
  for om in olist:
    string = "%s/%s_struct%s/struct%s" %(QM_direc,bb,om,om)
    debug(string)
    # create directory with nomenclature like alpha_struct1
    os.system("mkdir -p %s_struct%d" %(bb,om))

    # cd into the directory created
    os.chdir("./%s_struct%d" %(bb,om))

    # copy over rst7 from QM directory to this directory
    os.system("cp ../%s/filter/%s/struct%s.rst7 ./" %(QM_direc,bb,om))

    # cd out of the directory 
    os.chdir("../")

# first make a list of alpha and beta rst7 that are non suspect for making a trajectory 
frame_list = []
for bb in backbone:
  alist = open("%s_rst7_list" %(bb), "w")
  nonsuspect_list = np.genfromtxt('%s/filter/%s.nonSuspectlist' %(QM_direc,bb), usecols=0, dtype=str)
  rej_list=np.genfromtxt('%s/%s_reject.dat' %(QM_direc,bb), usecols=0, dtype=str)
  print(nonsuspect_list)
  print(rej_list)
  #reject_list= np.concatenate((nonsuspect_list,rej_list))
  for struct in nonsuspect_list:
    if struct not in rej_list: 
      frame_list.append("%s" %(struct))
      alist.write("%s_%s/%s\n" %(bb,struct[0:-5],struct[0:-5])) #0:-5 is to remove the .rst7 from the string
  alist.close()

# Make a trajectory of only alpha bb or only beta bb that are conformation  
for bb in backbone:
  gA.create_traj(top, "%s_nonsuspect" %(bb), nconf, 1, "./%s_rst7_list" %(bb))
 
# get the dihedral name and index for each dihedral  
dih_name = gA.get_all_dihedral(top, 'dih_name')
dih_index = gA.get_all_dihedral(top, 'dih_index')
dih_type = gA.get_all_dihedral(top, 'dih_type')
#print (dih_name)
#print (dih_index)
#print (dih_type)

# for each backbone
for bb in backbone:
  # read in suspect lists 
  suspect = np.genfromtxt('%s/filter/%s.suspect_list_sorted' %(QM_direc,bb), usecols=0, dtype=str)
  rej_list=np.genfromtxt('%s/%s_reject.dat' %(QM_direc,bb), usecols=0, dtype=str)
  debug(suspect)
  # energy file 
  eng = open("%s_energies.dat" %(bb),  "w")
  eng.write("%10s %17s %17s %17s %36s\n" %("#Structure", "QM", "MM0", "MM14SB", "dihedrals"))
  # create File that has the atom mask of all dihedrals
  fout = open("%sdihedrals_avg_info.dat" %(bb), "w")
  fout1 = open("%sdihedrals_all_info.dat" %(bb), "w")
  fout.write("#%-15s%-16s%-16s%-16s%-16s%-16s\n" %('Atom Name','Atom Number','Atom Type','Dihedral Avg','Dihedral min','Dihedral max'))
  fout1.write("#%-15s%-16s%-16s " %('Atom Name','Atom Number','Atom Type'))
  for frame in frame_list:
    fout1.write("Frame %s  " %(frame))
  fout1.write("\n")

  alpha_dihedrals=[] # arrays will hold all the dihedral for alpha and opt
  for atomname, atomnum, atomtype in zip(dih_name,dih_index,dih_type): 
    # write atom name and atom number index for a given dihedral,
    atna = '%s %s %s %s' %(atomname[0],atomname[1],atomname[2],atomname[3])
    atnu = '%s %s %s %s' %(atomnum[0],atomnum[1],atomnum[2],atomnum[3])
    aty  = '%s %s %s %s' %(atomtype[0],atomtype[1],atomtype[2],atomtype[3])
    fout.write("%-16s%-16s%-16s" %(atna,atnu,aty))
    fout1.write("%-16s%-16s%-16s" %(atna,atnu,aty))
    
    # use only the nonsuspect.nc    
    crd = '%s_nonsuspect.nc' %(bb) 
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
    alpha_dihedrals.append(alpha_data) # populate alpha or opt dihedrals 

    # calculate circular averages (we can use the avg or actual but they should be the same) 
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
  dih_data = np.genfromtxt("%sdihedrals_avg_info.dat" %(bb), delimiter=16, dtype=str, autostrip=True)
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
  outphi = open('phi.dat', "w")
  outpsi = open('psi.dat', "w")
  outom = open('omega.dat', "w")
  outchi = open('chi.dat', "w")

  # Now loop through all structures in nosuspect list and write chi.rst
  #nonsuspect_list = np.genfromtxt('./%s/filter/%s.nonSuspectlist' %(QM_direc,bb), usecols=0, dtype=str)
  #rej_list=np.genfromtxt('%s/%s_reject.dat' %(QM_direc,bb), usecols=0, dtype=str)
  #print (reject_list)
  #klist = np.arange(1, len(reject_list)+1, 1)
  nonsuspect_list = np.genfromtxt('%s/filter/%s.nonSuspectlist' %(QM_direc,bb), usecols=0, dtype=str)
  om = 0
  for struct in nonsuspect_list:
    #print (struct)
    if struct not in rej_list: 
      #print(struct)
      om = om + 1
      #print (om)
      # cd into the directory created
      stt = struct[0:-5]
      os.chdir("./%s_%s" %(bb,stt))
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
          alpha_r23.append("%.4f"%(act)) # if u want to use averge do avg, but I am using actual. They should be no difference 
          adih_index.append(ati_rst)
          outAme.write("Struct%s used %-19s %-16s %-16s\n" %(stt,atn_rst,ati_rst,act))

      # repeating the for loop to keep things sorted, need to figure out a better faster way 
      # let us do omega1 and omega2, we will use avg for all 
      for atn, ati, avg, act in zip(atom_name, atom_num, avg_dih, actual_alpha):
        atn = atn.split() 
        ati = ati.split()
        ati_rst = '%s,%s,%s,%s' %(ati[0],ati[1],ati[2],ati[3])
        atn_rst = '%s,%s,%s,%s' %(atn[0],atn[1],atn[2],atn[3])
        if (atn[1] == 'C' and atn[2] == 'N') or (atn[1] == 'N' and atn[2] == 'C'): 
          alpha_r23.append("%.4f"%(act)) # if u want to use averge do avg, but I am using actual. They should be no difference
          adih_index.append(ati_rst)
          outom.write("Struct%s used %-16s %-16s %-16s\n" %(stt,atn_rst,ati_rst,avg))
 
      # let us do phi 
      for atn, ati, avg, act in zip(atom_name, atom_num, avg_dih, actual_alpha):
        atn = atn.split() 
        ati = ati.split()
        ati_rst = '%s,%s,%s,%s' %(ati[0],ati[1],ati[2],ati[3])
        atn_rst = '%s,%s,%s,%s' %(atn[0],atn[1],atn[2],atn[3])
        if (atn[1] == 'N' and atn[2] == 'CA') or (atn[1] == 'CA' and atn[2] == 'N'): # phi
          alpha_r23.append("%.4f"%(act)) ## if u want to use averge do avg, but I am using actual. They should be no difference
          adih_index.append(ati_rst)
          outphi.write("Struct%s used %-16s %-16s %-16s\n" %(stt,atn_rst,ati_rst,avg))
   
      # let us do psi
      for atn, ati, avg, act in zip(atom_name, atom_num, avg_dih, actual_alpha):
        atn = atn.split() 
        ati = ati.split()
        ati_rst = '%s,%s,%s,%s' %(ati[0],ati[1],ati[2],ati[3])
        atn_rst = '%s,%s,%s,%s' %(atn[0],atn[1],atn[2],atn[3])
        if (atn[1] == 'CA' and atn[2] == 'C' ) or (atn[1] == 'C' and atn[2] == 'CA'): # psi central two
          alpha_r23.append("%.4f"%(act)) # if u want to use averge do avg, but I am using actual. They should be no difference
          adih_index.append(ati_rst)
          outpsi.write("Struct%s used %-16s %-16s %-16s\n" %(stt,atn_rst,ati_rst,avg))
 
      # let us do Nme 
      for atn, ati, avg, act in zip(atom_name, atom_num, avg_dih, actual_alpha):
        atn = atn.split() 
        ati = ati.split()
        ati_rst = '%s,%s,%s,%s' %(ati[0],ati[1],ati[2],ati[3])
        atn_rst = '%s,%s,%s,%s' %(atn[0],atn[1],atn[2],atn[3])
        if (atn[1] == 'N' and atn[2] == 'CH3') or (atn[1] =='CH3' and atn[2] == 'N'): # methyl rotation 
          alpha_r23.append("%.4f"%(act)) # if u want to use averge do avg, but I am using actual. They should be no difference
          adih_index.append(ati_rst)
          outNme.write("Struct%s used %-19s %-16s %-16s\n" %(stt,atn_rst,ati_rst,act))
 
      # Now side chains  
      chi_1st = ['N', 'CA', 'CB', 'SG', 'SD', 'SD'] # all the first atoms for each chi you are fitting
      chi_4th = ['SG', 'SD', 'CE', 'CZ', 'CH', 'CK'] # all the fourth atoms for each chi you are fitting
      for chi_1,chi_4 in zip(chi_1st,chi_4th):
        for atn, ati, avg, act in zip(atom_name, atom_num, avg_dih, actual_alpha):
          atn = atn.split() 
          ati = ati.split()
          ati_rst = '%s,%s,%s,%s' %(ati[0],ati[1],ati[2],ati[3])
          atn_rst = '%s,%s,%s,%s' %(atn[0],atn[1],atn[2],atn[3])
          # chis 
          if (atn[0] == chi_1 and atn[3] == chi_4) or (atn[0] == chi_4 and atn[3] == chi_1): # only so example N and SG, then CA and SD, 
            alpha_r23.append("%.4f"%(act)) # using actual since avg is already formatted
            adih_index.append(ati_rst)
            outchi.write("Struct%s used %-16s %-16s %-16.4f chi1\n" %(stt,atn_rst,ati_rst,act))
 
      # make r1 and r4 arrays, Here I am only using -360 and 360 but if you want to customize u can do it in the loop below
      array_length = len(atom_name)
      r1_array = np.full((array_length), -360)
      r4_array = np.full((array_length), 360)
      #rk23 is 2500 in this case since we are using the force constant in the chi.rst and not in the input file
      rk23_array = np.full((array_length), 25000.)
      debug(r1_array)
      debug(r4_array)
      debug(rk23_array)
      # make chi.rst
      gA.make_chirst(adih_index, r1_array, alpha_r23, r4_array, rk23_array, '%s.chi.rst' %(stt), 'M')
      # make input file for relaxation, 1st argument is filename and second arguments is
      # chirst name 

      # only do MM on structures that are not suspect, change here if you want to do suspects
      # do the minimization with zero dihedral MMzero
      
      gA.run_MMenergy("../%s"%(top_zero), maxcyc, '%s.chi.rst' %(stt), '%s' %(stt), '%s.zero' %(stt), "%s.in" %(stt)) # regular min
      # do the minimization with 14SB dihedral 
      gA.run_MMenergy("../%s"%(top), maxcyc, '%s.chi.rst' %(stt), '%s' %(stt), '%s.14SB' %(stt), "%s.in" %(stt)) # regular min
     
      # get QM energy for this
      qmeng = gA.extract_QMenergy("../%s/%s_%s/Energy.log" %(QM_direc,bb,stt), 1, 1)
      eng.write("%10s %18.4f" %("%s" %(stt), qmeng[0]))
      # get MM0 energy 
      mm0eng = gA.extract_MMenergy('%s.zero' %(stt), 1, 1, '%s.zero' %(stt), 'out')
      eng.write("%18.4f" %(mm0eng[-1])) # it will take the last energy in the output file which is maxcyc
      # get MMenergy
      mmeng = gA.extract_MMenergy('%s.14SB' %(stt), 1, 1, '%s.14SB' %(stt), 'out')
      eng.write("%18.4f" %(mmeng[-1])) # take the last energy
      # get dihedral 
      now_traj = pt.iterload('%s.zero.rst7' %(stt), "../%s" %(top_zero))
      for mask in dihatom_mask:
        dih = gA.cal_dih("../%s" %(top_zero), now_traj, mask, 1)
        eng.write("  %6.2f " %(dih[0]))
      eng.write("\n")

      # cd out of the directory 
      os.chdir("../")

    # let us make a special note for suspect directory 
    suspect_list = np.genfromtxt('./%s/filter/%s.suspect_list_sorted' %(QM_direc,bb), usecols=0, dtype=str)
    rej_list=np.genfromtxt('%s/%s_reject.dat' %(QM_direc,bb), usecols=0, dtype=str)
    reject_list= np.concatenate((suspect_list,rej_list))
    fk = open("%s_list_of_structures_not_used.dat" %(bb), "w")
    for sus_struct in reject_list: 
      # cd into the directory created
      os.chdir("./%s_%s" %(bb,sus_struct[0:-5]))
      fk1 = open("suspect_structure.dat", "w")
      fk1.write("This is a suspect structure")
      fk.write("%s_%s\n" %(bb,sus_struct[0:-5]))
      # cd out of the directory 
      os.chdir("../")

 


