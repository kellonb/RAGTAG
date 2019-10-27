import genA
import os
import sys
import numpy as np
import time

if len(sys.argv) < 2:
    error = '''Usage: 

if nbreaks is 1, which is 1 dataset
python create_genA_input.py -i input -at dihedral_info
  
if nbreaks is > 1, which is multiple datasets
python create_genA_input.py -i input -at dihedral_info -top top_info
'''
    print (error)
    sys.exit()

if sys.argv[1] == "-i":

    ##-----------------Get user inputs-------------------------##
    # This section read the input file to get the variables we need
    input_data = []
    # get all lines in the input array
    with open(sys.argv[2], 'r') as my_file:
        for line in my_file:
            line = line.strip()
            info = line.split("=")
            info = line.split(",")
            input_data.append(info)

    # turn array of arrays into one array,
    input_data = np.array(input_data)
    input_data = np.concatenate(input_data)
    input_data = [x for x in input_data if x]  # remove empty strings
    #print (input_data)

    input_arg = {}  # dictionary to hold all arguments
    arguments = ['nbreaks', 'mol_name', 'nconf', 'dihnum', 'flag_mm_read',
                 "mm_read_file", 'flag_qm_read', "qm_read_file", 'top']
    num = np.arange(0, len(input_data), 1)

    # populate the dictionary
    # input_arg[argumemts] = variable
    for i in num:
        a = str(input_data[i])
        b = a.replace(" ", "")  # remove spaces in the string
        #print (b)
        for arg in arguments:
            len_arg = len(arg)
            # this is the part of the string that is the arg
            arg_match = b[0:len_arg]
            if arg_match in arguments:
                #print (arg_match)
                c = b.replace("=", "")  # remove the equal sign
                len_c = len(c)
                input_arg[arg_match] = c[len_arg:len_c]
    # set the inputs
    nbreaks = int(input_arg['nbreaks'])  # number of different data set
    dihnum = int(input_arg['dihnum'])  # number of dihedral fitting
    # read parameter, if read = 1, then structures from file list
    mm_read = int(input_arg['flag_mm_read'])
    # read parameter, if read = 1, then structures from file list
    qm_read = int(input_arg['flag_qm_read'])
    # file with the MM list of structures
    mm_read_file = str(input_arg['mm_read_file'])
    # file with the qm eneries list
    qm_read_file = str(input_arg['qm_read_file'])

    if nbreaks == 1:  # only set top and mol_name and nconf if we have one dataset
        top = input_arg['top']  # toplogy name
        mol_name = input_arg['mol_name']  # molecule name
        nconf = int(input_arg['nconf'])  # number of conformation

    # get the atom mask information for calculating dihedrals etc, and the
    # fitting group and periodicity
    if sys.argv[3] == "-at":
        dih_name = np.genfromtxt(
            sys.argv[4], dtype='str', delimiter=",", skip_header=1, usecols=0)
        dih_mask = np.genfromtxt(
            sys.argv[4], dtype='str', delimiter=",", skip_header=1, usecols=1)
        dih_atomtype = np.genfromtxt(
            sys.argv[4], dtype='str', delimiter=",", skip_header=1, usecols=2)
        fit_group = np.genfromtxt(
            sys.argv[4], dtype='str', delimiter=",", skip_header=1, usecols=3)
        period = np.genfromtxt(
            sys.argv[4], dtype='str', delimiter=",", skip_header=1, usecols=4)
        dih_mask = np.char.strip(dih_mask)
        dih_atomtype = np.char.strip(dih_atomtype)
        fit_group = np.char.strip(fit_group)
        period = np.char.strip(period)
    print (dih_name)
    print (dih_mask)
    print (dih_atomtype)
    print(fit_group)
    print (period)
    if nbreaks > 1:
        if len(sys.argv) < 6:
            print("nbreaks is greater than 1, so use the -top option")
            sys.exit()
        elif sys.argv[5] == "-top":
            # , autostrip=True)
            top_data = np.genfromtxt(sys.argv[6], delimiter=",", dtype=str)
            print (len(top_data))
            top = top_data[:, 0]
            mol_name = top_data[:, 1]  # get the 1st and second column
            nconf = top_data[:, 2]
            #top = np.genfromtxt(sys.argv[6], dtype='str', delimiter=",", skip_header=1, usecols=0)
            #mol_name = np.genfromtxt(sys.argv[6], dtype='str', delimiter=",", skip_header=1, usecols=1)
            #nconf = np.genfromtxt(sys.argv[6], dtype='str', delimiter=",", skip_header=1, usecols=2)
            #weight = np.genfromtxt(sys.argv[6], dtype='str', delimiter=",", skip_header=1, usecols=3)
            if len(top_data) == 2:  # only have three columns
                weight = '1'
            else:
                weight = top_data[:, 3]
                weight = np.char.strip(weight)
            mol_name = np.char.strip(mol_name)
            nconf = np.char.strip(nconf)

#------------run genA functions--------------------------##

if nbreaks == 1:  # only have 1 dataset

    # read from a file with energy, keep energy in hatrees
    if qm_read == 0:  # orca log file was given
        QMeng = genA.extract_QMenergy(qm_read_file, nconf, 0)
    elif qm_read == 1:  # file with structure name was given but orca
        QMeng = genA.extract_QMenergy(qm_read_file, nconf, 1)
    elif qm_read == 2:  # file with energy was given, but hatress
        QMeng = genA.extract_QMenergy(qm_read_file, nconf, 2)
    elif qm_read == 3:  # file with energy was given, but in kcal/mol
        QMeng = genA.extract_QMenergy(qm_read_file, nconf, 3)

    # get MM energy and write trajectory
    if mm_read == 0:
        # write a trajectory based on the number of conformations
        genA.create_traj(top, mol_name, nconf, 0)
        # get MM energy
        MMeng = genA.extract_MMenergy(mol_name, nconf, 0)
    elif mm_read == 1:  # get energy from file
        # write a trajectory based on the number of conformations
        genA.create_traj(top, mol_name, nconf, 1, mm_read_file)
        # get MM energy
        MMeng = genA.extract_MMenergy(mol_name, nconf, 1, mm_read_file, 'out')
    elif mm_read == 2:  # get energy from file
        MMeng = genA.extract_MMenergy(mol_name, nconf, 2, mm_read_file, 'info')

    # calculate dihedral
    dihedral = genA.cal_dih(top, "%s.nc" % (mol_name), dih_mask, dihnum)

    # make genA input
    genA.make_genA_input(dihnum, mol_name, dih_name, dih_atomtype,
                         dihedral, QMeng, MMeng, nbreaks, fit_group, period)

if nbreaks > 1:

    # split up the qm and mm readfile into chunks based on the word break
    search = 'break'
    mmdata = []
    qmdata = []
    mmtemp = []
    qmtemp = []

    for line in open(mm_read_file):
        line = line.strip("\n")
        if line.startswith(search) and mmtemp:
            mmdata.append(mmtemp[:])
            mmtemp = []
        mmtemp.append(line)
    mmdata.append(mmtemp)
    for line in open(qm_read_file):
        line = line.strip("\n")
        if line.startswith(search) and qmtemp:
            qmdata.append(qmtemp[:])
            qmtemp = []
        qmtemp.append(line)
    qmdata.append(qmtemp)

    # remove break from the strings
    for i in np.arange(1, nbreaks, 1):
        qmdata[i].remove('break')
        mmdata[i].remove('break')

    # create temp files
    for i in range(nbreaks):
        t = open("temp%s" % (i), "w")
        for j in qmdata[i]:
            t.write(str(j) + "\n")
    t.close()

    QMeng = []
    for i in range(nbreaks):
        # read from a file with energy, keep energy in hatrees
        if qm_read == 0:  # `orca log file was given
            temp = genA.extract_QMenergy("temp%s" % (i), nconf[i], 0)
            QMeng.append(temp)
        elif qm_read == 1:  # file with structure name was given but orca
            temp = genA.extract_QMenergy("temp%s" % (i), nconf[i], 1)
            QMeng.append(temp)
        elif qm_read == 2:  # file with energy was given, but hatress
            temp = genA.extract_QMenergy("temp%s" % (i), nconf[i], 2)
            QMeng.append(temp)
        elif qm_read == 3:  # file with energy was given, but in kcal/mol
            temp = genA.extract_QMenergy("temp%s" % (i), nconf[i], 3)
            QMeng.append(temp)

    print (len(QMeng[0]))
    print (len(QMeng[1]))

    # delete the temp files
    for i in range(nbreaks):
        os.system("rm temp%s" % (i))

    # create temp files
    for i in range(nbreaks):
        t = open("temp%s" % (i), "w")
        for j in mmdata[i]:
            t.write(str(j) + "\n")
    t.close()

    # get MM energy and write trajectory, for nbreaks > 1 can only do mmread
    # =1 or 2
    MMeng = []
    dihedral = []
    for i in range(nbreaks):
        if mm_read == 1:  # get energy from file
            # write a trajectory based on the number of conformations
            genA.create_traj(top[i], mol_name[i], nconf[i], 1, "temp%s" % (i))
            # get MM energy
            temp = genA.extract_MMenergy(
                mol_name[i], nconf[i], 1, "temp%s" % (i), 'info')
            MMeng.append(temp)
        elif mm_read == 2:  # get energy from file
            temp = genA.extract_MMenergy(
                mol_name[i], nconf[i], 2, "temp%s" % (i), 'info')
            MMeng.append(temp)
    for i in range(nbreaks):
        if dihnum == 1:
            dih = genA.cal_dih(top[i], "./%s.nc" % (mol_name[i]),
                               dih_mask[i * dihnum:(i + 1) * dihnum][0], dihnum)
            dihedral.append(dih)
        else:
            dih = genA.cal_dih(top[i], "./%s.nc" % (mol_name[i]),
                               dih_mask[i * dihnum:(i + 1) * dihnum], dihnum)
            dihedral.append(dih)

    # delete the temp files
    for i in range(nbreaks):
        os.system("rm temp%s" % (i))

    # make genA input
    genA.make_genA_input(dihnum, mol_name, dih_name, dih_atomtype,
                         dihedral, QMeng, MMeng, nbreaks, fit_group, period, weight)
