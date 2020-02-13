import numpy as np 
import sys 
import itertools



# open the prod.E file from James calcEnergies and skip every stride 
with open('./prod.E') as Efile: 
   E = np.genfromtxt(itertools.islice(Efile, 0, None, 1))
print (len(E))
#print (E)

# create a dictionary to search for the Energy of a given frame 
frame_num = np.arange(1, len(E)+1, 1)
search_dict = dict(zip(frame_num, E))
#print (search_dict)

# number of dihedrals used to do the clustering 
num_dih=5
# total lines in the dihedral output file 
tot_lines = sum(1 for line in open('./prod.dihCpptraj'))
# nume of lines, excluding the first ndihedrals plus 2 lines 
num_lines = tot_lines - (num_dih+2)
frame_line = np.arange(1, num_lines,2)
'''
#Cluster          1         12 [   0  35  35   0   9  ]
#60753 72966 167026 176440 244210 247479 293530 312993 405705 408614 432235 457397
#Cluster          2         12 [   0  24  35  34  27  ]
#78842 175310 191609 229174 271836 300427 301868 330822 358703 395121 406582 498142
#Cluster          3         11 [  35  24   1  24   1  ]
#2204 43163 136212 169578 191693 315533 325892 395535 436611 454616 458803
'''

st = open('select_struct.dat', 'w')
st.write('CLuster Frame#   Energy \n')
# open the output file from clusterdihedral command in cpptraj 
with open('prod.dihCpptraj', "r") as f:
    # skip the first line, the next ndihedral line and then one more line
    for _ in range(num_dih+2):
        next(f)
    # for every line after this 
    for line in f:
        # split the lines
        l1 = line.strip().split()
        # if the line has the word Cluster get the total frame in that cluster
        if l1[0] == 'Cluster':
            tot_frame = l1[2]
            #print (tot_frame)
        # the next line is the frame numbers for that cluster
        else: 
            # create a list that has all the frame and Energies for that cluster 
            low_frame = []; low_E = []
            # for every frame in that cluster 
            for i in range(int(tot_frame)):
                # get the frame number
                target = l1[i]
                target = int(target)
                #print (target)
                # get the Energy for that frame number
                E = search_dict[target]
                #print (E) 
                # append E and frame to list
                low_frame.append(target); low_E.append(E) 
            # now get the lowest Energy structure in that cluster
            minE = min(low_E)
            # get the index of the lowerst Energy 
            index_min = min(range(len(low_E)), key=low_E.__getitem__) 
            # now get the frame of the lowest Energy 
            minFrame = low_frame[index_min]
            #print ("minE = %s, minFrame = %s" %(minE, minFrame))
            st.write("%s     %s\n" %(minFrame, minE))
            # clear the list for the next cluster
            low_frame.clear(); low_E.clear() 

