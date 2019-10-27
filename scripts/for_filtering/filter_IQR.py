import numpy as np 
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

'''
What are the criteria to identify an outlier?
    1) Data point that falls outside of 1.5 times of an interquartile range above the 3rd quartile and below the 1st quartile
    2) Data point that falls outside of 3 standard deviations. we can use a z score and if the z score falls outside of 2 standard deviati
'''

# get time
f = open("check.time.QM", "r")
time = []; struct=[]
for line in f:
   data = line.strip().split()
   structures= data[0]
   hrs = data[6] ; mins = data[8]; secs = data[10];
   ti = (float(hrs) * 60.0) + float(mins) + (float(secs) * 0.01666666667)
   time.append(ti)
   struct.append(structures)
   #print (ti)

# get energy 
energy = np.genfromtxt("E.MP2_kcalpmol", usecols=0, dtype=float)
mean = np.mean(energy)
std = np.std(energy)
print ("Avg Energy (kcal/mol) is %s +/- %s" %(mean,std))


# plot distribution
import seaborn as sns; sns.set(style="ticks")
fig = plt.figure(figsize=(8,6))
sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
gs = GridSpec(4, 4) #no of rows and no of columns in grid
# plot the scatter plot
ax1 = plt.subplot(gs[1:4, 0:3])
ax1.scatter(time, energy, s=2, c="red") #, label="%s" % (label1))
ax1.set_ylabel("Energy (kcal/mol)", fontsize=28)
ax1.set_xlabel("Time MP2 Calculation (mins)", fontsize=28)
ax1.tick_params(axis='both', which='both', labelsize=20, labelcolor='purple')

# plot x histogram   
ax2 = plt.subplot(gs[0,0:3])
ax2.hist(time, bins=np.arange(min(time), max(time), 1), color='blue')
ax2.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax2.yaxis.set_ticklabels([]) #hide the labels 

# plot y histogram 
ax3 = plt.subplot(gs[1:4, 3]) # sharey=ax1) # neeto share axis with 1
ax3.hist(energy, bins=np.arange(min(energy), max(energy), 10), color="purple", orientation='horizontal')
ax3.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) # turn off labels and ticks
ax3.yaxis.set_ticklabels([]) #hide the labels 
plt.tight_layout()
plt.savefig("energy_vs_time.png")


# plot box plot showing even more outliers that the 3 std deviations did not catch 
plt.figure()
sns.boxplot(data=energy, flierprops=dict(markerfacecolor='pink', marker='D'))
plt.savefig("boxplot_visualize_outlier_ALL.png")


#Calculate z scores 
zscore= open("zscore.dat", "w")
reja= open("alpha_reject.dat", "w")
rejo= open("opt_reject.dat", "w")
energy_zscore_pass=[]
struct_zscore_pass=[]
time_zscore_pass=[]
for e,s,ti in zip(energy,struct,time):
   threshold = 3 # 3 std deviation from the mean
   z = (e - mean) / std
   zscore.write("%s  %.4f\n" %(s,z))
   if abs(z) > threshold:
      if s[0:3] == 'alp': # if first three string is alp for alpha
         reja.write("%s.rst7  %.4f Zscore\n" %(s[6:len(s)],z))
      elif s[0:3] == 'opt':
         rejo.write("%s.rst7  %.4f Zscore\n" %(s[4:len(s)],z))
   else:
      energy_zscore_pass.append(e)
      struct_zscore_pass.append(s)
      time_zscore_pass.append(ti)

# plot box plot showing even more outliers that the 3 std deviations did not catch 
plt.figure()
sns.boxplot(data=energy_zscore_pass, flierprops=dict(markerfacecolor='pink', marker='D'))
plt.savefig("boxplot_visualize_outlier_afterZscore.png")

# remove more outliers 
# sort the data
sort = sorted(energy_zscore_pass)
# find 1st and 3rd quartile
q1, q3 = np.percentile(sort, [25,75])
# find the IQR
iqr = q3 - q1 
# find the lower and upper bound
lbound = q1 - (1.5 * iqr)
ubound = q3 + (1.5 * iqr)

energy_IQR_pass=[]
time_IQR_pass=[]
for en,st,tp in zip(energy_zscore_pass,struct_zscore_pass,time_zscore_pass):
   if (en <= lbound) or (en >= ubound):
      z = (en - mean) / std
      if st[0:3] == 'alp': # if first three string is alp for alpha
         reja.write("%s.rst7  %.4f IQR\n" %(st[6:len(st)],z))
      elif st[0:3] == 'opt':
         rejo.write("%s.rst7  %.4f IQR\n" %(st[4:len(st)],z))
   else: 
      energy_IQR_pass.append(en)
      time_IQR_pass.append(tp)

# plot box plot after
plt.figure()
sns.boxplot(data=energy_IQR_pass, flierprops=dict(markerfacecolor='pink', marker='D'))
plt.savefig("boxplot_visualize_outlier_afterIQR.png")

# plot distribution
import seaborn as sns; sns.set(style="ticks")
fig = plt.figure(figsize=(8,6))
sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
gs = GridSpec(4, 4) #no of rows and no of columns in grid
# plot the scatter plot
ax1 = plt.subplot(gs[1:4, 0:3])
ax1.scatter(time_IQR_pass, energy_IQR_pass, s=2, c="red") #, label="%s" % (label1))
ax1.set_ylabel("Energy (kcal/mol)", fontsize=28)
ax1.set_xlabel("Time MP2 Calculation (mins)", fontsize=28)
ax1.tick_params(axis='both', which='both', labelsize=20, labelcolor='purple')
#ax1.set_ylim(, 30)
#ax1.set_xlim(60, 100)

# plot x histogram   
ax2 = plt.subplot(gs[0,0:3])
ax2.hist(time_IQR_pass, bins=np.arange(min(time_IQR_pass), max(time_IQR_pass), 1), color='blue')
ax2.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax2.yaxis.set_ticklabels([]) #hide the labels 
#ax2.set_xticks(np.arange(60, 100, 5))
#ax2.set_xlim(60, 100)

# plot y histogram 
ax3 = plt.subplot(gs[1:4, 3]) # sharey=ax1) # neeto share axis with 1
ax3.hist(energy_IQR_pass, bins=np.arange(min(energy_IQR_pass), max(energy_IQR_pass), 1), color="purple", orientation='horizontal')
ax3.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) # turn off labels and ticks
ax3.yaxis.set_ticklabels([]) #hide the labels 
#ax3.set_yticks(np.arange(-40, 30, 10))
#ax3.set_ylim(-40, 30)
plt.tight_layout()
plt.savefig("energy_vs_time.afterfilter.png")


