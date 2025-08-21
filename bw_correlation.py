import os
import sys
import re
import bbi
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import gaussian_kde


bed_file=sys.argv[1]
bw_file1=sys.argv[2]
bw_file2=sys.argv[3]
out_file=sys.argv[4]

# bed_file="/lustre/fengyang/Xiaona/SC_G4/20220831/link_peaks/G4_48_rep1_2.narrowPeak"
# bw_file1="G4_48h_rep1.sorted.rmdup.filtered.RPGCnorm.bw"
# bw_file2="G4_48h_rep2.sorted.rmdup.filtered.RPGCnorm.bw"
# out_file=""
file_G4=open(bed_file)
bed_content=[]
for line in file_G4:
	line=line.strip("\n")
	array=line.split("\t")[0:3]
	bed_content.append(array)

file_G4.close()
filename_wt=bw_file1
filename_ko=bw_file2
bins_number=200

with bbi.open(filename_wt) as f1,bbi.open(filename_ko) as f2:
	array_wt=[]
	array_ko=[]
	for j in bed_content:
		matrix_wt=[]
		matrix_ko=[]
		try:
			x = f1.fetch(j[0], int(j[1]), int(j[2]), bins=bins_number)
			matrix_wt.append(x)
		except:
			continue
		try:
			x = f2.fetch(j[0], int(j[1]), int(j[2]), bins=bins_number)
			matrix_ko.append(x)
		except:
			continue
		matrix_wt=np.asarray(matrix_wt)
		matrix_ko=np.asarray(matrix_ko)
		matrix_wt=np.average(matrix_wt,axis=1)
		matrix_ko=np.average(matrix_ko,axis=1)
		# print(matrix_wt)
		# print(matrix_ko)
		array_wt.append(matrix_wt[0])
		array_ko.append(matrix_ko[0])

array_wt=np.asarray(array_wt)
array_wt=np.nan_to_num(array_wt)
array_ko=np.asarray(array_ko)
array_ko=np.nan_to_num(array_ko)


import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rcParams



x=np.log2(array_wt+1)
y=np.log2(array_ko+1)
result=stats.pearsonr(x, y)
print(result)
out_file=str(result[0])+"_"+str(result[1])+"_"+out_file
#pearsons_coefficient = np.corrcoef(x, y)
#print(pearsons_coefficient[0][1])
#out_file=str(pearsons_coefficient[0][1])+"_"+out_file
# array([[1.        , 0.95315636],
#        [0.95315636, 1.        ]])
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)
idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]
fig,ax=plt.subplots(figsize=(12,9),dpi=100)
scatter=ax.scatter(x,y,marker='o',c=z,s=15,label='LST',cmap='Spectral_r')
cbar=plt.colorbar(scatter,shrink=1,orientation='vertical',extend='both',pad=0.015,aspect=30,label='frequency') #orientation='horizontal'
font3={'family':'SimHei','size':16,'color':'k'}
plt.xlabel(bw_file1,fontdict=font3)
plt.ylabel(bw_file2,fontdict=font3)
#plt.xlim((0,6))
#plt.ylim((0,6))
plt.savefig(out_file,dpi=800,bbox_inches='tight',pad_inches=0)
#plt.show()
