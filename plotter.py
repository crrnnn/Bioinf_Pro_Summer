#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 16:23:52 2018

@author: ceromi
"""
import numpy as np
import matplotlib.pyplot as plt
import glob 


#def plot_hist (a)

files = glob.glob('svm_pssm_results_u2352/*.txt')
sorted_files = sorted(files)
print (sorted_files)


#%%

Single_Seq_CV = []
Multpl_Seq_Subs_Mtrx_CV =[]
Single_Seq_Freq_Mtrx_CV =[]
Single_Seq_TS = []
Multpl_Seq_Subs_Mtrx_TS =[]
Single_Seq_Freq_Mtrx_TS =[]

#data_cross_val = np.zeros(shape=(7, 3))
#data_test_score = np.zeros(shape=(7, 3))


for i in sorted_files:
    with open(i, 'r') as temp:
        tempfilehandle = temp.read().splitlines()
        #for j in range (len(tempfilehandle)):
            #data_cross_val.append([j+1],[j+6],[j+11])
            #data_test_score.append([j+2],[j+7],[j+12])
        Single_Seq_CV.append(float(tempfilehandle[1]))
        Multpl_Seq_Subs_Mtrx_CV.append(float(tempfilehandle[6]))
        Single_Seq_Freq_Mtrx_CV.append(float(tempfilehandle[11]))
        Single_Seq_TS.append(float(tempfilehandle[2]))
        Multpl_Seq_Subs_Mtrx_TS.append(float(tempfilehandle[7]))
        Single_Seq_Freq_Mtrx_TS.append(float(tempfilehandle[12]))

print (Single_Seq_CV)

#%%




fig, ax = plt.subplots()
index = np.arange(7)
bar_width = 0.30




bar_1 = ax.bar(index, Single_Seq_CV, bar_width, alpha= 0.70, 
               color='b', label='Single_Seq')
 
bar_2 = ax.bar(index + bar_width, Multpl_Seq_Subs_Mtrx_CV, bar_width, 
               alpha=0.70, color='g', label='Multpl_Seq_Subs_Mtrx')

bar_3 = ax.bar(index + bar_width*2, Single_Seq_Freq_Mtrx_CV, bar_width, 
               alpha=0.70, color='r', label='Single_Seq_Freq_Mtrx')
ax.set_ylim(0.5,1)
#ax.legend(loc='lower right')
ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.30),
          fancybox=True, shadow=True, ncol=5)
ax.set_xlabel('Window Size')
ax.set_ylabel('Score')
ax.set_title('Cross-Validation Scores')
ax.set_xticks(index + bar_width)
ax.set_xticklabels(('13','15', '17', '19', '21', '23', '25'))
#plt.show() 
plt.savefig('cross_validation_scores.png', bbox_inches='tight', tight_layout=True)

#%%

fig, ax = plt.subplots()
index = np.arange(7)
bar_width = 0.30




bar_4 = plt.bar(index, Single_Seq_TS, bar_width, alpha= 0.70, color='b', label='Single_Seq')
 
bar_5 = plt.bar(index + bar_width, Multpl_Seq_Subs_Mtrx_TS, bar_width, alpha=0.70, color='g', label='Multpl_Seq_Subs_Mtrx')

bar_6 = plt.bar(index + bar_width*2, Single_Seq_Freq_Mtrx_TS, bar_width, alpha=0.70, color='r', label='Single_Seq_Freq_Mtrx')

ax.set_ylim(0.4,0.8)
#plt.legend(loc='lower right')
plt.xlabel('Window Size')
plt.ylabel('Score')
plt.title('Test Scores')
plt.xticks(index + bar_width, ('13','15', '17', '19', '21', '23', '25'))
ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.30),
          fancybox=True, shadow=True, ncol=5)
#plt.show() 

plt.savefig('test_scores.png', bbox_inches='tight', tight_layout=True)




#%% 




