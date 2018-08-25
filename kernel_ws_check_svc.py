#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 14:04:59 2018

@author: ceromi
"""

import numpy as np
import more_itertools as mit
from sklearn.preprocessing import OneHotEncoder
import sys 
from sklearn.model_selection import cross_val_score
from sklearn.metrics import f1_score
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split


mydict={}



filereader = open('quarter_dataset.txt','r')
text = filereader.read().splitlines()


for i in range (0, len(text), 3):
    mydict[text[i].lstrip ('>')] = [list(text[i+1]), list(text[i+2])]

aa_list= list('ACDEFGHIKLMNPQRSTVWXY')
states_list= list ('HEC')

aa_dict = {}
states_dict= {}
for i,j in enumerate (aa_list):
    aa_dict[j]=i

#print (aa_dict)

for i,j in enumerate (states_list):
    states_dict[j]=i

#print (states_dict)

encode = OneHotEncoder(n_values=len(aa_list))

for i in mydict:
    for j in range (len(mydict [i][0])):
        mydict[i][0][j] = aa_dict[mydict[i][0][j]]

        mydict[i][0][j] = encode.fit_transform(mydict [i][0][j]).toarray().reshape(21)


for i in mydict:
    for j in range (len(mydict [i][1])):
        mydict[i][1][j] = states_dict[mydict[i][1][j]]


#%%
#print (mydict)

# sliding window
### window size of 19 or 21 gives the best results accourding to http://biomine.cs.vcu.edu/papers/CIBCB2006-1.pdf

#window_size = int(sys.argv[1])      
window_size = 19
pad = [np.zeros(shape=21)]
padding = pad * int((window_size-1)/2)

#%%

#padded_aa_dict= {}
for i in mydict:
    padded_temp = padding
    padded_temp = padded_temp + mydict[i][0]
    padded_temp = padded_temp + padding
    
    mydict[i][0] = padded_temp

#padded_aa_list = pad.extendlist(padded_aa_list).extendlist(pad)

#windows_dict = {}
for i in mydict:
    

    window_temp = list(mit.windowed(mydict[i][0], window_size))
    window_temp = [list(mit.flatten(i)) for i in window_temp]
    mydict[i][0] = window_temp
    
    
     
#%%


#############################################################################

kernella= sys.argv[2] 

X = []
Y = []

for i in mydict:
    X.extend(mydict[i][0])
    Y.extend(mydict[i][1])


train_X, test_X, train_Y, test_Y = train_test_split(X, Y, train_size = 0.9 , shuffle=True)
#%%

model = SVC(kernel= kernella, verbose=True)

model.fit(train_X, train_Y)

results = open ('kernel_ws_' + str(window_size) + '_'+ str(kernella) + '.txt', 'w') 

score_train_X_Y = model.score(train_X, train_Y)
print(score_train_X_Y)

score_test_X_Y = model.score(test_X, test_Y)
print(score_test_X_Y)

pred_Y = model.predict(test_X)
print(pred_Y)
print()

#############################
# Cross validation 5k
############################

scores= cross_val_score(model.fit(train_X, train_Y), X, Y, cv=5) 
print ('Cross validation scores: ')
print(scores)
print()

#
results.write('Single sequence' + '\n' + 'train score: ' + str(score_train_X_Y) + '\n' + 'test score: ' + str(score_test_X_Y) + 
               '\n' + '5k cross-validation score: ' + str(scores) + '\n' + 'prediction: ' + str(pred_Y) + '\n\n\n')

score = f1_score(test_Y, pred_Y, labels=None, pos_label=1, average='weighted', sample_weight=None)

results.write ('F1-score: ' + str(score) + '\n\n\n')



#%%
###################################################################################
# Add psi-blast information obtained from PSSM files #
###################################################################################

#sigmoid function to normalize PSSM frequency scores#

from math import e
def sigmoid(x):
    function = 1/(1+e**(-x))
    return function

# Read PSSM files and create separete lists for multiple sequence substituion matrix and single sequence frequency matrix #


def pssm (prot):
    
 
    filereader = open('dataset/PSSM/' + prot +'.fasta.pssm','r')
    text = filereader.read().splitlines()
    subs_list=[]
    freq_list=[]
    
    
    for i in range (3,len(text),1):
        if len(text[i]) == 0:
            break   
        
        int_subs_list = [int(j) for j in text[i].split()[2:42]]
    
        subs_list.append([sigmoid(k) for k in int_subs_list[0:20]])
        
        freq_list.append(int_subs_list[20:40])
        
        
    return [subs_list, freq_list]


# Create dictionary with PSSM information, takes protein Ids from mydict and values are multiple seq. substitution matrix and frequency matrix information

pssm_dict = {}

for i in mydict:
    pssm_dict[i] = pssm(i)

 
# Sliding windows for PSSM

#window_size_p = int(sys.argv[1])     
pad_pssm = [np.zeros(shape=20)]
padding_pssm = pad_pssm * int((window_size-1)/2)


#%%

# Add padding and create a dictionary with PSSM information

for i in pssm_dict:
    padded_pssm_temp = padding_pssm
    padded_pssm_temp = padded_pssm_temp + pssm_dict[i][0]
    padded_pssm_temp = padded_pssm_temp + padding_pssm
    
    pssm_dict[i][0] = padded_pssm_temp


for i in pssm_dict:
    padded_pssm_temp2 = padding_pssm
    padded_pssm_temp2 = padded_pssm_temp2 + pssm_dict[i][1]
    padded_pssm_temp2 = padded_pssm_temp2 + padding_pssm
    
    pssm_dict[i][1] = padded_pssm_temp2


for i in pssm_dict:
    

    window_pssm_temp = list(mit.windowed(pssm_dict[i][0], window_size))
    window_pssm_temp = [list(mit.flatten(i)) for i in window_pssm_temp]
    pssm_dict[i][0] = window_pssm_temp   
    
 
    window_pssm_temp2 = list(mit.windowed(pssm_dict[i][1], window_size))
    window_pssm_temp2 = [list(mit.flatten(i)) for i in window_pssm_temp2]
    pssm_dict[i][1] = window_pssm_temp2  
    
####################################################################################################   
# List for SVM, A for multiple sequence substutition matrix, B for frequency matrix, C for states    
#####################################################################################################        
A = []
B = []
C = []

for i in pssm_dict:
    
    A.extend(pssm_dict[i][0])
    B.extend(pssm_dict[i][1])
   
for i in mydict:
    C.extend(mydict[i][1])
    
    
#%% 

train_A, test_A, train_C, test_C = train_test_split(A, C, train_size = 0.9 , shuffle=True)

model = SVC(kernel= kernella, verbose=True)

model.fit(train_A, train_C)

score_train_A_C = model.score(train_A, train_C)
score_test_A_C = model.score(test_A, test_C)


pred_C = model.predict(test_A)


#############################
# Cross validation 5k
############################

scores= cross_val_score(model.fit(train_A, train_C), A, C, cv=5) 
print ('Cross validation scores: ')
print(scores)
print()


results.write('Multiple sequence substitution matrix' + '\n' + 'train score: ' + str(score_train_A_C) + '\n' + 
              'test score: ' + str(score_test_A_C) + 
               '\n' + '5k cross-validation score: ' + str(scores) + '\n' + 'prediction: ' + str(pred_C) + '\n\n\n')

score = f1_score(test_C, pred_C, labels=None, pos_label=1, average='weighted', sample_weight=None)

results.write ('F1-score: ' + str(score) + '\n\n\n')

  
#%%

train_B, test_B, train_D, test_D = train_test_split(B, C, train_size = 0.9 , shuffle=True)

model = SVC(kernel= kernella, verbose=True)

model.fit(train_B, train_D)

score_train_B_D = model.score(train_B, train_D)
score_test_B_D = model.score(test_B, test_D)

pred_D = model.predict(test_B)
    
#############################
# Cross validation 5k
############################

scores= cross_val_score(model.fit(train_B, train_D), B, C, cv=5) 
print ('Cross validation scores: ')
print(scores)
print()


results.write('Frequency matrix' + '\n' + 'train score: ' + str(score_train_B_D) + '\n' + 
              'test score: ' + str(score_test_B_D) + 
               '\n' + '5k cross-validation score: ' + str(scores) + '\n' + 'prediction: ' + str(pred_D) + '\n\n\n')

score = f1_score(test_D, pred_D, labels=None, pos_label=1, average='weighted', sample_weight=None)

results.write ('F1-score: ' + str(score) + '\n\n\n')




