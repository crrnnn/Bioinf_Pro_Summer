#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 00:10:32 2018

@author: ceromi
"""

import numpy as np
import more_itertools as mit
from sklearn.preprocessing import OneHotEncoder
import sys 
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report
from sklearn.model_selection import cross_val_score
from sklearn.metrics import f1_score

#########################################################################################
## Script that takes a sequence as an input and creates an input vector to sklearn ##
#########################################################################################
# Read the input file and split it into 3 lists (ids, sequences, and states) and put everyhting in a dictionary #
## keys = ids, values = sequence and structures ##

mydict={}
filereader = open('quarter_dataset.txt','r')
text = filereader.read().splitlines()


for i in range (0, len(text), 3):
    mydict[text[i].lstrip ('>')] = [list(text[i+1]), list(text[i+2])]
 
   
### One hot encoding ###
# Assign numbers to each amino acid and states #
aa_list= list('ACDEFGHIKLMNPQRSTVWXY')
states_list= list ('HEC')

aa_dict = {}
states_dict= {}
for i,j in enumerate (aa_list):
    aa_dict[j]=i



for i,j in enumerate (states_list):
    states_dict[j]=i


# Use OneHotEncoder to assign vectors to amino acids #
encode = OneHotEncoder(n_values=len(aa_list))

for i in mydict:
    for j in range (len(mydict [i][0])):
        mydict[i][0][j] = aa_dict[mydict[i][0][j]]

        mydict[i][0][j] = encode.fit_transform(mydict [i][0][j]).toarray().reshape(21)


for i in mydict:
    for j in range (len(mydict [i][1])):
        mydict[i][1][j] = states_dict[mydict[i][1][j]]


#%%


# Sliding windows for single sequence
### window size of 19 or 21 gives the best results accourding to http://biomine.cs.vcu.edu/papers/CIBCB2006-1.pdf ###
# Take the window size as a user input
        
window_size = int(sys.argv[1])      
pad = [np.zeros(shape=21)]
padding = pad * int((window_size-1)/2)


#%%

# Add padding to the beginning and ending
for i in mydict:
    padded_temp = padding
    padded_temp = padded_temp + mydict[i][0]
    padded_temp = padded_temp + padding
    
    mydict[i][0] = padded_temp

# Create windows
for i in mydict:
    

    window_temp = list(mit.windowed(mydict[i][0], window_size))
    window_temp = [list(mit.flatten(i)) for i in window_temp]
    mydict[i][0] = window_temp
    

##########################################################
# Create lists for SVM, X for amino acid, Y for states 
##########################################################


X = []
Y = []

for i in mydict:
    X.extend(mydict[i][0])
    Y.extend(mydict[i][1])

##############################################################################################
# Create train and test sets for svm, calculate the scores,  for single sequence information
##############################################################################################
print()
print ('Single sequence based train & test scores and prediction for window size: ' + str(window_size))    
train_X, test_X, train_Y, test_Y = train_test_split(X, Y, train_size = 0.9 , shuffle=True)

model = SVC(verbose=True)
model.fit(train_X, train_Y)

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
#%%
###################################################
# Grid search for single sequence predictor
###################################################

C_range = np.logspace(-2, 5, 8)
gamma_range = np.logspace(-5, 2, 8)

param_grid = [{'kernel': ['rbf'], 'gamma': gamma_range, 'C': C_range},
                    {'kernel': ['linear'], 'C': C_range},
                    {'kernel': ['poly'], 'gamma': gamma_range},
                    {'kernel': ['sigmoid'], 'gamma': gamma_range}
                    ]


#################################
# f1-score F1 = 2 * (precision * recall) / (precision + recall) #
#################################

score = f1_score(test_Y, pred_Y, labels=None, pos_label=1, average='weighted', sample_weight=None)
print ('F1-score: ', score)
print ()


print ('Grid-Search with 5k-fold Cross Validation // Single sequence matrix // Window size: ' + str(window_size) )
print("# Tuning hyper-parameters for f1_weighted")

#scoring = score + '_macro'
clf = GridSearchCV(SVC(), param_grid, cv=5,
                   scoring="f1_weighted", verbose=True,  n_jobs=-1)

clf.fit(train_X, train_Y)

print("Best parameters set found on development set:")
print()
print(clf.best_params_)
   
print("Grid scores on development set:")
print()
means = clf.cv_results_['mean_test_score']
stds = clf.cv_results_['std_test_score']
for mean, std, params in zip(means, stds, clf.cv_results_['params']):
    print("%0.3f (+/-%0.03f) for %r"
          % (mean, std * 2, params))
 

print("Detailed classification report:")
print()
print("The model is trained on the full development set.")
print("The scores are computed on the full evaluation set.")
print()
y_true, y_pred = test_Y, clf.predict(test_X)
print(classification_report(y_true, y_pred))
print('done \n') 
     
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
#################################################################################################################   
# Create train and test sets for svm, calculate the scores,  multiple sequence substutition matrix information 
#################################################################################################################    
print ('Multple sequence substitution matrix based train and test scores for window size: ' + str(window_size)) 
train_A, test_A, train_C, test_C = train_test_split(A, C, train_size = 0.9 , shuffle=True)

model = SVC(verbose=True)

model.fit(train_A, train_C)

score_train_A_C = model.score(train_A, train_C)
print(score_train_A_C)

score_test_A_C = model.score(test_A, test_C)
print(score_test_A_C)


pred_C = model.predict(test_A)
print(pred_C)
print()
#############################
# Cross validation 5k
############################


scores= cross_val_score(model.fit(train_A, train_C), A, C, cv=5) 
print ('Cross validation scores: ')
print(scores)
print()


######################################################################
#Grid search for multiple sequence substitution matrix predictor
#######################################################################
C_range = np.logspace(-2, 5, 8)
gamma_range = np.logspace(-5, 2, 8)

param_grid = [{'kernel': ['rbf'], 'gamma': gamma_range, 'C': C_range},
                    {'kernel': ['linear'], 'C': C_range},
                    {'kernel': ['poly'], 'gamma': gamma_range},
                    {'kernel': ['sigmoid'], 'gamma': gamma_range}
                    ]
#################################
# f1-score F1 = 2 * (precision * recall) / (precision + recall) #
#################################

score = f1_score(test_C, pred_C, labels=None, pos_label=1, average='weighted', sample_weight=None)

print ('F1-score: ', score)
print ()


print ('Grid-Search with 5k-fold Cross Validation // Multiple sequence substitution matrix // Window size: ' + str(window_size) )
print ()
print("# Tuning hyper-parameters for f1_weighted")


clf = GridSearchCV(SVC(), param_grid, cv=5,
                   scoring="f1_weighted", verbose=True,  n_jobs=-1)

clf.fit(train_A, train_C)

print("Best parameters set found on development set:")
print()
print(clf.best_params_)
   
print("Grid scores on development set:")
print()
means = clf.cv_results_['mean_test_score']
stds = clf.cv_results_['std_test_score']
for mean, std, params in zip(means, stds, clf.cv_results_['params']):
    print("%0.3f (+/-%0.03f) for %r"
          % (mean, std * 2, params))
 

print("Detailed classification report:")
print()
print("The model is trained on the full development set.")
print("The scores are computed on the full evaluation set.")
print()
y_true, y_pred = test_C, clf.predict(test_A)
print(classification_report(y_true, y_pred))
print('done \n')

    

#%% 
###########################################################################################    
# Create train and test sets for svm, calculate the scores,  frequency matrix information  
###########################################################################################
print ('Frequency matrix based train and test scores for window size: ' + str(window_size))    
train_B, test_B, train_D, test_D = train_test_split(B, C, train_size = 0.9 , shuffle=True)


model = SVC(verbose=True)

model.fit(train_B, train_D)

score_train_B_D = model.score(train_B, train_D)

print(score_train_B_D)

score_test_B_D = model.score(test_B, test_D)
print(score_test_B_D)

pred_D = model.predict(test_B)
print(pred_D)
print()

#############################
# Cross validation 5k
############################

scores= cross_val_score(model.fit(train_B, train_D), B, C, cv=5) 
print ('Cross validation scores: ')
print(scores)
print ()
#############################################
# Grid Search for frequency matrix predictor   
#############################################
    
C_range = np.logspace(-2, 5, 8)
gamma_range = np.logspace(-5, 2, 8)
 
param_grid = [{'kernel': ['rbf'], 'gamma': gamma_range, 'C': C_range},
                    {'kernel': ['linear'], 'C': C_range},
                    {'kernel': ['poly'], 'gamma': gamma_range},
                    {'kernel': ['sigmoid'], 'gamma': gamma_range}
                    ]

#################################
# f1-score F1 = 2 * (precision * recall) / (precision + recall) #
#################################

score = f1_score(test_D, pred_D, labels=None, pos_label=1, average='weighted', sample_weight=None)

print ('F1-score: ', score)
print ()


print ('Grid-Search with 5k-fold Cross Validation // Single sequence frequency matrix // Window size: ' + str(window_size) )
print ()
print("# Tuning hyper-parameters for f1_weighted")


clf = GridSearchCV(SVC(), param_grid, cv=5,
                   scoring="f1_weighted", verbose=True,  n_jobs=-1)

clf.fit(train_B, train_D)

print("Best parameters set found on development set:")
print()
print(clf.best_params_)
   
print("Grid scores on development set:")
print()
means = clf.cv_results_['mean_test_score']
stds = clf.cv_results_['std_test_score']
for mean, std, params in zip(means, stds, clf.cv_results_['params']):
    print("%0.3f (+/-%0.03f) for %r"
          % (mean, std * 2, params))
 

print("Detailed classification report:")
print()
print("The model is trained on the full development set.")
print("The scores are computed on the full evaluation set.")
print()
y_true, y_pred = test_D, clf.predict(test_B)
print(classification_report(y_true, y_pred))
print('done')

    
    
    
    
    
    
    


