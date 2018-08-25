# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

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
from sklearn.ensemble import RandomForestClassifier
from scipy.stats import randint as sp_randint
from sklearn.model_selection import RandomizedSearchCV

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
# Lists for classifier, X for amino acid, Y for states 
##########################################################


X = []
Y = []

for i in mydict:
    X.extend(mydict[i][0])
    Y.extend(mydict[i][1])

##############################################################################################
# Create train and test sets for random forest, calculate the scores,  for single sequence information
##############################################################################################
print()
print ('Single sequence based train & test scores and prediction for window size: ' + str(window_size)) 
print ()
print ('Random forest classifier')   
train_X, test_X, train_Y, test_Y = train_test_split(X, Y, train_size = 0.9 , shuffle=True)

############################################
# Random Forest #
############################################

clf = RandomForestClassifier()
clf.fit(train_X, train_Y)

score_train_X_Y = clf.score(train_X, train_Y)
print(score_train_X_Y)

score_test_X_Y = clf.score(test_X, test_Y)
print(score_test_X_Y)

pred_Y = clf.predict(test_X)
print(pred_Y)
print()
#############################
# Cross validation 5k
############################

scores= cross_val_score(clf.fit(train_X, train_Y), X, Y, cv=5) 
print ('Cross validation scores: ')
print(scores)
print() 


#################################
# f1-score F1 = 2 * (precision * recall) / (precision + recall) #
#################################

score = f1_score(test_Y, pred_Y, labels=None, pos_label=1, average='weighted', sample_weight=None)
print ('F1-score: ', score)
print ()



#%%
###################################################
# Randomized search for single sequence predictor
###################################################
print ('Randomized-Search with 5k-fold Cross Validation // Single sequence matrix // Window size: ' + str(window_size) )
print("# Tuning hyper-parameters for f1_weighted")

# specify parameters and distributions to sample from
param_dist = {"n_estimators": [50, 100, 150],
              "max_depth": [3, None],
              "max_features": sp_randint(1, 11),
              "min_samples_split": sp_randint(2, 11),
              "min_samples_leaf": sp_randint(1, 11),
              "bootstrap": [True, False],
              "criterion": ["gini", "entropy"]}

#%%
# Utility function to report best scores
#http://scikit-learn.org/stable/auto_examples/model_selection/plot_randomized_search.html
def report(results, n_top=3):
    for i in range(1, n_top + 1):
        candidates = np.flatnonzero(results['rank_test_score'] == i)
        for candidate in candidates:
            print("Model with rank: {0}".format(i))
            print("Mean validation score: {0:.3f} (std: {1:.3f})".format(
                  results['mean_test_score'][candidate],
                  results['std_test_score'][candidate]))
            print("Parameters: {0}".format(results['params'][candidate]))
            print("")

#%%
# run randomized search
n_iter_search = 10
random_search = RandomizedSearchCV(clf, param_distributions=param_dist, cv=5, 
                                   scoring="f1_weighted", verbose=True,  n_jobs=-1,
                                   n_iter=n_iter_search)

random_search.fit(train_X, train_Y)
print("RandomizedSearchCV for %d candidates parameter settings."  % (n_iter_search))
print(report(random_search.cv_results_))


print("Best parameters set found on development set:")
print()
print(random_search.best_params_)


print()

print("Grid scores on development set:")
print()
means = random_search.cv_results_['mean_test_score']
stds = random_search.cv_results_['std_test_score']
for mean, std, params in zip(means, stds, random_search.cv_results_['params']):
    print("%0.3f (+/-%0.03f) for %r"
          % (mean, std * 2, params))
print()

print("Detailed classification report:")
print()
print("The model is trained on the full development set.")
print("The scores are computed on the full evaluation set.")
print()
y_true, y_pred = test_Y, clf.predict(test_X)
print(classification_report(y_true, y_pred))
print('done \n') 

#%%
#####################################################
# Grid search for single sequence predictor
####################################################  

#print ('Grid-Search with 5k-fold Cross Validation // Single sequence matrix // Window size: ' + str(window_size) )
#print("# Tuning hyper-parameters for f1_weighted")
#      
## use a full grid over all parameters
#param_grid = {"n_estimators": [50, 100, 150],
#              "max_depth": [3, None],
#              "max_features": [1, 3, 10],
#              "min_samples_split": [2, 3, 10],
#              "min_samples_leaf": [1, 3, 10],
#              "bootstrap": [True, False],
#              "criterion": ["gini", "entropy"]}
#
## run grid search
#grid_search = GridSearchCV(clf, param_grid=param_grid, cv=5, scoring="f1_weighted", 
#                           verbose=True,  n_jobs=-1)
#
#grid_search.fit(train_X, train_Y)

#
#%%
#print("GridSearchCV took for %d candidate parameter settings."
#     % len(grid_search.cv_results_['params']))
#print(report(grid_search.cv_results_))      
#
#print("Best parameters set found on development set:")
#print()
#print(grid_search.best_params_)
#
#
#print()
#
#print("Grid scores on development set:")
#print()
#means = grid_search.cv_results_['mean_test_score']
#stds = grid_search.cv_results_['std_test_score']
#for mean, std, params in zip(means, stds, grid_search.cv_results_['params']):
#    print("%0.3f (+/-%0.03f) for %r"
#          % (mean, std * 2, params))
#print()
#
#print("Detailed classification report:")
#print()
#print("The model is trained on the full development set.")
#print("The scores are computed on the full evaluation set.")
#print()
#y_true, y_pred = test_Y, clf.predict(test_X)
#print(classification_report(y_true, y_pred))
#print()
##print('Best score: ')
##print (grid_search.best_score_)
#
#print('done \n') 
 
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
# Lists for classifier , A for multiple sequence substutition matrix, B for frequency matrix, C for states    
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
# Create train and test sets for random forest, calculate the scores,  multiple sequence substutition matrix information 
#################################################################################################################    
print()
print ('Multple sequence substitution matrix based train and test scores for window size: ' + str(window_size)) 
print ()
print ('Random forest classifier') 
train_A, test_A, train_C, test_C = train_test_split(A, C, train_size = 0.9 , shuffle=True)

clf = RandomForestClassifier()
clf.fit(train_A, train_C)

score_train_A_C = clf.score(train_A, train_C)
print(score_train_A_C)

score_test_A_C = clf.score(test_A, test_C)
print(score_test_A_C)

pred_C = clf.predict(test_A)
print(pred_C)
print()


#############################
# Cross validation 5k
############################


scores= cross_val_score(clf.fit(train_A, train_C), A, C, cv=5) 
print ('Cross validation scores: ')
print(scores)
print()


#################################
# f1-score F1 = 2 * (precision * recall) / (precision + recall) #
#################################

score = f1_score(test_C, pred_C, labels=None, pos_label=1, average='weighted', sample_weight=None)
print ('F1-score: ', score)
print ()


#%%
###################################################
# Randomized search for multiple sequence substitution matrix predictor
###################################################
print ('Randomized-Search with 5k-fold Cross Validation // Multiple sequence substitution matrix // Window size: ' + str(window_size) )
print ()
print("# Tuning hyper-parameters for f1_weighted")


# run randomized search
n_iter_search = 10
random_search = RandomizedSearchCV(clf, param_distributions=param_dist, cv=5, 
                                   scoring="f1_weighted", verbose=True,  n_jobs=-1,
                                   n_iter=n_iter_search)

random_search.fit(train_A, train_C)
print("RandomizedSearchCV for %d candidates parameter settings."  % (n_iter_search))
print(report(random_search.cv_results_))


print("Best parameters set found on development set:")
print()
print(random_search.best_params_)


print()

print("Grid scores on development set:")
print()
means = random_search.cv_results_['mean_test_score']
stds = random_search.cv_results_['std_test_score']
for mean, std, params in zip(means, stds, random_search.cv_results_['params']):
    print("%0.3f (+/-%0.03f) for %r"
          % (mean, std * 2, params))
print()

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
# Create train and test sets for random forest, calculate the scores,  frequency matrix information  
###########################################################################################
print()
print ('Frequency matrix based train and test scores for window size: ' + str(window_size))    
print()
print ('Random forest classifier') 
train_B, test_B, train_D, test_D = train_test_split(B, C, train_size = 0.9 , shuffle=True)

clf = RandomForestClassifier()
clf.fit(train_B, train_D)

score_train_B_D = clf.score(train_B, train_D)
print(score_train_B_D)

score_test_B_D = clf.score(test_B, test_D)
print(score_test_B_D)

pred_D = clf.predict(test_B)
print(pred_D)
print()

#############################
# Cross validation 5k
############################

scores= cross_val_score(clf.fit(train_B, train_D), B, C, cv=5) 
print ('Cross validation scores: ')
print(scores)
print ()


#################################
# f1-score F1 = 2 * (precision * recall) / (precision + recall) #
#################################

score = f1_score(test_D, pred_D, labels=None, pos_label=1, average='weighted', sample_weight=None)

print ('F1-score: ', score)
print ()


#%%
###################################################
# Randomized search for frequency matrix predictor
###################################################
print ('Grid-Search with 5k-fold Cross Validation // Single sequence frequency matrix // Window size: ' + str(window_size) )
print ()
print("# Tuning hyper-parameters for f1_weighted")


# run randomized search
n_iter_search = 10
random_search = RandomizedSearchCV(clf, param_distributions=param_dist, cv=5, 
                                   scoring="f1_weighted", verbose=True,  n_jobs=-1,
                                   n_iter=n_iter_search)

random_search.fit(train_B, train_D)
print("RandomizedSearchCV for %d candidates parameter settings."  % (n_iter_search))
print(report(random_search.cv_results_))


print("Best parameters set found on development set:")
print()
print(random_search.best_params_)


print()

print("Grid scores on development set:")
print()
means = random_search.cv_results_['mean_test_score']
stds = random_search.cv_results_['std_test_score']
for mean, std, params in zip(means, stds, random_search.cv_results_['params']):
    print("%0.3f (+/-%0.03f) for %r"
          % (mean, std * 2, params))
print()

print("Detailed classification report:")
print()
print("The model is trained on the full development set.")
print("The scores are computed on the full evaluation set.")
print()
y_true, y_pred = test_D, clf.predict(test_B)
print(classification_report(y_true, y_pred))
print('done \n') 
    

    
    
    
    
    
    
    


