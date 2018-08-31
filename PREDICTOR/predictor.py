#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 26 17:13:31 2018

@author: ceromi
"""


import numpy as np
import more_itertools as mit
import sys 
from sklearn.externals import joblib
from Bio import SeqIO
import os
#########################################################################################
## Script that takes a sequence as an input and creates an input vector to sklearn ##
#########################################################################################

print ("---Welcome to The PREDICTOR...---")
print ("Your prediction will be ready in seconds...")

mydict={}

if os.path.exists("Psi_blast") == False:
    os.mkdir("Psi_blast/")

for record in SeqIO.parse(sys.argv[1], "fasta"):
    mydict[str(record.id)] = str(record.seq)
    
    with open ("Psi_blast/" + record.id +'.fasta', 'w') as lines:    
        lines.write('>' + record.id +'\n' + str(record.seq))


############################################################################
## Run Psi-BLAST        
############################################################################
for i in mydict:
    Temp ="Psi_blast/"+i+".fasta"
    outfa = i + ".fasta"
    cmd = """#!/bin/bash
    if [ ! -f "{fa}.pssm" ] ; then
    echo "Running Psi-blast on {fa}..."
    psiblast -query "{fa}" -evalue 0.01 -db /local_uniref/uniref/uniref90/uniref90.db -num_iterations 3  -out Psi_blast/"{outfa}.psiblast" -out_ascii_pssm Psi_blast/"{outfa}.pssm"
    fi
    """.format(fa=Temp, outfa=outfa)
    
    os.system(cmd)


        

###################################################################################
# Add psi-blast information obtained from PSSM files #
###################################################################################

#sigmoid function to normalize PSSM frequency scores#

from math import e
def sigmoid(x):
    function = 1/(1+e**(-x))
    return function

###########################################################################################################################
# Read PSSM files and create separete lists for multiple sequence substituion matrix and single sequence frequency matrix #
###########################################################################################################################

def pssm (prot):
    
 
    filereader = open('Psi_blast/' + prot +'.fasta.pssm','r')
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


########################################################################################################################
# Create dictionary with PSSM information, takes protein Ids from mydict and values are multiple seq. substitution matrix and frequency matrix information
########################################################################################################################
pssm_dict = {}

for i in mydict:
    pssm_dict[i] = pssm(i)

############################## 
# Sliding windows for PSSM
##############################    
window_size = 7
pad_pssm = [np.zeros(shape=20)]
padding_pssm = pad_pssm * int((window_size-1)/2)



#############################################################
# Add padding and create a dictionary with PSSM information
#############################################################
for i in pssm_dict:
    padded_pssm_temp = padding_pssm
    padded_pssm_temp = padded_pssm_temp + pssm_dict[i][0]
    padded_pssm_temp = padded_pssm_temp + padding_pssm
    
    pssm_dict[i][0] = padded_pssm_temp



for i in pssm_dict:
    

    window_pssm_temp = list(mit.windowed(pssm_dict[i][0], window_size))
    window_pssm_temp = [list(mit.flatten(i)) for i in window_pssm_temp]
    pssm_dict[i][0] = window_pssm_temp   
    
   



#################################################################################################################   
# Import the model for substitution matrix and predict
#################################################################################################################    
clf = joblib.load('RFC_final_optimized_7.pkl') 
state_dict={0:'H', 1:'E', 2:'C'}


with open ('Prediction_for_' + str(sys.argv[1]) + '.txt', 'w') as output:
    for j in pssm_dict:
        prediction = clf.predict(pssm_dict[j][0])
        pred_list= ''.join([state_dict[int(i)] for i in prediction])
        
        output.write('>'+ j + '\n' + mydict[j] +'\n' + pred_list +'\n')
         
    
    
    
print ("Prediction complete, please check the predition file for results!")    
    
    
    
    
    
    