#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 22 22:44:32 2018

@author: ceromi
"""

filereader = open('dataset/jpred1.3line.txt','r')
text = filereader.read().splitlines()


for i in range(0, len(text), 3):
    with open ('dataset/'+text[i].lstrip('>')+'.fasta', 'w') as lines:    
        lines.write(text[i] + '\n' + text[i+1])
        
        
        
        
    
    