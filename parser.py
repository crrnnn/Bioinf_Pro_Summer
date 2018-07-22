import numpy as np
import more_itertools as mit
import itertools as it
from sklearn import svm
from sklearn.preprocessing import OneHotEncoder
#from sklearn.feature_extraction import DictVectorizer
## Function that takes a sequence as an input and creates an input vector to sklearn ##

# First, write a parser to read the input file and split it into 3 lists (ids, sequences, and states) and put everyhting in a dictionary #
## keys = ids, values = sequence and structures ##

#def parser(filename):

# ids=[]
# #seq=[]
# #states=[]
mydict={}
# mydict1 ={}
# mydict2= {}


filereader = open('small_test_data.txt','r')
text = filereader.read().splitlines()


for i in range (0, len(text), 3):
    mydict[text[i].lstrip ('>')] = [list(text[i+1]), list(text[i+2])]
    #print (list(text[i+1]))
    # mydict1 [text[i].lstrip ('>')] = list(text[i+1])
    # mydict2 [text[i].lstrip ('>')] = text[i+2]
    # seq.append(text[i+1])
    # states.append(text[i+2])

#mydict = dict(zip(ids, [seq,states]))  ### why does this take only the 1st prot name but no the rest?

    # mydict1 = dict(zip(ids, seq))
    # mydict2 = dict(zip(ids, states))
    #
    # for key in mydict1.keys() & mydict2.keys():
    #     fulldict [key] = (mydict1[key], mydict2[key])
    #text.close()

# print (mydict)

# set1 = []
# set2 = []
# set3 = []
# n = 0
# for i in mydict:
#     if 0<= n <= 269:
#         set1.append(i)
#         #n=+1
#
#     if 270<= n <=539:
#         set2.append(i)
#         #n=+1
#
#     if 540<= n <=829:
#         set3.append(i)
#     n+=1

#print (mydict)
# print ('set1:')
# print (set1)
# print ('set2:')
# print (set2)
# print ('set3:')
# print (set3)

### One hot encoding ###


## DictVectorizer didnt work, probably it takes only keys with a string, not values with a string ##

# v = DictVectorizer(sparse=False)
# D= [mydict1, mydict2]
# X = v.fit_transform(mydict1)
#
# print (X.shape)

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
#%%
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
       
window_size = 19
pad = [np.zeros(shape=21)]
padding = pad * int((window_size-1)/2)

# for i in mydict:
#     for j in range (len(mydict [i][0])):
#         sliding_window =list(mit.stagger(mydict[i][0][j], offsets=(-8, 0, 8), longest=True,  fillvalue=0 ))
#
#
# print (sliding_window)
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
    

    window_temp = list(mit.windowed(mydict[i][0], 19))
    window_temp = [list(mit.flatten(i)) for i in window_temp]
    mydict[i][0] = window_temp
    
    
     
#%%
print(len(mydict[i][0]))

#############################################################################

from sklearn.svm import SVC
from sklearn.model_selection import train_test_split

X = []
Y = []
for i in mydict:
    X.extend(mydict[i][0])
    Y.extend(mydict[i][1])


train_X, test_X, train_Y, test_Y = train_test_split(X, Y, train_size = 0.9 , shuffle=True)


model = SVC(verbose=True)

model.fit(train_X, train_Y)

model.score

model.score()

model.score(train_X, train_Y)

model.score(test_X, test_Y)





# if __name__=="__main__":
#  print (parser('jpred1.3line.txt'))



# if __name__=="__main__":
#  print (parser('jpred1.3line.txt'))
