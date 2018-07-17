import numpy as np
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

print (aa_dict)

for i,j in enumerate (states_list):
    states_dict[j]=i

print (states_dict)

encode = OneHotEncoder(n_values=len(aa_list))

for i in mydict:
    for j in range (len(mydict [i][0])):
        mydict[i][0][j] = aa_dict[mydict[i][0][j]]

        mydict[i][0][j] = encode.fit_transform(mydict [i][0][j]).toarray()


for i in mydict:
    for j in range (len(mydict [i][1])):
        mydict[i][1][j] = states_dict[mydict[i][1][j]]



print (mydict)








# sliding window

#https://stackoverflow.com/questions/8269916/what-is-sliding-window-algorithm-examples




# if __name__=="__main__":
#  print (parser('jpred1.3line.txt'))
