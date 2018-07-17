import numpy as np
from sklearn import svm
from sklearn.preprocessing import OneHotEncoder
from sklearn.feature_extraction import DictVectorizer
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

filereader = open('jpred1.3line.txt','r')
text = filereader.read().splitlines()

for i in range (0, len(text), 3):
    mydict[text[i].lstrip ('>')] = [text[i+1], text[i+2]]
    print (list(text[i+1]))
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

set1 = []
set2 = []
set3 = []
n = 0
for i in mydict:
    if 0<= n <= 269:
        set1.append(i)
        #n=+1

    if 270<= n <=539:
        set2.append(i)
        #n=+1

    if 540<= n <=829:
        set3.append(i)
    n+=1

print (mydict)
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




# sliding window

#https://stackoverflow.com/questions/8269916/what-is-sliding-window-algorithm-examples




# if __name__=="__main__":
#  print (parser('jpred1.3line.txt'))
