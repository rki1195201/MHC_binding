

import math
import random

def sigmoidal(x):
        num = 1/float(1 + math.exp(-x))
        return num
def g(x):
        num = sigmoidal(x)
        return num*(1-num)

Amino_bits={'A':'10000000000000000000','R':'01000000000000000000','N':'00100000000000000000','D':'00010000000000000000','C':'00001000000000000000','Q':'00000100000000000000','E':'00000010000000000000',
 'G':'00000001000000000000','H':'00000000100000000000','I':'00000000010000000000','L':'00000000001000000000','K':'00000000000100000000','M':'00000000000010000000','F':'00000000000001000000',
 'P':'00000000000000100000','S':'00000000000000010000','T':'00000000000000001000','W':'00000000000000000100','Y':'00000000000000000010','V':'00000000000000000001'}


fp3 = open("abc.txt","w")
e = 2.0

kk=0
storeb = [''] * 400
fp1 = open("data_b.dat","r")
for line in fp1:
        line = line.split('  ')
        storeb[kk] = line
        kk=kk+1
fp1.close()

kk=0
storenb = [''] * 400
fp2 = open("data_nb.dat","r")
for line in fp2:
        line = line.split(' ')
        storenb[kk] = line
        kk=kk+1
fp2.close()


hidden = 5
v = [0] * (hidden+1)*10
w = [0] * (hidden+1)


E = 0
bias_v = [0] * (hidden+1)
bias_w = 0
for j in range(1,hidden+1):
    bias_v[j] = random.uniform(0,1)
    for k in range(1,10):
        jk = int(str(j)+str(k))
        v[jk] = bias_v[j]

bias_w = random.uniform(0,1)
for i in range(1,hidden+1):
    w[i]  = bias_w



for ii in range(1,51):

    last = 0
    
    bp = 0
    nbp = 0
    predict = 0
    count_err = 0
    average_err = 0
    deltaV = [0] * (hidden+1)*10
    deltaW = [0] * (hidden+1)
    
    for all in range(400):
        
        random_num = random.uniform(0,1)
        if(random_num>0.5 and bp<200):
                line = storeb[bp]
                bp = bp + 1
        else:
                line = storenb[nbp]
                nbp = nbp + 1
        #print(line)
        seq = line[0]           #get peptide like: DLMGYIPAV
        label = line[1].strip('\n')      #get class label like: 1
        
        if(label == '1'):
                d = 1
        else:
                d = 0
        #print(seq,label,random_num)
        I = [0]*len(seq)
        
        h = [0] * (hidden+1)
        H = [0] * (hidden+1)
        for j in range(1,hidden+1):
            for length in range(len(seq)):
                I[length] = Amino_bits[seq[length]]  #get each amino acid of peptide
                jk = int(str(j)+str(length+1))
                h[j] = h[j] + (v[jk] * (1/(int(I[length],2))))
                
            h[j] = h[j] + bias_v[j]
            
            H[j] = sigmoidal(h[j]) #sigmoidal function
            #print(H[j])
        
        o = 0
        O = 0      
        for j in range(1,hidden+1):
            o = o + (w[j] * H[j])
        o = o + bias_w
        O = sigmoidal(o) #sigmoidal function
        
        E = (O - d)**2 / 2
        #print(E,all)
        if(O>=0.5):
            predict = 1
        else:
            predict = 0
                
        if(predict != d):
            count_err = count_err + E
        #print(E,O,predict,d,all)
        
        
        for j in range(1,hidden+1):
            for k in range(1,10):
                jk = int(str(j)+str(k))
                deltaV[jk] = -e * (O - d) *  g(o) * w[j] * g(h[j]) * (1/(int(I[length],2)))
                v[jk] = v[jk] + deltaV[jk]

        for i in range(1,hidden+1):
            deltaW[i] = -e * (O - d) * g(o) * H[i]
            w[i] = w[i] + deltaW[i]

    average_err = count_err /400
    print("average error: "+str(average_err))
    
    fp3.write(str(average_err)+"\n")
print("-- Training process --\nTrain data (b): 200\nTrain data (nb): 200\nHidden node: "+str(hidden))
print("Epsilon: 0.5\nTraining cycle: 50\nInitial bias w: %f"%(bias_w))
print("Inital bias v: %f %f %f %f %f"%(bias_v[1],bias_v[2],bias_v[3],bias_v[4],bias_v[5]))
print("Training error after the last cycle: " + str(average_err) )




bp = 20
nbp = 20
predict = 0
count = 0
average_err = 0
deltaV = [0] * (hidden+1)*10
deltaW = [0] * (hidden+1)
    
for all in range(400):
        
    random_num = random.uniform(0,1)
    if(random_num>0.5 and bp<400):
            line = storeb[bp]
            bp = bp + 1
    else:
            line = storenb[nbp]
            nbp = nbp + 1
        #print(line)
    seq = line[0]           #get peptide like: DLMGYIPAV
    label = line[1].strip('\n')      #get class label like: 1
        
    if(label == '1'):
        d = 1
    else:
        d = 0
        #print(seq,label,random_num)
    I = [0]*len(seq)
        
    h = [0] * (hidden+1)
    H = [0] * (hidden+1)
    for j in range(1,hidden+1):
        for length in range(len(seq)):
            I[length] = Amino_bits[seq[length]]  #get each amino acid of peptide
            jk = int(str(j)+str(length+1))
            h[j] = h[j] + (v[jk] * (1/(int(I[length],2))))
                
        h[j] = h[j] + bias_v[j]
            
        H[j] = sigmoidal(h[j]) #sigmoidal function
            #print(H[j])
        
    o = 0
    O = 0      
    for j in range(1,hidden+1):
        o = o + (w[j] * H[j])
    o = o + bias_w
    O = sigmoidal(o) #sigmoidal function
        
    E = (O - d)**2 / 2
        #print(E,all)
    if(O>=0.5):
        predict = 1
    else:
        predict = 0
                
    if(predict == d):
        count = count + 1
        #print(E,O,predict,d,all)
        
        
    for j in range(1,hidden+1):
        for k in range(1,10):
            jk = int(str(j)+str(k))
            deltaV[jk] = -e * (O - d) *  g(o) * w[j] * g(h[j]) * (1/(int(I[length],2)))
            v[jk] = v[jk] + deltaV[jk]

    for i in range(1,hidden+1):
        deltaW[i] = -e * (O - d) * g(o) * H[i]
        w[i] = w[i] + deltaW[i]

print("-- Validation process –-\nValidation data (b): 200\nValidation data (nb): 200")
print("Validation accuracy: "+str(count/400*100)+"%")