#!/usr/bin/python3
import numpy as np
import pickle
from sklearn.ensemble import RandomForestClassifier
from sklearn import svm
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import KFold
from sklearn import metrics
from matplotlib import pyplot as plt
from sklearn.metrics import confusion_matrix
from sklearn.metrics import matthews_corrcoef
from sklearn.preprocessing import StandardScaler

n = 8
n2 = 3
n3 = 6
n4 = 5
n5 = 6
data = open("RSC758.txt").readlines()
fasta = open("RSC758.fa").readlines()
residues = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","X","Y"]
residuenames = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]
F=dict()
for i in range(1,len(fasta),2):
    F[fasta[i-1].replace(">","").strip()] = fasta[i].strip()
data = [x.split() for x in data][1:]
SASA = pickle.load(open("FreeSASA/SASA.dat","rb"))
MD = pickle.load(open("ModellerRSC758/DIST.dat","rb"))
MDCARBON = pickle.load(open("ModellerRSC758/DISTCARBON.dat","rb"))
MDALPHA = pickle.load(open("ModellerRSC758/DISTRESIDUE.dat","rb"))
PROPKA = pickle.load(open("Propka/PropkaRSC758.dat","rb"))
PSSM = pickle.load(open("PSSM/PSSMRSC758.dat","rb"))
SS = pickle.load(open("SecondaryStructure/SSRSC758.dat","rb"))
y=np.array([float(x[2]) for x in data])
Xfa = [x[0] for x in data]
Xidx = [int(x[1]) for x in data]
X=np.zeros((len(data),n + n2 + n3 + 1 + n4 * 21 + n5 * 20 + 1 + 260 + 39))

def nearest(i):
    global Xfa
    global Xidx
    global n
    indices = np.array([i+1 for i,x in enumerate(list(F[Xfa[i]])) if x == "C"])
    distances = sorted(abs(indices-Xidx[i]))
    D = distances[1:n+1]
    if len(D) == 0:
        D = [25]*n
    avg = int(sum(D)/len(D))
    while len(D) < n:
        D.append(avg)
    return D

def modeller_cysteine_nearest(i):
    global Xfa
    global Xidx
    global n2
    distances = MD[(Xfa[i],Xidx[i])]
    distances = sorted(distances)
    D = distances[1:n2+1]
    if len(D) == 0:
        D = [25]*n2
    avg = int(sum(D)/len(D))
    while len(D) < n2:
        D.append(avg)
    return D

def modeller_carbon_nearest(i):
    global Xfa
    global Xidx
    global n3
    distances = MDCARBON[(Xfa[i],Xidx[i])]
    D = distances[0:n3]
    if len(D) == 0:
        D = [3]*n3
    avg = int(sum(D)/len(D))
    while len(D) < n3:
        D.append(avg)
    return D

def nearest_residue(i):
    global Xfa
    global Xidx
    global n4
    global residues
    global labels
    RD = list()
    for residue in residues:
        indices = np.array([i+1 for i,x in enumerate(list(F[Xfa[i]])) if x == residue])
        distances = sorted(abs(indices-Xidx[i]))
        if residue != "C":
            D = distances[0:n4]
        if residue == "C":
            D = distances[1:n4+1]
        if len(D) == 0:
            D = [25]*n4
        avg = int(sum(D)/len(D))
        while len(D) < n4:
            D.append(avg)
        RD = RD + D 
    return RD

def modeller_alpha_nearest(i):
    global Xfa
    global Xidx
    global n5
    distances = MDALPHA[(Xfa[i],Xidx[i])]
    RDM = list()
    for residue in residuenames:
        D = sorted(distances[residue])
        D = D[0:n5]
        if len(D) == 0:
            D = [25]*n5
        avg = int(sum(D)/len(D))
        while len(D) < n5:
            D.append(avg)
        RDM = RDM + D
    return RDM

labels=[]
for i in range(n):
    labels.append("Cysteine"+str(i))
for i in range(n2):
    labels.append("CysteineDist"+str(i))
for i in range(n3):
    labels.append("CarbonDist"+str(i))
labels.append("SASA")
for residue in residues:
    for i in range(n4):
        labels.append(residue+str(i))
for residue in residuenames:
    for i in range(n5):
        labels.append(residue+str(i))
labels.append("PROPKA")
for i in range(260):
    labels.append("PSSM"+str(i))
for i in range(39):
    labels.append("SS"+str(i))

labels = np.array(labels)

for i in range(len(Xfa)):
    #X[i,0:n] = nearest(i)
    #X[i,n:n+n2] = modeller_cysteine_nearest(i)
    #X[i,n+n2:n+n2+n3] = modeller_carbon_nearest(i)
    X[i,n+n2+n3] = SASA[(Xfa[i],Xidx[i])]
    X[i,n+n2+n3+1:n+n2+n3+1+n4*21] = nearest_residue(i)
    X[i,n+n2+n3+1+n4*21:n+n2+n3+1+n4*21+n5*20] = modeller_alpha_nearest(i)
    X[i,n+n2+n3+1+n4*21+n5*20] = PROPKA[Xfa[i]][Xidx[i]]
    X[i,n+n2+n3+1+n4*21+n5*20+1:n+n2+n3+1+n4*21+n5*20+1+260] = PSSM[Xfa[i]][Xidx[i]]
    X[i,n+n2+n3+1+n4*21+n5*20+1+260:n+n2+n3+1+n4*21+n5*20+1+260+39] = SS[Xfa[i]][Xidx[i]]

#Cysteine makes up 0.715 of AUC
#X2 = np.zeros((len(Xfa),18))
#for i in range(len(Xfa)):
#    k = 0
#    for j, label in enumerate(labels):
#        if "GLN" in label:
#            X2[i,k] = X[i,j]
#            k += 1
#X=X2
#X2 = np.append(X2,y[:,None],1)
#np.savetxt("Glycine.csv",X2,delimiter=",")
#print(X2)
#exit()

#Z-scale scaling for svm
#ss = StandardScaler()
#ss.fit(X[:,n+n2+n3+1:n+n2+n3+1+n4*21])
#X[:,n+n2+n3+1:n+n2+n3+1+n4*21] = ss.transform(X[:,n+n2+n3+1:n+n2+n3+1+n4*21])
#ss = StandardScaler()
#ss.fit(X[:,n+n2+n3+1+n4*21:n+n2+n3+1+n4*21+n5*20])
#X[:,n+n2+n3+1+n4*21:n+n2+n3+1+n4*21+n5*20] = ss.transform(X[:,n+n2+n3+1+n4*21:n+n2+n3+1+n4*21+n5*20])
#ss = StandardScaler()
#ss.fit(X[:,n+n2+n3+1+n4*21+n5*20+1:n+n2+n3+1+n4*21+n5*20+1+260])
#X[:,n+n2+n3+1+n4*21+n5*20+1:n+n2+n3+1+n4*21+n5*20+1+260] = ss.transform(X[:,n+n2+n3+1+n4*21+n5*20+1:n+n2+n3+1+n4*21+n5*20+1+260])

clf = RandomForestClassifier(n_estimators=500,max_depth=None,random_state=3,min_samples_leaf=1,n_jobs=-1)
#clf = svm.SVC(C=1,gamma=.005,probability=True,kernel="rbf",random_state=0,degree=3)
#clf = KNeighborsClassifier(n_neighbors=20,weights="distance",n_jobs=-1)
yhat = np.zeros((len(data)))
ypred = np.empty((len(data)))
clf.fit(X,y)
print(labels[clf.feature_importances_.argsort()[::-1][:20]])
print(sorted(clf.feature_importances_,reverse=True)[:20] )

data = open("BALOSCTdb.txt").readlines()
fasta = open("BALOSCTdb.fa").readlines()
F=dict()
for i in range(1,len(fasta),2):
    F[fasta[i-1].replace(">","").strip()] = fasta[i].strip()
data = [x.split() for x in data][1:]
SASA = pickle.load(open("FreeSASA/SASAOSCTdb.dat","rb"))
MD = pickle.load(open("ModellerOSCTdb/DIST.dat","rb"))
MDCARBON = pickle.load(open("ModellerOSCTdb/DISTCARBON.dat","rb"))
MDALPHA = pickle.load(open("ModellerOSCTdb/DISTRESIDUE.dat","rb"))
PROPKA = pickle.load(open("Propka/PropkaOSCTdb.dat","rb"))
PSSM = pickle.load(open("PSSM/PSSMOSCTdb.dat","rb"))
SS = pickle.load(open("SecondaryStructure/SSOSCTdb.dat","rb"))
y=[float(x[2]) for x in data]
Xfa = [x[0] for x in data]
Xidx = [int(x[1]) for x in data]

for i in range(1,len(fasta),2):
    name = fasta[i-1].replace(">","").strip()
    cysteines_in_name = [j+1 for j, x in enumerate(fasta[i]) if x == 'C']
    cysteines_labeled = [Xidx[j] for j, name1 in enumerate(Xfa) if name1 == name]
    difference = list(set(cysteines_in_name) - set(cysteines_labeled))
    for diff in difference:
        Xfa += [name]
        Xidx += [diff]
        y += [-1]
y = np.array(y)

X=np.zeros((len(Xfa),n + n2 + n3 + 1 + n4 * 21 + n5 * 20 + 1 + 260 + 39))
for i in range(len(Xfa)):
    #X[i,0:n] = nearest(i)
    #X[i,n:n+n2] = modeller_cysteine_nearest(i)
    #X[i,n+n2:n+n2+n3] = modeller_carbon_nearest(i)
    X[i,n+n2+n3] = SASA[(Xfa[i],Xidx[i])]
    X[i,n+n2+n3+1:n+n2+n3+1+n4*21] = nearest_residue(i)
    X[i,n+n2+n3+1+n4*21:n+n2+n3+1+n4*21+n5*20] = modeller_alpha_nearest(i)
    X[i,n+n2+n3+1+n4*21+n5*20] = PROPKA[Xfa[i]][Xidx[i]]
    X[i,n+n2+n3+1+n4*21+n5*20+1:n+n2+n3+1+n4*21+n5*20+1+260] = PSSM[Xfa[i]][Xidx[i]]
    X[i,n+n2+n3+1+n4*21+n5*20+1+260:n+n2+n3+1+n4*21+n5*20+1+260+39] = SS[Xfa[i]][Xidx[i]]


yhat = clf.predict_proba(X)[:,1]
ypred = clf.predict(X)


Ymcc = list()
Xthresh = list()
print(y,yhat)
maxmcc = float('-inf')
for i in np.arange(0,1,.001):
    threshold = i
    ypred[yhat >= threshold] = 1
    ypred[yhat < threshold] = -1
    mcc = matthews_corrcoef(y,ypred)
    Ymcc.append(mcc)
    Xthresh.append(threshold)
    if mcc > maxmcc:
        maxmcc = mcc
        maxthreshold = threshold
ypred[yhat >= maxthreshold] = 1
ypred[yhat < maxthreshold] = -1

np.savetxt("OSCTdbMCCThresh.csv",np.stack((Xthresh,Ymcc),axis=-1),delimiter=",")
    
print("MCC:", matthews_corrcoef(y,ypred))
print("Max Threshold:", maxthreshold)

fpr, tpr, thresholds = metrics.roc_curve(y, yhat)
auc = metrics.auc(fpr, tpr)
print("AUC:", auc)
tn, fp, fn, tp = confusion_matrix(y, ypred).ravel()
specificity = tn / (tn+fp)
sensitivity = tp / (tp+fn)
accuracy = (specificity + sensitivity)/2
print("Sensitivity:",sensitivity)
print("Specificity:",specificity)
print("Accuracy:", accuracy)
plt.plot(fpr, tpr)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.show(block=True)
