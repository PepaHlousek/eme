#  eme.py
#  
#  Photon statistics retrieval from measured multichannel coincidences using the EME method.

import math
import scipy.special as sc
import numpy as np
from scipy.misc import factorial

# response function (matrix) of multichannel detector
# eta - quantum efficiency
# mMax - nuber of detection channels
# nMax - nuber of elements of photon statistics
def DET(mMax,nMax,eta):
	det=np.zeros((mMax+1,nMax+1))
	for m in range(mMax+1):
		for n in range(nMax+1):
			if m>n:
				det[m][n]=0
			elif m<n:
				summary=[]
				for j in range(0,m+1):
					summary.append(((-1)**j)*sc.binom(m,j)*((1-eta)+((m-j)*eta)/mMax)**n)
				det[m][n]=sc.binom(mMax,m)*np.sum(summary)
			else:
				det[m][n]=(eta/mMax)**n*(factorial(mMax)/factorial(mMax-n))	
	return det

# photon statistics retrieval using EME algorithm				
# iterations - the maximum number of iterations (actual number is much lower, based on "epsilon")
# epsilon - the distance between n-th and (n+1)-th iterations when the process is stopped
# c - click statistics
# pn - photon statistics
def EME(mMax,nMax,eta,det,l,c):
	iterations = 10**10
	epsilon = 10**(-12)
	pn=np.array([1./(nMax+1)]*(nMax+1))
	iteration=0
	while (iteration<iterations):
		EM=np.dot(c/np.dot(det,pn),det)*pn
		E=l*(np.log(pn)-np.sum(pn*np.log(pn)))*pn
		E[np.isnan(E)]=0.0
		EME=EM-E
		dist = np.sqrt(np.sum((EME-pn)**2))
		if dist<=epsilon:
			break
		else:
			pn=EME
			iteration+=1
	return EME

# numerical values
mMax=10
nMax=50
eta=0.5
l=10**(-3)
det=DET(mMax,nMax,eta)

# measured click statistics
c=np.array([6.73794700e-03,4.37104954e-02,1.27601677e-01,2.20741125e-01\
,2.50599060e-01,1.95082729e-01,1.05461930e-01,3.90945126e-02\
,9.51054071e-03,1.37104223e-03,8.89424261e-05])

# retrieval
p=EME(mMax,nMax,eta,det,l,c)

print p
