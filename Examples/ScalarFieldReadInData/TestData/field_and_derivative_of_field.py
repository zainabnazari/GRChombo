#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 12:12:11 2020

@author: zainabnazari
"""

#Gaussian Random Field and field derivative with the power spectrum associated to Bunch-Davies.

import numpy as np

N = 16

x = np.array([a for a in range(N)])
y = np.array([a for a in range(N)])
z = np.array([a for a in range(N)])

import numpy as np
X, Y, Z = np.meshgrid(x, y, z)
data = np.random.normal(0,1,(N,N,N))


Fdata = np.fft.fftn(data);


k1 = np.array(2 * np.pi / N * np.arange(0, N/2))
k2 = np.array(2 * np.pi / N * np.arange(-N/2, 0))
k = np.asarray(np.concatenate((k1, k2), axis = None));


def PS(k, n, A):
    if k == 0:
        return 0;
    elif k != 0:
        return A * (k**(-n))
    
# here is the function for the field derivative
def DS(k, n, A):
     if k == 0:
        return 0;
     elif k != 0:
        return A * (k**(n))+(61.6955)/(k**(n)) 
    
    
a=np.meshgrid(k**2,k**2,k**2)[0]
b=np.meshgrid(k**2,k**2,k**2)[1]
c=np.meshgrid(k**2,k**2,k**2)[2]



ksqrd=np.asarray(np.sqrt(a+b+c))

  
  

ksqrdravel=np.asarray(ksqrd.ravel())



PSravel = np.asarray([PS(x,1,1/2) for x in ksqrdravel if x != 0])

PSravel0 = np.concatenate((np.array([0.]),PSravel))
RootPS = PSravel0.reshape(N,N,N)

# this part is for the field derivative

DSravel = np.asarray([DS(x,1,1/2) for x in ksqrdravel if x != 0])
#print(PSravel)
DSravel0 = np.concatenate((np.array([0.]),DSravel))
RootDS = DSravel0.reshape(N,N,N)

PSFourier = Fdata*RootPS # it multiplies element-wise
#print(PSFourier)#

# this part is for field derivative

DSFourier = Fdata*RootDS 


FieldRealization = np.fft.ifftn(PSFourier)
RealFieldRealization = np.real(FieldRealization)

#  for field derivative
DFieldRealization = np.fft.ifftn(DSFourier)
DRealFieldRealization = np.real(DFieldRealization)

#*********************************
file = open('field.txt', 'w')

for (idx1, val1) in enumerate(RealFieldRealization):
    for (idx2, val2) in enumerate(val1):
        for (idx3, val3) in enumerate(val2):
             file.write(str(idx1).ljust(20))
             file.write(str(idx2).ljust(20))
             file.write(str(idx3).ljust(20))
             file.write(str(val3).ljust(20)+"\n")
        
#*********************************

# this part is for the field derivative
file = open('dfield.txt', 'w')
for (idx1, val1) in enumerate(DRealFieldRealization):
    for (idx2, val2) in enumerate(val1):
        for (idx3, val3) in enumerate(val2):
             file.write(str(idx1).ljust(20))
             file.write(str(idx2).ljust(20))
             file.write(str(idx3).ljust(20))
             file.write(str(val3).ljust(20)+"\n")

