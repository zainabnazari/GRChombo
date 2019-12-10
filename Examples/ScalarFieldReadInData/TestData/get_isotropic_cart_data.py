import os, sys
import numpy as np
from scipy import interpolate, integrate, optimize
from scipy.integrate import odeint
from scipy.interpolate import interp1d

# function that returns dr/dR
# Es. D22 of https://arxiv.org/pdf/gr-qc/0410040.pdf
def drdR(r, R, a):
    out =  a * r / R
    return out

# Read in data in areal polar coordiantes
grrID    = np.loadtxt("grr001.csv")
PiID  = np.loadtxt("Pi001.csv")

# set up grid in radial direction in areal polar coordinates
R  = 0
dR = 0.01;
length = np.size(grrID)
R = np.linspace(0, 0.01*(length-1), num=length)

# invert the ordering since we start integrating and the outermost point
Rinv = R[::-1]
a    = np.sqrt(grrID[::-1])
r0   = Rinv[0]/a[0]*(0.5 + 0.5*np.sqrt(a[0]))**2
print("r0 is", r0)
print("R0 is", Rinv[0])

#storage for solution
r  = np.empty_like(Rinv)

#ODE integration
for i in range(1,length):
    Rspan  = [Rinv[i-1], Rinv[i]]
    r_out  = odeint(drdR, r0, Rspan, args=(a[i],)) 
    r[i] = r_out[1][0]
    r0 = r_out[1][0]

#revert ordering agian and compute psi
r = r[::-1]
psi  = np.sqrt(R/(r + 1e-100))

#set end points by hand, 0th order extrapolation for psi
r[0] = 0.0
psi[0] = psi[1]
r[length-1] = R[length-1]
psi[length-1] = psi[length-2]

dx = 1.0
N = 100
negative_offset = -20.0

#transform Pi data
Pi = PiID
f_Pi = interp1d(r, Pi)
file = open('Pi001_cartesian.txt', 'w')
for i in range(N) :
    for j in range(N) :
        for k in range(N) :
            x = i*dx + negative_offset
            y = j*dx + negative_offset
            z = k*dx + negative_offset
            R_new = np.sqrt(x*x + y*y + z*z)
            Pi_here = f_Pi(R_new)
            file.write(str(x).ljust(20))  
            file.write(str(y).ljust(20))  
            file.write(str(z).ljust(20))  
            file.write(str(Pi_here).ljust(20)+"\n")

#transform chi data
chi = psi**(-0.25)
f_chi = interp1d(r, chi)
file = open('chi001_cartesian.txt', 'w')
for i in range(N) :
    for j in range(N) :
        for k in range(N) :
            x = i*dx + negative_offset
            y = j*dx + negative_offset
            z = k*dx + negative_offset
            R_new = np.sqrt(x*x + y*y + z*z)
            chi_here = f_chi(R_new)
            file.write(str(x).ljust(20))  
            file.write(str(y).ljust(20))  
            file.write(str(z).ljust(20))  
            file.write(str(chi_here).ljust(20)+"\n")
