'''
Created on 8 Feb 2018

@author: jerald

Plots Energy band diagram using Kronig Penney model
'''
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
from scipy.constants.constants import hbar, electron_mass, pi


def func(x, RHS):
    return (P*np.sin(x)/(x)) + np.cos(x) - RHS
    
def funcE(x):
    return (1/(2*m)*(hbar*x/a)**2)
    
def jouletoRyd(x):
    return (x/(1.6e-19))/13.6

def solver(region, RHS):
    x0 = opt.fsolve(func, region, args=(RHS))
    return x0

############Constants#############        
m = electron_mass#Particle mass
V_0 = 1 #QW height in eV
#b = #QW barrier thickness
a = 10*5.2917721067e-11#QW period
hbar = hbar
# P = (m*V_0*b*a/(hbar**2))
P = 3*pi/2

#############Start################
plt.close("all")

#Solver works best for values slightly lower than integer * pi, weird fluctuations for factors>=1
initialGuessFactor = 0.9
for j in np.arange(1., 7., 1.): #Max order before fsolver breaks is 6
    qaRegion = j*initialGuessFactor*pi 
    
    y = []
    x = []
    kaRange = np.asarray(np.linspace(-pi,pi,1000))
    RHSRange = np.cos(kaRange)
    for i in range(len(RHSRange)):
        x0 = solver(qaRegion, RHSRange[i])
        E = jouletoRyd(funcE(x0))
        x.append(kaRange[i])
        y.append(E[0])
    plt.plot(x,y)
    
plt.show()
    