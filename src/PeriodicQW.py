'''
Created on 8 Feb 2018

@author: jerald
'''
import numpy as np
import scipy.optimize as opt
from scipy.constants.constants import hbar, electron_mass, pi

m = electron_mass#Particle mass
V_0 = 1 #QW height in eV
#b = #QW barrier thickness
a = 10*5.2917721067e-11#QW period
hbar = hbar

# P = (m*V_0*b*a/(hbar**2))
P = 3*pi/2
RHS = np.cos(pi)

def func(x):
    return (P*np.sin(x)/(x)) + np.cos(x) - RHS

def funcE(x):
    return (1/(2*m)*(hbar*x/a)**2)

def jouletoRyd(x):
    return (x/(1.6e-19))/13.6

x0 = opt.fsolve(func, 1.9*pi)
print x0

Energy = jouletoRyd(funcE(x0))
print Energy