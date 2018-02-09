'''
Created on 8 Feb 2018

@author: jerald

Plots Energy band diagram using Kronig Penney model
'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import scipy.optimize as opt
from scipy.constants.constants import hbar, electron_mass, pi, elementary_charge

plt.close("all")

def func(x, *data):
    P, RHS = data
    return (P*np.sin(x)/(x)) + np.cos(x) - RHS
    
def funcE(x,m):
    return (1/(2*m)*(hbar*x/a)**2)
    
def jouletoRyd(x):
    return (x/(1.6e-19))/13.6

def solver(region, data):
    x0 = opt.fsolve(func, region, args=(data))
    return x0

def plot(m, V0):
    P = (m*V0*b*a/(hbar**2))
    ys = []
    #Solver works best for values slightly lower than integer * pi, weird fluctuations for factors>=1
    initialGuessFactor = 0.8
    for j in np.arange(1., 3., 1.): #Max order before fsolver breaks is 6
        qaRegion = j*initialGuessFactor*pi 
        y = []
        x = []
        kaRange = np.asarray(np.linspace(-pi,pi,100))
        RHSRange = np.cos(kaRange)
        for i in range(len(RHSRange)):
            data = (P, RHSRange[i])
            x0 = solver(qaRegion, data)
            E = jouletoRyd(funcE(x0,m))
            x.append(kaRange[i])
            y.append(E[0])
        ys.append(y)
    return x, ys

############Constants############# 
#Fundamental Constants       
m_e = electron_mass#Particle mass
hbar = hbar

#Variables
V0i = 5. #QW height in eV
b = 1.0e-10 #QW barrier thickness
a = 1.5e-10 #QW period
mi = 1.0 * m_e #Mass factor * m_e

axcolor = 'lightgoldenrodyellow'
fig = plt.figure()
axdata = plt.subplot2grid((7,4),(0,0),colspan=4,rowspan=4)
axV0 = plt.subplot2grid((7,4),(-2,0),colspan=4, axisbg=axcolor)
axm = plt.subplot2grid((7,4),(-1,0),colspan=4, axisbg=axcolor)

#############Start################
x0points, y0s = plot(mi,V0i*elementary_charge)
l0, = axdata.plot(x0points, y0s[0])
l1, = axdata.plot(x0points, y0s[1])

sV0 = Slider(axV0, 'V0', 5., 20.0, valinit=V0i)
sm = Slider(axm, 'Mass', 0.1, 5.0, valinit=mi)

def update(val):
    V0 = sV0.val
    m = sm.val
    xpoints, ys = plot(m*m_e, V0*elementary_charge)
    l0.set_ydata(ys[0])
    l1.set_ydata(ys[1])
    fig.canvas.draw()
 
sV0.on_changed(update)
sm.on_changed(update)  
plt.show()

# button = Button(axreset, 'Reset', color=axcolor, hovercolor='0.975')
#  
# def reset(event):
#     sV0.reset()
#     sm.reset()
# button.on_clicked(reset)

# rax = plt.axes([0.025, 0.5, 0.15, 0.15], facecolor=axcolor)
# radio = RadioButtons(rax, ('red', 'blue', 'green'), active=0)
# 
# def colorfunc(label):
#     l0.set_color(label)
#     l1.set_color(label)
#     fig.canvas.draw_idle()
# radio.on_clicked(colorfunc)