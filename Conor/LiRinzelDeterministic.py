# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 15:23:29 2020

@author: Conor
"""
import numpy as np
from scipy.integrate import odeint
from matplotlib import pyplot as plt


c1 = 0.185
v1 = 6.0
v2 = 0.11
v3 = 0.9
k3 = 0.1
d1 = 0.13
d2 = 1.049
d3 = 0.9434
d5 = 0.08234
a2 = 0.2
c0 = 2.0
I = 0.3
def a(I):
    return (a2*d2)*((I+d1)/(I+d3))

def b(c):
    return a2*c

def minf(I):
    return I/(I+d1)


def ninf(c):
    return c/(c+d5)

def Jchannel(c,h,I,cER):
    n = ninf(c)
    m = minf(I)
    return c1*v1*(m**3)*(n**3)*h**3*(c - cER)


def Jpump(c):
    return (v3*(c**2))/((k3**2)+(c**2))

def Jleak(c,cER):
    return c1*v2*(c-cER)

def sys(x,t):
    
    c,h = x[0],x[1]
    cER = (c0 - c)/c1
    dc = -Jchannel(c,h,I,cER) - Jpump(c) - Jleak(c,cER)
    dh = a(I)*(1-h) - b(c)*h
    
    return np.array([dc,dh])

t = np.arange(0,100,0.0001)
y0 = np.array([0.2,0.5])

y = odeint(sys,y0,t)

plt.plot(t,y[:,0])
#
#
#plt.figure()
#t = np.arange(0,80,0.001)
#y0 = np.array([0.2,0.5])
#
#Irng = np.arange(0,0.9,0.001)
#
#
#for i in np.arange(0,len(Irng),1):
#    I = Irng[i]
#    y = odeint(sys,y0,t)
#  
#    lasts = y[-20000:,0]
#    
#    
#    x = np.repeat(Irng[i],1)
#    
#    plt.plot(x,np.amax(lasts),'rx')
#    plt.plot(x,np.amin(lasts),'rx')
#    plt.xlabel('[IP3]'), plt.ylabel('[Ca2+]')
