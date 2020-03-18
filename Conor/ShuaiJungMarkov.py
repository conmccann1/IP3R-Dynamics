# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 13:45:41 2020

@author: Conor
"""
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
plt.close('all')

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
I = 0.8
N = 20

tmax = 500

def ah(I):
    return a2*d2*((I + d1)/(I + d3))

def bh(c):
    return a2*c


def minf(I):
    return I/(I+d1)


def ninf(c):
    return c/(c+d5)

def Jchannel(c,I,cER,NopenT):
    n = ninf(c)
    m = minf(I)
    return c1*v1*(m**3)*(n**3)*(NopenT/N)*(c - cER)


def Jpump(c):
    return (v3*(c**2))/((k3**2)+(c**2))

def Jleak(c,cER):
    return c1*v2*(c-cER)

dt = 0.01
                     # each column represents a gate
t = np.arange(0,tmax,dt)              
h = np.random.randint(2,size=(N,3)) # each row represents a channel
c = np.zeros(len(t))
Nopen=0
for j in np.arange(0,len(h),1):
    
    if sum(h[j,:]) == 3:
        Nopen += 1
NopenT = np.zeros(len(t))
NopenT[0] = Nopen

for k in np.arange(0,len(t)-1,1):
    Nopen = 0
    print(k)
#    print('\n')
    
    for j in np.arange(0,len(h),1):     
        
        for i in np.array([0,1,2]):
            
            
            mewT = ah(I)+bh(c[k])
#            print(mewT)
            u1,u2 = sp.rand(2)
            
            if h[j,i] == 0 and u1 < np.exp(-ah(I)*dt):
#                print(u1)
                h[j,i] = 0
                
#                print('same')
            
                
#                
            elif h[j,i] == 1 and u2 < np.exp(-bh(c[k])*dt):
                h[j,i] = 1
#                print('same')
                
            else:
                h[j,i] = 1 - h[j,i]
#                print('my guy')
                
#            print(h)
#            if h[j,i] == 0:
#                p = np.random.exponential(ah(I)*dt)
#                
#                if p <= 1 - np.exp(-ah(I)*dt):
#                    h[j,i] = 1
#                else:
#                    h[j,i] = 0
#                    
#            elif h[j,i] == 1:
#                p = np.random.exponential(bh(c[k])*dt)
#                
#                if p <= 1 - np.exp(-bh(c[k])*dt):
#                    h[j,i] = 0
#            else:
#                continue
#                print('error')
                
        if sum(h[j,:]) == 3:
            Nopen += 1
#            print(Nopen)
#    print(h) 
    NopenT[k] = Nopen
#    print(NopenT[k])
    cER = (c0 - c[k])/c1

    c[k+1] = c[k] + dt*(-Jchannel(c[k],I,cER,NopenT[k]) - Jpump(c[k]) - Jleak(c[k],cER))

Nfrac = NopenT/N

    
#    print(Nopen)
#    print(h)


    
plt.figure()
plt.plot(t,c)
plt.figure()
plt.plot(t,NopenT/N)





    
    
    
    
    
#    j = np.zeros(3)
#    tot = ah(I) + bh(I)
#    tau = np.random.exponential(1/tot)
#    t.append(tau)
#    
#    for i in np.array([0,1,2]):        
#        u1 = np.random.uniform(0,1)
#        u2 = np.random.uniform(0,1)
#    
#        if u2*tot <= bh(I):
#            j[i] = 0 # closed
#        else:
#            j[i] = 1 # open
#    
#    if sum(j) == 3:
#        channel = 'open'
#    else:
#        channel = 'closed'
        
    



















