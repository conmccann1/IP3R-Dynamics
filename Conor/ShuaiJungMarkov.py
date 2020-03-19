# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 13:45:41 2020

@author: Conor
"""
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
plt.close('all')

""" parameters """
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

""" parameters you might want to play with """
I = 0.3
N = 20
dt = 0.01
tmax = 500

""" Functions from Paper """
def Jchannel(c,I,cER,NopenT): #(3)
    n = ninf(c)
    m = minf(I)
    return c1*v1*(m**3)*(n**3)*(NopenT/N)*(c - cER)

def Jpump(c): #(4)
    return (v3*(c**2))/((k3**2)+(c**2))

def Jleak(c,cER): #(5)
    return c1*v2*(c-cER)

def ah(I): #(6)
    return a2*d2*((I + d1)/(I + d3))

def bh(c): #(6)
    return a2*c

def minf(I): #(6)
    return I/(I+d1)

def ninf(c): #(6)
    return c/(c+d5)


""" Set up for simulation """
# set up time vector
t = np.arange(0,tmax,dt)


"""     IMPORTANT     """            
# 'h' is a matrix with 20 rows, representing channels
# and 3 columns representing gates. 0 represents a closed gate and
# 1 represents an open gate. A channel (row) is open if all 
# 3 gates (columns) are open (=1).
h = np.random.randint(2,size=(N,3)) #random initial config
""" ^^^ IMPORTANT ^^^ """


# initialse [Ca2+] vector
c = np.zeros(len(t))

# we want to keep track of number of gates (columns) open (=1) so...
Nopen=0

# check initial h config for how many gates (columns) are open (equal to 1)
for j in np.arange(0,len(h),1):
    if sum(h[j,:]) == 3:
        Nopen += 1
        
# we want to keep track of the total number of open channels so...
NopenT = np.zeros(len(t))
NopenT[0] = Nopen


""" Begin Simulation """
# for every time point we want to see if a gate opens, closes or stays the same
for k in np.arange(1,len(t)-1,1):
    # as above
    Nopen = 0
    
    # print some values of k to keep track of how long is lef in the console
    if np.mod(k,1000) == 0:    
        print(k)
        
    # for every ROW in h
    for j in np.arange(0,len(h),1):      # for every roq
        
        # for every COLUMN in h
        for i in np.array([0,1,2]):
            
            # generate 2 random numbers (0,1) for open and closed gates             
            u1,u2 = sp.rand(2)
            
            # if a gate is closed AND remains closed (u1<np.exp(-ah...))
            if h[j,i] == 0 and u1 < np.exp(-ah(I)*dt):
                # gate remains closed
                h[j,i] = 0
            
            # if gate is open AND remains open (u2<...)
            elif h[j,i] == 1 and u2 < np.exp(-bh(c[k])*dt):
                # gate remains open
                h[j,i] = 1
            #otherwise gate switcehs configuration (O->C or C->O)
            else:
                h[j,i] = 1 - h[j,i]
        
        # find out what channels are open and update no. of gates open param
        if sum(h[j,:]) == 3:
            Nopen += 1
    
    # record total open receptors for time k
    NopenT[k] = Nopen
    
    # Use Eulers method to find out calcium conc using NopenT and N (see paper)
    cER = (c0 - c[k])/c1
    c[k+1] = c[k] + dt*(-Jchannel(c[k],I,cER,NopenT[k]) - Jpump(c[k]) - Jleak(c[k],cER))

# fraction of gates open at each time point
Nfrac = NopenT/N

# plot calcium conc and Nfrac over time
plt.figure()
plt.subplot(211)
plt.plot(t,c)
plt.subplot(212)
plt.plot(t,NopenT/N)





    
    
    
    

        
    



















