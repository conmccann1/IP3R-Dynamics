# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 15:43:03 2020

DeYoung + Keizer (1992)

@author: Conor
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import fsolve
plt.close('all')

""" Parameters """
#Binding affinities
a1 = 400
a2 = 0.2
a3 = 400
a4 = 0.2
a5 = 20
a = np.array([a1,a2,a3,a4,a5])

# Kd/Dissociation constants
d1 = 0.130
d2 = 1.049
d3 = 0.9434
d4 = 0.1445
d5 = 0.08234
d = np.array([d1,d2,d3,d4,d5])

# Dissociation rates
b = np.zeros(5)
for i in np.arange(0,5,1):
    b[i] = a[i]*d[i]
        
c0 = 2
c1 = 0.185
u1 = 6
u2 = 0.11
u3 = 0.9
k3 = 0.1
v = 0  
cER = 0.125

Istar = 0.25
Ir = 1
Ip = 8
""" Steady state open probabilities """

# Steady state open probability function

#def pO(c,i):
#    return ((c*i*d2) / ( (c*i + i*d2 + d1*d2 +c*d3)*(c+d5) ))**3
#
## Plot pO for range of Ca2+ concs
#for i in np.array([2,1,0.5,0.25]):
#    c = np.arange(0.01,10,0.001)
#
#    plt.plot(np.log10(c),pO(c,i),label='IP3 = %2f' %i)
#    plt.xlabel('Log10([Ca2+])')
#    plt.ylabel('Probability of open channel')
#plt.legend()
#
## Plot pO for range of IP3 concs
#plt.figure()
#c = 0.1
#i = np.arange(0.015,20,0.001)
#
#plt.plot(np.log10(i),pO(c,i))
#plt.xlabel('Log10([Ca2+])')    
#plt.ylabel('Probability of open channel')
#plt.ylim([0,0.15])


""" Ca2+ Oscillations """
# Functions

#Outward Flux
def J1(x110,cER,c):
    return c1*(u1*(x110**3) + u2)*(cER - c)

#Inward Flux
def J2(c):
    return u3*((c**2) / (c**2 + k3**2))

#Impulse
def f(t):
    if v == 1:
        f = 0
    else:
        f = 1
    return f
        

def sys(x,t,p):
    s0,s1,s2,s3,v0,v1,v2,v3,c,I = x[0:10] 
    
    if t > 110 and t < 110.5:
        f = 1
    elif t > 60 and t < 60.5:
        f = 1
    elif t > 40 and t < 40.5:
        f=1
    else:
        f = 0


        
    ds0 = -a[0]*I*s0 + b[0]*s2 - a[3]*c*s0 + b[3]*s1 - a[4]*c*s0 + b[4]*v0
    ds1 =  a[3]*c*s0 - b[3]*s1 - a[2]*I*s1 + b[2]*s3 - a[4]*c*s1 + b[4]*v1
    ds2 =  a[0]*I*s0 - b[0]*s2 - a[1]*c*s2 + b[1]*s3 - a[4]*c*s2 + b[4]*v2
    ds3 =  a[1]*c*s2 - b[1]*s3 + a[2]*I*s1 - b[2]*s3 - a[4]*c*s3 + b[4]*v3

    dv0 = -a[0]*I*v0 + b[0]*v2 - a[3]*c*v0 + b[3]*v1 + a[4]*c*s0 - b[4]*v0
    dv1 =  a[3]*c*v0 - b[3]*v1 - a[2]*I*v1 + b[2]*v3 + a[4]*c*s1 - b[4]*v1
    dv2 =  a[0]*I*v0 - b[0]*v2 - a[1]*c*v2 + b[1]*v3 + a[4]*c*s2 - b[4]*v2
    dv3 =  a[1]*c*v2 - b[1]*v3 + a[2]*I*v1 - b[2]*v3 + a[4]*c*s3 - b[4]*v3
    
    cER = (c0-c)/c1
    dc = J1(v2,cER,c) - J2(c) 
    
    if p == 1:    
        dI = Ir*(Istar-I) + Ip*f 
    else:
        dI = 0
    
    return [ds0,ds1,ds2,ds3,dv0,dv1,dv2,dv3,dc,dI]

""" IP3 Pulse """
###parameter changes
Istar = 0.24

t = np.arange(0,250,0.0001)
           #ds0, ds1,ds2,ds3,dv0,dv1,dv2,dv3,dc,dI
y0 = np.array([0.5,0,0,0,0.5,0,0,0,0.125,0.24])
p = 1
y = odeint(sys,y0,t,args=(p,))
# plt.figure()
plt.plot(t,y[:,-1])
plt.ylabel('[IP3]')
plt.figure()
plt.plot(t,y[:,-2])
plt.ylabel('[Ca2+]')
plt.xlabel('time')
plt.figure()
plt.plot(t,y[:,-4])
plt.ylabel('x110')
plt.xlabel('time')



""" Hopf Bifurcation """
#Istar = 0.25
#plt.figure()
#t = np.arange(0,80,0.001)
#y0 = np.array([0.5,0,0,0,0.5,0,0,0,0.125,0])
#
#Irng = np.arange(0,0.9,0.001)
#p = 0
#for i in np.arange(0,len(Irng),1):
#    y0[-1] = Irng[i]
#    y = odeint(sys,y0,t,args=(p,))
#  
#    lasts = y[-20000:,-2]
#    
#    print(i)
#    x = np.repeat(Irng[i],len(lasts))
#    
#    plt.plot(x,lasts,'x')
#    plt.xlabel('[IP3]'), plt.ylabel('[Ca2+]')

""" x0ik vs x1ik """
#plt.close('all')
#Istar = 0.5
#
#t = np.arange(0,100,0.001)
#y0 = np.array([0.5,0,0,0,0.5,0,0,0,0.125,0.5])
#p = 0
#y = odeint(sys,y0,t,args=(p,))
#plt.figure()
#plt.plot(t,y[:,0])
#plt.figure()
#plt.plot(t,y[:,2])
#plt.figure()
#plt.plot(t,y[:,-2])
#plt.figure()
#plt.plot(y[:,0],y[:,2])





