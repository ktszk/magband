#!/usr/bin/env python
#-*- coding:utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

fname="fort.60"
#fname="fort.72"
#fname2="fort.73"
sw=True
temp=1.e-2
if sw:
    klist=[0,16,32]
    #klabel=['(1,0,0)','(1,0,1)']
    #klist=[0,32,64]
    #klabel=['(0,0)','(1,0)','(1,1)','(0,0)']
    klabel=['(0,0)','(1,0)','(1,1)']
else:
    klist=[0,16]
    klabel=['0','$\pi$']

ncnt=100
ic=2

data=np.loadtxt(fname)
x=data[:,0]
nE=(np.where(x==0.)[0]).size
nk=x.size//nE
x=x.reshape(nk,nE)
y=data[:,1].reshape(nk,nE)*1.0e3

try:
    fname2
except NameError:
    z=(data[:,2] if sw else data[:,3]).reshape(nk,nE)
else:
    z1=(data[:,2] if sw else data[:,3]).reshape(nk,nE)
    data2=np.loadtxt(fname2)
    z2=(data2[:,2] if sw else data2[:,3]).reshape(nk,nE)
    z=z1+z2

#fb=np.array([[1. if j==0 else (.5+.5/np.tanh(.5*j/temp))/(.5+.5/np.tanh(0.5*y[0][1]/temp)) for j in yy] for yy in y])
#z=z*fb

colors=['hot','spectral','jet','summer','winter','hsv']

plt.rcParams['font.family']='serif'
plt.rcParams['text.usetex']=True
plt.contourf(x,y,z,ncnt)
if sw:
    klist=klist+[np.max(x)]
    plt.ylabel('Energy (meV)')
    plt.xlabel('(1,0,L)')
else:
    plt.yticks(klist,klabel)
exec('plt.%s()'%colors[ic])
plt.xticks(klist,klabel)
plt.colorbar()
plt.show()
