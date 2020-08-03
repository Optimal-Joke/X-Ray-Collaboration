#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 18:05:13 2020

@author: susannahabrams
"""

import math
pi = math.pi
import matplotlib.pyplot as plt
import numpy as np
from random import randrange
import random

"""randomly select a grain of dust (ra, dec, d). Each  has an associated 
angle through which a scattered photon would reach us.
Randomly select a scattering angle (from the distribution), 
and plot all of the photons that reach us for blocks of time.""" 

#distribution of scattering angles:
x = np.arange(0,0.4, step = 0.0001)
y =[]
def equation(x):
    for n in x:
        c = 22.523*np.sqrt(2*np.pi)/0.0292
        f= c*np.sin(n+0.0292)
        exp_num = (n+0.0292)**2
        g = np.exp(exp_num/-0.0017059)
        m = f*g
        y.append(m)        
equation(x)       
p_total = sum(y)
p_list = [n / p_total for n in y]
a = random.choices(x, weights = p_list, k =10000)
np.histogram(a, bins = 100) 

none = []
times = []
Decs = []
RAs = []
radii = []

#making the layer of dust:
event = (180.0,45.0)
photons = np.arange(0,10**8)
for photon in photons:
    ra = randrange(0, 3600, step=2)/10
    dec = randrange(0, 900, step =2)/10
    r = np.sqrt((event[0]-ra)**2+(event[1]-dec)**2)
    #d = randrange(1.0,100.0)
    d= 50
    ang = np.arctan(r/d)
    ang = round(ang, 3)
    
#scattering the photons:
    sa = round(random.choice(a), 3)
    
#collecting all the points in the layer of dust where the photons 
#scattered to us, along with how long it took those photons to reach us:
    if ang == sa:
        dist = d/np.cos(ang)
        c = 3
        t = dist/c
        times.append(t)
        Decs.append(dec)
        RAs.append(ra)
        radii.append(r)
    else:
        none.append(sa)
        
#sort the lists by time:
times, Decs, RAs,radii = zip(*sorted(zip(times,Decs,RAs,radii)))

#change the range of ra and dec to get different halos:
plt.scatter(Decs[100:200],RAs[100:200], color = "blue")



     