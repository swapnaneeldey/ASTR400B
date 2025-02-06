#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import pandas as pd
import astropy.units as u


# In[26]:


def Read(filename):
    file = open(filename, "r")
    line1 = file.readline()
    label, value = line1.split()
    time = float(value)*u.Myr
    line2 = file.readline()
    label2, value2 = line2.split()
    total_particles = float(value2)
    file.close()
    data = np.genfromtxt(filename, dtype = None, names = True, skip_header=3)
    return time, total_particles, data





