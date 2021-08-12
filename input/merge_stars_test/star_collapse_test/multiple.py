#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 14:12:49 2021

@author: Aoife Flood
Description:
           This program takes files named "merge_data_XX.txt" output by data.py and uses the to create plots of mass, no. of particles and momentum for multiple merge radii.
"""
import os.path
import matplotlib.pyplot as plt


def getmomentum(filename):
    """
    This function opens a file and puts the data in it into arrays
    
    Parameters
    ----------
    filename: string
    
    Returns
    -------
    p         : array
        array for the magnitude of the momentum
    particles : array
        array for the number of particles
    t         : array
        array for time/collapse time
    mass      : array
        array for total mass at each time interval
    """
    in_file = open(filename, 'r')#Open file to read
    t_collapse=4.8987e16
    G=6.67e-8  
    R=3.086e24
    mass=[]
    time=[]
    particles=[]
    vx=[]
    vy=[]
    vz=[]
    p=[]
    for line in in_file:
        values=line.split()#split columns
        vx.append(float(values[0]))
        vy.append(float(values[1]))
        vz.append(float(values[2]))
        mass.append(float(values[3]))
        particles.append(float(values[4]))
        p.append(float(values[5]))
        time.append(float(values[6]))

    in_file.close()

    return p, particles, time, mass




rad=0.1
data=[]
radii=[]
for i in range(20):
    
    radstr=str(rad)
    filename="merge_data_"+radstr+".txt"
    if os.path.isfile(filename) == True:
        data.append(getmomentum(filename))
        radii.append(rad)
    rad=round(rad+0.1,2)

fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=[10,8], sharex=True)
fig.suptitle("Momentum and No. of Particles VS Time")
plt.xlabel("Time/Collapse Time")
ax1.set(ylabel="Mass (g)")
ax2.set(ylabel="Momentumsqrt(GM^3R)")
ax3.set(ylabel="No. of Particles")

for i in range(len(data)):
    ax1.plot(data[i][2],data[i][3], label= str(radii[i]))
    ax2.plot(data[i][2],data[i][0], label= str(radii[i]))
    ax3.plot(data[i][2],data[i][1], label= str(radii[i]))
    
handles, labels = ax1.get_legend_handles_labels()
fig.legend(handles, labels, loc='upper left')
fig.savefig("Mass_momentum_particles_graph_many_radii.png")
