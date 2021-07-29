#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 14:12:49 2021

@author: aoife
"""
import os.path
import matplotlib.pyplot as plt


def getmomentum(filename):
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
    for line in in_file:
        values=line.split()#split columns
        vx.append(float(values[0]))
        vy.append(float(values[1]))
        vz.append(float(values[2]))
        mass.append(float(values[3]))
        particles.append(float(values[4]))
        time.append(float(values[5]))

    in_file.close()

    M=mass[0]
    p=[]
    t=[]
    GMR=((G*(M)**3)/R)**(1/2)
    for i in range(len(mass)):
        p.append(((vx[i]**2+vy[i]**2+vz[i]**2)**(1/2)*mass[i])/GMR)
        t.append(time[i]/t_collapse)
    return p, particles, t, mass




rad=0.1
data=[]
radii=[]
for i in range(10):
    
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
