import yt
import numpy as np
import matplotlib.pyplot as plt
import cello_parse as cp
import io
import argparse as ap
import matplotlib
import os.path

G=6.67e-8
matplotlib.use("Agg")

parser = ap.ArgumentParser(
    description="""                                                                                                                                                                                          
    Reads a data directory containing a series of snapshots of the Shu Collapse                                                                                                                              
    problem, and makes density slice images through the centre of collapse.                                                                                                                                  
    """
    )

parser.add_argument(
    "-d",
    "--data-directory",
    help="The directory we are going to be reading. Default: ./",
    required=False,
    default=".",
    type=str,
)

parser.add_argument(
    "-p",
    "--parameter-file",
    help="""                                                                                                                                                                                                 
    The name of the Enzo-E parameter file corresponding to the data. Must be in                                                                                                                              
    the data directory. Default = parameters.out                                                                                                                                                             
    """,
    required=False,
    type=str,
    default="parameters.out",
)
args = vars(parser.parse_args())
with io.open(args["data_directory"] + "/" + args["parameter_file"],
             encoding = "utf-8") as f:
    params = cp.load(f)

#Get merge radius to use for file name                                                                                                                                                                       
Merge_radius = str(params.Method.merge_stars.merging_radius_cells)
filename="merge_data_"+Merge_radius+".txt"

#Get parameters from param file                                                                                                                                                                              
infall_speed=params.Initial.collapse_stars.infall_speed
t_collapse=params.Stopping.time
method_list=params.Method.list
radius=params.Initial.collapse_stars.truncation_radius

#Check if file already exists for this test and delete it                                                                                                                                                    
if os.path.isfile(filename)== True:
    os.remove(filename)

yt.enable_parallelism()
ts =  yt.DatasetSeries("%s/Dir_Collapse-STARS_????/Dir_Collapse-STARS_????.block_list" %args["data_directory"])

storage = {}

#arrays for total velocities, masses and number of particles at each time step                                                                                                                               
t_mass=[]
time=[]
t_particles=[]
p=[]
t=[]


#Getting total velocities and mass for each timestep from directories                                                                                                                                        
for ds in ts.piter():

    box = ds.box(left_edge = -ds.domain_width/2.0,right_edge = ds.domain_width/2.0)
    vx = box["star","vx"]
    vy = box["star","vy"]
    vz = box["star","vz"]
    mass = box["star","mass"]
    particles=len(vx)
    vvx = float(np.sum(vx))
    vvy = float(np.sum(vy))
    vvz = float(np.sum(vz))
    m = float(np.sum(mass))

    
    t_mass.append(m)
    t_particles.append(particles)
    M=t_mass[0]
    
    #Check if gravity is on to determine how to normalise momentum                                                                                                                                               
    if 'gravity' in method_list:
        N=((G*(M**3))/radius)**(1/2)
        
    else:
        N=M*infall_speed
    
    momentum=(float(((np.sum(vx))**2+(np.sum(vy))**2+(np.sum(vz))**2)**(1/2)*(np.sum(mass))/N))
    time=(float((float(ds.current_time))/t_collapse))

    p.append(momentum)
    t.append(time)
    filename="merge_data_"+Merge_radius+".txt"
    out_file=open(filename,'a')
    print(vvx,vvy,vvz,m,particles,momentum,time,file=out_file)
    out_file.close()


#Plot Mass, No. of Particles and Momentum VS Time on 3 different subplots                                                                                                                                    
fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=[10,8], sharex=True)

fig.suptitle("Mass, Momentum and No. of Particles VS Time")
plt.xlabel("Time/Collapse Time")


ax1.set(ylabel="Mass (g)")
ax1.plot(t,t_mass)
handles, labels = ax1.get_legend_handles_labels()

fig.legend(handles, labels, loc='upper left')
ax2.set(ylabel="No. of Particles")
ax2.plot(t,t_particles)

ax3.set(ylabel="Momentum")
ax3.plot(t,p)
fig.savefig("Mass_momentum_particles_graph2.png")


