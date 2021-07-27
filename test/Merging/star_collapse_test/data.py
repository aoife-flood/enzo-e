import yt
import numpy as np
import matplotlib.pyplot as plt
import cello_parse as cp
import io
import argparse as ap
import matplotlib
import os.path
 
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

Merge_radius = str(params.Method.merge_stars.merging_radius_cells)
filename="merge_data_"+Merge_radius+".txt"

infall_speed=params.Initial.collapse_stars.infall_speed
t_collapse=params.Stopping.time
method_list=params.Method.list
radius=params.Initial.collapse_stars.truncation_radius

if os.path.isfile(filename)== True:
    os.remove(filename)

yt.enable_parallelism()
ts =  yt.DatasetSeries("%s/Dir_Collapse-STARS_????/Dir_Collapse-STARS_????.block_list" %args["data_directory"])

storage = {}
vvx=[]
vvy=[]
vvz=[]
t_mass=[]
time=[]
t_particles=[]

for ds in ts.piter():

    box = ds.box(left_edge = -ds.domain_width/2.0,right_edge = ds.domain_width/2.0)
    vx = box["star","vx"]
    vy = box["star","vy"]
    vz = box["star","vz"]
    mass = box["star","mass"]
    particles=len(vx)

    vvx.append(float(np.sum(vx)))
    vvy.append(float(np.sum(vy)))
    vvz.append(float(np.sum(vz)))
    time.append(float(ds.current_time))
    t_mass.append(float(np.sum(mass)))
    t_particles.append(particles)
    filename="merge_data_"+Merge_radius+".txt"
    out_file=open(filename,'a')
    print(float(np.sum(vx)),float(np.sum(vy)),float(np.sum(vz)),float(np.sum(mass)),particles,float(ds.current_time),file=out_file)
    out_file.close()

G=6.67e-8
M=t_mass[0]
p=[]
t=[]

#Check if gravity is on                                                                                                              
                                                                                                                                      
if 'gravity' in method_list:
    N=((G*(M**3))/radius)**(1/2)

else:
    N=M*infall_speed

for i in range(len(t_mass)):
    p.append(((vvx[i]**2+vvy[i]**2+vvz[i]**2)**(1/2)*t_mass[i])/N)
    t.append(time[i]/t_collapse)



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
fig.savefig("Mass_momentum_particles_graph.png")

