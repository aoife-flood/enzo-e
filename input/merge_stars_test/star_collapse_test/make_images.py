#65;6003;1cimport numpy as np
import yt
import matplotlib.pyplot as plt
import cello_parse as cp
import io
import argparse as ap
import matplotlib
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

collapse_centre = params.Initial.collapse_stars.centre

yt.enable_parallelism()
ts =  yt.DatasetSeries("%s/Dir_Collapse-STARS_????/Dir_Collapse-STARS_????.block_list" %args["data_directory"])

    
momentum=[]
storage = {}
for ds in ts.piter():

    box = ds.box(left_edge = -ds.domain_width/2.0,right_edge = ds.domain_width/2.0)
    x = box["star","x"]
    y = box["star","y"]
    z = box["star","z"]
    
    fig,ax = plt.subplots()
    ax.plot(x,y,marker = "x",linestyle = "None")
    ax.set_xlim(-ds.domain_width[0]/2.0,ds.domain_width[0]/2.0)
    ax.set_ylim(-ds.domain_width[0]/2.0,ds.domain_width[0]/2.0)
    ax.set_xlabel("x (cm)")
    ax.set_ylabel("y (cm)")
    collapse_time=4.89e16
    time=float(ds.current_time)/collapse_time
    ax.set_title(time)
    fig.savefig("image_%d.png" %ds.parameters["current_cycle"])
