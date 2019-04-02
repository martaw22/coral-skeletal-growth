# coral-skeletal-growth

##How to use the code
The code is set up to use the inorganic rate laws for nucleation and bulk mineral growth for aragonite to grow a skeleton on a flat surface.  Nuclei are randomly distributed and each grows at the same rate.  

The inputs to the code are:

Omega - this will impact both rates laws
Max Time - how long do you want the code to run
Timestep - how fast do you want to step through time? - there is a chance a nucleus is deposited each timestep, and that chance is determined by omega - if this timestep is too fast, then you will not deposit as many nuclei as you should because the chance a nucleus is deposited as at least 100% - but too slow and it will take longer to run

Optional inputs:
X_length, Y_length, and Z_length - the 3D space that you are growing the skeleton in
max_height - if you want to know how long it takes to grow to a paricular height






##Running this code on the computer clusters
Use ssh to connect to the clusters, called Mox.

ssh -X martaw@mox.hyak.uw.edu

You are automatically in your home directory.  I have created a folder there called coral-skeletal-growth that contains my code (spheremaster) and contains a folder called results.  

The version of the code present currently is from the branch cluster-version.

The default anaconda3 environment in the computer clusters is not compatible with my code - specifically, the version of matplotlib present makes the code crash immediately.  It is also not possible to update or change any of the versions of packages for the clusters without permission.

To be able to run the code, create a local, permanent environment in your home directory with python3.  To do this, first load anaconda in Mox.  

To find which versions of anaconda are present:
module avail anaconda

To load the version of anaconda you want:
module load anaconda3_3.5
-make sure that you are loading anaconda3

To create a local virtual environment:
conda create --name myenv pip

Then load the already created local environment:
source activate martaenv

Notes on the local environment:
-it does not have anaconda, but it does have python3
-you can install whatever version of packages you want using:
pip install numpy (for example)
- to exit:
source deactivate
-to delete:
remove --martaenv --all

Once you are in the local environment, then you can enter the coenv nodes and run a job.
-do not run a job on Mox - it is for scheduling only
Enter the coenv nodes using:
srun -p coenv -A coenv --time=2:00:00 --mem=20G --pty /bin/bash

-the time needs to be set so you have enough of it to run your job - it will automatically exit when you run out of time
-the memory can be set as well, but there's a limit of ~120G, don't use more than you need

To check what's in the queue for the coenv nodes:
squeue -p coenv

There are two nodes.

Once inside the coenv nodes, you can run python3 as you normally would.