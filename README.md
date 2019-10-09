# coral-skeletal-growth

## Prerequisites

-Python 3

## Running the code

-Edit the constants in [spheresmaster.py](spheresmaster.py) (see below for documentation)

-'python3 spheresmaster.py'

-the script will generate three files in a 'results' folder in the your current working directory: a .csv file with the output results and two .png files plotting the results, all named with the constants that you chose

## How to use the constants
The code is set up to use the inorganic rate laws for nucleation and bulk mineral growth for aragonite to grow a skeleton on a flat surface.  Nuclei are randomly distributed and each grows at the same rate.  

The inputs to the code are:

Omega - this will impact both rate laws

Max Time - how long do you want the code to run

Time step - how fast do you want to step through time? There is a probability a nucleus is deposited each time step, and that probability is determined by omega.  If this time step is too fast, then you will not deposit as many nuclei as you should because the chance a nucleus is deposited is at least 100%, but if it is too slow it will take longer to run

Optional inputs:

X_length, Y_length, and Z_length - the 3D space that you are growing the skeleton in

max_height - if you want to know how long it takes to grow to a particular height

alpha_multiplier - you can change the ratio of nucleation rate to growth rate by changing alpha in the nucleation rate equation, representing a change in interfacial energy
