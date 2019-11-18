'''Contains global variables.'''

#Define the dimensions of a control volume in microns(µm):
X_LENGTH = 100 
Y_LENGTH = 100 
Z_LENGTH = 100 

#The time step in seconds
DELTA_T = 100 

#The amount of time that the simulation will take
MAX_T = 1000

#MAX_HEIGHT is chosen for each OMEGA
#It is the bar to reach for vertical extension
#A height that is high enough to not be influence by the first layer of nuclei
MAX_HEIGHT = 15

#All nuclei start at a specified size = seed size, determined by the radius in µm
SEED_RADIUS = 0.005

#OMEGA is the primary driver of growth rate and nucleation rate
OMEGA = 10

#molar volume of aragonite in µm3/mol = MW (g/mol) / density (g/cm3) * 1E12
MOLARV_ARAG = 100.09/2.93*1E12

#These are the constants in the nucleation rate equation, determined experimentally
#J = Aexp(-balpha^3 / (ln(OMEGA))^2) in nuclei  µm-2 sec-1
A = 1858506.76
ALPHA_MULTIPLIER = 1
BALPHA3 = -19.44553408

#The discrete surface grid is made up of boxes with sides of length GRIDSIZEINPUT_FORSURFACE
#The GRIDSIZE_MULTIPLIER is a ratio that makes up for the sides not being 1
#Have to pick a multiple of 0.5 for the code to work
GRIDSIZEINPUT_FORSURFACE = 1
GRIDSIZE_MULTIPLIER = 1/GRIDSIZEINPUT_FORSURFACE
