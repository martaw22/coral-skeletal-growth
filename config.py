'''Contains global variables.'''

#Define the dimensions of a control volume in microns:
X_LENGTH = 100 #µm
Y_LENGTH = 100 #µm
Z_LENGTH = 100 #µm
DELTA_T = 100 #time step in seconds

MAX_T = 1000

#MAX_HEIGHT is chosen for each OMEGA
#It is the bar to reach for vertical extension
#A height that is high enough to not be influence by the first layer of nuclei
MAX_HEIGHT = 15

#All nuclei start at a specified size = seed size
SEED_RADIUS = 0.005  #µm radius

#OMEGA is the primary driver of growth rate and nucleation rate
OMEGA = 10

#molar volume of aragonite in µm3/mol = MW (g/mol) / density (g/cm3) * 1E12
MOLARV_ARAG = 100.09/2.93*1E12

#Nucleation rate is also constant at constant OMEGA
#J = Aexp(-balpha^3 / (ln(OMEGA))^2) in nuclei  µm-2 sec-1
A = 1858506.76
ALPHA_MULTIPLIER = 1
BALPHA3 = -19.44553408

#The discrete surface grid is made up of boxes with sides of length GRIDSIZEINPUT_FORSURFACE
#The GRIDSIZE_MULTIPLIER is a ratio that makes up for the sides not being 1
#Right now - have to pick a multiple of 0.5 for the code to work
GRIDSIZEINPUT_FORSURFACE = 1
GRIDSIZE_MULTIPLIER = 1/GRIDSIZEINPUT_FORSURFACE
