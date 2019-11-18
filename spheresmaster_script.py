'''This script deposits and grows nuclei on a surface. It contains all of the constants and the input variables that can be changed.'''

from timeit import default_timer as timer
import random
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import itertools
import skeletal_module as sk
import config

start = timer()
###############################################################################################
##Simulation

#The growth rate is dependent on OMEGA, in units of mol m-2 sec-1
Growth_Rate = ((config.OMEGA-1)**1.7)*(11*1E-9)
#change units to mol/um/s
Growth_Rate = (Growth_Rate/1000/1000/1000/1000)
#Nucleation rate is constant at constant OMEGA
#nuclei/m2/s converted to nuclei/um2
J_rate = config.A * np.exp(config.BALPHA3/np.log(config.OMEGA)**2)/1000/1000/1000/1000
#Start with one nuclei randomly distrubuted on the XY plane.
##[x, y, z, radius]
nuclei = np.array([[random.random()*config.X_LENGTH, random.random()*config.Y_LENGTH, 0, config.SEED_RADIUS]])
#recording the time step of each nuclei that is formed and the time that each percentage of the 2d xy grid is covered
nuclei_timeofdeposition = [0]
ttime = np.arange(0,config.MAX_T,config.DELTA_T)        
for t in ttime:
    print(t)
    ## Grow Spheres
    #Each sphere is described by a radius and an origin. 
    #This information saved in an array nuclei that grows in length with more nuclei.            
    #After each timestep grow each sphere according to the inorganic rate law. 
    sk.growEachNucleus(nuclei, Growth_Rate)           
    #The surface is made up of boxes and the vertices of the boxes are saved in the array surface_points
    surface_points = sk.Discrete3dSurface(config.GRIDSIZEINPUT_FORSURFACE, nuclei)            
    areas = sk.areasofEachCellinGrid(surface_points)                        
    #now want to be able to pick a point within each small cell, once that cell is picked to nucleate within
    #also want to be able to take the sum of all areas, then pick a number randomly from 0 to the sum, then find out which cell that number is referring to
    sum_areas = sum(areas)
    #pick an number randomly from 0 through sum_areas
    random_location = random.random()*sum_areas           
    ## Nucleation        
    #Based on probability of a nucleus being deposited (based on OMEGA and area)
    new_nuclei = None
    Nuclei_flag = 0       
    prob_nuc = config.X_LENGTH*config.Y_LENGTH*J_rate*config.DELTA_T
    if prob_nuc > 1:
        print('ERROR: reduce time step size - nucleation too fast')            
    randomnumber = random.random() 
    new_coordinates = sk.findNewXYZ(areas, random_location, nuclei)
    if randomnumber<=prob_nuc:
        #generate a new nuclei on x,y,z surface
        new_nuclei = np.array([[ new_coordinates[0],  new_coordinates[1], new_coordinates[2], config.SEED_RADIUS]])
        Nuclei_flag = 1              
    #if flag = 0, don't add the nuclei, but if it's 1, then it's ok to add
    if Nuclei_flag == 1:                         
        #adding the nuclei continuously so that they can be added on top of each other        
        nuclei = np.append(nuclei, [new_nuclei[0,:]], axis=0)
        print(nuclei)
        nuclei_timeofdeposition = np.append(nuclei_timeofdeposition,t)              
    
    percentcoverage_firstlayer = sk.getSampledPercentageAreaOccupiedByNuclei(nuclei[:,0],nuclei[:,1],nuclei[:,3])*100
    #Record when ground is covered to certain percentages
    time_groundcoverage = sk.findTimetoCoverGround(percentcoverage_firstlayer)
    volume = sk.totalVolume(surface_points)
    timed_output = sk.outputAtCertainTimes(t, volume, surface_points, nuclei, percentcoverage_firstlayer)    
    
#plotting spheres in 3d  
sk.plot3DSpheres(nuclei, config.OMEGA)
#plotting the surface of the nuclei     
sk.plotSurfaceGrid(nuclei, config.OMEGA, surface_points)

#Output parameters
print('Number of nuclei:', np.size(nuclei[:,0]))


#save the following information in a file: nuclei, time of deposition, OMEGA, and time step
cwd = os.getcwd()
filename = sk.uniqueFileName(str(cwd) + '/results/'+str(config.OMEGA)+'_'+str(config.MAX_T)+ '_' + str(config.SEED_RADIUS) + 'r_' + str(config.ALPHA_MULTIPLIER) + 'alpha', 'txt')

file = open(filename, 'a')
file.write('\n' + 'Nucleus X, Nucleus Y, Nucleus Z, Nucleus R, Time of Deposition, Timestep' + '\n')
time_between_dep = []
for nucleus in range(len(nuclei)):
    file.write(str(nuclei[nucleus]) + ', ')
    time_dep = nuclei_timeofdeposition[nucleus]
    file.write(str(time_dep) + ', ')
    file.write(str(config.DELTA_T) + '\n')
    if nucleus != 0:
        time_between_dep = np.append(time_between_dep, time_dep - nuclei_timeofdeposition[nucleus-1])
file.write('\n' + 'OMEGA:' + str(config.OMEGA) + '\n')
file.write('\n' + 'Total Number Nuclei:' +str(np.size(nuclei[:,0])) + '\n')    
file.write('\n' + 'Total Calculated Volume:' +str(volume) + ' um3' + '\n')
file.write('\n' + 'Actual Volume with Overlaps:' + str(sk.growEachNucleus(nuclei, Growth_Rate)) + ' um3' + '\n')
file.write('\n' + 'Time to Cover Ground 10%:' + str(time_groundcoverage[0]) + 's' + '\n')
file.write('\n' + 'Time to Cover Ground 25%:' + str(time_groundcoverage[1]) + 's' + '\n')
file.write('\n' + 'Time to Cover Ground 50%:' + str(time_groundcoverage[2]) + 's' + '\n')
file.write('\n' + 'Time to Cover Ground 75%:' + str(time_groundcoverage[3]) + 's' + '\n')
file.write('\n' + 'Time to Cover Ground 90%:' + str(time_groundcoverage[4]) + 's' + '\n')
file.write('\n' + 'Total Time:' + str(config.MAX_T) + ' s' + '\n')
file.write('\n' + 'Average Time Step Between Depositions:' + str(np.average(time_between_dep)) + 's' + '\n')
file.write('\n' + 'Number Nuclei on Ground Level:' + str(sk.nucleiGroundDensity(nuclei)[0]) + '\n')
file.write('\n' + 'Nuclei Ground Density:' + str(sk.nucleiGroundDensity(nuclei)[1]) + ' nuclei/um2' + '\n')
file.write('\n' + 'Ratio of Nuclei to Growth:' + str(np.size(nuclei[:,0])/volume) + '\n')
for key in sorted(timed_output):
    file.write('\n' + str(key) + ',' + str(timed_output[key]) + '\n')

file.close()
    
end = timer()
#time in seconds
print('time elapsed:',end - start)
