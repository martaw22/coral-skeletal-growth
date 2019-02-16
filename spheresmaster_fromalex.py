# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 11:29:09 2018

@author: Marta
"""
from timeit import default_timer as timer
import random
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pickle
import math
#from mayavi import mlab
import matplotlib.cm


'''main things to work on:
1. making the surface plots easier to look at
2. figure out what i want to do about edges in the plot - do I want to have them appear on other sides?? boundary conditions'''


#import time
#from matplotlib.animation import FuncAnimation

#SphereMaster - from Alex in matlab - trying to recreate in python

##assumes union of independent spheres is same as growth of surface made by
#spherulites. Need to make sure this is true Note ths may not be easliy
#adapted to other shapes in the way that a polygon approximation of a
#surface would be. 

#This makes it so that the decimal points are not cut off: format long
#I think pretty much same thing as a float
#this measures the time it takes to run the program: tic

#python version of timing
#from: https://stackoverflow.com/questions/7370801/measure-time-elapsed-in-python
start = timer()

##USER DEFINED STUFF
#Define the dimensions of a control volume in microns:
x_length = 100 #µm
y_length = 100 #µm
z_length = 100 #µm
delta_t = 10 #time step in seconds
max_t = 200

#All nuclei start at a specified size = seed size 
seed_radius = 0.5  #µm radius

#Growth only depends on radius and is independent of
#origin. Use rate law Rate = k(Omega-1)^n where k = 11 nmol m-2 s-1 and
#n=1.7. If Omega is constant then this Growth_Rate is always the same
#It is not clear that this bulk growth rate scales down to this scale
omega = 20
Growth_Rate = ((omega-1)**1.7)*11/1e9 #mol m-2 sec-1
Growth_Rate = Growth_Rate*1e8*2

#molar volume of aragonite in µm/mol = MW (g/mol) / density (g/ml) *1e3*1e3
molarV_arag = 100/2.93*1e3*1e3

#Nucleation rate is also constant at constant Omega. 
#J = Aexp(-balpha^3 / (ln(omega))^2) in nuclei  µm-2 sec-1
J_rate = 10/60/60/100/100*100*10

## Simulation
#Start with one nuclei randomly distrubuted on the XY plane. 
NUCLEI = np.array([[ random.random()*x_length,  random.random()*y_length, 0, seed_radius]]); #made a 2D array because can't concatenate it correctly eitherwise - it will only append to the end of the same array, it will not add a new array to the array of arrays
##                       [x                                  y            z     radius]

#Created an array for the second row of nuclei that I will store separately than the first row - had to have values in it to concatenate with new values, but the values are outside of the grid
NUCLEI_2ROW = np.array([[-100,-100,0,0.000000000001]])
NUCLEI_3ROW = np.array([[-100,-100,0,0.000000000001]])



ttime = np.arange(0,max_t,delta_t)

for t in ttime:
## Grow Spheres
#each sphere is described by a radius and an origin. Can keep this
#information in an array NUCLEI that grows in length with more nuclei.

#After each timestep grow each sphere according to the inorganic rate law. 
#Have a function do this in the future and have it operate on the 
#whole NUCLEI array. Growth only depends on radius and is independent of
#origin. Use rate law Rate = k(Omega-1)^n where k = 11 nmol m-2 s-1 and
#n=1.7
    print(t)    
    SA = NUCLEI[:,3]*NUCLEI[:,3] * np.pi*4/1E6/1e6  #surface area of each sphere in m2  
    Delta_moles = Growth_Rate*SA*delta_t #mol per sphere
    Delta_volume = Delta_moles*molarV_arag
    VOL = NUCLEI[:,3]*NUCLEI[:,3]*NUCLEI[:,3] * np.pi*4/3 #surface area of each sphere in µm2
    NEW_VOL = VOL + Delta_volume;
    NEW_R = (NEW_VOL*3/4/np.pi)**(1/3)  
    NUCLEI[:,3] = NEW_R             #new radius based on growth in this step is input back into the array

    #Grow 2nd layer of nuclei at the same rate as the first layer
    SA_2ROW = NUCLEI_2ROW[:,3]*NUCLEI_2ROW[:,3] * np.pi*4/1E6/1e6  #surface area of each sphere in m2  
    Delta_moles_2ROW = Growth_Rate*SA_2ROW*delta_t #mol per sphere
    Delta_volume_2ROW = Delta_moles_2ROW*molarV_arag
    VOL_2ROW = NUCLEI_2ROW[:,3]*NUCLEI_2ROW[:,3]*NUCLEI_2ROW[:,3] * np.pi*4/3 #surface area of each sphere in µm2
    NEW_VOL_2ROW = VOL_2ROW + Delta_volume_2ROW;
    NEW_R_2ROW = (NEW_VOL_2ROW*3/4/np.pi)**(1/3)  
    NUCLEI_2ROW[:,3] = NEW_R_2ROW             #new radius based on growth in this step is input back into the array
    
    #Grow 3nd layer of nuclei at the same rate as the first layer
    SA_3ROW = NUCLEI_3ROW[:,3]*NUCLEI_3ROW[:,3] * np.pi*4/1E6/1e6  #surface area of each sphere in m2  
    Delta_moles_3ROW = Growth_Rate*SA_3ROW*delta_t #mol per sphere
    Delta_volume_3ROW = Delta_moles_3ROW*molarV_arag
    VOL_3ROW = NUCLEI_3ROW[:,3]*NUCLEI_3ROW[:,3]*NUCLEI_3ROW[:,3] * np.pi*4/3 #surface area of each sphere in µm2
    NEW_VOL_3ROW = VOL_3ROW + Delta_volume_3ROW;
    NEW_R_3ROW = (NEW_VOL_3ROW*3/4/np.pi)**(1/3)  
    NUCLEI_3ROW[:,3] = NEW_R_3ROW     
    
    
    #master list of all the nuclei in any row
    NUCLEI_MASTER = np.append(NUCLEI, NUCLEI_2ROW, axis=0)
    NUCLEI_MASTER = np.append(NUCLEI_MASTER, NUCLEI_3ROW, axis=0)
    
#Function for determining if a given point intersects with any nucleus in NUCLEI - right now only in NUCLEI and not in NUCLEI_2ROW    
    def doesPointIntersectAnyNucleus(x,y,z, nuclei_x, nuclei_y, nuclei_z, nuclei_r):
        "does a given x,y point intersect a nucleus?"
        "compute distance from x,y to nucleus center and see if it's less than nucleus radius"    

        distance = np.sqrt((x-nuclei_x)**2 + (y-nuclei_y)**2 + (z-nuclei_z)**2)
        if np.any(distance <= nuclei_r):
            return True
        else:
            return "no intersection"
            
            
    ##These functions are used to find the area of the grid occupied by nuclei    
    
    def getSampledPercentageAreaOccupiedByNuclei(nuclei_xlayer, nuclei_ylayer, nuclei_rlayer):
        "What is the area of the xy plane occupied by nuclei at each timepoint?"
        numIntersecting = 0
        numPointsTested = 0
        for xx in range(x_length):
            
            for yy in range(y_length):
                
                if doesPointIntersectAnyNucleus(xx,yy,0, nuclei_xlayer, nuclei_ylayer, 0, nuclei_rlayer) == True:
                                        
                    numIntersecting = numIntersecting + 1
                numPointsTested = numPointsTested + 1
            yy = yy + 1
            xx = xx + 1
        percentcoverage = numIntersecting/numPointsTested
        return percentcoverage
            


    ##Question: what is the z elevation at a given point of the potential nuclei there?

    def getElevationonNucleus(x,y,x0,y0,r):
        "want to measure the height of a nucleus from any point within the radius"
        "Calculated z based on trig of point on edge of sphere, which is related to radius, x, y, and two angles theta and phi"
        "x and y are the point that we are testing the height of, while x0 and y0 are the center of the sphere"    
        #z = r * np.sqrt((r**2 - np.absolute((y-y0)**2)-np.absolute((x-x0)**2))/(r**2))        
        inside_equation = r**2-(x-x0)**2-(y-y0)**2         
        if inside_equation <= 0:
            z = 0
        else:
            z = np.sqrt(r**2-(x-x0)**2-(y-y0)**2)
        return z
        
           
        
    def getZElevation(x,y):
        "want to know what the elevation of a nucleus is at a given point"  
        "the given point is x,y, and it loops through the NUCLEI matrix, checking to see if it overlaps with any nuclei - if it does, the height at that spot should be >0"
        highestZSeen = 0
        nucleisize = np.array(np.size(NUCLEI_MASTER, axis=0))
        for row in range(nucleisize):
            height = getElevationonNucleus(x, y, NUCLEI_MASTER[row,0], NUCLEI_MASTER[row,1], NUCLEI_MASTER[row, 3])
            
            if height > highestZSeen:
                highestZSeen = height

        return highestZSeen    

    
#Need to update the z of nuclei in the 2nd row as the nuclei beneath them grow higher
    numberrows2 = np.size(NUCLEI_2ROW[:,0])   
    for row in range(numberrows2):
        newz2 = getZElevation(NUCLEI_2ROW[row,0],NUCLEI_2ROW[row,1])
        NUCLEI_2ROW[row,2] = newz2
        
    numberrows3 = np.size(NUCLEI_3ROW[:,0])
    for row in range(numberrows3):
        newz3 = getZElevation(NUCLEI_3ROW[row,0],NUCLEI_3ROW[row,1])
        NUCLEI_3ROW[row,2] = newz3

    ##doing the weighting of areas that are within nuclei to account for surface area
    def distanceFormula(x,y,z,x0,y0,z0):
        '''distance between two points'''        
        distance = np.sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
        return distance
    
    def volumeofOverlappingSphereSections(x,y,z,x0,y0,z0,r0,r1):
        '''calculate the volume of the overlapping parts of two spheres'''
        h0 = r0*(1-(r0**2+distanceFormula(x,y,z,x0,y0,z0)**2-r1**2)/(2*r0*distanceFormula(x,y,z,x0,y0,z0)))
        h1 = r1*(1-(r1**2+distanceFormula(x,y,z,x0,y0,z0)**2-r0**2)/(2*r1*distanceFormula(x,y,z,x0,y0,z0)))         
               
        volume = np.pi()/3*(3*r0*h0**2 - h0**3 + 3*r1*h1**2 - h1**3)
        return volume
    
    #Find which nuclei overlap with others
    #r0 is from nuclei in each, the one being tested against all others, r1 is from nuclei being tested against
    #was initially using this as a way to do weighting of the surface area of spheres, but am not using it currently    
    '''nextrow = 0
    for each in np.size(NUCLEI[:,0]):
        distance_between_nuclei = distanceFormula(NUCLEI[each,0],NUCLEI[each,1],NUCLEI[each,2],NUCLEI[each+1,0],NUCLEI[each+1,1],NUCLEI[each+1,2])
        if distance_between_nuclei < NUCLEI[each,3]:
            volume_overlap = volumeofOverlappingSphereSections(NUCLEI[each,0],NUCLEI[each,1],NUCLEI[each,2],NUCLEI[each+1,0],NUCLEI[each+1,1],NUCLEI[each+1,2],NUCLEI[each,3],NUCLEI[each+1,3])    
            percentvolume_notoverlapping_each = 1 - volume_overlap/VOL[each]
            percentvolume_notoverlapping_nextnucleus = 1 - volume_overlap/VOL[each+1]
            #running into a problem here: this only works for two overlapping spheres, but as soon as you have more than two, you start to discount things multiple times
        nextrow = nextrow+1'''
        
    
    ###Testing making a surface with discrete points in x,y,z space
    #change this so that it finds the z at each point, so that we're not rounding z
    def Discrete3dSurface():
        '''querying every integer point in the x,y,z box to determine where the surface is'''
        numPointsTested = 0
        surface_points = [[-100,-100,-100]]
        for xx in range(x_length+1):
            
            for yy in range(y_length+1):
                
                
                zz = getZElevation(xx,yy)               
                surface_points = np.append(surface_points,[[xx,yy,zz]], axis=0)
                    #if doesPointIntersectAnyNucleus(xx,yy,zz,nuclei_xlayer, nuclei_ylayer, nuclei_zlayer, nuclei_rlayer) == True:
                    
                numPointsTested = numPointsTested + 1
            yy = yy + 1
            xx = xx + 1
            
        
        return surface_points
    
    #making a function to make arrays of the walls of the grid that will grow with the nuclei that touch the edge of the grid
    def wallsofGrid():
        '''if x or y equals 0 or 100, then find z at those points and store those points in an array, so that they can be plotted'''
        wall_points_front = [[-100,-100,-100]]
        wall_points_back = [[-100,-100,-100]]
        for xx in range(x_length+1):
            if xx == 0:
                for yy in range(y_length+1):
                    zz = getZElevation(xx,yy)               
                    wall_points_back = np.append(wall_points_back,[[xx,yy,zz]], axis=0)
                    for zz in np.arange(zz):
                        wall_points_back = np.append(wall_points_back, [[xx,yy,zz]], axis=0)
        for xx in range(x_length+1):
            if xx == 100:
                for yy in range(y_length+1):
                    zz = getZElevation(xx,yy)               
                    wall_points_front = np.append(wall_points_front,[[xx,yy,zz]], axis=0)
                    for zz in np.arange(zz):
                        wall_points_front = np.append(wall_points_front, [[xx,yy,zz]], axis=0)
        for yy in range(y_length+1):
            if yy == 0:
                for xx in range(y_length+1):
                    zz = getZElevation(xx,yy)               
                    wall_points_front = np.append(wall_points_front,[[xx,yy,zz]], axis=0)
                    for zz in np.arange(zz):
                        wall_points_front = np.append(wall_points_front, [[xx,yy,zz]], axis=0)
        for yy in range(y_length+1):
            if yy == 100:
                for xx in range(y_length+1):
                    zz = getZElevation(xx,yy)               
                    wall_points_back = np.append(wall_points_back,[[xx,yy,zz]], axis=0)
                    for zz in np.arange(zz):
                        wall_points_back = np.append(wall_points_back, [[xx,yy,zz]], axis=0)  
        wall_points_front = np.delete(wall_points_front,0,0)
        wall_points_back = np.delete(wall_points_back,0,0)
        return wall_points_front, wall_points_back

    
    
    #after making a discrete grid of surface points, the next step is finding the area of boxes in that grid so that I can weight 
    #those areas for depositing nuclei
    def areaofIrregularQuad(l1,l2,l3,l4,p1x,p1y,p1z,p3x,p3y,p3z):
        '''area of a quadrilateral that is irregular, i.e. not a square/rectangle'''
        '''get all four lengths of the sides
        divide the quad into two triangles with a diagonal down the middle
        find the length of the diagonal using the distance formula between points 1 and 3 of the quad
        find areas of both triangles using heron's formula
        add the areas of the two triangles together'''
        diagonal_len = distanceFormula(p1x,p1y,p1z,p3x,p3y,p3z)
        s1 = (l1+l2+diagonal_len)/2
        s2 = (l3+l4+diagonal_len)/2
        area_tri_1 = np.sqrt(s1*(s1-l1)*(s1-l2)*(s1-diagonal_len))
        area_tri_2 = np.sqrt(s2*(s2-l3)*(s2-l4)*(s2-diagonal_len))
        area_total = area_tri_1 + area_tri_2
        return area_total
    
    
    
    surface_points = Discrete3dSurface()
    #surface_points isn't returning x,y = 0,0 right now - not sure why, so i'm initializing it to have a 0,0 and then keeping that point even though the height isn't right
    #so that means I shouldn't delete the first point anymore    
    surface_points = np.delete(surface_points, 0, 0)
    #use distance formula to get lenths of all four sides of each box
    #delete the very first point because it is the initializing point for the array and is not actually representing a nucleus     
    surface_points_2 = surface_points[0:202]     
    range_surfacepoints_2 = np.size(surface_points_2, axis=0)-102
    ##one thing to be careful of - at the end of each row, I don't want it to wrap around to the next row
    #think I solved that issue by calculating the area row by row - from 0 to 99 for each x row
    #calculating the box areas moving from one side of the grid to the other, so calculating the area from 0 to 1 in x lengths, then from 1 to 2, and so on, so from 99 to 100 is the last one calculating and we don't want to calculate from 100 to anything, but simply move on to the next row
    #for this reason, we can just use arange ending in 100 because we don't want to actually use that number or the box that starts with that number
    #this loop provides the area of 10000 boxes in the grid with x lenght and y length equal to 100 in the array all_areas    
    all_areas = []
    for ones in range(y_length):
       
        grid_range = np.array(np.arange(ones*100+ones,x_length + ones*100+ones))
    
        for point in grid_range:
            l1 = distanceFormula(surface_points[point,0],surface_points[point,1],surface_points[point,2],surface_points[point+1,0], surface_points[point+1,1], surface_points[point+1,2])
            l2 = distanceFormula(surface_points[point+1,0], surface_points[point+1,1], surface_points[point+1,2], surface_points[point+102,0],surface_points[point+102,1],surface_points[point+102,2])
            l3 = distanceFormula(surface_points[point+102,0],surface_points[point+102,1],surface_points[point+102,2],surface_points[point+101,0], surface_points[point+101,1], surface_points[point+101,2])
            l4 = distanceFormula(surface_points[point+101,0], surface_points[point+101,1], surface_points[point+101,2], surface_points[point,0],surface_points[point,1],surface_points[point,2])
            area_cell = areaofIrregularQuad(l1,l2,l3,l4,surface_points[point,0],surface_points[point,1],surface_points[point,2],surface_points[point+102,0],surface_points[point+102,1],surface_points[point+102,2])
            all_areas = np.append(all_areas,area_cell)
 
    #now want to be able to pick a point within each small cell, once that cell is picked to nucleate within
    #also want to be able to take the sum of all areas, then pick a number randomly from 0 to the sum, then find out which cell that number is referring to
    sum_areas = sum(all_areas)
    #pick an number randomly from 0 through sum_areas
    random_location = random.random()*sum_areas
    #go through the cell areas, adding them up until you get to the value you just generated in random_location
     
    #go through the cell areas and add them up until you get to the value you generated in random_location
    which_cell_area = 0
    area_count = 0
    for area in all_areas:
            
        if which_cell_area <= random_location:
            which_cell_area += area
            area_count += 1
    #finding the point1 of the square that is going to have the nucleus in it - point1 is determined by the  lowest x and y value pair of the four corners
    point1_y = area_count%x_length
    point1_x = (area_count - point1_y)/x_length
    
    
    #pick a point randomly in the (1,1) x,y box that has a vertice of point1 - basically making an assumption that weighting of area across a box this small is negligibly different
    #this is the new nucleation location on a surface!    
    new_weighted_x = random.random() + point1_x
    new_weighted_y = random.random() + point1_y    
    new_weighted_z = getZElevation(new_weighted_x, new_weighted_y)    
            
    
    ## Nucleation
#dont forget to nucleate on XY plane too (then check later if within the volume
#of a sphere)

#get a likleyhood of nucleation per SA per time step. 
#Need to select from the right probability distribution -> think about
#this! For now, assume in steady-state regiem so a uniform distribution
#seems reasonable
    NEW_NUCLEI = None
    Nuclei_flag = 0

    prob_nuc = x_length*y_length*J_rate*delta_t
    print(prob_nuc)
    if prob_nuc > 1:
        print('ERROR: reduce time step size - nucleation too fast')
    
    randomnumber = random.random()
    
    #this is the way to add nuclei previously, with some incorrect weighting and mostly ignoring the z factor
      
    '''print(randomnumber)
    new_x = random.random()*x_length
    new_y = random.random()*y_length
    
    if doesPointIntersectAnyNucleus(new_x,new_y, 0, NUCLEI[:,0],NUCLEI[:,1],0,NUCLEI[:,3])==True:
        weightedprob_nuc = prob_nuc*2
        if randomnumber<=weightedprob_nuc:
            NEW_NUCLEI = np.array([[ new_x,  new_y, 0, seed_radius]]);  #put an extra set of brackets around this to make it a 2D array (even though it only has 1 row and is a 1d array)
            Nuclei_flag = 1
    
    elif randomnumber<=prob_nuc:     
        #generate an new Nuclei on XY surface
        NEW_NUCLEI = np.array([[ new_x,  new_y, 0, seed_radius]]);  #put an extra set of brackets around this to make it a 2D array (even though it only has 1 row and is a 1d array)
        Nuclei_flag = 1        '''                                                                              #did this so that I could in the next if section call the number of rows of this array, even if the number is 1 (which I think it might always be...) because if you have a 1D array, then the number of rows is not callable because it's just 1

    if randomnumber<=prob_nuc:
        #generate a new nuclei on x,y,z surface
        NEW_NUCLEI = np.array([[ new_weighted_x,  new_weighted_y, new_weighted_z, seed_radius]]);  #put an extra set of brackets around this to make it a 2D array (even though it only has 1 row and is a 1d array)
        Nuclei_flag = 1  
    
    
    #Nucleate on surface of sphere based on whole surface area of sphere and do
    #so randomly on sphere

    #for periodic boundary conditions, clone every nuclei into adjacent copies
    #of the control volume (adjacent in X&Y, not Z) 
    
    #Check if Nuclei are within the volume of other spheres. Just check
    #distance from center of each other sphere and compare to radius
    percentcoverage_firstlayer = getSampledPercentageAreaOccupiedByNuclei(NUCLEI[:,0],NUCLEI[:,1],NUCLEI[:,3])*100
    percentcoverage_secondlayer = getSampledPercentageAreaOccupiedByNuclei(NUCLEI_2ROW[:,0],NUCLEI_2ROW[:,1],NUCLEI_2ROW[:,3])*100    
    
    
    if Nuclei_flag == 1:  #added the =1 part, assume this is what was meant, because if flag = 0, then an error occurs and don't add the nuclei, but if it's 1, then it's ok to add
        #v = np.arange(1,NEW_NUCLEI.shape[0]+1)#NN_shape = NEW_NUCLEI.shape
        if percentcoverage_firstlayer < 10:
             #it doesn't seem to remember that i made NEW_NUCLEI a 2D array up above - so it doesn't think you can call the number of rows from it if it's a 1D array
            nuc_xyz = NEW_NUCLEI[0,0:3]*np.ones((NUCLEI.shape[0],3)) #extract coordinates of the origin #want to get the first three columns of new nuclei in the first row - again i think there should only be one row always 
            dist_matrix = np.sqrt(np.sum((NUCLEI[:,0:3]-nuc_xyz)**2,axis=1)) #calcualte distnace from every other nuclei  #the first three columns in NUCLEI - the xyz coordinates of the NEWest NUCLEI, squared, then sum up the rows of this new array - so distance between new point and old points squared and then x,y,andz components summed
            if np.all(dist_matrix > NUCLEI[:,3]):   #add new nuclei to the NUCLEI array, but remove if any distance are shorter than the corresponding radius
                NUCLEI = np.append(NUCLEI, [NEW_NUCLEI[0,:]], axis=0)
                print(NUCLEI)
        
            else:
                print('New nuclei too close to existing nuclei')
        '''##Once coverage on the grid is over a certain percent, new nuclei can be placed on top of old nuclei at the height of the nuclei - however, they still can't be placed on top of each other in this layer, and they will be saved in their own matrix
        if percentcoverage_firstlayer >= 10 and percentcoverage_secondlayer <10:
             #it doesn't seem to remember that i made NEW_NUCLEI a 2D array up above - so it doesn't think you can call the number of rows from it if it's a 1D array         
            #determine height to be the top of the nucleus already present - the radius of the new sphere, because we want the middle of the sphere to be at the point            
            height = getZElevation(NEW_NUCLEI[0,0], NEW_NUCLEI[0,1])
            NEW_NUCLEI[0,2] = height
            nuc2_xyz = NEW_NUCLEI[0,0:3]*np.ones((NUCLEI_2ROW.shape[0],3)) #extract coordinates of the origin #want to get the first three columns of new nuclei in the first row - again i think there should only be one row always 
            dist_matrix2 = np.sqrt(np.sum((NUCLEI_2ROW[:,0:3]-nuc2_xyz)**2,axis=1)) #calcualte distnace from every other nuclei  #the first three columns in NUCLEI - the xyz coordinates of the NEWest NUCLEI, squared, then sum up the rows of this new array - so distance between new point and old points squared and then x,y,andz components summed
            if np.all(dist_matrix2 > NUCLEI_2ROW[:,3]):   #add new nuclei to the NUCLEI array, but remove if any distance are shorter than the corresponding radius
                NUCLEI_2ROW = np.append(NUCLEI_2ROW, [NEW_NUCLEI[0,:]], axis=0) #I think the brackets around NEW_NUCLEI are needed when NEW_NUCLEI is set to none for each loop
                print(NUCLEI_2ROW)
   
        
        if percentcoverage_secondlayer >= 10:
            height = getZElevation(NEW_NUCLEI[0,0], NEW_NUCLEI[0,1])
            NEW_NUCLEI[0,2] = height
            nuc2_xyz = NEW_NUCLEI[0,0:3]*np.ones((NUCLEI_3ROW.shape[0],3)) #extract coordinates of the origin #want to get the first three columns of new nuclei in the first row - again i think there should only be one row always 
            dist_matrix2 = np.sqrt(np.sum((NUCLEI_3ROW[:,0:3]-nuc2_xyz)**2,axis=1)) #calcualte distnace from every other nuclei  #the first three columns in NUCLEI - the xyz coordinates of the NEWest NUCLEI, squared, then sum up the rows of this new array - so distance between new point and old points squared and then x,y,andz components summed
            if np.all(dist_matrix2 > NUCLEI_3ROW[:,3]):   #add new nuclei to the NUCLEI array, but remove if any distance are shorter than the corresponding radius
                NUCLEI_3ROW = np.append(NUCLEI_3ROW, [NEW_NUCLEI[0,:]], axis=0) #I think the brackets around NEW_NUCLEI are needed when NEW_NUCLEI is set to none for each loop
                print(NUCLEI_3ROW)'''
        
        #doing things without layers, just continuously adding nuclei
        #adding an initial layer where nuclei cannot be added on top of each other (that's the if statement above), then just adding the nuclei continuously so that they can be added on top of each other        
        if percentcoverage_firstlayer >= 10:
            NUCLEI = np.append(NUCLEI, [NEW_NUCLEI[0,:]], axis=0)
            print(NUCLEI)
        
        
        
        
#Generate new sphere objects for each new nuclei

#calculating volume of union of a family of spheres could be complicated.
#Could do this seqentially with some recursion, get volume of first sphere,
#then check if second sphere overlaps, if so get only the non overlapping
#volume, then check if third sphere overlaps with first sphere or second
#sphere, only get region not overlapping either, etc. Is there  numerical
#integration we can do based on the union of the equations? Counting when IF
#ANY?


##plotting circles: https://stackoverflow.com/questions/9215658/plot-a-circle-with-pyplot
    

    '''fig, ax = plt.subplots()
    ax.set_xlim((0,x_length))
    ax.set_ylim((0,y_length))
    plt.hold(True)                          #keeps previous plot open to addd to - not sure if this and the next step are actually necessary or not
    plt.cla()                               #clears axis of the plot that already exists
   
    
    shape = np.arange(0,NUCLEI.shape[0])    ##making an array from 0 to the number of rows of NUCLEI, increasing by 1
    for row in np.nditer(shape):            #for each nucleus, or 'row' in the in the NUCLEUS array, plot as a circle
        circle_row = plt.Circle((NUCLEI[row-1,0:2]), NUCLEI[row-1,3], color='r', fill=False)
        ax = plt.gca()
        ax.add_artist(circle_row)     
   
    plt.pause(.001)                         #this pause step enables the program to plot in real time
    plt.show()

    if t == max_t-delta_t:                                                  #save the last timestep image, and this also keeps that plot open
        fig.savefig('/Users/Marta/Documents/Python/Plots/modelnuclei_om10.png')
    else:                                                                   #close all the other plots from the other timepoints - definitely do not keep them all open!
        plt.close()'''
    

###plotting spheres in 3d
#
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim((0,x_length))
    ax.set_ylim((0,y_length))
    ax.set_zlim((0,z_length))
    # Make data
    #useful website for understanding spherical coordinates: https://mathinsight.org/spherical_coordinates
    #so radius determines the size of the spheres, and the angles determine how much of a sphere they are (ie half, full, other shapes, etc)
    shape = np.arange(0,NUCLEI.shape[0])    ##making an array from 0 to the number of rows of NUCLEI, increasing by 1
    #shape2 = np.arange(0,NUCLEI_2ROW.shape[0])     
    #shape3 = np.arange(0,NUCLEI_3ROW.shape[0])           
            
    '''for row2 in np.nditer(shape2):   
            r2 = NUCLEI_2ROW[row2-1,3]
            u2 = np.linspace(0,2*np.pi,10)
            v2 = np.linspace(0,np.pi,10)
            x2 = r2 * np.outer(np.cos(u2),np.sin(v2))+NUCLEI_2ROW[row2-1,0]
            y2 = r2*np.outer(np.sin(u2),np.sin(v2))+NUCLEI_2ROW[row2-1,1]
            z2 = r2*np.outer(np.ones(np.size(u2)),np.cos(v2))+NUCLEI_2ROW[row2-1,2]
            
        # Plot the surface    
            
            ax.plot_surface(x2, y2, z2, color='r', alpha=0.5, clip_on = True)'''
    #adding everything as hemispheres - think about if we want to add things as spheres on top of bottom layer
    for row in np.nditer(shape):       
        
            r=NUCLEI[row-1,3]                                 ##radius
    #the 100 in u and v just make it so that each line drawn between points around the edge of the sphere is very small (1/100 of a circle) so that it looks like a smooth circle at the end, can make this number smaller to get a circle that has edges
            u = np.linspace(0, 2*np.pi, 10)     ##theta in spherical coordinates - (how far around the x and y axes in a full circle we go)
            v = np.linspace(0, np.pi/2, 10)    ##phi in spherical coordinates - determines how far along z axis sphere is (so forms a full half circle from fully positive to fully negative z when at full pi units)
            x = r * np.outer(np.cos(u), np.sin(v))+NUCLEI[row-1,0]
            y = r * np.outer(np.sin(u), np.sin(v))+NUCLEI[row-1,1] 
            z = r * np.outer(np.ones(np.size(u)), np.cos(v))
            ax.plot_surface(x, y, z, color='b', alpha=0.5, clip_on = True)    
   
    '''for row3 in np.nditer(shape3):       
        
            r3 = NUCLEI_3ROW[row3-1,3]
            u3 = np.linspace(0,2*np.pi,10)
            v3 = np.linspace(0,np.pi,10)
            x3 = r3 * np.outer(np.cos(u3),np.sin(v3))+NUCLEI_3ROW[row3-1,0]
            y3 = r3 * np.outer(np.sin(u3),np.sin(v3))+NUCLEI_3ROW[row3-1,1]
            z3 = r3 * np.outer(np.ones(np.size(u3)),np.cos(v3))+NUCLEI_3ROW[row3-1,2]
            ax.plot_surface(x3, y3, z3, color='g', alpha=0.5, clip_on = True)'''
   
    plt.pause(.001) 
    
    if t == max_t-delta_t:                                                  #save the last timestep image, and this also keeps that plot open
        fig.savefig('/Users/Marta/Documents/Python/Plots/modelnucleisphere_om10_3layers.png')
        
    else:                                                                   #close all the other plots from the other timepoints - definitely do not keep them all open!
        plt.close()
    
    plt.show()
    
    
###plotting the surface of the nuclei    
         
    plt.pause(.001)    
    if t == max_t-delta_t:                                                  #save the last timestep image, and this also keeps that plot open
        fig.savefig('/Users/Marta/Documents/Python/Plots/modelnucleisurface_om10_3layers.png')
        
    else:                                                                   #close all the other plots from the other timepoints - definitely do not keep them all open!
        plt.close()    
    
    plt.show()
####trying to save this 3d image - don't know if there's a file format that will let me, but can do one of two things: save data to a format and have separate file that just plots from saved data, makes it go faster to play with plots, or can use pickling, which creates a form of the plot that's interactive    
#   https://stackoverflow.com/questions/7290370/store-and-reload-matplotlib-pyplot-object    
#   https://stackoverflow.com/questions/29127593/trying-to-write-a-cpickle-object-but-get-a-write-attribute-type-error
#    with open('/Users/Marta/Documents/Python/Plots/modelnucleisphere_om10.png', 'wb') as pickle_file:    
#        pickle.dump(fig, pickle_file.pickle)
#    
#    engine = Engine()
#    engine.start()
#    fig = mlab.figure(figure=None, engine=engine)
#    contour3d = mlab.contour3d(surface_points[:,0], surface_points[:,1], surface_points[:,2],figure=fig)
#
#      
      

    ##So now we know the percent coverage by nuclei at a given time, and we can find the elevation at a given spot
    ##Now: above a certain percent coverage, next time you generate a spot to nucleate on, you can put it at that spot with the z at the height of the nucleus in the first row        
    percentfirstlayer = getSampledPercentageAreaOccupiedByNuclei(NUCLEI[:,0],NUCLEI[:,1],NUCLEI[:,3])*100
    percentsecondlayer = getSampledPercentageAreaOccupiedByNuclei(NUCLEI_2ROW[:,0],NUCLEI_2ROW[:,1],NUCLEI_2ROW[:,3])*100
    
    
    
print("percent of first layer covered by nuclei:", percentfirstlayer)   
print('percent of second layer covered by nuclei:', percentsecondlayer)
#end of the timer function, at the very end of all the code
end = timer()
#time in seconds
print('time elapsed:',end - start)

'''#make a function to find all of the points that are within any circles present, excluding overlapping area

    
    def getcircleareaoverlap():
        "What is the area of the circles present minus any overlaps?"
        circlepoints = [[0,0]]        
        numIntersecting = 0
        numPointsTested = 0
        for xx in range(x_length):
            
            for yy in range(y_length):
                
                if doesPointIntersectAnyNucleus(xx,yy) == True:
                    numIntersecting = numIntersecting + 1
                    if np.any([[xx,yy]]==circlepoints):
                    #if np.any(x == circlepoints[:,0] and y == circlepoints[:,1]):
                    #if isPointMemberOfArrayOfPoints(xx,yy,circlepoints):
                        print('overlapping area')
                    else:
                        circlepoints = np.append(circlepoints,[[xx,yy]], axis=0)               
                        print(circlepoints)
                numPointsTested = numPointsTested + 1
            yy = yy + 1
            xx = xx + 1
        
        return circlepoints
        ##right now this code seems to be returning only points that do not overlap... - so that might be useful actually because then I will know what percentage of points
        ##in circles overlaps other points in circles, and I will know to only weight them once instead of twice'''
    