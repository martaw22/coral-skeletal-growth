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
#import pickle
#import math
#from mayavi import mlab
import matplotlib.cm
import sys


'''main things to work on:
1. Saving data from each run - save x, y, z, r and time of each nuclei - should be all you need - maybe 
save each plot? - did this
2. Clean up all of the code and annote better
3. Add in any things that Alex and I talked about in terms of outputs

Next things:
1. Call this code from a different code and run it so that the time is not the input - it just runs until 
you tell it to stop because the floor is covered to a certain extent
2. Fix the aragonite growth rate
3. Start working on the part of random distribution of nuclei/clustering'''

'''Things to change or consider every round:
1. Omega
2. Timestep
3. Ground coverage percent
4. Change name for saving plots
5. Maybe size of x,y,z'''

#from matplotlib.animation import FuncAnimation


##assumes union of independent spheres is same as growth of surface made by
#spherulites. Need to make sure this is true Note ths may not be easliy
#adapted to other shapes in the way that a polygon approximation of a
#surface would be. 



#python version of timing
#from: https://stackoverflow.com/questions/7370801/measure-time-elapsed-in-python
start = timer()

##USER DEFINED STUFF
#Define the dimensions of a control volume in microns:
x_length = 100 #µm
y_length = 100 #µm
z_length = 100 #µm
delta_t = 10 #time step in seconds
max_t = 100

#All nuclei start at a specified size = seed size 
seed_radius = 0.5  #µm radius

#Growth only depends on radius and is independent of
#origin. Use rate law Rate = k(Omega-1)^n where k = 11 nmol m-2 s-1 and
#n=1.7 (from Alex's summary figure that he sent me). If Omega is constant then this Growth_Rate is always the same
#It is not clear that this bulk growth rate scales down to this scale
omega = 25
Growth_Rate = ((omega-1)**1.7)*(11*10E-9) #mol m-2 sec-1
Growth_Rate = Growth_Rate/1000/1000/1000/1000    #change units to mol/um/s

#molar volume of aragonite in µm3/mol = MW (g/mol) / density (g/cm3) * 1E12
molarV_arag = 100.09/2.93*1E12

#Nucleation rate is also constant at constant Omega. 
#J = Aexp(-balpha^3 / (ln(omega))^2) in nuclei  µm-2 sec-1
A = 1858506.76
Balpha3 = -19.44553408
#nuclei/m2/s converted to nuclei/um2 for every 10 seconds
J_rate = A * np.exp(Balpha3/np.log(omega)**2)/1000/1000/1000/1000*10

## Simulation
#Start with one nuclei randomly distrubuted on the XY plane. 
NUCLEI = np.array([[ random.random()*x_length,  random.random()*y_length, 0, seed_radius]]); #made a 2D array because can't concatenate it correctly eitherwise - it will only append to the end of the same array, it will not add a new array to the array of arrays
##                       [x                                  y            z     radius]

#recording the time step of each nuclei that is formed
nuclei_timeofdeposition = [0]
time_groundcover = [0]

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
    SA = NUCLEI[:,3]*NUCLEI[:,3] * np.pi*4/2  #surface area of each hemisphere in um2  
    Delta_moles = Growth_Rate*SA*delta_t #mol per sphere
    Delta_volume = Delta_moles*molarV_arag  #um of growth to add
    VOL = (NUCLEI[:,3]*NUCLEI[:,3]*NUCLEI[:,3] * np.pi*4/3)/2 #volume of each hemisphere in um3
    NEW_VOL = VOL + Delta_volume;
    NEW_R = (NEW_VOL*3/4*2/np.pi)**(1/3)  
    NUCLEI[:,3] = NEW_R             #new radius based on growth in this step is input back into the array
  
     
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
        nucleisize = np.array(np.size(NUCLEI, axis=0))
        for row in range(nucleisize):
            height = getElevationonNucleus(x, y, NUCLEI[row,0], NUCLEI[row,1], NUCLEI[row, 3])
            
            if height > highestZSeen:
                highestZSeen = height

        return highestZSeen    

    
#Need to update the z of nuclei that are on top of other nuclei as the nuclei beneath them grow higher
    numberrows = np.size(NUCLEI[:,0])    
    for row in range(numberrows):
        if NUCLEI[row,2] > 0:
            newz = getZElevation(NUCLEI[row,0], NUCLEI[row,1])
            NUCLEI[row,2] = newz
    
    
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
    def Discrete3dSurface(gridsize):
        '''querying every integer point in the x,y,z box to determine where the surface is
        - put in the size of the walls of each grid in the surface as the input'''
        numPointsTested = 0
        surface_points = [[-100,-100,-100]]
        for xx in np.arange(0,x_length+gridsize,gridsize):
            
            for yy in np.arange(0,y_length+gridsize,gridsize):
                
                
                zz = getZElevation(xx,yy)               
                surface_points = np.append(surface_points,[[xx,yy,zz]], axis=0)
                    #if doesPointIntersectAnyNucleus(xx,yy,zz,nuclei_xlayer, nuclei_ylayer, nuclei_zlayer, nuclei_rlayer) == True:
                    
                numPointsTested = numPointsTested + 1
            yy = yy + 1
            xx = xx + 1
            
        surface_points = np.delete(surface_points, 0, 0)
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
    def areaofIrregularQuad(l1,l2,l3,l4):
        '''area of a quadrilateral that is irregular, i.e. not a square/rectangle'''
        '''get all four lengths of the sides
        divide the quad into two triangles with a diagonal down the middle
        find the length of the diagonal using the distance formula between points 1 and 3 of the quad
        find areas of both triangles using heron's formula
        add the areas of the two triangles together'''
        diagonal_len = distanceFormula(surface_points[point,0],surface_points[point,1],surface_points[point,2],surface_points[point+(gridsize_multiplier*x_length+2),0],surface_points[point+(gridsize_multiplier*x_length+2),1],surface_points[point+(gridsize_multiplier*x_length+2),2])
        s1 = (l1+l2+diagonal_len)/2
        s2 = (l3+l4+diagonal_len)/2
        area_tri_1 = np.sqrt(s1*(s1-l1)*(s1-l2)*(s1-diagonal_len))
        area_tri_2 = np.sqrt(s2*(s2-l3)*(s2-l4)*(s2-diagonal_len))
        area_total = area_tri_1 + area_tri_2
        return area_total
    
    
    gridsizeinput_forsurface = .75
    gridsize_multiplier = 1/gridsizeinput_forsurface
    surface_points = Discrete3dSurface(gridsizeinput_forsurface)
    
    #use distance formula to get lenths of all four sides of each box
       
    ##at the end of each row, I don't want it to wrap around to the next row
    #think I solved that issue by calculating the area row by row - from 0 to 99 for each x row
    #calculating the box areas moving from one side of the grid to the other, so calculating the area from 0 to 1 in x lengths, then from 1 to 2, and so on, so from 99 to 100 is the last one calculating and we don't want to calculate from 100 to anything, but simply move on to the next row
    #for this reason, we can just use arange ending in 100 because we don't want to actually use that number or the box that starts with that number
    #this loop provides the area of 10000 boxes in the grid with x lenght and y length equal to 100 in the array all_areas        
    all_areas = []
    for ones in np.arange(0,y_length*gridsize_multiplier):
      
        grid_range = np.array(np.arange(ones*x_length*gridsize_multiplier+ones,x_length*gridsize_multiplier + ones*x_length*gridsize_multiplier+ones))
    
        for point in grid_range:
            l1 = distanceFormula(surface_points[point,0],surface_points[point,1],surface_points[point,2], surface_points[point+1,0], surface_points[point+1,1], surface_points[point+1,2])
            l2 = distanceFormula(surface_points[point+1,0], surface_points[point+1,1], surface_points[point+1,2], surface_points[point+(gridsize_multiplier*x_length+2),0],surface_points[point+(gridsize_multiplier*x_length+2),1],surface_points[point+(gridsize_multiplier*x_length+2),2])
            l3 = distanceFormula(surface_points[point+(gridsize_multiplier*x_length+2),0],surface_points[point+(gridsize_multiplier*x_length+2),1],surface_points[point+(gridsize_multiplier*x_length+2),2], surface_points[point+(gridsize_multiplier*x_length+1),0], surface_points[point+(gridsize_multiplier*x_length+1),1], surface_points[point+(gridsize_multiplier*x_length+1),2])
            l4 = distanceFormula(surface_points[point+(gridsize_multiplier*x_length+1),0], surface_points[point+(gridsize_multiplier*x_length+1),1], surface_points[point+(gridsize_multiplier*x_length+1),2], surface_points[point,0],surface_points[point,1],surface_points[point,2])            
                      
            area_cell = areaofIrregularQuad(l1,l2,l3,l4)
            all_areas = np.append(all_areas,area_cell)
 
    #now want to be able to pick a point within each small cell, once that cell is picked to nucleate within
    #also want to be able to take the sum of all areas, then pick a number randomly from 0 to the sum, then find out which cell that number is referring to
    sum_areas = sum(all_areas)
    #pick an number randomly from 0 through sum_areas
    random_location = random.random()*sum_areas
     
    #go through the cell areas and add them up until you get to the value you generated in random_location
    which_cell_area = 0
    area_count = 0
    for area in all_areas:
            
        if which_cell_area <= random_location:
            which_cell_area += area
            area_count += gridsizeinput_forsurface
    #finding the point1 of the square that is going to have the nucleus in it - point1 is determined by the  lowest 
    #x and y value pair of the four corners
    point1_y = (area_count%(x_length*gridsize_multiplier))/gridsize_multiplier
    point1_x = (area_count - point1_y)/(x_length*gridsize_multiplier)
    
    
    #pick a point randomly in the (gridsize,gridsize) x,y box that has a vertice of point1 - basically making 
    #an assumption that weighting of area across a box this small is negligibly different
    #this is the new nucleation location on a surface! 

    new_weighted_x = random.random()-gridsizeinput_forsurface + point1_x
    new_weighted_y = random.random()-gridsizeinput_forsurface + point1_y    
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
    
    if percentcoverage_firstlayer > 90 and percentcoverage_firstlayer <= 100:
        time_groundcover = np.append(time_groundcover, t)
    
    #pick the percent coverage that you want to seed the ground with - nuclei won't be able to build on top of each other until this coverage percent of first layer is passed
    coverage_percent = 90
   
    
    if Nuclei_flag == 1:  #if flag = 0, don't add the nuclei, but if it's 1, then it's ok to add
        if percentcoverage_firstlayer < coverage_percent:
            nuc_xyz = NEW_NUCLEI[0,0:3]*np.ones((NUCLEI.shape[0],3)) #extract coordinates of the origin #want to get the first three columns of new nuclei in the first row 
            dist_matrix = np.sqrt(np.sum((NUCLEI[:,0:3]-nuc_xyz)**2,axis=1)) #calcualte distnace from every other nuclei  #the first three columns in NUCLEI - the xyz coordinates of the NEWest NUCLEI, squared, then sum up the rows of this new array - so distance between new point and old points squared and then x,y,andz components summed
            if np.all(dist_matrix > NUCLEI[:,3]):   #add new nuclei to the NUCLEI array, but remove if any distance are shorter than the corresponding radius
                NUCLEI = np.append(NUCLEI, [NEW_NUCLEI[0,:]], axis=0)
                print(NUCLEI)
                nuclei_timeofdeposition = np.append(nuclei_timeofdeposition,t)                
                
            else:
                print('New nuclei too close to existing nuclei')                        
        #after seeding of ground, adding the nuclei continuously so that they can be added on top of each other        
        if percentcoverage_firstlayer >= coverage_percent:
            NUCLEI = np.append(NUCLEI, [NEW_NUCLEI[0,:]], axis=0)
            print(NUCLEI)
            nuclei_timeofdeposition = np.append(nuclei_timeofdeposition,t)
        
        
        

    def totalVolume():
        '''find the total volume under the surface
        can do this by dividing each of the squares in the grid into two triangles
        if the points of each square are 1-4, with the top left being 1 and then proceeding clockwise,
        draw the line for the triangle from point 1 to point 3. Then the formula for area of each truncated
        triangular prism is V = A 1/3 (z1 + z2 + z3) where A is area of base of prism'''
        #Area of the base of the triangular prism is half the area of one flat grid cell (z=0), with an area of 0.5*0.5 (grid wall lengths)
        A = gridsizeinput_forsurface**2/2
        Vol_total = 0
        number_squares = 0
        for ones in np.arange(0,y_length*gridsize_multiplier):
       
            grid_range = np.array(np.arange(ones*x_length*gridsize_multiplier+ones,x_length*gridsize_multiplier + ones*x_length*gridsize_multiplier+ones))
    
            for point in grid_range:
                point1 = [surface_points[point,0],surface_points[point,1],surface_points[point,2]]
                point4 = [surface_points[point+(gridsize_multiplier*x_length+1),0], surface_points[point+(gridsize_multiplier*x_length+1),1], surface_points[point+(gridsize_multiplier*x_length+1),2]] 
                point3 = [surface_points[point+(gridsize_multiplier*x_length+2),0],surface_points[point+(gridsize_multiplier*x_length+2),1],surface_points[point+(gridsize_multiplier*x_length+2),2]]
                point2 = [surface_points[point+1,0], surface_points[point+1,1], surface_points[point+1,2]]
                
                V1_tri = (point1[2]+point2[2]+point3[2])*(1/3)*A
                V2_tri = (point1[2]+point3[2]+point4[2])*(1/3)*A
                Vol_square = V1_tri + V2_tri
                Vol_total += Vol_square
                number_squares += 1
        return Vol_total
        
    
    def totalVolume_anotherWay():
        '''find the total volume under the surface by using the signed volume of a tetrahedron
        this link is sort of useful: https://stackoverflow.com/questions/1406029/how-to-calculate-the-volume-of-a-3d-mesh-object-the-surface-of-which-is-made-up
        again split each cell into two triangles
        then you can find the volume of a part of one triangle by finding the volume of all the combinations of smaller tetrahedrons that can be made from the vertices'''
        vol_total = 0        
        for ones in np.arange(0,y_length*gridsize_multiplier):
       
            grid_range = np.array(np.arange(ones*x_length*gridsize_multiplier+ones,x_length*gridsize_multiplier + ones*x_length*gridsize_multiplier+ones))
    
            for point in grid_range:
                point1 = [[surface_points[point,0],surface_points[point,1],surface_points[point,2]]]
                point4 = [[surface_points[point+(gridsize_multiplier*x_length+1),0], surface_points[point+(gridsize_multiplier*x_length+1),1], surface_points[point+(gridsize_multiplier*x_length+1),2]]]
                point3 = [[surface_points[point+(gridsize_multiplier*x_length+2),0],surface_points[point+(gridsize_multiplier*x_length+2),1],surface_points[point+(gridsize_multiplier*x_length+2),2]]]
                point2 = [[surface_points[point+1,0], surface_points[point+1,1], surface_points[point+1,2]]]
                
                tri1 = np.append(point1, point2, axis=0)
                tri1 = np.array(np.append(tri1, point3, axis=0))
                vol_tri1 = (1/6)*np.linalg.det(tri1)

                tri2 = np.append(point3, point4, axis=0)
                tri2 = np.array(np.append(tri2, point1, axis=0))
                vol_tri2 = (1/6)*np.linalg.det(tri2)
                vol_cell = vol_tri1+vol_tri2
                vol_total += vol_cell
        return vol_total
        
   
    
    #function to determine density of nuclei in total skeleton
#    def nucleiDensity():
#        '''we have total growth overall from the above function.  If we say that nuclei deposition is one process
#        that has it's own equation that governs it and that nuclei growth is a different process governed by a 
#        different equation, then we can say that the ratio of nuclei volume with initial radius to overall
#        growth volume - this link is useful https://math.stackexchange.com/questions/2371139/volume-of-truncated-prism'''
#        total_nuclei = np.size(NUCLEI[:,0])
#        #divide by 2 because they're hemispheres
#        total_volume = 4/3*np.pi*seed_radius**3/2*total_nuclei
#        total_growth = totalGrowth()
#        arag_growth = total_growth[0] - total_volume
#        ratio = total_volume/arag_growth
#        return ratio
        
        
    
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
            
    #adding everything as hemispheres - think about if we want to add things as spheres on top of bottom layer
for row in np.nditer(shape):       
        
    r=NUCLEI[row-1,3]                                 ##radius
    #the 10 in u and v just determines the number of lines drawn between points around the edge of the sphere (1/10 of a circle) - a higher number looks more like a circle but takes more time to draw
    u = np.linspace(0, 2*np.pi, 10)     ##theta in spherical coordinates - (how far around the x and y axes in a full circle we go)
    v = np.linspace(0, np.pi/2, 10)    ##phi in spherical coordinates - determines how far along z axis sphere is (so forms a full half circle from fully positive to fully negative z when at full pi units)
    x = r * np.outer(np.cos(u), np.sin(v))+NUCLEI[row-1,0]
    y = r * np.outer(np.sin(u), np.sin(v))+NUCLEI[row-1,1] 
    z = r * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, color='b', alpha=0.5, clip_on = True)    
   
    plt.pause(.001) 
    
    if t == max_t-delta_t:                                                  #save the last timestep image, and this also keeps that plot open
        fig.savefig('/Users/Marta/Documents/Python/nucleation_model_output/plots/modelnucleisphere_om' + str(omega) + '_' + str(max_t) + '.png')
        
                                                                     #close all the other plots from the other timepoints - definitely do not keep them all open!
plt.close()
    
    
    
    
###plotting the surface of the nuclei   
    
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim((0,x_length))
ax.set_ylim((0,y_length))
ax.set_zlim((0,z_length))
    
Z = surface_points[:,2]
walls = wallsofGrid()
walls_front = walls[0]
walls_back = walls[1]
    
Z_front = walls_front[:,2]
Z_back = walls_back[:,2]
    
ax.scatter(walls_back[:,0], walls_back[:,1], walls_back[:,2], c=Z_back, alpha=0.5, clip_on = True, s=20, lw=1)
ax.scatter(walls_front[:,0], walls_front[:,1], walls_front[:,2], c=Z_front, alpha=0.5, clip_on = True, s=20, lw=1)
ax.scatter(surface_points[:,0], surface_points[:,1], surface_points[:,2], c=Z, alpha=0.5, clip_on = True, s=20, lw=1)  
         
#    plt.pause(.001)    
#    if t == max_t-delta_t:                                                  #save the last timestep image, and this also keeps that plot open
fig.savefig('/Users/Marta/Documents/Python/nucleation_model_output/plots/modelnucleisurface_om' + str(omega) + '_' + str(max_t) +  '.png')
        
                                                                       #close all the other plots from the other timepoints - definitely do not keep them all open!
plt.close()    
    
    #plt.show()
####trying to save this 3d image - don't know if there's a file format that will let me, but can do one of two things: save data to a format and have separate file that just plots from saved data, makes it go faster to play with plots, or can use pickling, which creates a form of the plot that's interactive    
#   https://stackoverflow.com/questions/7290370/store-and-reload-matplotlib-pyplot-object    
#   https://stackoverflow.com/questions/29127593/trying-to-write-a-cpickle-object-but-get-a-write-attribute-type-error
#    with open('/Users/Marta/Documents/Python/Plots/modelnucleisphere_om10.png', 'wb') as pickle_file:    
#        pickle.dump(fig, pickle_file.pickle)
#    
#
#      
      

    ##So now we know the percent coverage by nuclei at a given time, and we can find the elevation at a given spot
    ##Now: above a certain percent coverage, next time you generate a spot to nucleate on, you can put it at that spot with the z at the height of the nucleus in the first row        
percentfirstlayer = getSampledPercentageAreaOccupiedByNuclei(NUCLEI[:,0],NUCLEI[:,1],NUCLEI[:,3])*100
    
    
#Output parameters
print('Number of nuclei:', np.size(NUCLEI[:,0]))
total_volume = totalVolume()
total_vol_2 = totalVolume_anotherWay()
print('Total vol way 1:', total_volume)
print('Total vol way 2:', total_vol_2)
print('Total sum of nuclei vol:', sum(VOL))
#density of nuclei on the ground
nuclei_ground_count = 0
for nucleus in range(len(NUCLEI)):
    if NUCLEI[nucleus,2] == 0:
        nuclei_ground_count += 1
nuclei_density = nuclei_ground_count/(x_length*y_length)

print("percent of surface covered by nuclei:", percentfirstlayer)   
#print('percent of second layer covered by nuclei:', percentsecondlayer)
#end of the timer function, at the very end of all the code

#save the following information in a file: NUCLEI, time of deposition, omega, and time step
filename = '/Users/Marta/Documents/Python/nucleation_model_output/'+str(omega)+'_'+str(max_t)+'_nucleationmodelrun'

file = open(filename, 'a')
file.write('\n' + 'Nucleus X, Nucleus Y, Nucleus Z, Nucleus R, Time of Deposition, Timestep' + '\n')
for nucleus in range(len(NUCLEI)):
    file.write(str(NUCLEI[nucleus]) + ', ')
    file.write(str(nuclei_timeofdeposition[nucleus]) + ', ')
    file.write(str(delta_t) + '\n')
file.write('\n' + 'Total Number Nuclei:' +str(np.size(NUCLEI[:,0])) + '\n')    
file.write('\n' + 'Total Calculated Volume:' +str(total_volume) + ' um3' + '\n')
file.write('\n' + 'Actual Volume with Overlaps:' + str(sum(VOL)) + ' um3' + '\n')
file.write('\n' + 'Percent of Ground Covered:' +str(percentfirstlayer) + '%' + '\n')
file.write('\n' + 'Time to Cover Ground:' + str(time_groundcover) + 's' + '\n')
file.write('\n' + 'Total Time:' + str(max_t) + ' s' + '\n')
file.write('\n' + 'Number Nuclei on Ground Level:' + str(nuclei_ground_count) + '\n')
file.write('\n' + 'Nuclei Ground Density:' + str(nuclei_density) + ' nuclei/um2' + '\n')
file.write('\n' + 'Ratio of Nuclei to Growth:' + str(np.size(NUCLEI[:,0])/total_volume) + '\n')

file.close()
    

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
    