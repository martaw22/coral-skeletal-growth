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
import os
import itertools


'''Next things:
1. Call this code from a different code and run it so that the time is not the input - it just runs until 
you tell it to stop because the floor is covered to a certain extent

3. Start working on the part of random distribution of nuclei/clustering
4. Think about shape of nuclei
5. Check and make sure that you can't ever make a negative deposition'''

#Function for determining if a given point intersects with any nucleus in nuclei    
def doesPointIntersectAnyNucleus(x,y,z, nuclei_x, nuclei_y, nuclei_z, nuclei_r):
    '''does a given x,y point intersect a nucleus?
    compute distance from x,y to nucleus center and see if it's less than nucleus radius'''    

    distance = np.sqrt((x-nuclei_x)**2 + (y-nuclei_y)**2 + (z-nuclei_z)**2)
    if np.any(distance <= nuclei_r):
        return True
    else:
        return "no intersection"
    
    
##These functions are used to find the area of the 2D xy grid occupied by nuclei    

def getSampledPercentageAreaOccupiedByNuclei(nuclei_xlayer, nuclei_ylayer, nuclei_rlayer):
    "What is the area of the xy plane occupied by nuclei at each timepoint?"
    numIntersecting = 0
    numPointsTested = 0
    for xx in range(X_LENGTH+1):

        for yy in range(Y_LENGTH+1):

            if doesPointIntersectAnyNucleus(xx,yy,0, nuclei_xlayer, nuclei_ylayer, 0, nuclei_rlayer) == True:
                                    
                numIntersecting = numIntersecting + 1
            numPointsTested = numPointsTested + 1

    percentcoverage = numIntersecting/numPointsTested
    return percentcoverage

##Question: what is the z elevation at a given point of the potential nuclei there?

def getElevationonNucleus(x,y,x0,y0,z0,r):
    "want to measure the height of a nucleus from any point within the radius"
    "Calculated z based on trig of point on edge of sphere, which is related to radius, x, y, and two angles theta and phi"
    "x and y are the point that we are testing the height of, while x0 and y0 are the center of the sphere"   
    '''This checks the relative height of the nucleus, meaning it adds the z value to the radius, essentially, to get the overall height'''       
    inside_equation = r**2-(x-x0)**2-(y-y0)**2             
    if inside_equation <= 0:
        z = 0
    else:
        if z0 == 0:   
            z = np.sqrt(r**2-(x-x0)**2-(y-y0)**2)
        else:
            z = np.sqrt(r**2-(x-x0)**2-(y-y0)**2) + z0
    return z

def getElevationUnderNucleus(x,y,x0,y0,r):
    "want to measure the height of a nucleus from any point within the radius"
    "Calculated z based on trig of point on edge of sphere, which is related to radius, x, y, and two angles theta and phi"
    "x and y are the point that we are testing the height of, while x0 and y0 are the center of the sphere"   
    '''This checks the relative height of the nucleus, meaning it adds the z value to the radius, essentially, to get the overall height'''       
    inside_equation = r**2-(x-x0)**2-(y-y0)**2             
    if inside_equation <= 0:
        z = 0
    else:
        z = np.sqrt(r**2-(x-x0)**2-(y-y0)**2)
    return z

    
def getZElevation(x,y):
    "want to know what the elevation of a nucleus is at a given point"  
    "the given point is x,y, and it loops through the nuclei matrix, checking to see if it overlaps with any nuclei - if it does, the height at that spot should be >0"
    highestZSeen = 0
    nucleisize = np.array(np.size(nuclei, axis=0))
    nuclei_x = nuclei[:,0]
    nuclei_y = nuclei[:,1]
    nuclei_z = nuclei[:,2]
    nuclei_r = nuclei[:,3]
    for row in range(nucleisize):       
        height = getElevationonNucleus(x, y, nuclei_x[row], nuclei_y[row], nuclei_z[row], nuclei_r[row])
        if height > highestZSeen:
            highestZSeen = height
    return highestZSeen     

def getZElevationUnderNucleus(x,y):
    "want to know what the elevation of a nucleus is at a given point"  
    "the given point is x,y, and it loops through the nuclei matrix, checking to see if it overlaps with any nuclei - if it does, the height at that spot should be >0"
    highestZSeen = 0
    nucleisize = np.array(np.size(nuclei, axis=0))
    nuclei_x = nuclei[:,0]
    nuclei_y = nuclei[:,1]    
    nuclei_r = nuclei[:,3]
    for row in range(nucleisize):       
        height = getElevationUnderNucleus(x, y, nuclei_x[row], nuclei_y[row], nuclei_r[row])
        if height > highestZSeen:
            highestZSeen = height
    return highestZSeen        

##doing the weighting of areas that are within nuclei to account for surface area
def distanceFormula(x,y,z,x0,y0,z0):
    '''distance between two points'''        
    distance = np.sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
    return distance


##Making a surface with discrete points in x,y,z space
def Discrete3dSurface(gridsize):
    '''querying every integer point in the x,y,z box to determine where the surface is
    - put in the size of the walls of each grid in the surface as the input'''
    numPointsTested = 0    
    surface_points = np.zeros(((X_LENGTH/gridsize + 1) * (Y_LENGTH/gridsize + 1), 3))
    for xx in np.arange(0,X_LENGTH+gridsize,gridsize):
        for yy in np.arange(0,Y_LENGTH+gridsize,gridsize):
            zz = getZElevation(xx,yy) 
            surface_points[numPointsTested][0] = xx
            surface_points[numPointsTested][1] = yy
            surface_points[numPointsTested][2] = zz
            numPointsTested = numPointsTested + 1    
    return surface_points

#making a function to make arrays of the walls of the grid that will grow with the nuclei that touch the edge of the grid
def wallsofGrid(gridsize):
    '''if x or y equals 0 or 100, then find z at those points and store those points in an array, so that they can be plotted'''
    wall_points_front = [[-100,-100,-100]]
    wall_points_back = [[-100,-100,-100]]
    
    xx_back = 0
    for yy in np.arange(0,Y_LENGTH+gridsize,gridsize):
        zz = getZElevation(xx_back,yy)               
        wall_points_back = np.append(wall_points_back,[[xx_back,yy,zz]], axis=0)
        for possible_zz in np.arange(0,zz,0.1):
            wall_points_back = np.append(wall_points_back, [[xx_back,yy,possible_zz]], axis=0)
    possible_x = np.arange(0,X_LENGTH+gridsize,gridsize)
    xx_front = possible_x[-1]
    for yy in np.arange(0,Y_LENGTH+gridsize,gridsize):
        zz = getZElevation(xx_front,yy)               
        wall_points_front = np.append(wall_points_front,[[xx_front,yy,zz]], axis=0)
        for possible_zz in np.arange(0,zz,0.1):
            wall_points_front = np.append(wall_points_front, [[xx_front,yy,possible_zz]], axis=0)
    
    yy_front = 0
    for xx in np.arange(0,X_LENGTH+gridsize,gridsize):
        zz = getZElevation(xx,yy_front)               
        wall_points_front = np.append(wall_points_front,[[xx,yy_front,zz]], axis=0)
        for possible_zz in np.arange(0,zz,0.1):
            wall_points_front = np.append(wall_points_front, [[xx,yy_front,possible_zz]], axis=0)
    possible_y = np.arange(0,Y_LENGTH+gridsize,gridsize)
    yy_back = possible_y[-1]
    for xx in np.arange(0,X_LENGTH+gridsize,gridsize):
        zz = getZElevation(xx,yy_back)               
        wall_points_back = np.append(wall_points_back,[[xx,yy_back,zz]], axis=0)
        for possible_zz in np.arange(0,zz,0.1):
            wall_points_back = np.append(wall_points_back, [[xx,yy_back,possible_zz]], axis=0)  
    wall_points_front = np.delete(wall_points_front,0,0)
    wall_points_back = np.delete(wall_points_back,0,0)
    return wall_points_front, wall_points_back

#used this thread to figure out how to save each file with a new unique name: https://stackoverflow.com/questions/33691187/how-to-save-the-file-with-different-name-and-not-overwriting-existing-one        
def uniqueFileName(basename, ext):
    '''Each time the code runs, it is saved with a unique file name'''
    actualname = "%s.%s" % (basename, ext)
    c = itertools.count()
    while os.path.exists(actualname):
        actualname = "%s (%d).%s" % (basename, next(c), ext)
    return actualname

#after making a discrete grid of surface points, the next step is finding the area of boxes in that grid so that I can weight 
#those areas for depositing nuclei
def areaofIrregularQuad(point,l1,l2,l3,l4):
    '''area of a quadrilateral that is irregular, i.e. not a square/rectangle'''
    '''get all four lengths of the sides
    divide the quad into two triangles with a diagonal down the middle
    find the length of the diagonal using the distance formula between points 1 and 3 of the quad
    find areas of both triangles using heron's formula
    add the areas of the two triangles together'''
    diagonal_len = distanceFormula(surface_points[point,0],surface_points[point,1],surface_points[point,2],surface_points[point+(GRIDSIZE_MULTIPLIER*X_LENGTH+2),0],surface_points[point+(GRIDSIZE_MULTIPLIER*X_LENGTH+2),1],surface_points[point+(GRIDSIZE_MULTIPLIER*X_LENGTH+2),2])
    s1 = (l1+l2+diagonal_len)/2
    s2 = (l3+l4+diagonal_len)/2
    area_tri_1 = np.sqrt(s1*(s1-l1)*(s1-l2)*(s1-diagonal_len))
    area_tri_2 = np.sqrt(s2*(s2-l3)*(s2-l4)*(s2-diagonal_len))
    area_total = area_tri_1 + area_tri_2
    return area_total

def totalVolume(surface_points):
    '''find the total volume under the surface
    can do this by dividing each of the squares in the grid into two triangles
    if the points of each square are 1-4, with the top left being 1 and then proceeding clockwise,
    draw the line for the triangle from point 1 to point 3. Then the formula for area of each truncated
    triangular prism is V = A 1/3 (z1 + z2 + z3) where A is area of base of prism'''
    #Area of the base of the triangular prism is half the area of one flat grid cell (z=0), with an area of 0.5*0.5 (grid wall lengths)
    A = GRIDSIZEINPUT_FORSURFACE**2/2
    Vol_total = 0
    number_squares = 0
    for ones in np.arange(0,Y_LENGTH*GRIDSIZE_MULTIPLIER):
   
        grid_range = np.array(np.arange(ones*X_LENGTH*GRIDSIZE_MULTIPLIER+ones,X_LENGTH*GRIDSIZE_MULTIPLIER + ones*X_LENGTH*GRIDSIZE_MULTIPLIER+ones))

        for point in grid_range:
            point1 = [surface_points[point,0],surface_points[point,1],surface_points[point,2]]
            point4 = [surface_points[point+(GRIDSIZE_MULTIPLIER*X_LENGTH+1),0], surface_points[point+(GRIDSIZE_MULTIPLIER*X_LENGTH+1),1], surface_points[point+(GRIDSIZE_MULTIPLIER*X_LENGTH+1),2]] 
            point3 = [surface_points[point+(GRIDSIZE_MULTIPLIER*X_LENGTH+2),0],surface_points[point+(GRIDSIZE_MULTIPLIER*X_LENGTH+2),1],surface_points[point+(GRIDSIZE_MULTIPLIER*X_LENGTH+2),2]]
            point2 = [surface_points[point+1,0], surface_points[point+1,1], surface_points[point+1,2]]
            
            V1_tri = (point1[2]+point2[2]+point3[2])*(1/3)*A
            V2_tri = (point1[2]+point3[2]+point4[2])*(1/3)*A
            Vol_square = V1_tri + V2_tri
            Vol_total += Vol_square
            number_squares += 1
    return Vol_total
 

#Grows each nucleus according to inorganic rate law - calculates new r
def growEachNucleus(nucleus_array):
    '''After each timestep, grow each nucleus according to inorganic rate law defined above
    Updates r with every iteration of time, and updates z if the nucleus is not located on the ground'''
    #surface area of each hemisphere in um2 
    SA = nucleus_array[:,3]*nucleus_array[:,3] * np.pi*4/2 
    #mol per sphere
    Delta_moles = Growth_Rate*SA*DELTA_T 
    #um of growth to add
    Delta_volume = Delta_moles*MOLARV_ARAG 
    #volume of each hemisphere in um2
    VOL = (nucleus_array[:,3]*nucleus_array[:,3]*nucleus_array[:,3] * np.pi*4/3)/2 
    NEW_VOL = VOL + Delta_volume;
    NEW_R = (NEW_VOL*3/4*2/np.pi)**(1/3)  
    #new radius based on growth in this step is input back into the array
    nucleus_array[:,3] = NEW_R  
    return sum(VOL)

#use distance formula to get lenths of all four sides of each box               
##at the end of each row, I don't want it to wrap around to the next row
#think I solved that issue by calculating the area row by row - from 0 to 99 for each x row
#calculating the box areas moving from one side of the grid to the other, so calculating the area from 0 to 1 in x lengths, then from 1 to 2, and so on, so from 99 to 100 is the last one calculating and we don't want to calculate from 100 to anything, but simply move on to the next row
#for this reason, we can just use arange ending in 100 because we don't want to actually use that number or the box that starts with that number
#this loop provides the area of 10000 boxes in the grid with x lenght and y length equal to 100 in the array all_areas        
def areasofEachCellinGrid(surface_points):
    '''Uses the distance formula between the vertices of each box of the grid and calculates the 
    area of each triangular half of each box, then adds them together for each box, and saves each area in an array called all_areas'''
    all_areas = np.zeros(X_LENGTH/GRIDSIZEINPUT_FORSURFACE * Y_LENGTH/GRIDSIZEINPUT_FORSURFACE)
    area_count = 0
    for ones in np.arange(0,Y_LENGTH*GRIDSIZE_MULTIPLIER):
        
        grid_range = np.array(np.arange(ones*X_LENGTH*GRIDSIZE_MULTIPLIER+ones,X_LENGTH*GRIDSIZE_MULTIPLIER + ones*X_LENGTH*GRIDSIZE_MULTIPLIER+ones))

        for point in grid_range:
            l1 = distanceFormula(surface_points[point,0],surface_points[point,1],surface_points[point,2], surface_points[point+1,0], surface_points[point+1,1], surface_points[point+1,2])
            l2 = distanceFormula(surface_points[point+1,0], surface_points[point+1,1], surface_points[point+1,2], surface_points[point+(GRIDSIZE_MULTIPLIER*X_LENGTH+2),0],surface_points[point+(GRIDSIZE_MULTIPLIER*X_LENGTH+2),1],surface_points[point+(GRIDSIZE_MULTIPLIER*X_LENGTH+2),2])
            l3 = distanceFormula(surface_points[point+(GRIDSIZE_MULTIPLIER*X_LENGTH+2),0],surface_points[point+(GRIDSIZE_MULTIPLIER*X_LENGTH+2),1],surface_points[point+(GRIDSIZE_MULTIPLIER*X_LENGTH+2),2], surface_points[point+(GRIDSIZE_MULTIPLIER*X_LENGTH+1),0], surface_points[point+(GRIDSIZE_MULTIPLIER*X_LENGTH+1),1], surface_points[point+(GRIDSIZE_MULTIPLIER*X_LENGTH+1),2])
            l4 = distanceFormula(surface_points[point+(GRIDSIZE_MULTIPLIER*X_LENGTH+1),0], surface_points[point+(GRIDSIZE_MULTIPLIER*X_LENGTH+1),1], surface_points[point+(GRIDSIZE_MULTIPLIER*X_LENGTH+1),2], surface_points[point,0],surface_points[point,1],surface_points[point,2])                              
            area_cell = areaofIrregularQuad(point,l1,l2,l3,l4)
            
            all_areas[area_count] = area_cell
            area_count += 1
    return all_areas

def findNewXYZ(areas_gridcells, random_location):
    '''Find which cell grid you should put a new nucleus in (the grid is area weighted) and 
    then place the nucleus randomly within that grid cell'''
    #go through the cell areas and add them up until you get to the value you generated in random_location
    which_cell_area = 0
    area_count = 0
    for area in areas_gridcells:
            if which_cell_area <= random_location:
                which_cell_area += area
                area_count += GRIDSIZEINPUT_FORSURFACE
    #finding the point1 of the square that is going to have the nucleus in it - point1 is determined by the  lowest 
    #x and y value pair of the four corners
    area_count = area_count - GRIDSIZEINPUT_FORSURFACE
    point1_y = (area_count%(Y_LENGTH*GRIDSIZE_MULTIPLIER))/GRIDSIZE_MULTIPLIER
    point1_x = (area_count - point1_y)/(X_LENGTH*GRIDSIZE_MULTIPLIER)
    #pick a point randomly in the (gridsize,gridsize) x,y box that has a vertice of point1 - basically making 
    #an assumption that weighting of area across a box this small is negligibly different
    #this is the new nucleation location on a surface! 
    new_weighted_x = random.random() + point1_x
    new_weighted_y = random.random() + point1_y    
    new_weighted_z = getZElevation(new_weighted_x, new_weighted_y)  
    return new_weighted_x, new_weighted_y, new_weighted_z

time_groundcover_10 = []
time_groundcover_25 = []
time_groundcover_50 = []
time_groundcover_75 = []
time_groundcover_final = []

def findTimetoCoverGround(percentcoverage_firstlayer):
    '''amount of time it takes to cover the ground to different percentages'''  
    if percentcoverage_firstlayer > 10 and percentcoverage_firstlayer < 15:
        time_groundcover_10.append(t)
    
    if percentcoverage_firstlayer > 25 and percentcoverage_firstlayer < 30:
        time_groundcover_25.append(t) 
    
    if percentcoverage_firstlayer > 50 and percentcoverage_firstlayer < 55:
        time_groundcover_50.append(t)    
    
    if percentcoverage_firstlayer > 75 and percentcoverage_firstlayer < 80:
        time_groundcover_75.append(t)     

    if percentcoverage_firstlayer > 90 and percentcoverage_firstlayer < 100:
        time_groundcover_final.append(t)
    return time_groundcover_10, time_groundcover_25, time_groundcover_50, time_groundcover_75, time_groundcover_final 



def porosityofSkeleton(nuclei, volume, surface_points):
    '''find the porosity by finding all of the places in the skeleton where there are no nuclei, 
    then subtracting that number from total volume'''
    numIntersecting = 0
    numPointsTested = 0
    nuclei_x = nuclei[:,0]
    nuclei_y = nuclei[:,1]
    nuclei_z = nuclei[:,2]
    nuclei_r = nuclei[:,3]
    
    for gridpoint in range(np.size(surface_points[:,0])):
        z = surface_points[gridpoint, 2]
        if z > 0:
            for zz in np.arange(0,z,0.1):
                if doesPointIntersectAnyNucleus(surface_points[gridpoint,0],surface_points[gridpoint,1],zz, nuclei_x, nuclei_y, nuclei_z, nuclei_r) ==True:
                    numIntersecting = numIntersecting + 1
                numPointsTested = numPointsTested + 1
    if numPointsTested > 0:            
        porosity_percent = 1 - numIntersecting/numPointsTested
    else:
        porosity_percent = 0
    return numIntersecting, numPointsTested, porosity_percent, volume*porosity_percent

#density of nuclei on the ground
def nucleiGroundDensity(nuclei):
    '''total number of nuclei on the ground divided by ground area'''
    nuclei_ground_count = 0
    for nucleus in range(len(nuclei)):
        if nuclei[nucleus,2] == 0:
            nuclei_ground_count += 1
    nuclei_density = nuclei_ground_count/(X_LENGTH*Y_LENGTH)
    return nuclei_ground_count, nuclei_density

def verticalExtension(max_height, surface_points, t):
    '''Counts the zs at every iteration that are above or equal to a certain arbitrary height that is chosen for each omega - if
    90% of the zs are higher than or equal to that height, it divides by amount of time it took to get that high, resultin gin vertical 
    extension'''  
    z_count = 0        
    for z in surface_points[:,2]:
        if z >= max_height:
            z_count += 1
        if z_count >= 0.9 * len(surface_points[:,2]):
            vertical_extension = max_height/t
        else:
            vertical_extension = 0
    return vertical_extension

dict_times_output = {}
def outputAtCertainTimes(t, volume):
    '''Defines a dictionary of each time as the key and the number of nuclei, amount of floor covered, total growth at that time as the values, ratio of nuclei/growth, calcification, nuclei ground count, vertical extension rate, average height across the surface, std deviation of heights'''  
    vertical_extension = verticalExtension(max_height, surface_points, t)
    z_sum = sum(surface_points[:,2])
    number_z = len(surface_points[:,2])
    ave_z = z_sum/number_z
    std_dev_z = np.std(surface_points[:,2])
    dict_times_output[t] = np.size(nuclei[:,0]), percentcoverage_firstlayer, volume, np.size(nuclei[:,0])/volume, volume/(X_LENGTH*Y_LENGTH*t), nucleiGroundDensity(nuclei)[0], vertical_extension, ave_z, std_dev_z
    return dict_times_output


            

def plot3DSpheres(nuclei, omega, max_t):
    '''plotting hemispheres for each nucleus'''            
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim((0,X_LENGTH))
    ax.set_ylim((0,Y_LENGTH))
    ax.set_zlim((0,Z_LENGTH))
    #useful website for understanding spherical coordinates: https://mathinsight.org/spherical_coordinates
    #so radius determines the size of the spheres, and the angles determine how much of a sphere they are (ie half, full, other shapes, etc)
    shape = np.arange(0,nuclei.shape[0])    ##making an array from 0 to the number of rows of nuclei, increasing by 1          
            
    #adding everything as hemispheres - think about if we want to add things as spheres on top of bottom layer
    for row in np.nditer(shape):       
        
        r=nuclei[row-1,3]                                 ##radius
        #the 10 in u and v just determines the number of lines drawn between points around the edge of the sphere (1/10 of a circle) - a higher number looks more like a circle but takes more time to draw
        u = np.linspace(0, 2*np.pi, 10)     ##theta in spherical coordinates - (how far around the x and y axes in a full circle we go)
        v = np.linspace(0, np.pi/2, 10)    ##phi in spherical coordinates - determines how far along z axis sphere is (so forms a full half circle from fully positive to fully negative z when at full pi units)
        x = r * np.outer(np.cos(u), np.sin(v))+nuclei[row-1,0]
        y = r * np.outer(np.sin(u), np.sin(v))+nuclei[row-1,1] 
        z = r * np.outer(np.ones(np.size(u)), np.cos(v))
        ax.plot_surface(x, y, z, color='b', alpha=0.5, clip_on = True)    

    plt1_name = uniqueFileName('/Users/Marta/Documents/Python/nucleation_model_output/plots_VE/modelnucleisphere_om' + str(omega) + '_' +  str(max_t) + '_',  'png')
    fig.savefig(plt1_name)
    plt.close()

def plotSurfaceGrid(nuclei, omega, max_t):  
    '''plot the 3D grid surface of the nuclei and the walls'''  
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim((0,X_LENGTH+GRIDSIZEINPUT_FORSURFACE))
    ax.set_ylim((0,Y_LENGTH+GRIDSIZEINPUT_FORSURFACE))
    ax.set_zlim((0,Z_LENGTH))
    
    Z = surface_points[:,2]
    walls = wallsofGrid(GRIDSIZEINPUT_FORSURFACE)
    walls_front = walls[0]
    walls_back = walls[1]
    
    Z_front = walls_front[:,2]
    Z_back = walls_back[:,2]
    
    ax.scatter(walls_back[:,0], walls_back[:,1], walls_back[:,2], c=Z_back, alpha=0.5, clip_on = True, s=20, lw=1)
    ax.scatter(walls_front[:,0], walls_front[:,1], walls_front[:,2], c=Z_front, alpha=0.5, clip_on = True, s=20, lw=1)
    ax.scatter(surface_points[:,0], surface_points[:,1], surface_points[:,2], c=Z, alpha=0.5, clip_on = True, s=20, lw=1)  
         
    plt2_name = uniqueFileName('/Users/Marta/Documents/Python/nucleation_model_output/plots_VE/modelnucleisurface_om' + str(omega) + '_' + str(max_t) + '_',  'png')
    fig.savefig(plt2_name)                                                                      
    plt.close() 



#python version of timing
#from: https://stackoverflow.com/questions/7370801/measure-time-elapsed-in-python
start = timer()

#####################################################################################
########################USER DEFINED STUFF
#Define the dimensions of a control volume in microns:
X_LENGTH = 100 #µm
Y_LENGTH = 100 #µm
Z_LENGTH = 100 #µm
DELTA_T = 10 #time step in seconds

maximum_t = [1000]

#max_height is chosen for each omega to be the bar to reach for vertical extension - a height that is high enough to not be influence
#by the first layer of nuclei on the ground
max_height = 15

#All nuclei start at a specified size = seed size 
SEED_RADIUS = 0.5  #µm radius

#Growth only depends on radius and is independent of
#origin. Use rate law Rate = k(Omega-1)^n where k = 11 nmol  m-2 s-1 and
#n=1.7 (from Alex's summary figure that he sent me). If Omega is constant then this Growth_Rate is always the same
#It is not clear that this bulk growth rate scales down to this scale
omega_values = [90]

#molar volume of aragonite in µm3/mol = MW (g/mol) / density (g/cm3) * 1E12
MOLARV_ARAG = 100.09/2.93*1E12

#Nucleation rate is also constant at constant Omega. 
#J = Aexp(-balpha^3 / (ln(omega))^2) in nuclei  µm-2 sec-1
A = 1858506.76
BALPHA3 = -19.44553408

#The discrete surface grid is made up of boxes with sides of length GRIDSIZEINPUT_FORSURFACE
#The GRIDSIZE_MULTIPLIER is a ratio that makes up for the sides not being 1
#Right now - have to pick a multiple of 0.5 for the code to work
GRIDSIZEINPUT_FORSURFACE = 1
GRIDSIZE_MULTIPLIER = 1/GRIDSIZEINPUT_FORSURFACE

###############################################################################################
################## Simulation
for omega in omega_values:    
    Growth_Rate = ((omega-1)**1.7)*(11*10E-9) #mol m-2 sec-1
    Growth_Rate = (Growth_Rate/1000/1000/1000/1000)      #change units to mol/um/s
    #nuclei/m2/s converted to nuclei/um2 for every 10 seconds
    J_rate = A * np.exp(BALPHA3/np.log(omega)**2)/1000/1000/1000/1000*10 
    for max_t in maximum_t:
        #Start with one nuclei randomly distrubuted on the XY plane. 
        nuclei = np.array([[ random.random()*X_LENGTH,  random.random()*Y_LENGTH, 0, SEED_RADIUS]]); #made a 2D array because can't concatenate it correctly eitherwise - it will only append to the end of the same array, it will not add a new array to the array of arrays
        ##                       [x                                  y            z     radius]
        #recording the time step of each nuclei that is formed and the time that each percentage of the 2d yx grid is covered
        nuclei_timeofdeposition = [0]
        ttime = np.arange(0,max_t,DELTA_T)        
        for t in ttime:
            print(t)
            ## Grow Spheres
            #Each sphere is described by a radius and an origin. 
            #This information saved in an array nuclei that grows in length with more nuclei.            
            #After each timestep grow each sphere according to the inorganic rate law. 
            growEachNucleus(nuclei)           
            #The surface is made up of boxes and the vertices of the boxes are saved in the array surface_points
            surface_points = Discrete3dSurface(GRIDSIZEINPUT_FORSURFACE)            
            areas = areasofEachCellinGrid(surface_points)                        
            #now want to be able to pick a point within each small cell, once that cell is picked to nucleate within
            #also want to be able to take the sum of all areas, then pick a number randomly from 0 to the sum, then find out which cell that number is referring to
            sum_areas = sum(areas)
            #pick an number randomly from 0 through sum_areas
            random_location = random.random()*sum_areas           
            ## Nucleation        
            #Based on probability of a nucleus being deposited (based on omega and area)
            #Need to select from the right probability distribution -> think about
            #this! For now, assume in steady-state regime so a uniform distribution
            #seems reasonable
            new_nuclei = None
            Nuclei_flag = 0       
            prob_nuc = X_LENGTH*Y_LENGTH*J_rate*DELTA_T
            if prob_nuc > 1:
                print('ERROR: reduce time step size - nucleation too fast')            
            randomnumber = random.random() 
            new_coordinates = findNewXYZ(areas, random_location)
            if randomnumber<=prob_nuc:
                #generate a new nuclei on x,y,z surface
                new_nuclei = np.array([[ new_coordinates[0],  new_coordinates[1], new_coordinates[2], SEED_RADIUS]]);  #put an extra set of brackets around this to make it a 2D array (even though it only has 1 row and is a 1d array)
                Nuclei_flag = 1              
            if Nuclei_flag == 1:  #if flag = 0, don't add the nuclei, but if it's 1, then it's ok to add                       
                #adding the nuclei continuously so that they can be added on top of each other        
                nuclei = np.append(nuclei, [new_nuclei[0,:]], axis=0)
                print(nuclei)
                nuclei_timeofdeposition = np.append(nuclei_timeofdeposition,t)              
            
            percentcoverage_firstlayer = getSampledPercentageAreaOccupiedByNuclei(nuclei[:,0],nuclei[:,1],nuclei[:,3])*100
            #Record when ground is covered to certain percentages
            time_groundcoverage = findTimetoCoverGround(percentcoverage_firstlayer)
            volume = totalVolume(surface_points)
            timed_output = outputAtCertainTimes(t, volume)    
            
        #plotting spheres in 3d            
        plot3DSpheres(nuclei, omega, max_t)
        #plotting the surface of the nuclei     
        plotSurfaceGrid(nuclei, omega, max_t)
   
        
    
    #Output parameters
    print('Number of nuclei:', np.size(nuclei[:,0]))
    
    
    #save the following information in a file: nuclei, time of deposition, omega, and time step
    filename = uniqueFileName('/Users/Marta/Documents/Python/nucleation_model_output/text_files_VE/'+str(omega)+'_'+str(max_t)+'_', 'txt')
    
    file = open(filename, 'a')
    file.write('\n' + 'Nucleus X, Nucleus Y, Nucleus Z, Nucleus R, Time of Deposition, Timestep' + '\n')
    time_between_dep = []
    for nucleus in range(len(nuclei)):
        file.write(str(nuclei[nucleus]) + ', ')
        time_dep = nuclei_timeofdeposition[nucleus]
        file.write(str(time_dep) + ', ')
        file.write(str(DELTA_T) + '\n')
        if nucleus is not 0:
            time_between_dep = np.append(time_between_dep, time_dep - nuclei_timeofdeposition[nucleus-1])
    file.write('\n' + 'Omega:' + str(omega) + '\n')
    #file.write('\n' + '.5J' + '\n' )
    file.write('\n' + 'Total Number Nuclei:' +str(np.size(nuclei[:,0])) + '\n')    
    file.write('\n' + 'Total Calculated Volume:' +str(volume) + ' um3' + '\n')
    file.write('\n' + 'Actual Volume with Overlaps:' + str(growEachNucleus(nuclei)) + ' um3' + '\n')
    file.write('\n' + 'Time to Cover Ground 10%:' + str(time_groundcoverage[0]) + 's' + '\n')
    file.write('\n' + 'Time to Cover Ground 25%:' + str(time_groundcoverage[1]) + 's' + '\n')
    file.write('\n' + 'Time to Cover Ground 50%:' + str(time_groundcoverage[2]) + 's' + '\n')
    file.write('\n' + 'Time to Cover Ground 75%:' + str(time_groundcoverage[3]) + 's' + '\n')
    file.write('\n' + 'Time to Cover Ground 90%:' + str(time_groundcoverage[4]) + 's' + '\n')
    file.write('\n' + 'Total Time:' + str(max_t) + ' s' + '\n')
    file.write('\n' + 'Average Time Step Between Depositions:' + str(np.average(time_between_dep)) + 's' + '\n')
    file.write('\n' + 'Number Nuclei on Ground Level:' + str(nucleiGroundDensity(nuclei)[0]) + '\n')
    file.write('\n' + 'Nuclei Ground Density:' + str(nucleiGroundDensity(nuclei)[1]) + ' nuclei/um2' + '\n')
    file.write('\n' + 'Ratio of Nuclei to Growth:' + str(np.size(nuclei[:,0])/volume) + '\n')
    for key in sorted(timed_output):
        file.write('\n' + str(key) + ',' + str(timed_output[key]) + '\n')

file.close()
    

end = timer()
#time in seconds
print('time elapsed:',end - start)
        