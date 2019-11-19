'''This module contains all of the functions that are used in the spheresmaster.py script.'''
import random
import os
import itertools
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import config

cwd = os.getcwd()

#Function for determining if a given point intersects with any nucleus in nuclei
def doesPointIntersectAnyNucleus(x, y, z, nuclei_x, nuclei_y, nuclei_z, nuclei_r):
    '''This function computes distance from x,y to nucleus center and see if it's less than nucleus radius.

    The inputs are a coordinates of a point on the grid of the surface of the skeleton and the coordinates and radius of a nucleus.
    It returns True if the point does intersect with the nucleus.
    '''
    distance = np.sqrt((x-nuclei_x)**2 + (y-nuclei_y)**2 + (z-nuclei_z)**2)
    if np.any(distance <= nuclei_r):
        return True
    else:
        return "no intersection"

##These functions are used to find the area of the 2D xy grid occupied by nuclei
def getSampledPercentageAreaOccupiedByNuclei(nuclei_xlayer, nuclei_ylayer, nuclei_rlayer):
    '''This function calculates the area of the xy plane occupied by nuclei at each timepoint.
    
    The inputs are arrays of the x and y coordinates and radii of the nuclei that are present at that timepoint.
    It returns a float.'''
    numIntersecting = 0
    numPointsTested = 0
    for xx in range(config.X_LENGTH+1):
        for yy in range(config.Y_LENGTH+1):
            if doesPointIntersectAnyNucleus(xx, yy, 0, nuclei_xlayer, nuclei_ylayer, 0, nuclei_rlayer) == True:                                    
                numIntersecting = numIntersecting + 1
            numPointsTested = numPointsTested + 1
    percentcoverage = numIntersecting/numPointsTested
    return percentcoverage

##Determines the z elevation at a given point of the potential nuclei there
def getElevationonNucleus(x, y, x0, y0, z0, r):
    '''This function determines the height of a nucleus if it is present at a specific point on the surface grid.

    x and y are the point location that we are testing the height of, while x0 and y0 are the center of the nucleus.
    The output is z, the height of the nucleus from the ground at that point on the surface grid.
    If the nucleus is not present at that point, the height is 0.
    '''
    inside_equation = r**2-(x-x0)**2-(y-y0)**2             
    if inside_equation <= 0:
        z = 0
    else:
        if z0 == 0:
            z = np.sqrt(r**2-(x-x0)**2-(y-y0)**2)
        else:
            z = np.sqrt(r**2-(x-x0)**2-(y-y0)**2) + z0
    return z

def getElevationUnderNucleus(x, y, x0, y0, r):
    '''This function calculates the height of a nucleus from any point within the radius of that nucleus from the bottom of the hemisphere.
    
    x and y are the point that we are testing the height of, while x0 and y0 are the center of the nucleus.
    The output is z, the height of hemisphere at that point within the nucleus' radius.
    It is not calculating the overall height of that point on the nucleus from the ground.
    '''
    inside_equation = r**2-(x-x0)**2-(y-y0)**2             
    if inside_equation <= 0:
        z = 0
    else:
        z = np.sqrt(r**2-(x-x0)**2-(y-y0)**2)
    return z
    
def getZElevation(x, y, nuclei):
    '''Determines the elevation of a specific nucleus from the ground at a specific point in the xy grid.
    
    The given point's coordinate is the input x, y, and the code loops through the existing nuclei coordinates, input as an array containing the coordinates and radii of all the nuclei.
    It checks to see if any nuclei overlap with the given point, and stores the highest height recorded at that spot from any nuclei that are present.
    The output is a float.
    '''
    highestZSeen = 0
    nucleisize = np.array(np.size(nuclei, axis=0))
    nuclei_x = nuclei[:, 0]
    nuclei_y = nuclei[:, 1]
    nuclei_z = nuclei[:, 2]
    nuclei_r = nuclei[:, 3]
    for row in range(nucleisize):
        height = getElevationonNucleus(x, y, nuclei_x[row], nuclei_y[row], nuclei_z[row], nuclei_r[row])
        if height > highestZSeen:
            highestZSeen = height
    return highestZSeen

def getZElevationUnderNucleus(x, y, nuclei):
    '''Determines the height of a specific nucleus from the bottom of the hemisphere at a specific point on the grid.
    
    The coordinates of the given point are the inputs x, y, and the nuclei coordinates and radii are input as an array.
    The output is a float.
    '''
    highestZSeen = 0
    nucleisize = np.array(np.size(nuclei, axis=0))
    nuclei_x = nuclei[:, 0]
    nuclei_y = nuclei[:, 1]
    nuclei_r = nuclei[:, 3]
    for row in range(nucleisize):
        height = getElevationUnderNucleus(x, y, nuclei_x[row], nuclei_y[row], nuclei_r[row])
        if height > highestZSeen:
            highestZSeen = height
    return highestZSeen

##doing the weighting of areas that are within nuclei to account for surface area
def distanceFormula(x, y, z, x0, y0, z0):
    '''Calculates the distance between two points in a 3D grid.
    
    The inputs are the coordinates of the two points, and they must each have an x, y, and z value.'''
    distance = np.sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
    return distance

##Making a surface with discrete points in x,y,z space
def Discrete3dSurface(gridsize, nuclei):
    '''This function creates an array that represents the surface grid of the skeleton that has formed so far.

    It takes as inputs the xy length of the size of one box in the surface grid, and the array of the nuclei coordinates and radii.
    The function queries each integer point in the x, y, z box to determine where the surface is.
    It returns an array of arrays, where each internal array is the vertice (x,y,z) of one of the boxes of the surface grid.
    '''
    numPointsTested = 0
    surface_points = np.zeros(((int(config.X_LENGTH/gridsize) + 1) * (int(config.Y_LENGTH/gridsize) + 1), 3))
    for xx in np.arange(0, config.X_LENGTH+gridsize, gridsize):
        for yy in np.arange(0, config.Y_LENGTH+gridsize, gridsize):
            zz = getZElevation(xx, yy, nuclei)
            surface_points[numPointsTested][0] = xx
            surface_points[numPointsTested][1] = yy
            surface_points[numPointsTested][2] = zz
            numPointsTested = numPointsTested + 1
    return surface_points

#making a function to make arrays of the walls of the grid that will grow with the nuclei that touch the edge of the grid
def wallsofGrid(gridsize, nuclei):
    '''This function creates two arrays that represent vertices of the 3D grid that forms at the edges of the box we are building the skeleton in.

    The inputs are the length of the xy grid cell that makes up the grid and the array of coordinates and radii of the nuclei.
    We are interested in the edges of the xy grid of the box and finding the z values of the grid.
    The output arrays will be used with the surface grid to plot the overall surface of the skeleton.
    '''
    wall_points_front = [[-100, -100, -100]]
    wall_points_back = [[-100, -100, -100]]
    xx_back = 0
    for yy in np.arange(0, config.Y_LENGTH+gridsize, gridsize):
        zz = getZElevation(xx_back, yy, nuclei)       
        wall_points_back = np.append(wall_points_back,[[xx_back, yy, zz]], axis=0)
        for possible_zz in np.arange(0, zz, 0.1):
            wall_points_back = np.append(wall_points_back, [[xx_back, yy, possible_zz]], axis=0)
    possible_x = np.arange(0, config.X_LENGTH+gridsize, gridsize)
    xx_front = possible_x[-1]
    for yy in np.arange(0, config.Y_LENGTH+gridsize, gridsize):
        zz = getZElevation(xx_front, yy, nuclei)               
        wall_points_front = np.append(wall_points_front, [[xx_front, yy, zz]], axis=0)
        for possible_zz in np.arange(0, zz, 0.1):
            wall_points_front = np.append(wall_points_front, [[xx_front, yy, possible_zz]], axis=0)

    yy_front = 0
    for xx in np.arange(0, config.X_LENGTH+gridsize, gridsize):
        zz = getZElevation(xx, yy_front, nuclei)
        wall_points_front = np.append(wall_points_front, [[xx, yy_front, zz]], axis=0)
        for possible_zz in np.arange(0, zz, 0.1):
            wall_points_front = np.append(wall_points_front, [[xx, yy_front, possible_zz]], axis=0)
    possible_y = np.arange(0, config.Y_LENGTH+gridsize, gridsize)
    yy_back = possible_y[-1]
    for xx in np.arange(0, config.X_LENGTH+gridsize, gridsize):
        zz = getZElevation(xx, yy_back, nuclei)
        wall_points_back = np.append(wall_points_back, [[xx, yy_back, zz]], axis=0)
        for possible_zz in np.arange(0, zz, 0.1):
            wall_points_back = np.append(wall_points_back, [[xx, yy_back, possible_zz]], axis=0)
    wall_points_front = np.delete(wall_points_front, 0, 0)
    wall_points_back = np.delete(wall_points_back, 0, 0)
    return wall_points_front, wall_points_back

#used this thread to figure out how to save each file with a new unique name: https://stackoverflow.com/questions/33691187/how-to-save-the-file-with-different-name-and-not-overwriting-existing-one        
def uniqueFileName(basename, ext):
    '''Each time the code runs, the outputs are saved with a unique file name.

    The inputs are the the string format of the cwd and the constants that were chosen, as well as the format type.
    '''
    actualname = "%s.%s" % (basename, ext)
    c = itertools.count()
    while os.path.exists(actualname):
        actualname = "%s (%d).%s" % (basename, next(c), ext)
    return actualname

#after making a discrete grid of surface points, the next step is finding the area of boxes in that grid so that I can weight 
#those areas for depositing nuclei
def areaofIrregularQuad(point, l1, l2, l3, l4, surface_points):
    '''Finding the area of a quadrilateral that is irregular, i.e. not a square/rectangle
    
    The inputs are the x,y,z of the vertice of one of the grid boxes, the lengths of each of the sides of that grid box, 
    and the array of all of the vertices of the surface grid.
    It works by dividing the quad into two triangles with a diagonal down the middle.
    Then it finds the length of the diagonal using the distance formula between points 1 and 3 of the quad.
    It finds areas of both triangles using heron's formula, adn then adds the areas of the two triangles together.
    The output is a float.
    '''
    diagonal_len = distanceFormula(surface_points[point, 0], surface_points[point, 1], surface_points[point, 2], surface_points[point+(int(config.GRIDSIZE_MULTIPLIER*config.X_LENGTH)+2), 0], surface_points[point+(int(config.GRIDSIZE_MULTIPLIER*config.X_LENGTH)+2), 1],surface_points[point+(int(config.GRIDSIZE_MULTIPLIER*config.X_LENGTH)+2), 2])
    s1 = (l1+l2+diagonal_len)/2
    s2 = (l3+l4+diagonal_len)/2
    area_tri_1 = np.sqrt(s1*(s1-l1)*(s1-l2)*(s1-diagonal_len))
    area_tri_2 = np.sqrt(s2*(s2-l3)*(s2-l4)*(s2-diagonal_len))
    area_total = area_tri_1 + area_tri_2
    return area_total

def totalVolume(surface_points):
    '''This function finds the total volume under the surface grid on top of the already formed skeleton.
    
    It does this by dividing each of the squares in the grid into two triangles.
    If the points of each square are 1-4, with the top left being 1 and then proceeding clockwise,
    draw the line for the triangle from point 1 to point 3. Then the formula for area of each truncated
    triangular prism is V = A 1/3 (z1 + z2 + z3) where A is area of base of prism.
    The output is a float.
    '''
    #Area of the base of the triangular prism is half the area of one flat grid cell (z=0), with an area of 0.5*0.5 (grid wall lengths)
    A = config.GRIDSIZEINPUT_FORSURFACE**2/2
    Vol_total = 0
    number_squares = 0
    for ones in np.arange(0, config.Y_LENGTH*config.GRIDSIZE_MULTIPLIER):
   
        grid_range = np.array(np.arange(ones*config.X_LENGTH*config.GRIDSIZE_MULTIPLIER+ones, config.X_LENGTH*config.GRIDSIZE_MULTIPLIER + ones*config.X_LENGTH*config.GRIDSIZE_MULTIPLIER+ones))

        for point in grid_range:
            point = int(point)
            point1 = [surface_points[point, 0], surface_points[point, 1],surface_points[point, 2]]
            point4 = [surface_points[point+(int(config.GRIDSIZE_MULTIPLIER*config.X_LENGTH)+1), 0], surface_points[point+(int(config.GRIDSIZE_MULTIPLIER*config.X_LENGTH)+1), 1], surface_points[point+(int(config.GRIDSIZE_MULTIPLIER*config.X_LENGTH)+1), 2]] 
            point3 = [surface_points[point+(int(config.GRIDSIZE_MULTIPLIER*config.X_LENGTH)+2), 0], surface_points[point+(int(config.GRIDSIZE_MULTIPLIER*config.X_LENGTH)+2), 1], surface_points[point+(int(config.GRIDSIZE_MULTIPLIER*config.X_LENGTH)+2), 2]]
            point2 = [surface_points[point+1,0], surface_points[point+1,1], surface_points[point+1,2]]
            V1_tri = (point1[2]+point2[2]+point3[2])*(1/3)*A
            V2_tri = (point1[2]+point3[2]+point4[2])*(1/3)*A
            Vol_square = V1_tri + V2_tri
            Vol_total += Vol_square
            number_squares += 1
    return Vol_total

#Grows each nucleus according to inorganic rate law - calculates new r
def growEachNucleus(nucleus_array, Growth_Rate):
    '''After each timestep, grow each nucleus according to inorganic rate law defined above.

    Updates the nuclei array: r with every iteration of time, and updates z if the nucleus is not located on the ground.
    Returns a float of the sum of all of the volumes of nuclei present so far.
    '''
    #surface area of each hemisphere in um2 
    SA = nucleus_array[:, 3]*nucleus_array[:, 3] * np.pi*4/2 
    #mol per sphere
    Delta_moles = Growth_Rate*SA*config.DELTA_T 
    #um of growth to add
    Delta_volume = Delta_moles*config.MOLARV_ARAG 
    #volume of each hemisphere in um2
    VOL = (nucleus_array[:, 3]*nucleus_array[:, 3]*nucleus_array[:, 3] * np.pi*4/3)/2 
    NEW_VOL = VOL + Delta_volume
    NEW_R = (NEW_VOL*3/4*2/np.pi)**(1/3)  
    #new radius based on growth in this step is input back into the array
    nucleus_array[:, 3] = NEW_R
    #Need to update the z of nuclei that are on top of other nuclei as the nuclei beneath them grow higher
    numberrows = np.size(nucleus_array[:, 0])    
    for row in range(numberrows):
        if nucleus_array[row,2] > 0:
            newz = getZElevationUnderNucleus(nucleus_array[row, 0], nucleus_array[row, 1], nuclei)
            nucleus_array[row,2] = newz
    return sum(VOL)

def areasofEachCellinGrid(surface_points):
    '''Uses the distance formula between the vertices of each box of the grid and calculates the 
    area of each triangular half of each box, then adds them together for each box, and saves each area in an array called all_areas.
    
    '''
    all_areas = np.zeros(int(config.X_LENGTH/config.GRIDSIZEINPUT_FORSURFACE) * int(config.Y_LENGTH/config.GRIDSIZEINPUT_FORSURFACE))
    area_count = 0
    for ones in np.arange(0,config.Y_LENGTH*config.GRIDSIZE_MULTIPLIER):
        grid_range = np.array(np.arange(ones*config.X_LENGTH*config.GRIDSIZE_MULTIPLIER+ones, config.X_LENGTH*config.GRIDSIZE_MULTIPLIER + ones*config.X_LENGTH*config.GRIDSIZE_MULTIPLIER+ones))
        for point in grid_range:
            point = int(point)
            l1 = distanceFormula(surface_points[point, 0], surface_points[point, 1], surface_points[point, 2], surface_points[point+1, 0], surface_points[point+1, 1], surface_points[point+1, 2])
            l2 = distanceFormula(surface_points[point+1,0], surface_points[point+1,1], surface_points[point+1,2], surface_points[point+(int(config.GRIDSIZE_MULTIPLIER*config.X_LENGTH)+2),0],surface_points[point+(int(config.GRIDSIZE_MULTIPLIER*config.X_LENGTH)+2),1],surface_points[point+(int(config.GRIDSIZE_MULTIPLIER*config.X_LENGTH)+2),2])
            l3 = distanceFormula(surface_points[point+(int(config.GRIDSIZE_MULTIPLIER*config.X_LENGTH)+2),0],surface_points[point+(int(config.GRIDSIZE_MULTIPLIER*config.X_LENGTH)+2),1],surface_points[point+(int(config.GRIDSIZE_MULTIPLIER*config.X_LENGTH)+2),2], surface_points[point+(int(config.GRIDSIZE_MULTIPLIER*config.X_LENGTH)+1),0], surface_points[point+(int(config.GRIDSIZE_MULTIPLIER*config.X_LENGTH)+1),1], surface_points[point+(int(config.GRIDSIZE_MULTIPLIER*config.X_LENGTH)+1),2])
            l4 = distanceFormula(surface_points[point+(int(config.GRIDSIZE_MULTIPLIER*config.X_LENGTH)+1),0], surface_points[point+(int(config.GRIDSIZE_MULTIPLIER*config.X_LENGTH)+1),1], surface_points[point+(int(config.GRIDSIZE_MULTIPLIER*config.X_LENGTH)+1),2], surface_points[point,0],surface_points[point,1],surface_points[point,2])                              
            area_cell = areaofIrregularQuad(point,l1,l2,l3,l4, surface_points)
            
            all_areas[area_count] = area_cell
            area_count += 1
    return all_areas

def findNewXYZ(areas_gridcells, random_location, nuclei):
    '''Find which cell grid you should put a new nucleus in (the grid is area weighted) and 
    then place the nucleus randomly within that grid cell.
    
    '''
    #go through the cell areas and add them up until you get to the value you generated in random_location
    which_cell_area = 0
    area_count = 0
    for area in areas_gridcells:
            if which_cell_area <= random_location:
                which_cell_area += area
                area_count += config.GRIDSIZEINPUT_FORSURFACE
    #finding the point1 of the square that is going to have the nucleus in it - point1 is determined by the  lowest 
    #x and y value pair of the four corners
    area_count = area_count - config.GRIDSIZEINPUT_FORSURFACE
    point1_y = (area_count%(config.Y_LENGTH*config.GRIDSIZE_MULTIPLIER))/config.GRIDSIZE_MULTIPLIER
    point1_x = (area_count - point1_y)/(config.X_LENGTH*config.GRIDSIZE_MULTIPLIER)
    #pick a point randomly in the (gridsize,gridsize) x,y box that has a vertice of point1 - basically making 
    #an assumption that weighting of area across a box this small is negligibly different
    #this is the new nucleation location on a surface! 
    new_weighted_x = random.random() + point1_x
    new_weighted_y = random.random() + point1_y
    new_weighted_z = getZElevation(new_weighted_x, new_weighted_y, nuclei)
    return new_weighted_x, new_weighted_y, new_weighted_z

#Used in the findTimetoCoverGround function
time_groundcover_10 = []
time_groundcover_25 = []
time_groundcover_50 = []
time_groundcover_75 = []
time_groundcover_final = []

def findTimetoCoverGround(percentcoverage_firstlayer):
    '''Calculates the amount of time it takes to cover the ground to different percentages.
    
    '''  
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
    '''Calculates the porosity by finding all of the places in the skeleton where there are no nuclei, 
    then subtracting that number from total volume.
    
    This is not turned on in the spheresmaster_script.py because it makes the code take a lot longer to run.'''
    numIntersecting = 0
    numPointsTested = 0
    nuclei_x = nuclei[:, 0]
    nuclei_y = nuclei[:, 1]
    nuclei_z = nuclei[:, 2]
    nuclei_r = nuclei[:, 3]
    for gridpoint in range(np.size(surface_points[:, 0])):
        z = surface_points[gridpoint, 2]
        if z > 0:
            for zz in np.arange(0, z, 0.1):
                if doesPointIntersectAnyNucleus(surface_points[gridpoint, 0],surface_points[gridpoint, 1], zz, nuclei_x, nuclei_y, nuclei_z, nuclei_r) ==True:
                    numIntersecting = numIntersecting + 1
                numPointsTested = numPointsTested + 1
    if numPointsTested > 0:
        porosity_percent = 1 - numIntersecting/numPointsTested
    else:
        porosity_percent = 0
    return numIntersecting, numPointsTested, porosity_percent, volume*porosity_percent

#density of nuclei on the ground
def nucleiGroundDensity(nuclei):
    '''Calculates the total number of nuclei on the ground divided by ground area.
    
    '''
    nuclei_ground_count = 0
    for nucleus in range(len(nuclei)):
        if nuclei[nucleus, 2] == 0:
            nuclei_ground_count += 1
    nuclei_density = nuclei_ground_count/(config.X_LENGTH*config.Y_LENGTH)
    return nuclei_ground_count, nuclei_density

def verticalExtension(MAX_HEIGHT, surface_points):
    '''Counts the zs at every iteration that are above or equal to a certain arbitrary height that is chosen for each OMEGA 
    
    If 90% of the zs are higher than or equal to that height, it divides by amount of time it took to get that high, resulting in vertical 
    extension calculation.
    '''  
    z_count = 0
    for z in surface_points[:, 2]:
        if z >= MAX_HEIGHT:
            z_count += 1
        if z_count >= 0.9 * len(surface_points[:, 2]):
            vertical_extension = MAX_HEIGHT/config.MAX_T
        else:
            vertical_extension = 0
    return vertical_extension

#Used in the outputAtCertainTimes function
dict_times_output = {}
def outputAtCertainTimes(t, volume, surface_points, nuclei, percentcoverage_firstlayer):
    '''Outputs a dictionary of results for each timepoint.
    
    Each time as the key. 
    The number of nuclei, amount of floor covered, total growth at that time as the values, ratio of nuclei/growth, calcification, nuclei ground count, vertical extension rate are the values.
    '''  
    vertical_extension = verticalExtension(config.MAX_HEIGHT, surface_points)
    z_sum = sum(surface_points[:, 2])
    number_z = len(surface_points[:, 2])
    ave_z = z_sum/number_z
    std_dev_z = np.std(surface_points[:, 2])
    dict_times_output[t] = np.size(nuclei[:, 0]), percentcoverage_firstlayer, volume, np.size(nuclei[:, 0])/volume, volume/(config.X_LENGTH*config.Y_LENGTH*t), nucleiGroundDensity(nuclei)[0], vertical_extension, ave_z, std_dev_z
    return dict_times_output

def plot3DSpheres(nuclei, OMEGA):
    '''Plots hemispheres to represent each nucleus.
    
    The figure is plotted, saved, and closed.
    It is given a unique filename and saved in the 'results' folder in the cwd.
    '''            
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim((0, config.X_LENGTH))
    ax.set_ylim((0, config.Y_LENGTH))
    ax.set_zlim((0, config.Z_LENGTH))
    #so radius determines the size of the spheres, and the angles determine how much of a sphere they are (ie half, full, other shapes, etc)
    shape = np.arange(0, nuclei.shape[0])    ##making an array from 0 to the number of rows of nuclei, increasing by 1          
    #adding everything as hemispheres - think about if we want to add things as spheres on top of bottom layer
    for row in np.nditer(shape):
        r=nuclei[row-1, 3]
        #the 10 in u and v just determines the number of lines drawn between points around the edge of the sphere (1/10 of a circle) - a higher number looks more like a circle but takes more time to draw
        u = np.linspace(0, 2*np.pi, 10)     ##theta in spherical coordinates - (how far around the x and y axes in a full circle we go)
        v = np.linspace(0, np.pi/2, 10)    ##phi in spherical coordinates - determines how far along z axis sphere is (so forms a full half circle from fully positive to fully negative z when at full pi units)
        x = r * np.outer(np.cos(u), np.sin(v))+nuclei[row-1, 0]
        y = r * np.outer(np.sin(u), np.sin(v))+nuclei[row-1, 1]
        z = r * np.outer(np.ones(np.size(u)), np.cos(v))
        ax.plot_surface(x, y, z, color='b', alpha=0.5, clip_on =True)

    plt1_name = uniqueFileName(str(cwd) + '/results/om' + str(config.OMEGA) + '_' +  str(config.MAX_T) + '_' + str(config.SEED_RADIUS) + 'r_' + str(config.ALPHA_MULTIPLIER) + 'alpha',  'png')
    fig.savefig(plt1_name)
    plt.close()

def plotSurfaceGrid(nuclei, OMEGA, surface_points):  
    '''Plots the 3D grid surface and walls of the nuclei present.
    
    The figure is plotted, saved, and closed.
    It is given a unique filename and saved in the 'results' folder in the cwd.
    '''  
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim((0, config.X_LENGTH+config.GRIDSIZEINPUT_FORSURFACE))
    ax.set_ylim((0, config.Y_LENGTH+config.GRIDSIZEINPUT_FORSURFACE))
    ax.set_zlim((0, config.Z_LENGTH))
    
    Z = surface_points[:, 2]
    walls = wallsofGrid(config.GRIDSIZEINPUT_FORSURFACE, nuclei)
    walls_front = walls[0]
    walls_back = walls[1]
    
    Z_front = walls_front[:, 2]
    Z_back = walls_back[:, 2]
    
    ax.scatter(walls_back[:, 0], walls_back[:, 1], walls_back[:, 2], c=Z_back, alpha=0.5, clip_on =True, s=20, lw=1)
    ax.scatter(walls_front[:, 0], walls_front[:, 1], walls_front[:, 2], c=Z_front, alpha=0.5, clip_on =True, s=20, lw=1)
    ax.scatter(surface_points[:, 0], surface_points[:, 1], surface_points[:, 2], c=Z, alpha=0.5, clip_on =True, s=20, lw=1)

    plt2_name = uniqueFileName(str(cwd) + '/results/om' + str(config.OMEGA) + '_' + str(config.MAX_T) + '_' + str(config.SEED_RADIUS) + 'r_' + str(config.ALPHA_MULTIPLIER) + 'alpha',  'png')
    fig.savefig(plt2_name)                       
    plt.close()
