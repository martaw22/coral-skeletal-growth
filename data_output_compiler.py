# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 14:45:43 2019

@author: Marta
"""

'''This program will open all the text files that are in a specific folder that are the outputs of the model 
nucleation code.  It will open them one by one, read out all the outputs, and then write them to new text files
that will each contain only one parameter'''

#import numpy as np
import os

#helpful thread about opening files in a directory: https://stackoverflow.com/questions/18262293/how-to-open-every-file-in-a-folder


filename = '/Users/Marta/Documents/Python/nucleation_model_output/masterlist_nucleationmodelrun_2'

masterfile = open(filename, 'w')
masterfile.write('File Name, Omega, Total Number Nuclei, Total Calculated Volume (um3), Actual Volume with Overlaps (um3),Percent of Ground Covered (%), Percent of Ground Covered at 500 s (%), Percent of Ground Covered at 1000 s (%), Total Calculated Volume at 500 s (um3), Total Calculated Volume at 1000 s (um3), Time to Cover Ground 10% (s), Time to Cover Ground 25% (s), Time to Cover Ground 50% (s), Time to Cover Ground 75% (s), Time to Cover Ground 90% (s), Total Time (s), Average Time Step Between Depositions (s), Number Nuclei on Ground Level, Nuclei Ground Density (nuclei/um2), Ratio of Nuclei to Growth' + '\n')

path = '/Users/Marta/Documents/Python/nucleation_model_output/text_files_2/'
#can use os.getcwd() for current directory or can use this path to get there from anywhere

test = '/Users/Marta/Documents/Python/nucleation_model_output/text_files_2/60_3000_(0).txt'
testopen = open(test)
for line in testopen:
    print(line)
testopen.close()


number_files = 0

for file in os.listdir(path):
    number_files += 1
    file_working = open(path + file)
    masterfile.write(str(file) + ', ')
    for line_of_text in file_working:
        line_of_text = line_of_text.replace('\n', ' ')
        if 'Omega' in line_of_text:
            masterfile.write(line_of_text[6:] + ', ')
        
        if 'Total Number Nuclei' in line_of_text:
    
            masterfile.write(line_of_text[20:] + ', ')
        if 'Total Calculated' in line_of_text:
            
            masterfile.write(line_of_text[24:-5] + ', ')
        if 'Actual Volume Volume (um3)' in line_of_text:
            
            masterfile.write(line_of_text[28:-5] + ', ')
        if 'Percent of Ground Covered (%)' in line_of_text: 
            
            masterfile.write(line_of_text[26:-2] + ', ')
            
        if 'Percent of Ground Covered at 500' in line_of_text: 
            
            masterfile.write(line_of_text[26:-2] + ', ')  
        if 'Percent of Ground Covered at 1000' in line_of_text: 
            
            masterfile.write(line_of_text[26:-2] + ', ')
        if 'Total Calculated Volume at 500 s' in line_of_text: 
            
            masterfile.write(line_of_text[26:-2] + ', ')
        if 'Total Calculated Volume at 1000 s' in line_of_text: 
            
            masterfile.write(line_of_text[26:-2] + ', ')
        if 'Time to Cover Ground 10' in line_of_text:
            
            masterfile.write(line_of_text[26:-1] + ', ')    
        if 'Time to Cover Ground 25' in line_of_text:
            
            masterfile.write(line_of_text[26:-1] + ', ')    
        
        if 'Time to Cover Ground 50' in line_of_text:
            #print(line_of_text[25:-1])
            masterfile.write(line_of_text[26:-1] + ', ')    
        
        if 'Time to Cover Ground 75' in line_of_text:
            
            masterfile.write(line_of_text[26:-1] + ', ')    
    
        if 'Time to Cover Ground 90' in line_of_text:
            
            masterfile.write(line_of_text[26:-1] + ', ')
        if 'Total Time' in line_of_text:
            
            masterfile.write(line_of_text[11:-3] + ', ')
        if 'Average Time' in line_of_text:
            
            masterfile.write(line_of_text[38:-2] + ', ')
        if 'Number Nuclei on Ground' in line_of_text:
            
            masterfile.write(line_of_text[30:] + ', ')
        if 'Nuclei Ground Density' in line_of_text:
            
            masterfile.write(line_of_text[22:-12] + ', ')
        if 'Ratio of Nuclei' in line_of_text:
            
            masterfile.write(line_of_text[26:] + ' ' + '\n')
    file_working.close()
            
masterfile.close()        