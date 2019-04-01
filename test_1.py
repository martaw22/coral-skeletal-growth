# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 11:18:23 2019

@author: Marta
"""

import os

cwd = os.getcwd()

filename = str(cwd) + '/results/test_2.csv'
    
file = open(filename, 'a')
file.write('this is the second test')

file.close()