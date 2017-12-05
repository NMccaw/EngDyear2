#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 15:47:50 2017

@author: nicholas
"""
import sys
import csv
import numpy as np

path2points = '/home/nicholas/Documents/Summer_Proj/Python3_codes/'
sys.path.insert(0,path2points)

up_points = []
low_points = []
keep_flag_up = 0
keep_flag_low = 0

N_probes = 300
points_x_up = np.zeros(N_probes)
points_y_up = np.zeros_like(points_x_up)

points_x_low = np.zeros(N_probes)
points_y_low = np.zeros_like(points_x_up)

f = open(path2points+'sampleDict') 
lines = f.readlines()
f.close()
w = open(path2points+'sampleDict','w')    
w.writelines([item for item in lines[:-2]])
w.close()





with open(path2points+'sampleDict') as f:
    for i,line in enumerate(f):
        if keep_flag_up == 1 and '(' and ')' in line:
            up_points.append(line.strip())
        if keep_flag_low == 1 and '(' and ')' in line:
            low_points.append(line.strip())
        if 'surfaceProbesUpper23' in line:
            keep_flag_up = 1
            keep_flag_low = 0
        if 'surfaceProbesLower23' in line:
            keep_flag_low = 1
            keep_flag_up = 0
        if 'surfaceProbesUpper' in line:
            keep_flag_low = 0

for i,line in enumerate(up_points):
    newline = ''
    char_keep = ''
    for letter in line:
        if letter == '(' or letter == ')' or letter == ';' or letter == ');':
            pass
        else:
            char_keep = letter
        newline += str(char_keep) 
    if newline.split() != []:
        i -= 1
        points_x_up[i] = float(newline.split()[0])
        points_y_up[i] = float(newline.split()[1])

for i,line in enumerate(low_points):
    newline = ''
    char_keep = ''
    for letter in line:
        if letter == '(' or letter == ')' or letter == ';' or letter == ');':
            pass
        else:
            char_keep = letter
        newline += str(char_keep) 
    if newline.split() != []:
        i -= 1
        points_x_low[i] = float(newline.split()[0])
        points_y_low[i] = float(newline.split()[1])

percent_span_list = np.linspace(80,100,21)

with open('sampleDict','a') as f:
    for j,percent_span in enumerate(percent_span_list):
        f.write('\n' +'surfaceProbesUpper' + str(int(percent_span)) + '\n' + '{' +'\n')
        f.write('type patchCloud;' + '\n'+ 'axis xyz;'+'\n' + 'maxDistance 0.1;'+'\n' +'patches' +'\n' +'(')
        f.write('\n' + 'naca_ascii'+'\n'+');'+ '\n'+'points'+'\n'+'('+'\n' )
        for i in range(len(points_x_up)):
            f.write('(' +str(points_x_up[i]) + ' ' + str(points_y_up[i]) + ' ' + str(percent_span/100)+')'+'\n' )
        f.write('\n'+');'+'\n'+'}' + '\n')
            
            
with open('sampleDict','a') as f:
    for j,percent_span in enumerate(percent_span_list):
        f.write('\n' +'surfaceProbesLower' + str(int(percent_span)) + '\n' + '{' +'\n')
        f.write('type patchCloud;' + '\n'+ 'axis xyz;'+'\n' + 'maxDistance 0.1;'+'\n' +'patches' +'\n' +'(')
        f.write('\n' + 'naca_ascii'+'\n'+');'+ '\n'+'points'+'\n'+'('+'\n' )
        for i in range(len(points_x_up)):
            f.write('(' +str(points_x_low[i]) + ' ' + str(points_y_low[i]) + ' ' + str(percent_span/100)+')'+'\n' )
        f.write('\n'+');'+'\n'+'}' + '\n')
    f.write(');')
    f.close()