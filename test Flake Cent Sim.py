# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 14:55:47 2016

@author: Administrator
"""
from scipy import linspace
from scipy import pi,sqrt,exp
from scipy.special import erf
#import math
import numpy as np
import matplotlib.pyplot as plt

def pdf(x):
    return 1/sqrt(2*pi) * exp(-x**2/2)
def cdf(x):
    return (1 + erf(x/sqrt(2))) / 2
def skew(x,e=0,w=1,a=0):
    t = (x-e) / w
    return 2 / w * pdf(t) * cdf(a*t) 

# make population distribution for length and aspect ratio len/thickness not N(layers)
x_len = linspace(20,990,98)
p_len = skew(x_len, 100.0, 200.0, 4) # location, scale, skew
norm_p_len = p_len/sum(p_len)


height_low_lim, height_top_lim, height_spacing = 1, 31, 1  ### height in mm
length_low_lim, length_top_lim, length_spacing = 20, 1000, 10  ### length in nm
thickness_low_lim, thickness_top_lim, thickness_spacing = 1, 100, 1  ### layer number

mat_height, mat_length, mat_thickness = (height_top_lim-height_low_lim), ((length_top_lim-length_low_lim)/length_spacing), ((thickness_top_lim-thickness_low_lim)/thickness_spacing)
flake_array = np.ones((mat_length,mat_thickness))

# make heights
height_list = []
for height in range(height_low_lim,height_top_lim,height_spacing):
    height = float(height)/1000
    height_list.append(height)
    
# make flakes
flake_thickness_list = []
flake_length_list = []
for thickness in range(thickness_low_lim,thickness_top_lim,thickness_spacing):
    thickness = float(thickness)
    flake_thickness_list.append(thickness)
for length in range(length_low_lim,length_top_lim,length_spacing):
    length = float(length)/1e9
    flake_length_list.append(length)
    
# make fraction for each flake length & thickness combination
len_counter = 0
for i in flake_length_list:
    thi_counter = 0
    for j in flake_thickness_list:
        aspect_ratio = i/(j*1e-9)
        if aspect_ratio < 1:
            p_ar = 0
        else:
            p_ar = skew(aspect_ratio, 5.0, 11.0, 4) # location, scale, skew
        flake_array[len_counter][thi_counter] = p_ar*p_len[len_counter]
        thi_counter += 1
    len_counter += 1

### don't forget to normalise the frxn for each length across the thicknesses so they sum to p_len[x]
### then normalise so that the total flake array sums to 1

total_flake_array = np.tile(flake_array, (mat_height,1,1))
norm_flake_array = total_flake_array/(sum(sum(sum(total_flake_array))))  # Normalised array of frxn of all length and thickness and height in vial
test_sum = sum(sum(sum(norm_flake_array)))

#dist_list_list = []
#dist_list_list.append(flake_sizes)
  
# Experimental Parameters   
impulse_list = [45]#, 45, 90, 150, 200, 270, 350, 400, 500, 600, 750, 1000]
rpm_list = [500]#, 500, 1000, 2000, 3000, 4000, 5000, 7500, 10000]
exponent_list = [2]#, 2.025, 2.05, 2.075, 2.1, 2.125, 2.15, 2.175, 2.2]
rho_water = 1000
rho_flake = 2200
rho_surf = 1000
eta_water = 0.00089
d_surf_layer = float(2)/1e9  #thickness of surfactant layer
layer_thickness = 6e-10  #thickness of single layer/sheet in nm
L_t_ratio = 1.2  #Length/Nlayers ratio
width_due_to_back_diff_etc = 30
length_of_rotor = 0.1

flake_avg_length_array = np.zeros((len(exponent_list),len(impulse_list),len(rpm_list)))
flake_avg_thickness_array = np.zeros((len(exponent_list),len(impulse_list),len(rpm_list)))
frxn_rem_length_array = np.zeros((len(exponent_list),len(impulse_list),len(rpm_list),len(flake_length_list)))
frxn_rem_thickness_array = np.zeros((len(exponent_list),len(impulse_list),len(rpm_list),len(flake_thickness_list)))
 
# Start of Simulation
exponent_counter = 0
for exponent in exponent_list:
    impulse_counter = 0
    for impulse in impulse_list:
        rpm_counter = 0
        for rpm in rpm_list:
            total_frxn_remain_matrix = norm_flake_array
            cent_time = (float(impulse)/(0.000001*rpm*rpm))*60 # seconds
            
            # Finding Position for each flake dimension at each height
            position_array = np.zeros((mat_height,mat_length,mat_thickness))
            mat_height_element = 0
            for height in height_list:
                mat_length_element = 0
                for length in flake_length_list:
                    mat_thickness_element = 0
                    for thickness in flake_thickness_list:
                        position = (length_of_rotor + (height_list[-1]-height))*(exp((((thickness*layer_thickness*(rho_flake-rho_water))+(2*d_surf_layer*(rho_surf-rho_water)))*(length*(rpm**exponent)*cent_time))/((10.5*eta_water)*(sqrt(L_t_ratio)))))
                        position_array[mat_height_element][mat_length_element][mat_thickness_element] = position
                        mat_thickness_element += 1
                    mat_length_element += 1
                mat_height_element += 1
            
            # Fraction Left of each Flake Dimension at each height
            frxn_remaining_array = np.zeros((mat_height,mat_length,mat_thickness))
            pos_diff_array = np.zeros((mat_height,mat_length,mat_thickness))
    
            for mat_h in range(mat_height):
                for mat_l in range(mat_length):
                    for mat_t in range(mat_thickness):
                        position_diff = (length_of_rotor + height_list[-1]) - position_array[mat_h][mat_l][mat_t]
                        pos_diff_array[mat_h][mat_l][mat_t] = position_diff
                        
                        ### ignores any diffusion effects
                        if position_diff <= 0:
                            frxn_remaining_array[mat_h][mat_l][mat_t] = 0
                        else:
                            frxn_remaining_array[mat_h][mat_l][mat_t] = total_frxn_remain_matrix[mat_h][mat_l][mat_t]
                        
                        ### includes diffusion
                        #try:
                        #    frxn_remaining = (1/(1+(math.exp(-width_due_to_back_diff_etc*position_diff))))
                        #except:
                        #    frxn_remaining = 0
                        #frxn_remaining_array[mat_h][mat_l][mat_t] = frxn_remaining*total_frxn_remain_matrix[mat_h][mat_l][mat_t]

            
            bottom_line = sum(sum(sum(frxn_remaining_array[:,:,:])))
            counter = 0
            top_line = 0
            top_line_length = 0
            counter_length = 0
            
            # Average Flake Thickness
            for thickness in flake_thickness_list:
                f = sum(sum(frxn_remaining_array[:,:,counter]))
                f_thick = f*thickness ### multiplying total frxn remaining by thickness
                top_line += f_thick
                frxn_rem_thickness_array[exponent_counter][impulse_counter][rpm_counter][counter] = f
                counter += 1
            mean_flake_thickness = top_line/bottom_line
            flake_avg_thickness_array[exponent_counter][impulse_counter][rpm_counter] = mean_flake_thickness
            
            
            # Average Flake Length
            for length in flake_length_list:
                h = (sum(sum(frxn_remaining_array[:,counter_length,:])))
                h_length = h*length
                top_line_length += h_length
                frxn_rem_length_array[exponent_counter][impulse_counter][rpm_counter][counter_length] = h
                counter_length += 1
            mean_flake_length = top_line_length/bottom_line
            flake_avg_length_array[exponent_counter][impulse_counter][rpm_counter] = mean_flake_length
            
            rpm_counter += 1
        impulse_counter += 1
    exponent_counter += 1
    
#Aspect Ratio Array
aspect_ratio_array = flake_avg_length_array/flake_avg_thickness_array

#np.savetxt("flake_avg_length_array.csv", flake_avg_length_array, delimiter=",")
#np.savetxt("flake_avg_thickness_array.csv", flake_avg_thickness_array, delimiter=",")
#np.savetxt("aspect_ratio_array.csv", aspect_ratio_array, delimiter=",")

        
#plt.ylabel("flake length")
#plt.xlabel("impulse")
plt.subplot(211)
plt.plot(x_len,norm_p_len)
x_ar = linspace(2,98,97)
p_ar_list = skew(x_ar, 5.0, 11.0, 4) # location, scale, skew
plt.subplot(212)
plt.plot(x_ar,p_ar_list)
