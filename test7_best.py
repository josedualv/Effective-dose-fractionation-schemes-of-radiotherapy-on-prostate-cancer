# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 16:23:37 2020

@author: josed
"""


import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def tumor_PC3( y, t, cprolif_c, cprolif_r, cCapacity_c, cCapacity_r, lambda_c, lambda_r):
    Vc, Vr = y

    dVcdt = cprolif_c * Vc * (1 - (Vc/cCapacity_c) - (lambda_r*(Vr/cCapacity_c)))
    dVrdt = cprolif_r * Vr * (1 - (Vr/cCapacity_r) - (lambda_c*(Vc/cCapacity_r)))
    
    return dVcdt, dVrdt

# Given parameters
    
cprolif_c = 0.015; cprolif_r = 0.02;

cCapacity_c = 0.85; cCapacity_r = 2;

lambda_c = 0.2; lambda_r = 0;


def surviving_fraction( alpha, beta, dosage):
    
    sigma = math.exp( -alpha * dosage - beta * dosage * dosage )

    return sigma

t = np.linspace(0, 336, 337)

# Initial Conditions (1):
y0 = [0.5, 0.5] 

sol = odeint( tumor_PC3, y0, t, args = ( cprolif_c, cprolif_r, cCapacity_c, cCapacity_r, lambda_c, lambda_r) )

plt.title("PC3")
plt.figure(1)
plt.plot(t, sol[:, 0], 'b', label='sensitive')
plt.plot(t, sol[:, 1], 'g', label='resistant')
plt.xlabel('time (hours)')
plt.ylabel('volume')
plt.legend(loc='best')
plt.grid()

# Given Parameters:

tin = 336

a_res = 0.300
b_res = 0.0402

a_sen = 0.430
b_sen = 0.0407

area_res = 0
area_sen = 0

### Resistant Population ###

# Week 1

for dos_count in range (0, 12, 3):

    y0_ = sol[-1, :]
    y0 = y0_ * surviving_fraction(a_res, b_res, 3)
    t = np.linspace(tin, tin + 24, tin + 25)
    sol = odeint( tumor_PC3, y0, t, args = ( cprolif_c, cprolif_r, cCapacity_c, cCapacity_r, lambda_c, lambda_r) )
    t = np.append( tin, t )
    sol = np.append( [y0_], sol, 0 )        
    plt.plot(t, sol[:, 1], 'g', label = 'resistant')
    tin += 24

# Trapezoidal Rule: 
    
    if tin >= 361:
        long = sol[0, 1]
        area_res += 24/2 * (short + long)
    short = sol[1,1]

# First Friday:

y0_ = sol[-1, :]
y0 = y0_ * surviving_fraction(a_res, b_res, 4)
t = np.linspace(tin, tin + 72, tin + 73)
sol = odeint( tumor_PC3, y0, t, args = ( cprolif_c, cprolif_r, cCapacity_c, cCapacity_r, lambda_c, lambda_r) )

# First Weekend:

t = np.append( tin, t )
sol = np.append( [y0_], sol, 0 )        
plt.plot(t, sol[:, 1], 'g', label='resistant')
tin += 72

# Trapezoidal Rule: 

long = sol[0,1]
#the number below (24) might be 72 (due to "height" og the trapezoid - fix later)
area_res += 24/2 * (short + long)
short = sol[1,1]

# Weeks 2-6: 

tin = 504
x = 0

for tin in range (tin, 1344, 168):
   
    for tin in range (tin, tin + (5 * 24), 48):
        
        if tin == 1272:
            dos = 2
        else:
            dos = 3
            
        y0_ = sol[-1, :]
        y0 = y0_ * surviving_fraction(a_res, b_res, dos)
        t = np.linspace(tin, tin + 48, tin + 49)
        sol = odeint( tumor_PC3, y0, t, args = ( cprolif_c, cprolif_r, cCapacity_c, cCapacity_r, lambda_c, lambda_r) )
        t = np.append( tin, t )
        sol = np.append( [y0_], sol, 0 )        
        plt.plot(t, sol[:, 1], 'g', label='resistant')
        
# Big vs. Small Trapezoid:
        
        if tin == 504 + x:
            length = 72
        else:
            length = 48
        
        long = sol[0, 1]
        area_res += length/2 * (short + long)
        short = sol [1,1]
    
    t = np.linspace(tin + 48, tin + 72, tin + 73)
    sol = odeint( tumor_PC3, sol[-1,:], t, args = ( cprolif_c, cprolif_r, cCapacity_c, cCapacity_r, lambda_c, lambda_r) )
    plt.plot(t, sol[:, 1], 'g', label='resistant')
    x += 168

final_vol_res = sol[1344, 1]
area_res += 72/2 * (short + final_vol_res)

# Results:

print('Area under resistant tumor curve = ', area_res)
print('Final volume of resistant tumor = ', final_vol_res)

# Reset the graph:

tin = 336
t = np.linspace(0, 336, 337)

# Initial Conditions (2):
y0 = [0.5, 0.5] 
sol = odeint( tumor_PC3, y0, t, args = ( cprolif_c, cprolif_r, cCapacity_c, cCapacity_r, lambda_c, lambda_r) )

### Sensitive Population ###

# Week 1

for dos_count in range (0, 12, 3):

    y0_ = sol[-1, :]
    y0 = y0_ * surviving_fraction(a_sen, b_sen, 3)
    t = np.linspace(tin, tin + 24, tin + 25)
    sol = odeint( tumor_PC3, y0, t, args = ( cprolif_c, cprolif_r, cCapacity_c, cCapacity_r, lambda_c, lambda_r) )
    t = np.append( tin, t )
    sol = np.append( [y0_], sol, 0 )        
    plt.plot(t, sol[:, 0], 'b', label='sensitive')
    tin += 24

# Trapezoidal Rule: 
    
    if tin >= 361:
        long = sol[0, 0]
        area_sen += 24/2 * (short + long)
    short = sol[1,0]
    
# First Weekend:

y0_ = sol[-1, :]
y0 = y0_ * surviving_fraction(a_sen, b_sen, 4)
t = np.linspace(tin, tin + 72, tin + 73)
sol = odeint( tumor_PC3, y0, t, args = ( cprolif_c, cprolif_r, cCapacity_c, cCapacity_r, lambda_c, lambda_r) )
t = np.append( tin, t )
sol = np.append( [y0_], sol, 0 )        
plt.plot(t, sol[:, 0], 'b', label='sensitive')
tin += 72

# Trapezoidal Rule: 

long = sol[0,0]
area_sen += 24/2 * (short + long)
short = sol[1,0]

# Weeks 2-6: 

tin = 504
x = 0

for tin in range (tin, 1344, 168):
   
    for tin in range (tin, tin + (5 * 24), 48):

# **dosage = 2 on the last friday
        
        if tin == 1272:
            dos = 2
        else:
            dos = 3
            
        y0_ = sol[-1, :]
        y0 = y0_ * surviving_fraction(a_sen, b_sen, dos)
        t = np.linspace(tin, tin + 48, tin + 49)
        sol = odeint( tumor_PC3, y0, t, args = ( cprolif_c, cprolif_r, cCapacity_c, cCapacity_r, lambda_c, lambda_r) )
        t = np.append( tin, t )
        sol = np.append( [y0_], sol, 0 )        
        plt.plot(t, sol[:, 0], 'b', label='sensitive')   
        
# Big vs. Small Trapezoid:
    
        if tin == 504 + x:
            length = 72
        else:
            length = 48
        
        long = sol[0, 0]
        area_sen += length/2 * (short + long)
        short = sol [1,0]
    
    t = np.linspace(tin + 48, tin + 72, tin + 73)
    sol = odeint( tumor_PC3, sol[-1, :], t, args = ( cprolif_c, cprolif_r, cCapacity_c, cCapacity_r, lambda_c, lambda_r) )
    plt.plot(t, sol[:, 0], 'b', label='sensitive')
    x += 168

final_vol_sen = sol[1344, 0]
area_sen += 72/2 * (short + final_vol_sen)

# Results:

print('Area under sensitive curve = ', area_sen)
print('Final volume of sensitive tumor = ', final_vol_sen)





