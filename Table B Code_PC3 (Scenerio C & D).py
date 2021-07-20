# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 19:56:01 2020

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
y0 = [0.5, 0.5] 
# y0 = [0.1, 0.9] 

sol = odeint( tumor_PC3, y0, t, args = ( cprolif_c, cprolif_r, cCapacity_c, cCapacity_r, lambda_c, lambda_r) )

plt.title("PC3")
plt.figure(1)
plt.plot(t, sol[:, 0], 'b', label='parental')
plt.plot(t, sol[:, 1], 'g', label='resistant')
plt.xlabel('time (hours)')
plt.ylabel('volume')
plt.legend(loc='best')
plt.grid()

tin = 336
total_hours = 1344
dos_count = 0


# dos_weekday = 2.5
# dos_friday = 4
dos_weekday = float(input("Enter weekday dose: "))
dos_friday = float(input("Enter friday dose: "))



a_res = 0.300
b_res = 0.0402

a_sen = 0.430
b_sen = 0.0407

tinit = 336

area_res = 0
area_sen = 0

count = 0 

short = 0 
long = 0 

# for twk in range (tinit, 1344, 168):
for twk in range (tinit, 672, 168):
    for tin in range (twk, twk + (3 * 24), 48):
      
        y0_ = sol[-1, :]
        y0 = y0_ * surviving_fraction(a_res, b_res, dos_weekday)
        t = np.linspace(tin, tin + 48, tin + 49)
        sol = odeint( tumor_PC3, y0, t, args = ( cprolif_c, cprolif_r, cCapacity_c, cCapacity_r, lambda_c, lambda_r) )
        t = np.append( tin, t )
        sol = np.append( [y0_], sol, 0 )        
        plt.plot(t, sol[:, 1], 'g', label='resistant')

        
        if tin >= twk + 48:
            long = sol[0, 1]
            # print ("long = ", long)
            area_res += 24/2 * (short + long)
            count += 1 
            # print( 'growth after', count, '-th treatment added' )
         
        short = sol[1, 1]
        # print ("short = ", short)
  
#Friday Dose:
            
    y0_ = sol[-1, :]
    y0 = y0_ * surviving_fraction(a_res, b_res, dos_friday)
    t = np.linspace(tin + 48, tin + 48, tin + 49)
    sol = odeint( tumor_PC3, y0, t, args = ( cprolif_c, cprolif_r, cCapacity_c, cCapacity_r, lambda_c, lambda_r) )
    t = np.append( tin + 48, t )
    sol = np.append( [y0_], sol, 0 )        
    plt.plot(t, sol[:, 1], 'g', label='resistant')
    
    long = sol[0,1]
    # print ("long = ", long)
    area_res += 48/2 * (short + long)
    count += 1
    # print( 'growth after', count, '-th treatment added' )
    

### Weekend Period
    
    t = np.linspace(tin + 48, tin + 120, tin + 121)
    sol = odeint( tumor_PC3, sol[-1,:], t, args = ( cprolif_c, cprolif_r, cCapacity_c, cCapacity_r, lambda_c, lambda_r) )
    plt.plot(t, sol[:, 1], 'g', label='resistant')

    count += 1
    short = sol[0,1]
    print(short)
    # print("short = ", short)
    long = sol[-1, 1]
    # print ("long = ", long)
    area_res += 72/2 * (short + long)
    # print( 'growth after', count, '-th treatment added' )
        
print(area_res)
# print(sol[1343, 1])

count = 0
tinit = 336

print("break")

t = np.linspace(0, 336, 337)
y0 = [0.5, 0.5] 
sol = odeint( tumor_PC3, y0, t, args = ( cprolif_c, cprolif_r, cCapacity_c, cCapacity_r, lambda_c, lambda_r) )

# for twk in range (tinit, 1344, 168):
for twk in range (tinit, 672, 168):
    
    for tin in range (twk, twk + (3 * 24), 48):
      
        y0_ = sol[-1, :]
        y0 = y0_ * surviving_fraction(a_sen, b_sen, dos_weekday)
        t = np.linspace(tin, tin + 48, tin + 49)
        sol = odeint( tumor_PC3, y0, t, args = ( cprolif_c, cprolif_r, cCapacity_c, cCapacity_r, lambda_c, lambda_r) )
        t = np.append( tin, t )
        sol = np.append( [y0_], sol, 0 )        
        plt.plot(t, sol[:, 0], 'b', label='sensitive')
        
        if tin >= twk + 48:
            long = sol[0, 0]
            area_sen += 24/2 * (short + long)
            count += 1 
         
        short = sol[1, 0]
  
#Friday Dose:
            
    y0_ = sol[-1, :]
    y0 = y0_ * surviving_fraction(a_sen, b_sen, dos_friday)
    t = np.linspace(tin + 48, tin + 48, tin + 49)
    sol = odeint( tumor_PC3, y0, t, args = ( cprolif_c, cprolif_r, cCapacity_c, cCapacity_r, lambda_c, lambda_r) )
    t = np.append( tin + 48, t )
    sol = np.append( [y0_], sol, 0 )        
    plt.plot(t, sol[:, 0], 'b', label='sensitive')
    
    long = sol[0,0]
    area_sen += 48/2 * (short + long)
    count += 1
    

### Weekend Period
    
    t = np.linspace(tin + 48, tin + 120, tin + 121)
    sol = odeint( tumor_PC3, sol[-1,:], t, args = ( cprolif_c, cprolif_r, cCapacity_c, cCapacity_r, lambda_c, lambda_r) )
    plt.plot(t, sol[:, 0], 'b', label='sensitive')

    count += 1
    short = sol[0,0]
    print(short)
    long = sol[-1, 0]
    area_sen += 72/2 * (short + long)

print(area_sen)
# print(sol[1343, 0])

     
    

    