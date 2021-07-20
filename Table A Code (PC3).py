# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 23:12:59 2020

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
# y0 = [0.9, 0.1]

sol = odeint( tumor_PC3, y0, t, args = ( cprolif_c, cprolif_r, cCapacity_c, cCapacity_r, lambda_c, lambda_r) )

plt.title("PC3")
plt.figure(1)
plt.plot(t, sol[:, 0], 'b', label='parental')
plt.plot(t, sol[:, 1], 'g', label='resistant')
plt.xlabel('Time (hours)')
plt.ylabel('Volume')
plt.legend(loc='best')
plt.grid()

a_res = 0.300
b_res = 0.0402

a_sen = 0.430
b_sen = 0.0407

tin = 0

tin = 336
val = 24 * int(input("Enter time interval (days): "))
# val = 24 * 6
dosage = 0
dos = int(input("Enter dosage: "))
area_sen = 0
area_res = 0
var = int(60/dos)

curve_sen = (336/2) * (0.5 + sol[336,0])

for dos_count in range (0, var, 1):
        
    y0_ = sol[-1, :]
    y0 = y0_ * surviving_fraction(a_sen, b_sen, dos)
    t = np.linspace(tin, tin + val, tin + val + 1)
    sol = odeint( tumor_PC3, y0, t, args = ( cprolif_c, cprolif_r, cCapacity_c, cCapacity_r, lambda_c, lambda_r) )
    t = np.append( tin, t )
    sol = np.append( [y0_], sol, 0 )        
    plt.plot(t, sol[:, 0], 'b', label='sensitive')
    
    dosage += dos
    
    if tin >= 336 + val:
        long = sol[0, 0]
        area_sen += (val/2) * (short + long)
        if dosage == 60:
            final_vol = sol[tin + val, 0]
        elif tin >= 1344:
            final_vol = sol[1344,0]
            break
    
    short = sol[1,0]
    
    tin += val
    
area_sen += curve_sen

# print('Final Volume of Sensitive Tumor = ', final_vol)
print(final_vol)
# print('Total Area Under Sensitive Tumor Curve = ', area_sen)
print(area_sen)
    
tin = 336
dosage = 0

t = np.linspace(0, 336, 337)
y0 = [0.5, 0.5] 
sol = odeint( tumor_PC3, y0, t, args = ( cprolif_c, cprolif_r, cCapacity_c, cCapacity_r, lambda_c, lambda_r) )

curve_res = (336/2) * (0.5 + sol[336,1])

for dos_count in range (0, 60, dos):
        
    y0_ = sol[-1, :]
    y0 = y0_ * surviving_fraction(a_res, b_res, dos)
    t = np.linspace(tin, tin + val, tin + val + 1)
    sol = odeint( tumor_PC3, y0, t, args = ( cprolif_c, cprolif_r, cCapacity_c, cCapacity_r, lambda_c, lambda_r) )
    t = np.append( tin, t )
    sol = np.append( [y0_], sol, 0 )        
    plt.plot(t, sol[:, 1], 'g', label='resistant')
    
    dosage += dos
    
    if tin >= 336 + val:
        long = sol[0, 1]
        area_res += (val/2) * (short + long)
        if dosage == 60:
            final_vol = sol[tin + val, 1]
        elif tin >= 1344:
            final_vol = sol[1344, 1]
            break
        
    short = sol[1, 1]
    
    tin += val
    
area_res += curve_res

# print('Final Volume of Resistant Tumor = ', final_vol)
print(final_vol)
# print('Total Area Under Resistant Tumor Curve = ', area_res)
print(area_res)

    