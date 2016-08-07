#!/usr/bin/python
# coding: utf-8

import locale
locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

import sys
from string import *
import math
import numpy as np

N = 11			# wielkosc siatki
z = 3				# wielkosc rozpatrywanej siatki
d = 0.001			# wielkosc komorki
theta = 0
deltami = 0.2
mi0 = 0.01
ms = 4
L = 1930000000			# cieplo topnienia 
Teq = 1811			# temperatura rownowagowa 
gamma = 0.00000019
ro_L = 7000.0			# gestosc fazy cieklej
ro_S = 7250.0			# gestosc fazy stalej
cw_L = 5740000.0		# cieplo wlasciwe fazy cieklej poprzednia wartosc - 574
cw_S = 5730000.0		# cieplo wlasciwe fazy stalej poprzednia wartosc - 573
lambdaL = 35.0			# lambda fazy cieklej
lambdaS = 33.0			# lambda fazy stalej
max_fs_change = 15
omega = 0.8
X1 = 0.01			#  wspolczynniki do kroku czasowego
X2 = 0.1	

first = N / 2 - z / 2 
last = N / 2 + z / 2 
a_Solid = a_Liquid = 0 
current_step = eksp = 0
iks = once = iks_2 = 1, 
con = wait = 0
fs_change = 1 
nr_of_iter = 0
#  u - predkosc, t - czas, x - wielkosc komorki, theta - kat 
# // K - krzywa, R - wspl. przewodzenia *					
czas = 0
u = t = dx = 0	
mi = Tmi = Tk = Tuc =  K = R =  hp = cp = ro =  dt = dt_0 =  u1 =  dTeq =  dTk =  max_current_fs = nb_of_steps = max_T_change = p2 = save_step = dT = neighbourhood = p1 =  duration = count_dt = 0

'''
start_time 
current_time;
'''
a_coeff_struct = [('a_P','f4'),('a_W_ave','f4'),('a_E_ave','f4'),('a_N_ave','f4'),('a_S_ave','f4')]
cell_struct = [ ('state','f4'), ('temp','f4'), ('ent','f4'), ('fs','f4'), ('a_coeff',a_coeff_struct), ('number','i2'), ('K','f4')]

current = np.zeros((N,N),dtype= cell_struct)
previous = np.zeros((N,N),dtype= cell_struct)


save_file = open('temp.txt','w')
np.savetxt(save_file,current['temp'],fmt='%s',delimiter='\t',newline='\n')

num = 0


zapis_state = open("stan.txt","w")
zapis_temp = open("temperatura.txt","w")
zapis_dane = open("data.txt","w")
zapis_dane_2 = open("data_2.txt","w")
zapis_K = open("wartosci_K.txt", "w")
input_ = open("input.txt", "r")

input_var = []
for line in open('input.txt','r'):
    input_var.append(line.strip().split('\t')[1])
nb_of_steps, max_T_change, dx, save_step, p1, p2, dT, neighbourhood = [ float(x) for x in input_var ]

dt = dx * max_fs_change / mi0 / dT ;
dt_0 = dt;

mi = mi0/(1 - deltami * math.cos( ms * theta));		

def set_values(some_array, start, end, value):
    some_array[ start : end + 1, start : end + 1 ] = value
    
def seeding(array,method):
    if method == 1:
        array[ N/2, N/2 ] = 100
        array[ N/2, N/2 - 1 ] = 100
        array[ N/2, N/2 + 1 ] = 100
        array[ N/2 - 1, N/2 ] = 100
        array[ N/2 + 1, N/2 ] = 100

def przyrost_fazy_stalej(cos):



#boundary conditions - niepotrzebne, bo i tak pozniej wszystko ustawiam na ta sama wartosc
set_values(current['temp'],first-z/2,last+z/2,Teq-dT)
set_values(current['temp'],first,last,0)
print current['temp']

#initial conditions
set_values(previous['temp'],first,last,Teq-dT)
set_values(previous['state'],first,last,0)

#seeding
seeding(previous['state'],1)
print previous['temp']
print previous['state']

#zerowanie - niepotrzebne
set_values(previous['fs'],first,last,0)
print previous['fs']a

ma[1:4,1:4] = ma [1:4,1:4] + ma[0:3,1:4] + ma[2:5,1:4] + ma[1:4,0:3] + ma[1:4,2:5]
