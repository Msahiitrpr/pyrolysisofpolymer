# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 21:54:11 2020

@author: Manmohan
"""

#Simulation time
t=15 #min
dt=0.01 #s
#Kinetic Parameters A0=[1/s] Ea= [kJ/mol K]
A1=4.21E+16/60.0
A2=3.34E+13/60.0
A3=1.36E+18/60.0
A4=2.36E+11/60.0
A5=1.457244265/60.0
A6=1.32E-64/60.0
EA1=214.5900427
EA2=184.4210617
EA3=244.778692
EA4=160.4760931
EA5=35.17583039
EA6=-622.4955284
#Heat of reaction for the different reactions (Are not required, but can be use for total heat of reaction estimations)
H1=0
H2=0
H3=0
H4=0
H5=0
H6=0
#Temperature options
T= 400.0 #Start temperature
Tr= 1 # Temperature regime 1=Isothermal 2=Ramped
Rate=15 #Temperature rate if ramped Temperature.