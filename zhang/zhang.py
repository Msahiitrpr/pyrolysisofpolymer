# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 23:34:25 2020
 
@author: Avinash
"""
 
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import parameters
#Time Properties
t=60*parameters.t
dt=0.01
n=int(round(t/dt))
R=8.314/1000.
#=========================
#Kinetic properties
A1=parameters.A1
A2=parameters.A2
A3=parameters.A3
A4=parameters.A4
A5=parameters.A5
A6=parameters.A5
EA1=parameters.EA1
EA2=parameters.EA2
EA3=parameters.EA3
EA4=parameters.EA4
EA5=parameters.EA5
EA6=parameters.EA5
#===========================
#Generating arrays
xp=np.zeros(n+1)
xw=np.zeros(n+1)
xg=np.zeros(n+1)
xl=np.zeros(n+1)
time=np.zeros(n+1)
T=np.zeros(n+1)
hrx=np.zeros(n+1)
q=np.zeros(n+1)
hrx_1=np.zeros(n+1)
hrx_2=np.zeros(n+1)
hrx_3=np.zeros(n+1)
hrx_4=np.zeros(n+1)
hrx_5=np.zeros(n+1)
hrx_6=np.zeros(n+1)
Temp=np.zeros(n+1)
qacc=0
#============================
#Reaction heats
H1=parameters.H1
H2=parameters.H2
H3=parameters.H3
H4=parameters.H4
H5=parameters.H5
H6=parameters.H6
#=========================
T[0]=parameters.T+273.15
mass=1 #kg
xp[0]=1*mass
#Temperature regime
Tmp=parameters.Tr # Temperature regime 1=Isothermal 2=Temperature rate change
Tstep=(parameters.Rate*dt)/60.0 #C/min
for i in range (0,n):
    K1=A1*np.exp(-EA1/(T[i]*R))
    K2=A2*np.exp(-EA2/(T[i]*R))
    K3=A3*np.exp(-EA3/(T[i]*R))
    K4=A4*np.exp(-EA4/(T[i]*R))
    K5=A5*np.exp(-EA5/(T[i]*R))
    K6=A6*np.exp(-EA6/(T[i]*R))
    xp[i+1]=(-xp[i]*(K1+K2+K3))*dt+xp[i]
    xw[i+1]=(xp[i]*K1-xw[i]*(K4+K5))*dt +xw[i]
    xl[i+1]=(xp[i]*K2 + xw[i]*K4 -xl[i]*K6)*dt +xl[i]
    xg[i+1]=(xp[i]*K3 + xw[i]*K5 +xl[i]*K6)*dt + xg[i]
    time[i+1]=i*dt
    hrx_1[i+1]=H1*K1*xp[i]
    hrx_2[i+1]=H2*K2*xp[i]
    hrx_3[i+1]=H3*K3*xp[i]
    hrx_4[i+1]=H4*K4*xw[i]
    hrx_5[i+1]=H5*K5*xw[i]
    hrx_6[i+1]=H6*K6*xl[i]
    hrx[i+1]=hrx_1[i+1]+hrx_2[i+1]+hrx_3[i+1]+hrx_4[i+1]+hrx_5[i+1]+hrx_6[i+1]
    if Tmp==1: # Constant temperature
            T[i+1]=T[i]
            q[i+1]=hrx[i]
    if Tmp==2: # Temperature step
            T[i+1]=T[i]+Tstep
            q[i+1]=hrx[i]
time=time/60
#time=time/3
plt.plot(time,xp, label="Polymer")
plt.plot(time,xw, label="Wax")
plt.plot(time,xl, label="Liquid")
plt.plot(time,xg, label="Gas")
plt.xlabel("Time (min)")
plt.ylabel("Xi")
plt.legend(loc="upper left")
plt.show()
plt.plot(time,hrx_1, label="1")
plt.plot(time,hrx_2, label="2")
plt.plot(time,hrx_3, label="3")
plt.plot(time,hrx_4, label="4")
plt.plot(time,hrx_5, label="5")
plt.plot(time,hrx_6, label="6")
plt.plot(time,q, label="q(total)")
plt.legend(loc="upper left")
plt.xlabel("Time (min)")
plt.ylabel("kJ")
p3=plt.show()
plt.plot(time,T-273.15, label="Temperature")
plt.xlabel("Time (min)")
plt.ylabel("T (C)")
plt.show()
totkwh=(np.sum(q)/(3600))*dt
rxkwh= (np.sum(hrx_1)/(3600)+ np.sum(hrx_2)/(3600) +
np.sum(hrx_3)/(3600)+np.sum(hrx_4)/(3600) + np.sum(hrx_5)/(3600))*dt
pr="total energy usage " +str(totkwh)+ "kWh, used by reactions " + str(rxkwh) + " kWh, for heating" +str(totkwh-rxkwh) +"kWh"
print(pr)