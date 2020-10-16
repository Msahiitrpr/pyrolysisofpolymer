# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import parameters
#Time Properties
t=60*parameters.t
dt=parameters.dt
n=int(round(t/dt))
R=8.314/1000.
#=========================
#Kinetic properties
A1=parameters.A1
A2=parameters.A2
A3=parameters.A3
A4=parameters.A4
A5=parameters.A5
EA1=parameters.EA1
EA2=parameters.EA2
EA3=parameters.EA3
EA4=parameters.EA4
EA5=parameters.EA5
#===========================
#Generating arrays
xp=np.zeros(n+1)
xh=np.zeros(n+1)
xm=np.zeros(n+1)
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
Temp=np.zeros(n+1)
qacc=0
#============================
#Reaction heats
H1=parameters.H1
H2=parameters.H2
H3=parameters.H3
H4=parameters.H4
H5=parameters.H5
#========================= fit data
#Polymer	
xpol = [0.425532,	59.5745,	90.2128,	120,	150.213,	200]
ypol =[	0.997199,	0.806723,	0.742297,	0.593838,	0.568627,	0.470588]
						
#Heavy Fraction	
xhfn= [ 0,	        60,	        90.2128,	120,	    150.638,	200]
yhfn=[	0.00287264,	0.112957,	0.18901,	0.29587,	0.366327,	0.473461]
						
#Middle Distillates	
xmd =[0         	,60	       ,89.7872,	120,	     149.787,	200]
ymd	=[0.00287264,	0.0569343,	0.0181358,	0.0689791,	0.0497884,	0.0532928]
						
#Light fraction	
xlf =[0	            ,59.5745	  ,89.7872  	,120	 ,149.787  	    ,200]
ylf	=[0.00287264	,0.026116	,0.0517492	,0.0465701	,0.0273795, 	0.0252816]



#==========================
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
    xp[i+1]=(-xp[i]*(K1+K2+K3))*dt+xp[i]
    xh[i+1]=(xp[i]*K1-xh[i]*(K4+K5))*dt +xh[i]
    xm[i+1]=(xp[i]*K2 + xh[i]*K4)*dt +xm[i]
    xl[i+1]=(xp[i]*K3 + xh[i]*K5)*dt + xl[i]
    time[i+1]=i*dt
    hrx_1[i]=H1*K1*xp[i]
    hrx_2[i]=H2*K2*xp[i]
    hrx_3[i]=H3*K3*xp[i]
    hrx_4[i]=H4*K4*xh[i]
    hrx_5[i]=H5*K5*xh[i]
    hrx[i]=hrx_1[i]+hrx_2[i]+hrx_3[i]+hrx_4[i]+hrx_5[i]
    if Tmp==1:# Constant temperature
        T[i+1]=T[i]
        q[i+1]=hrx[i]
    if Tmp==2:# Temperature step
            T[i+1]=T[i]+Tstep
            q[i+1]=hrx[i]
time=time/60
plt.plot(time,xp, label="Polymer")
plt.plot(time,xh, label="Heavy product")
plt.plot(time,xm, label="Medium product")
plt.plot(time,xl, label="Light product")
plt.title(T[0]-273.15)
plt.xlabel("Time (min)")
plt.ylabel("Xi")
plt.legend(loc="upper left")
plt.scatter(xpol,ypol) 
plt.scatter(xhfn ,yhfn)
plt.scatter(xmd,ymd)
plt.scatter( xlf , ylf )

p1=plt.show()
plt.plot(time,hrx_1, label="1")
plt.plot(time,hrx_2, label="2")
plt.plot(time,hrx_3, label="3")
plt.plot(time,hrx_4, label="4")
plt.plot(time,hrx_5, label="5")
plt.plot(time,q, label="q(total)")
plt.legend(loc="upper left")
plt.xlabel("Time (min)")
plt.ylabel("kJ")
plt.show()
rxkwh= (np.sum(hrx_1)/(3600)+ np.sum(hrx_2)/(3600) +np.sum(hrx_3)/(3600)+np.sum(hrx_4)/(3600) + np.sum(hrx_5)/(3600))*dt
pr="total energy usage " +str(rxkwh) +"kWh"
print(pr)
