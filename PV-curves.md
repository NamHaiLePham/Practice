# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 10:28:08 2021

@author: Hai
"""

import matplotlib.pyplot as plt
import numpy as np
import math

#This is code for two bus system in assign 4
#First draw the PV curve and find the peak of curve.
#Vr can be found by geting value from i to Vs and step=0.01
#This peak is Pcrit 


Z=133.1
#the equivalent reactance, in assignment 4 
#Xeq= Xt+ Xl1//Xl2, Xeq=133.1ohm
Vs=220
#the voltage after transformer, in assignment 4 220kV
PF=float(input('Input PF:'))
k=int(input('Inductive PF:0, Capacitive PF:1, choose:'))
i=1
m=[]
n=[]
#Plot PV curve
while i<Vs/10-1:
    i=i+0.01 #step
    Vr=i*10 #set random Vr from i to Vs/10
    cosphi=Vr*PF/Vs
    phi1=math.acos(cosphi)
    if k==1:#capacitive
     phi2=phi1+math.acos(PF)
    else:#inductive
     phi2=phi1-math.acos(PF)
     
    t=math.sin(phi2)
    P=Vs*Vr*t/Z
    Vr=round(Vr,2)
    P=round(P,2)
    m.append(Vr)
    n.append(P)
    q=max(n)#Pcrit
    for s in range(len(n)):
        if n[s]-q==0:
            a=s
    

#Calculate Qsending
v=q*Z/(Vs*m[a])
phicri=math.asin(v)
Qs=Vs**2/Z-(Vs*m[a]*math.cos(phicri))/Z
Qs=round(Qs,2)

#Calculate Qreceiving
Qr=-m[a]**2/Z+(Vs*m[a]*math.cos(phicri))/Z
Qr=round(Qr,2)

print('the critical voltage is:',m[a])
print('the critical power is',q) 
print('the Q sending is',Qs)  
print('the Q receiving is',Qr) 
plt.plot(n,m)
plt.xlabel('P[kW]')
plt.ylabel('Vr[kV]')
plt.grid()
plt.title(label="PV curve", fontsize=15, color="red")



 

    
    
    


