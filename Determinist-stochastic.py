# -*- coding: utf-8 -*-
"""
Created on Sat Mar  4 13:14:37 2023

@author: Isabel
"""
import numpy as np
import random
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import odeint



#tiempo final
tf=35

#Parameter values
fm=1.0
fp=1.0
k=1000
K=602
nequ=2
alpham=10100*fm
deltam=1.0*fm
alphap=10.0*fp
deltap=0.1*fp
V=0.6022


#Ecuaciones deterministas
#initial conditions
y0 =  [] 
y0.append(0.0)
y0.append(0.0)
#y0.append(alphap*alpham/(deltap*deltam))


#Derivatives
def dy_dt(y,t):    
    dy = []
    dy.append(alpham*(k**2)/(k**2+y[1]**2)-deltam*y[0])
    dy.append(alphap*y[0]-deltap*y[1])
   
    return dy    

#Integration of ODEs
t = np.linspace(start=0, stop=50,num=1000)
prot = np.linspace(start=0, stop=15000,num=1000)
sol = odeint(func=dy_dt, y0=y0, t=t)
sol = pd.DataFrame(sol, columns=["mRNA","protein"])

sol.index = t
#sol.index = sol.index/60


N=np.linspace(0, 1, 2)



for i in N:
    M=0
    P=0

    t=0
    M_array=[]
    P_array=[]
    t_array=[]
    while t<tf:
        
        
        a1=V*alpham*(K**2)/(K**2+P**2)
        a2=deltam*M
        a3=alphap*M
        a4=deltap*P
        
        a0=a1+a2+a3+a4
        
        r1=np.random.uniform(0,1)
        
        
        M_array.append(M/V)
        P_array.append(P/V)
        t_array.append(t)
            
        tau=-np.log(r1)/a0
        
        t=t+tau
        
        r2=np.random.uniform(0,a0)
        
        if r2<=a1:
            M=M+1
        elif a1<r2<=(a1+a2):
            M=M-1
        elif (a1+a2)<r2<=(a1+a2+a3):
            P=P+1
        elif (a1+a2+a3)<r2<=(a1+a2+a3+a4):
            P=P-1
        
    if i==0:
        M_array1=M_array
        P_array1=P_array
        t_array1=t_array
    elif i==1:
        M_array2=M_array
        P_array2=P_array
        t_array2=t_array


print(sol.mRNA)
print(sol.protein)       
        
#Plotting
fig, ax = plt.subplots(3,figsize=(8,6))


ax[0].plot(sol.index,sol.mRNA, label="mRNA (nM)")
ax[0].set_xlabel("time (min)")
ax[0].set_ylabel("concentration (nM)")
ax[0].set_title("Deterministic dynamics")
ax[0].legend()

ax[1].plot(t_array1,M_array1, label="mRNA (nM)")
ax[1].set_xlabel("time (min)")
ax[1].set_ylabel("concentration (nM)")
ax[1].set_title("Stochastic dynamics")
ax[1].legend()

ax[2].plot(t_array2,M_array2, label="mRNA (nM)")
ax[2].set_xlabel("time (min)")
ax[2].set_ylabel("concentration (nM)")
ax[2].set_title("Stochastic dynamics")
ax[2].legend()


# ax[0].plot(sol.index,sol.protein, label="P (nM)")
# ax[0].set_xlabel("time (min)")
# ax[0].set_ylabel("concentration (nM)")
# ax[0].set_title("Deterministic dynamics")
# ax[0].legend()


# ax[1].plot(t_array1,P_array1, label="P (nM)")
# ax[1].set_xlabel("time (min)")
# ax[1].set_ylabel("concentration (nM)")
# ax[1].set_title("Stochastic dynamics")
# ax[1].legend()

# ax[2].plot(t_array2,P_array2, label="P (nM)")
# ax[2].set_xlabel("time (min)")
# ax[2].set_ylabel("concentration (nM)")
# ax[2].set_title("Stochastic dynamics")
# ax[2].legend()

# fig, ax = plt.subplots(2,figsize=(8,6))

# ax[0].plot(sol.protein,sol.mRNA, c="violet", label="protein(nM)")
# ax[0].set_xlabel("protein (nM)")
# ax[0].set_ylabel("mRNA (nM)")
# ax[0].set_title("Deterministic phase portrait")

# ax[1].plot(P_array1,M_array1, c="violet", label="protein(nM)")
# ax[1].set_xlabel("protein (nM)")
# ax[1].set_ylabel("mRNA (nM)")
# ax[1].set_title("Stochastic phase portrait")


fig.subplots_adjust(top=0.8)
fig.align_ylabels()
plt.tight_layout()
plt.show()




