# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 12:54:08 2023

@author: Isabel
"""
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import math
R=np.loadtxt('Proteinas.txt')
#print(R)



mu=10000
sigma=math.sqrt(mu)
a=np.linspace(mu-3*sigma,mu+3*sigma,100)
plt.plot(a,stats.norm.pdf(a,mu,sigma),color="b")
plt.hist(R,bins=10,density=True,color="lightblue")
plt.xlabel("p(nM)")
plt.ylabel("Probab. Density")
plt.show()

