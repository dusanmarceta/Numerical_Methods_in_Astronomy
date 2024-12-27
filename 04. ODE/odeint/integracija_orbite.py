from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

mu=1.32712440018e20
AJ=1.5e11


def sistem(u,t):
    r=np.sqrt(u[0]**2+u[1]**2)
    return [u[2],
            u[3],
            -mu*u[0]/r**3,
            -mu*u[1]/r**3]
    
u0=[AJ,0,0,4.2e4]
t=np.linspace(0,1e10,10000) 
u = odeint(sistem, u0, t)
u=np.transpose(u)

x=u[0]
y=u[1]

plt.plot(x/AJ,y/AJ)
plt.axis('equal')
plt.grid()

r=np.sqrt(u0[0]**2+u0[1]**2)
v=np.sqrt(u0[2]**2+u0[3]**2)

a=r*mu/(2*mu-r*v**2)/AJ
#plt.ylim([0,300])


#plt.axis([0, 600, 0, max(y)]) 