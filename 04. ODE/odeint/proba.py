from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt



def sistem(u,t):

    return [u[2],
            u[3],
            0,
            -9.81]
    
u0=[0,0,50,50]
t=np.linspace(0,100,1001) 
u = odeint(sistem, u0, t)
u=np.transpose(u)

x=u[0]
y=u[1]

plt.plot(x,y)
#plt.ylim([0,300])


plt.axis([0, 600, 0, max(y)]) 