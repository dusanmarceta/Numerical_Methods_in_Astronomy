from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

M=5.972e24 #masa Zemlje
G=6.67e-11 #gravitaciona konstanta


# sistem diferencijalnih jednacina
def sistem(u,t):
    r=np.sqrt(u[0]**2+u[1]**2)
    return [u[2],
            u[3],
            -M*G*u[0]/r**3,
            -M*G*u[1]/r**3]

#pocetni uslovi
x0=7e6
y0=0
vx0=0
vy0=8000
    
u0=[x0,y0,vx0,vy0]

#vreme integracije
t=np.linspace(0,1e5,10000) 

# resavanje sistema
u = odeint(sistem, u0, t)
u=np.transpose(u)

# koordinate satelita
x=u[0]
y=u[1]
#brzine satelita
vx=u[2]
vy=u[3]

plt.plot(x,y)
plt.axis('equal')
plt.grid()
