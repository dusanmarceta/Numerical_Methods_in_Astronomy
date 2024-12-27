import numpy as np
import matplotlib.pyplot as plt

def year2second(x):
    return x*31557600.0

def day2second(x):
    return x*86400

def au2m(x):
    return x*au

def m2au(x):
    return x/au

gm=1.32712440018e20
au=149597870700.

#pocetni uslovi
x=[au2m(1)]
y=[0]
z=[0]

vx=[0.]
vy=[4e4]
vz=[0.]

# parametri integracije
t=year2second(20)
dt=day2second(1)

tt=0

while tt<t:
    
    r=np.sqrt(x[-1]**2+y[-1]**2+z[-1]**2) #trenutno heliocentrično rastojanje
    
    # komponente ubrzanja
    ax=-gm/r**3*x[-1]
    ay=-gm/r**3*y[-1]
    az=-gm/r**3*z[-1]
    
    # koordinate
    x.append(x[-1]+vx[-1]*dt)
    y.append(y[-1]+vy[-1]*dt)
    z.append(z[-1]+vz[-1]*dt)
    
    # brzine
    vx.append(vx[-1]+ax*dt)
    vy.append(vy[-1]+ay*dt)
    vz.append(vz[-1]+az*dt)
    
    # vreme
    tt=tt+dt
    
plt.plot(m2au(np.array(x)),m2au(np.array(y)),'.')
plt.axis('equal')
    