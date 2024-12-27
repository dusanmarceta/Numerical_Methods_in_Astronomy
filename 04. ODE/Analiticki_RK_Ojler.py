import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time

def year2second(x):
    return x*31557600.0

def day2second(x):
    return x*86400

def au2m(x):
    return x*au

def m2au(x):
    return x/au

def kepler(e,M,tacnost):
    delta=2*tacnost
    E=M
    while np.abs(delta)>tacnost:
        f=E-e*np.sin(E)-M
        fprim=1-e*np.cos(E)
        delta=f/fprim
        E=E-delta
    return(E)

gm=1.32712440018e20
au=149597870700.


Ms = 1.989e30 # masa Sunca
gama = 6.67e-11 # gravitaciona konstanta

x=au
y=0
z=0
vx=0
vy=4e4


a=1/(2/x-vy**2/Ms/gama)
e=1-x/a

f=np.linspace(0,2*np.pi,1000)

r=a*(1-e**2)/(1+e*np.cos(f))

X=r*np.cos(f)
Y=r*np.sin(f)

plt.plot(X/au,Y/au,label='analiticki')
plt.grid()
plt.axis('equal')

T=np.sqrt(4*np.pi**2*a**3/Ms/gama)

t_lim=5*T
dt = day2second(10.) # korak integracije




def ubrzanje(x,y, Ms, gama):
    #x,y,z - helicentricne koordinate objekta

    # Ms - masa Sunca
    # gama - gravitaciona konstanta 
    r=(x ** 2 + y ** 2)**(1/2) # helicenticni radijus vektor objekta
   
    ax = -gama * (Ms / r**3 * x)
    ay = -gama * (Ms / r**3 * y)
 
    return ax, ay



# Runge-Kutta
pocetak=time.time()


x=au
y=0
z=0
vx=0
vy=4e4

t = 0. # početno vreme

x_rk=[]
y_rk=[]


while t < t_lim:

    xx=x
    vxx=vx
    yy=y
    vyy=vy

    
    # 1. korak
    Kx_1=dt*vxx
    Kvx_1=dt*ubrzanje(xx, yy, Ms, gama)[0]
    Ky_1=dt*vyy
    Kvy_1=dt*ubrzanje(xx, yy, Ms, gama)[1]

    
    xx=x+Kx_1/2
    vxx=vx+Kvx_1/2
    yy=y+Ky_1/2
    vyy=vy+Kvy_1/2

    
    # 2. korak
    Kx_2=dt*vxx
    Kvx_2=dt*ubrzanje(xx, yy, Ms, gama)[0]
    Ky_2=dt*vyy
    Kvy_2=dt*ubrzanje(xx, yy, Ms, gama)[1]

    
    xx=x+Kx_2/2
    vxx=vx+Kvx_2/2
    yy=y+Ky_2/2
    vyy=vy+Kvy_2/2

    
    # 3. korak
    Kx_3=dt*vxx
    Kvx_3=dt*ubrzanje(xx, yy, Ms, gama)[0]
    Ky_3=dt*vyy
    Kvy_3=dt*ubrzanje(xx, yy, Ms, gama)[1]

    
    xx=x+Kx_3
    vxx=vx+Kvx_3
    yy=y+Ky_3
    vyy=vy+Kvy_3

    
    # 4. korak
    Kx_4=dt*vxx
    Kvx_4=dt*ubrzanje(xx, yy, Ms, gama)[0]
    Ky_4=dt*vyy
    Kvy_4=dt*ubrzanje(xx, yy, Ms, gama)[1]


    # racunanje koeficijenata K
    Kx=1/6*(Kx_1+2*Kx_2+2*Kx_3+Kx_4)
    Kvx=1/6*(Kvx_1+2*Kvx_2+2*Kvx_3+Kvx_4)
    
    Ky=1/6*(Ky_1+2*Ky_2+2*Ky_3+Ky_4)
    Kvy=1/6*(Kvy_1+2*Kvy_2+2*Kvy_3+Kvy_4)
  
    x=x+Kx
    vx=vx+Kvx
    
    y=y+Ky
    vy=vy+Kvy
    
    x_rk.append(x)
    y_rk.append(y)
        
    t+=dt

vreme_rk=time.time()-pocetak
    
plt.plot(np.array(x_rk)/au,np.array(y_rk)/au,label='RK') 


# Ojler
dt=dt/100
pocetak=time.time()
x=au
y=0
vx=0
vy=4e4


x_o=[]
y_o=[]


t=0
while t < t_lim:
    
  
    x+=dt*vx
    y+=dt*vy
    vx+=dt*ubrzanje(x, y, Ms, gama)[0]
    vy+=dt*ubrzanje(x, y, Ms, gama)[1]
   
    x_o.append(x)
    y_o.append(y)

        
    t+=dt
 
vreme_o=time.time()-pocetak
plt.plot(np.array(x_o)/au,np.array(y_o)/au,label='Ojler') 


# Ojler-Kromer
pocetak=time.time()
x=au
y=0
vx=0
vy=4e4


x_ok=[]
y_ok=[]
z_ok=[]

t=0
while t < t_lim:
    
   
    vx+=dt*ubrzanje(x, y, Ms, gama)[0]
    vy+=dt*ubrzanje(x, y, Ms, gama)[1]
    x+=dt*vx
    y+=dt*vy

    x_ok.append(x)
    y_ok.append(y)
    z_ok.append(z)
        
    t+=dt

vreme_ok=time.time()-pocetak
 
plt.plot(np.array(x_ok)/au,np.array(y_ok)/au,label='Ojler-Kromer') 
plt.legend()  
 

n=np.sqrt(Ms*gama/a**3)
x_a=[]
y_a=[]
t=0
while t < t_lim:
    
    M=t*n
    E=kepler(e,M,1e-6)
    x_a.append(a*(np.cos(E)-e))
    y_a.append(a*np.sqrt(1-e**2)*np.sin(E))  
    t+=dt
    
x_rk=np.array(x_rk)
y_rk=np.array(y_rk)
x_o=np.array(x_o)
y_o=np.array(y_o)
x_ok=np.array(x_ok)
y_ok=np.array(y_ok)
x_a=np.array(x_a)
y_a=np.array(y_a)

delta_rk=np.sqrt((x_a-x_rk)**2+(y_a-y_rk)**2)
delta_o=np.sqrt((x_a-x_o)**2+(y_a-y_o)**2)
delta_ok=np.sqrt((x_a-x_ok)**2+(y_a-y_ok)**2)

plt.figure()
plt.plot(delta_rk/au,label='RK')
plt.plot(delta_o/au,label='Ojler')
plt.plot(delta_ok/au,label='Ojler-Kromer')
plt.legend()

plt.figure()
plt.plot(delta_o/delta_rk, label='Ojler/RK')
plt.plot(delta_ok/delta_rk, label='Ojler-Kromer/RK')
plt.legend()