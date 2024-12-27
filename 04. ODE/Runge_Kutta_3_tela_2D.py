import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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


Ms = 1.989e30 # masa Sunca
Mz=5.972e24 # masa Zemlje
#Mz=0 # masa Zemlje
gama = 6.67e-11 # gravitaciona konstanta

r0z = 6.378e6  # poluprecnik Zemlje
v_zemlja = 2 * np.pi * au2m(1.) / (365.25 * 86400)

v0 = 1e3 # pocetna brzina u odnosu na Zemlju





def ubrzanje(x,y,z,xp,yp,zp,Ms,Mp, gama):
    #x,y,z - helicentricne koordinate objekta
    #xp,yp,zp - heliocentricne koordinate planete
    # Ms - masa Sunca
    # Mp - masa planete
    # gama - gravitaciona konstanta
    
    x_pc=x-xp
    y_pc=y-yp
    z_pc=z-zp
    
    r=(x ** 2 + y ** 2 + z ** 2)**(1/2) # helicenticni radijus vektor objekta
    r_pc=(x_pc ** 2 + y_pc ** 2 + z_pc** 2)**(1/2) # planetocentricni radijus vektor objekta


    ax = -gama * (Ms / r**3 * x + Mp / r_pc**3 * x_pc)
    ay = -gama * (Ms / r**3 * y + Mp / r_pc**3 * y_pc)
    az = -gama * (Ms / r**3 * z + Mp / r_pc**3 * z_pc)
    
    return ax, ay, az


# pocetni uslovi
x=au2m(1.) + 50 * r0z
y=0
z=0
vx=0
vy=v_zemlja+v0
vz=v0/2


# koordinate Zemlje
x_zemlja = au2m(1.)
y_zemlja = 0
z_zemlja = 0

# ugaona brzina kretanja Zemlje oko Sunca
sk = 2 * np.pi / year2second(1.) # obiđe 2*pi za jednu godinu


dt = 60. # korak integracije
t = 0. # početno vreme

X=[]
Y=[]
Z=[]

while t < year2second(1.):

    xx=x
    vxx=vx
    yy=y
    vyy=vy
    zz=z
    vzz=vz
    
    # 1. korak
    Kx_1=dt*vxx
    Kvx_1=dt*ubrzanje(xx,yy,zz,x_zemlja, y_zemlja, 0, Ms, Mz, gama)[0]
    Ky_1=dt*vyy
    Kvy_1=dt*ubrzanje(xx,yy,zz,x_zemlja, y_zemlja, 0, Ms, Mz, gama)[1]
    Kz_1=dt*vzz
    Kvz_1=dt*ubrzanje(xx,yy,zz,x_zemlja, y_zemlja, 0, Ms, Mz, gama)[2]
    
    xx=x+Kx_1/2
    vxx=vx+Kvx_1/2
    yy=y+Ky_1/2
    vyy=vy+Kvy_1/2
    zz=z+Kz_1/2
    vzz=vz+Kvz_1/2
    
    # 2. korak
    Kx_2=dt*vxx
    Kvx_2=dt*ubrzanje(xx,yy,zz,x_zemlja, y_zemlja, 0, Ms, Mz, gama)[0]
    Ky_2=dt*vyy
    Kvy_2=dt*ubrzanje(xx,yy,zz,x_zemlja, y_zemlja, 0, Ms, Mz, gama)[1]
    Kz_2=dt*vzz
    Kvz_2=dt*ubrzanje(xx,yy,zz,x_zemlja, y_zemlja, 0, Ms, Mz, gama)[2]
    
    xx=x+Kx_2/2
    vxx=vx+Kvx_2/2
    yy=y+Ky_2/2
    vyy=vy+Kvy_2/2
    zz=z+Kz_2/2
    vzz=vz+Kvz_2/2
    
    # 3. korak
    Kx_3=dt*vxx
    Kvx_3=dt*ubrzanje(xx,yy,zz,x_zemlja, y_zemlja, 0, Ms, Mz, gama)[0]
    Ky_3=dt*vyy
    Kvy_3=dt*ubrzanje(xx,yy,zz,x_zemlja, y_zemlja, 0, Ms, Mz, gama)[1]
    Kz_3=dt*vzz
    Kvz_3=dt*ubrzanje(xx,yy,zz,x_zemlja, y_zemlja, 0, Ms, Mz, gama)[2]
    
    xx=x+Kx_3
    vxx=vx+Kvx_3
    yy=y+Ky_3
    vyy=vy+Kvy_3
    zz=z+Kz_3
    vzz=vz+Kvz_3
    
    # 4. korak
    Kx_4=dt*vxx
    Kvx_4=dt*ubrzanje(xx,yy,zz,x_zemlja, y_zemlja, 0, Ms, Mz, gama)[0]
    Ky_4=dt*vyy
    Kvy_4=dt*ubrzanje(xx,yy,zz,x_zemlja, y_zemlja, 0, Ms, Mz, gama)[1]
    Kz_4=dt*vzz
    Kvz_4=dt*ubrzanje(xx,yy,zz,x_zemlja, y_zemlja, 0, Ms, Mz, gama)[2]

    # racunanje koeficijenata K
    Kx=1/6*(Kx_1+2*Kx_2+2*Kx_3+Kx_4)
    Kvx=1/6*(Kvx_1+2*Kvx_2+2*Kvx_3+Kvx_4)
    
    Ky=1/6*(Ky_1+2*Ky_2+2*Ky_3+Ky_4)
    Kvy=1/6*(Kvy_1+2*Kvy_2+2*Kvy_3+Kvy_4)
    
    Kz=1/6*(Kz_1+2*Kz_2+2*Kz_3+Kz_4)
    Kvz=1/6*(Kvz_1+2*Kvz_2+2*Kvz_3+Kvz_4)
    
    x=x+Kx
    vx=vx+Kvx
    
    y=y+Ky
    vy=vy+Kvy
    
    z=z+Kz
    vz=vz+Kvz
    
    X.append(x-x_zemlja)
    Y.append(y-y_zemlja)
    Z.append(z-z_zemlja)

    
    lz = t * sk  # longituda Zemlje 

    # pravougle koordinate Zemlje
    x_zemlja = au2m(1.) * np.cos(lz)
    y_zemlja = au2m(1.) * np.sin(lz)
    
    t+=dt
    
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(np.array(X)/r0z,np.array(Y)/r0z,np.array(Z)/r0z)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
    

    