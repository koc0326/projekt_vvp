import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.animation import FuncAnimation
from scipy.sparse import csr_matrix, csc_matrix, coo_matrix, lil_matrix

zadani = "cihly"

if (zadani == "spirala"):
    lam0 = np.load("spiralaRHO_c.npy")
    rho0 = np.load("spiralaLAM.npy")
    u0 = np.load("spiralaU_initial.npy")
elif (zadani == "cihly"):
    rho0 = np.load("cihlyRHO_c.npy")
    lam0 = np.load("cihlyLAM.npy")
    u0 = np.load("cihlyU_initial.npy")



#lam0 = 3*np.ones((3,5))
#rho0 = 2*np.ones((3,5))
#u0 = 10*np.ones((3,5))

[a,b] = lam0.shape
print("a=",a,"b=",b)


Nx = b          # pocet dx
Ny = a          # pocet dy
print("Nx=",Nx,"Ny=",Ny)


dx = 1 #int(a/(Nx-1))       # vyska dilku = 1
dy = 1 #int(b/(Ny-1))       # sirka dilku dy = 1
print("dx=",dx,"dy=",dy)


start = 0
stop = 100
dt = 5
time = np.linspace(start,stop,int(stop/dt))
time = (np.ceil(time)).astype(int)


print(time)

U0 = np.pad(u0, pad_width=1, mode='constant', constant_values=0)
lam = np.pad(lam0, pad_width=1, mode='constant', constant_values=1)
lam_matice = np.pad(lam0, pad_width=1, mode='constant', constant_values=1)
rho = np.pad(rho0, pad_width=1, mode='constant', constant_values=1)
rho_matice = np.pad(rho0, pad_width=1, mode='constant', constant_values=1)



#print(U0)
#print(lam)
#print(rho)

[rows,cols] = lam.shape
U = np.zeros((rows,cols))

c = np.full((rows,cols),0.88103)         # merna tepelna kapacita betonu, predpokladame betonovou cihlu i spiralu
c_matice = np.full((a,b),0.88103)

time = 0
stop = 100
dt = 1

def Ucykly(U0,lam,rho,rows,cols,c,dt):
        for x in range(1,rows-1):
                for y in range(1,cols-1):
                        U[x,y] = U0[x,y] + (2/(rho[x,y]*c[x,y])) * dt * (((U0[x-dx,y] - U0[x,y]) / (1/lam[x,y] + 1/lam[x-dx,y])* dx**2)) + (((U0[x,y+dy] - U0[x,y]) / (1/lam[x,y] + 1/lam[x,y+dy])* dy**2)) + (((U0[x+dx,y] - U0[x,y]) / (1/lam[x,y] + 1/lam[x+dx,y])* dx**2)) + (((U0[x,y+dy] - U0[x,y]) / (1/lam[x,y] + 1/lam[x,y+dy])* dy**2))
        Unew = U[1:rows-1,1:cols-1]
        return [U,Unew]

def VypoctiUCylky(U0,time,stop,dt):
        while time < stop:
                [U,Unew] = Ucykly(U0,lam,rho,rows,cols,c,time)
                U0 = U
                time += dt
                #nazev = str(time) + ".png"
                ax = sns.heatmap(Unew,cmap = "viridis", xticklabels = False, yticklabels=False, linewidth=0.5)
                if time < stop:
                        #plt.savefig(str(time)+".png", dpi=400)
                        plt.show(block=False)
                        #plt.pause(0.2)  
                        plt.close()
                elif time == stop:
                        plt.savefig("cykly.png", dpi=400)
                        plt.show(block=True)
                        #plt.pause(0.5)  
                        #plt.close()
        
#VypoctiUCylky(U0,100,2000,100)

def Umatice(U0,lam_matice,rho,c,dt):
        # x-dx
        Uleft = ((U0[1:a+1,0:b]) - u0)
        lamleft = (1/lam0 + 1/(lam_matice[1:a+1,0:b])) * dx**2
        # y+dy
        Uup = ((U0[0:a,1:b+1]) - u0)
        lamup = (1/lam0 + 1/(lam_matice[0:a,1:b+1])) * dy**2
        # x+dx
        Uright = ((U0[1:a+1,-b:]) - u0)
        lamright = (1/lam0 + 1/(lam_matice[1:a+1,-b:])) * dx**2
        # y-dy
        Udown = ((U0[-a:,1:b+1]) - u0)
        lamdown = (1/lam0 + 1/(lam_matice[-a:,1:b+1])) * dy**2

        U = (2/(np.dot(rho,c))) * dt * ((Uleft / lamleft) + (Uup / lamup) + (Uright / lamright) + (Udown / lamdown))
        ax = sns.heatmap(U,cmap = "viridis", xticklabels = False, yticklabels=False, linewidth=0.5)
        plt.show()
        return print(U)

#Umatice(U0,lam_matice,rho0,c_matice,time)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import seaborn as sns

u0 = np.full((2,5),5)
u = np.tile(u0, (4,1,1))
u[1] = u[0]+5
u[3] = u[1]-3
print(u)

fig = plt.figure()
def init():
    sns.heatmap(u0, vmax=.8, square=True, cbar=False)

def animate(i):
    sns.heatmap(u[i], vmax=.8, square=True, cbar=False)


anim = animation.FuncAnimation(fig, animate, init_func=init, frames=10, repeat = False)

savefile = r"test3.gif"
pillowwriter = animation.PillowWriter(fps=20)
anim.save(savefile, writer=pillowwriter)

plt.show()