import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import seaborn as sns
from scipy.sparse import csr_matrix, csc_matrix, coo_matrix, lil_matrix
import scipy
from time import time as t


import warnings

# potlaceni warningu
warnings.filterwarnings('ignore')


class Obrazec:
    def __init__(self,u_0,lam_0,rho_0,c):
        self.lam = lam_0
        self.rho = rho_0
        self.u_0 = u_0
        self.c = c

    def vypocti_u(self,start,stop,dt):
        """
        Vypocet zmen teplot pomoci vzorce U 
        Parametry:
            start - pocatecni cas
            stop - konecny cas
            dt - casovy krok
            
        Return:
            U - pole matic

        Pouziti:
            vypocti_u(start,stop,dt)
        """
        [a,b] = self.u_0.shape
        [dx,dy] = 1,1
        U=[]
        u = self.u_0.copy()
        U_new = u
        r_x_c = (2/(self.rho*self.c))* dt
        for time in range(start, stop):
            for x in range(0, a):
                for y in range(0, b):
                    if(x > 0):
                        U_new += r_x_c * ((u[x-1,y] - u[x,y]) / ((1/self.lam[x,y] + 1/self.lam[x-1,y])* dx**2))
                    if(y > 0):
                        U_new += r_x_c * ((u[x,y-1] - u[x,y]) / ((1/self.lam[x,y] + 1/self.lam[x,y-1])* dy**2))
                    if(x < a-1):
                        U_new += r_x_c * ((u[x+1,y] - u[x,y]) / ((1/self.lam[x,y] + 1/self.lam[x+1,y]) * dx**2)) 
                    if(y < b-1):
                        U_new += r_x_c * ((u[x,y+1] - u[x,y]) / ((1/self.lam[x,y] + 1/self.lam[x,y+1])* dy**2))    
            u = U_new
            U.append(u)         
            time += 1
            print(np.unique(u))
        return U
    
    def vykresli_cas(self,start,stop,dt):
        u = self.vypocti_u(start,stop,dt)
        fig = plt.figure()  
        ax = sns.heatmap(u[stop-1])
        plt.show()
    

def run_u(S,typ,start,stop,dt):
    """
    Vypocet matic U a ulozeni animace
    Parametry:
        S - trida Object
        typ - cihla/spirala dle zadani
        start - pocatecni cas
        stop - konecny cas
        dt - casovy krok
                
    Return:
        None

    Pouziti:
        run_u(S,typ,pstart,stop,dt)
    """
    fig = plt.figure()
    u = S.vypocti_u(start,stop,dt)
    def init():
        plt.clf()
        ax = sns.heatmap(u[0])

    def animate(i):
        plt.clf()
        ax = sns.heatmap(u[i])

    anim = animation.FuncAnimation(fig, animate, init_func=init, interval=500,frames = stop-1, cache_frame_data=False)
    #HTML(anim.to_html5_video())
    pillowwriter = animation.PillowWriter(fps=10)
    if (typ == "spirala"):
        name = "spirala.gif"
    elif (typ == "cihla"):
        name = "cihla.gif"
    anim.save(name, writer=pillowwriter)

    plt.show()


def porovnej(S,start,stop,dt):
    begin = t()
    S.vypocti_u(start,stop,dt)
    end = t()
    print("U:", end - begin)
    #
    #begin = t()
    #S.vypocti_sparse(start,stop,dt)
    #end = t()
    #print("Sparse:", end - begin)


"""
typ = "spirala"

3if (typ == "spirala"):
    lam_0 = np.load("spiralaRHO_c.npy")
    rho_0 = np.load("spiralaLAM.npy")
    u_0 = np.load("spiralaU_initial.npy")
elif (typ == "cihla"):
    rho_0 = np.load("cihlyRHO_c.npy")
    lam_0 = np.load("cihlyLAM.npy")
    u_0 = np.load("cihlyU_initial.npy")

#c = 1   # merna tepelna kapacita
#S = Obrazec(u_0,lam_0,rho_0,c)
#start = 0
#stop = 2
#dt =   0.001

#run_u(S,typ,start,stop,dt)

#print(lam_0.shape)
"""
