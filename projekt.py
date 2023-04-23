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
    def __init__(self,u0,lam0,rho0):
        self.lam = lam0
        self.rho = rho0
        self.u0 = u0

    def vypocti_u(self,start,stop,dt,dx,dy):
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
        [a,b] = self.u0.shape
        U=[]
        U_new = self.u0.copy()
        for time in range(start, stop+1):
            u = U_new.copy()
            for y in range(0, b):
                for x in range(0, a):
                    r_x_c = (2/self.rho[x,y])* dt
                    #sum = 0
                    if(x > 0):
                        U_new[x,y] += r_x_c * ((u[x-1,y] - u[x,y]) / ((1/self.lam[x,y] + 1/self.lam[x-1,y])* (dx**2)))
                        #sum += (u[x-1,y] - u[x, y])/(((1/(self.lam[x,y])) + (1/(self.lam[x - 1, y])))*(dx**2))
                    if(y > 0):
                        U_new[x,y] += r_x_c * ((u[x,y-1] - u[x,y]) / ((1/self.lam[x,y] + 1/self.lam[x,y-1])* (dy**2)))
                        #sum += (u[x, y - 1] - u[x, y])/(((1/(self.lam[x, y])) + (1/(self.lam[x, y - 1])))*(dy**2))
                    if(x < a-1):
                        U_new[x,y] += r_x_c * ((u[x+1,y] - u[x,y]) / ((1/self.lam[x,y] + 1/self.lam[x+1,y]) * (dx**2))) 
                        #sum += (u[x + 1, y] - u[x, y])/(((1/(self.lam[x, y])) + (1/(self.lam[x + 1, y])))*(dx**2))
                    if(y < b-1):
                        U_new[x,y] += r_x_c * ((u[x,y+1] - u[x,y]) / ((1/self.lam[x,y] + 1/self.lam[x,y+1])* (dy**2)))     
                        #sum += (u[x, y + 1] - u[x, y])/(((1/(self.lam[x, y])) + (1/(self.lam[x, y + 1])))*(dy**2))
                    #U_new[x, y] = u[x, y] + r_x_c*sum
            U.append(U_new) 
            print(np.unique(U_new)) 
        return U
    

    def vypocti_u_sparse(self,start,stop,dt,dx,dy):
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
        [a,b] = self.u0.shape
        lam = lil_matrix(self.lam)
        u0 = lil_matrix(self.u0)
        U=[]
        U_new = u0.copy()
        soucet = lil_matrix(np.empty(shape = (a, b), dtype = float))
        u = lil_matrix(np.empty(shape = (a, b), dtype = float))
        rho = lil_matrix(self.rho)
        for time in range(start, stop+1):
            for y in range(0, b):
                for x in range(0, a):
                    r_x_c = (2/rho[x,y])* dt
                    if(x > 0):
                        U_new[x,y] += r_x_c * ((u[x-1,y] - u[x,y]) / ((1/self.lam[x,y] + 1/self.lam[x-1,y])* (dx**2)))
                    if(y > 0):
                        U_new[x,y] += r_x_c * ((u[x,y-1] - u[x,y]) / ((1/self.lam[x,y] + 1/self.lam[x,y-1])* (dy**2)))
                    if(x < a-1):
                        U_new[x,y] += r_x_c * ((u[x+1,y] - u[x,y]) / ((1/self.lam[x,y] + 1/self.lam[x+1,y]) * (dx**2))) 
                    if(y < b-1):
                        U_new[x,y] += r_x_c * ((u[x,y+1] - u[x,y]) / ((1/self.lam[x,y] + 1/self.lam[x,y+1])* (dy**2)))     
            U.append(U_new.todense()) 
            print(np.unique(U_new)) 
        return U

    def vykresli_cas(self,start,stop,dt):
        u = self.vypocti_u(start,stop,dt)
        fig = plt.figure()  
        ax = sns.heatmap(u[stop])
        plt.show()
    

def run_u(S,typ,postup,start,stop,dt,dx,dy):
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
    if (postup == "sparse"):
        u = S.vypocti_u_sparse(start,stop,dt,dx,dy)
    else:
        u = S.vypocti_u(start,stop,dt,dx,dy)
    def init():
        plt.clf()
        plt.title(f"Teplota v t = {0:.3f}")
        plt.xlabel("x")
        plt.ylabel("y")
        ax = sns.heatmap(u[0])

    def animate(i):
        plt.clf()
        plt.title(f"Teplota v t = {int(i):.3f}")
        plt.xlabel("x")
        plt.ylabel("y")
        ax = sns.heatmap(u[i])

    anim = animation.FuncAnimation(fig, animate, init_func=init, interval=500,frames = stop-1, cache_frame_data=False)
    pillowwriter = animation.PillowWriter(fps=10)
    if (typ == "spirala"):
        if (postup == "sparse"):
            name = "spirala_sparse.gif"
        else:
            name = "spirala.gif"
    else:
        if (postup == "sparse"):
            name = "cihla_sparse.gif"
        else:
            name = "cihla.gif"
    anim.save(name, writer=pillowwriter)

    plt.show()



def porovnej(S,start,stop,dt,dx,dy):
    begin = t()
    S.vypocti_u(start,stop,dt,dx,dy)
    end = t()
    begin_S = t()
    S.vypocti_u_sparse(start,stop,dt,dx,dy)
    end_S = t()
    print("Vypocet U:", end - begin)
    print("Vypocet sparse U:", end_S - begin_S)
