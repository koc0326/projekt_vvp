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

    def vypocti_u(self,stop,dt,dx,dy):
        """
        Vypocet zmen teplot pomoci vzorce U 
        Parametry:
            stop - konecny cas
            dt - casovy krok
            dx, dy - rozmery podoblasti
            
        Return:
            U - pole matic

        Pouziti:
            vypocti_u(stop,dt,dx,dy)
        """
        [a,b] = self.u0.shape
        U=[]
        U_new = self.u0.copy()
        for time in range(0, stop+1):
            u = U_new.copy()
            U.append(u) 
            for y in range(0, b):
                for x in range(0, a):
                    r_x_c = (2/self.rho[x,y])* dt
                    sum = 0
                    if(x > 0):
                        sum += (u[x-1,y] - u[x, y])/(((1/(self.lam[x,y])) + (1/(self.lam[x - 1, y])))*(dx**2))
                    if(y > 0):
                        sum += (u[x, y - 1] - u[x, y])/(((1/(self.lam[x, y])) + (1/(self.lam[x, y - 1])))*(dy**2))
                    if(x < a-1):
                        sum += (u[x + 1, y] - u[x, y])/(((1/(self.lam[x, y])) + (1/(self.lam[x + 1, y])))*(dx**2))
                    if(y < b-1):
                        sum += (u[x, y + 1] - u[x, y])/(((1/(self.lam[x, y])) + (1/(self.lam[x, y + 1])))*(dy**2))
                    U_new[x, y] = u[x, y] + r_x_c*sum
        return U
    

    def vypocti_u_sparse(self,stop,dt,dx,dy):
        """
        Vypocet zmen teplot pomoci vzorce U 
        Parametry:
            stop - konecny cas
            dt - casovy krok
            dx, dy - rozmery podoblasti
            
        Return:
            U - pole matic

        Pouziti:
            vypocti_u_sparse(stop,dt,dx,dy)
        """
        [a,b] = self.u0.shape
        U=[]
        U_new = csc_matrix(self.u0.copy())
        soucet = csc_matrix(np.zeros((a,b)))
        r_x_c = csc_matrix(np.empty((a,b)))
        for i in range(0,a):
            for j in range(0,b):
                r_x_c[i,j] =  (2/self.rho[i,j]) * dt 
        #r_x_c = csc_matrix((2/self.rho)*dt)
        for time in range(0, stop+1):
            mat_u_prev = U_new.copy() 
            U.append(U_new) 
            soucet[1:a, :] += (mat_u_prev[:a-1, :] - mat_u_prev[1:a, :])/(((1/(self.lam[1:a, :])) + (1/(self.lam[:a-1, :])))*(dx**2))
            soucet[:, 1:b] += (mat_u_prev[:, :b-1] - mat_u_prev[:, 1:b])/(((1/(self.lam[:, 1:b])) + (1/(self.lam[:, :b-1])))*(dy**2))
            soucet[:a-1, :] += (mat_u_prev[1:a, :] - mat_u_prev[:a-1, :])/(((1/(self.lam[:a-1, :])) + (1/(self.lam[1:a, :])))*(dx**2))
            soucet[:, :b-1] += (mat_u_prev[:, 1:b] - mat_u_prev[:, :b-1])/(((1/(self.lam[:, :b-1])) + (1/(self.lam[:, 1:b])))*(dy**2))
            U_new = mat_u_prev + r_x_c.multiply(soucet)
            U_new = csc_matrix(U_new)
        return U

    def vykresli_cas(self,postup,stop,dt,dx,dy):
        """
        Vykresli teplotu v danem case stop pomoci vypocti_u
        Parametry:
            stop - konecny cas
            dt - casovy krok
            dx, dy - rozmery podoblasti
            
        Return:
            None

        Pouziti:
            vykresli_cas(stop,dt,dx,dy)
        """
        U = self.vypocti_u(stop,dt,dx,dy)
        color =sns.color_palette("coolwarm", as_cmap=True)
        fig = plt.figure()  
        plt.title(f"Teplota v čase t = {int(stop)}")
        ax = sns.heatmap(U[stop],cmap=color)
        plt.show()
    

def run_u(S,typ,stop,dt,dx,dy):
    """
    Vypocet matic U pomoci vypocti_u a ulozeni animace
    Parametry:
        S - trida Object
        typ - cihla/spirala dle zadani
        stop - konecny cas
        dt - casovy krok
        dx,dy - rozmery podoblasti   
    Return:
        None

    Pouziti:
        run_u(S,typ,stop,dt)
    """
    color =sns.color_palette("coolwarm", as_cmap=True)

    fig = plt.figure()
    U = S.vypocti_u(stop,dt,dx,dy)
    
    def init():
        plt.clf()
        plt.title(f"Teplota v t = 0")
        ax = sns.heatmap(U[0], cmap=color)

    def animate(i):
        plt.clf()
        plt.title(f"Teplota v t = {int(i)}")
        ax = sns.heatmap(U[i], cmap = color)

    anim = animation.FuncAnimation(fig, animate, init_func=init, interval=100, frames = stop+1, cache_frame_data=False, repeat = False)
    pillowwriter = animation.PillowWriter(fps=10)
    if (typ == "spirala"):
            name = "spirala.gif"
    else:
            name = "cihla.gif"
    anim.save(name, writer=pillowwriter)
    plt.show()

def run_u_sparse(S,typ,stop,dt,dx,dy):
    """
    Vypocet matic U pomoci vypocti_u_sparse a ulozeni animace
    Parametry:
        S - trida Object
        typ - cihla/spirala dle zadani
        stop - konecny cas
        dt - casovy krok
        dx,dy - rozmery podoblasti
                
    Return:
        None

    Pouziti:
        run_u_sparse(S,typ,stop,dt)
    """
    color =sns.color_palette("coolwarm", as_cmap=True)

    fig = plt.figure()
    U = S.vypocti_u_sparse(stop,dt,dx,dy)
    
    def init():
        plt.clf()
        plt.title(f"Teplota v t = 0")
        ax = sns.heatmap(U[0].todense(),cmap=color)
        #plt.imshow(U[0].todense())

    def animate(i):
        plt.clf()
        plt.title(f"Teplota v t = {int(i)}")
        ax = sns.heatmap(U[i].todense(), cmap = color)
        #plt.imshow(U[i].todense())

    anim = animation.FuncAnimation(fig, animate, init_func=init, interval=100, frames = stop+1, cache_frame_data=False, repeat = False)
    pillowwriter = animation.PillowWriter(fps=10)
    if (typ == "spirala"):
            name = "spirala_sparse.gif"
    else:
            name = "cihla_sparse.gif"
    anim.save(name, writer=pillowwriter)
    plt.show()

def porovnej(S,stop,dt,dx,dy):
    """
    Porovnani metod vypocti_u  a vypocti_u_sparse
    Parametry:
        S - trida Object
        stop - konecny cas
        dt - casovy krok
        dx,dy - rozmery podoblasti
                
    Return:
        None

    Pouziti:
        porovnej(S,typ,stop,dt,dx,dy)
    """
    begin = t()
    S.vypocti_u(stop,dt,dx,dy)
    end = t()
    begin_S = t()
    S.vypocti_u_sparse(stop,dt,dx,dy)
    end_S = t()
    print("Vypocet U:", end - begin)
    print("Vypocet sparse U:", end_S - begin_S)


