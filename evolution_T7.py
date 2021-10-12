import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import kde
import random as rd

###############################
#######INITIALISATION
###############################

class User:
    
    def __init__(self) :
        self.pool = [0]
        self.taille = []
        self.facteurs = []
        self.seuil = 5
    
    def autoremplissage(self) :
        self.taille.append(300)
        self.facteurs.append(0.001)
        self.double = False
        self.g = 20
    
    def remplissage(self) :
        tmp = str(input("Mutation double (d) ou simple (s) ? "))
        if tmp == 'd' :
            self.double = True
        else : 
            self.double = False

        if double :
            self.taille.append(int(input("nombre de sites pour le premier facteur de mutation : ")))
            self.taille.append(int(input("nombre de sites pour le deuxième facteur de mutation : ")))
        else : 
            self.taille.append(int(input("nombre de sites pour le facteur de mutation : ")))

        if double :
            tmp = float(input("taux de mutation du premier facteur : "))
            self.facteurs.append(tmp)
            tmp = float(input("taux de mutation du deuxième facteur : "))
            self.facteurs.append(tmp)
        else :
            tmp = float(input("taux de mutation du facteur : "))
            self.facteurs.append(tmp)

        self.g = int(input("nombre de generations (<20) : "))

#valeurs = User()
#valeurs.autoremplissage()

###############################
#######EXECUTION_NB
###############################



class Simulation :
    """pour lancer une nouvelle simulation du nombre de mutation : instancier la classe puis choisir une simulation puis un ou plusieurs types de graphes simultanément"""

    def __init__(self):
        self.pool = valeurs.pool.copy()
        self.facteurs = valeurs.facteurs.copy()
        self.taille = valeurs.taille.copy()
        self.g = valeurs.g
        self.double = valeurs.double
        self.moy = []
        self.var = []
        self.seuil = 4


    def nb_mut_force(self) :
        gg = 0
        
        while gg < self.g :
            i = 0
            while i < len(self.pool) :
                tmp = self.pool[i]
                if self.double :
                    self.pool[i] = self.pool[i]+(self.taille[1]-self.pool[i])*self.facteurs[1]
                self.pool.append(tmp+(self.taille[0]-tmp)*self.facteurs[0])
                i += 2
            print(len(self.pool))
            self.moy.append(np.mean(self.pool))
            self.var.append(np.var(self.pool))
            gg += 1
    
    def nb_mut_force_seuil(self) :
        gg = 0
        
        while gg < self.g :
            i = 0
            while i < len(self.pool) :
                tmp = self.pool[i]
                if tmp < self.seuil :
                    if self.double :
                        self.pool[i] = self.pool[i]+(self.taille[1]-self.pool[i])*self.facteurs[1]
                    self.pool.append(tmp+(self.taille[0]-tmp)*self.facteurs[0])
                    i += 2
                else :
                    rado = rd.random()
                    if (rado < (1/2)) :
                        self.pool.pop(i)
            print(len(self.pool))
            self.moy.append(np.mean(self.pool))
            self.var.append(np.var(self.pool))
            gg += 1
    
    def graph_mean(self) :
         x = [i for i in range(self.g)]
         plt.plot(x, self.moy)
         plt.show()
    
    def graph_var(self) :
         x = [i for i in range(self.g)]
         plt.plot(x, self.var)
         plt.show()
            
    def graph_density(self) :#type de graph à revoir
        density = kde.gaussian_kde(self.pool)
        #density.covariance_factor = lambda : .75
        #density._compute_covariance()
        x = np.linspace(0,self.taille[0]*self.facteurs[0]*self.g, 70)
        y = density(x)
        
        plt.plot(x, y)
        plt.show()
    
#sim = Simulation()
#sim.nb_mut_force_seuil()
#sim.graph_density()
#sim.graph_var()



        
        
        
        
        

    
