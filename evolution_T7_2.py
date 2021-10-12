#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 16:39:05 2021

@author: adrien
"""

import numpy as np
import matplotlib.pyplot as plt
import random as rd
import scipy.special as sc
import csv

class user :
    
    def __init__(self):
        self.seq = ""
        self.facteurs = []
        self.g = 3
        self.mu0 = 100
    
    def ordi(self, alpha, g, mu0) :
        self.facteurs.append(alpha)
        self.g = g
        self.mu0 = mu0
    
    def humain(self) :
        self.seq = str(input("séquence à muter : "))
        self.facteurs.append(float("taux de mutation du domaine utilisé : "))
        self.g = int(input("Nombre de generations avant arret de la mutation : "))
        self.mu0 = int(input("nombre de cellules à t0 : "))
    
    def default(self) :
        self.seq = "ATGTCGAAAGCTACATATAAGGAACGTGCTGCTACTCATCCTAGTCCTGTTGCTGCCAAGCTATTTAATATCATGCACGAAAAGCAAACAAACTTGTGTGCTTCATTGGATGTTCGTACCACCAAGGAATTACTGGAGTTAGTTGAAGCATTAGGTCCCAAAATTTGTTTACTAAAAACACATGTGGATATCTTGACTGATTTTTCCATGGAGGGCACAGTTAAGCCGCTAAAGGCATTATCCGCCAAGTACAATTTTTTACTCTTCGAAGACAGAAAATTTGCTGACATTGGTAATACAGTCAAATTGCAGTACTCTGCGGGTGTATACAGAATAGCAGAATGGGCAGACATTACGAATGCACACGGTGTGGTGGGCCCAGGTATTGTTAGCGGTTTGAAGCAGGCGGCGGAAGAAGTAACAAAGGAACCTAGAGGCCTTTTGATGTTAGCAGAATTGTCATGCAAGGGCTCCCTAGCTACTGGAGAATATACTAAGGGTACTGTTGACATTGCGAAGAGCGACAAAGATTTTGTTATCGGCTTTATTGCTCAAAGAGACATGGGTGGAAGAGATGAAGGTTACGATTGGTTGATTATGACACCCGGTGTGGGTTTAGATGACAAGGGAGACGCATTGGGTCAACAGTATAGAACCGTGGATGATGTGGTCTCTACAGGATCTGACATTATTATTGTTGGAAGAGGACTATTTGCAAAGGGAAGGGATGCTAAGGTAGAGGGTGAACGTTACAGAAAAGCAGGCTGGGAAGCATATTTGAGAAGATGCGGCCAGCAAAACTAA"
        self.facteurs.append(0.00354)
  
class simulation :

    def __init__(self) :
        us = user()
        us.default()
        self.g = us.g
        self.facteurs = us.facteurs
        self.seq = us.seq
        self.freq = {'A' : 0, 'T' : 0, 'C' : 0, 'G' : 0}
        self.mu0 = us.mu0
    
    def set(self, g, facteur, mu0) :
        self.g = g
        self.facteurs = facteur
        self.mu0 = mu0
    
    def traitement(self) :
        for i in self.seq :
            self.freq[i] += 1
    
    def simulation(self) :
        g = self.g
        pool = []
        #nb_cell = self.mu0*(2**g)
        x = 0
        while x <= g :
            p = sc.binom(g,x)*self.mu0
            j = 1
            while j <= p :
                
                z = np.random.geometric(p=self.facteurs[0], size=self.freq['C'])
                nb_mut = (z <= x).sum()
                pool.append([])
                tmp = [f for f in range(self.freq['C'])]

                k = 0
                while k < nb_mut :
                    indice = rd.randint(0, len(tmp)-1)
                    pool[-1].append(tmp[indice])
                    tmp.pop(indice)
                    k += 1
                j += 1
            x += 1
        self.pool = pool
        print("simulation faite, taille finale de la population : ", len(pool))
     
    def diversityStricte(self) :
        """à lancer après une simulation, diversité comptant [1,2] et [1,3] comme different"""
        tmp_pool = self.pool.copy()
        i = 0
        while i < len(tmp_pool) : 
            if tmp_pool.count(tmp_pool[i]) > 1 :
                tmp_pool.pop(i)
            else :
                i += 1
        print("Diversité stricte obtenue : ", (len(tmp_pool)/len(self.pool)))
        return((len(tmp_pool)/len(self.pool)))
     
    def diversity2(self) :
        """à lancer après une simulation, diversité comptant [1,2] et [1,3] comme 1/2 different"""
        tmp_pool = self.pool
        

#sim = simulation()
#sim.traitement()
#sim.simulation()


def statistiques(sim) :
    print("A partir de ", sim.mu0, " cellules, nous sommes arrivées en ", sim.g, "générations à une population de ", len(sim.pool) , " individus")
    repartition = [len(k) for k in sim.pool]
    print("Nombre max de mutation sur un même individu : ", max(repartition))
    print("Nombre min de mutation sur un même individu : ", min(repartition))
    print("Nombre d'individus sans mutations : ", (repartition.count(0)), " (", (repartition.count(0)/len(sim.pool))*100 , "%)")
    print("Moyenne : ", np.mean(repartition), " Variance : ", np.var(repartition))
    print("---------------------------------------------------------------------")
    print(sim.diversityStricte())
    print(sim.diversity2())

def graph_div1(gg) :
    g = 1
    x= [i for i in range(gg)]
    y = []
    while g <= gg :
        sim = simulation()
        sim.set(g, [0.00354], 4)
        sim.traitement()
        sim.simulation()
        y.append(sim.diversityStricte())
        g += 1
    plt.plot(x,y)
    plt.show()

def graph_div_k(gg, kk) :
    x= [i for i in range(gg)]
    y = []
    k = 0
    while k < kk :
        g = 1
        while g <= gg :
            sim = simulation()
            sim.set(g, [0.00354], 4)
            sim.traitement()
            sim.simulation()
            y.append(sim.diversityStricte())
            g += 1
        k += 1
    y1 = []
    print(y, len(y))
    for i in range(gg) :
        yy = 0
        for j in range(kk) :
            yy += y[j+(i*kk)]
        y1.append(yy/kk)
    print(y1)
    plt.plot(x,y1)
    plt.show()

def w_csv(yy, gg, alpha) :
    g = 1
    y = []
    while g <= gg :
        sim = simulation()
        sim.set(g, [alpha], 4)
        sim.traitement()
        sim.simulation()
        y.append(sim.diversityStricte())
        g += 1
    yy.append(y)

    
#statistiques(sim)
yy = []
#graph_div1(10)
liste_alpha = [0.00714, 0.06304, 0.01525, 0.005229]
for a in liste_alpha :  
    for i in range(10) :
        w_csv(yy, 11, a)

kwargs = {'newline': ''}
with open('all_test.csv', 'w', **kwargs) as fp:
    writer = csv.writer(fp, delimiter=';')
    # writer.writerow(["your", "header", "foo"])  # write header
    writer.writerows(yy)

#graph_div_k(7, 1)




