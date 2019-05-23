#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 20 15:26:15 2019

@author: JKong
"""

import numpy as np
import random
from scipy.optimize import minimize
from sympy import diff



class producer(object):
    def __init__ (self,pr,A,K,L,omega,alpha,theta):
        #2 Producers in economy 
        #Producer 1 produces good 1, producer 2 producers good 2
        #pr :producers
        #K: capital
        #L:labor
        #A: techonology
        
        self.pr = pr
        self.A = A
        self.K = K
        self.L = L
        self.alpha = alpha
        self.theta = theta
        self.omega = omega
        production_array = np.array(self.pr)
        self.pa = production_array
            
#Production function
    def production_ne(self,w,r):
        #Production function for firm1 and firm2 (no reaction when no externality)
        for i in range(0,self.pa):
            self.pa[i] = self.A * (self.K**self.alpha) * (self.L**(1-self.alpha))
            
            return self.pa - self.pa  
    
    
    def production_e(self,w,r):
        #Production function for firm1
        self.pa[0] = self.A * (self.K**self.alpha) * (self.L**(1-self.alpha))
        #Define beta
        beta = self.omega * (1/abs(1-self.pa[0])) + (1 - self.omega)
        #Production function for firm2
        self.pa[1] = (self.A * (self.K**self.theta) * (self.L**(1-self.theta))) * beta
        
        return self.pa - self.pa
    
    def profit(self,p,w,r):
        goods_earnings = np.empty(2)
        for i in range(0,2):
            goods_earnings[i] = self.pa[i]*p
            
        factors_earnings = np.empty(2)
        for j in range(0,2):
            factors_earnings[j] = r*self.K + w* self.L
            
        profit_array = np.empty(2)
        for p in range(0,2):
            profit_array[p] = goods_earnings[p] - factors_earnings[p] 
            return (-1)* profit_array
           
    #Maximized under no externality
    def maximize_ne(self,p,w,r):
        
        guess = np.empty(4)
        for i in range(0,4):
            guess[i] = 1
        productionCon = {'type' : 'eq', 'fun' : self.production_ne, 'args' : (r,w)}
        constraint = [productionCon]
        solution = minimize(self.profit, guess, args = (p,r,w), method = 'SLSQP', constraints = constraint)
        
        print solution
    
    #Maximized with externality 
    def maximize_e(self,p,w,r):
        
        guess = np.empty(4)
        for i in range(0,4):
            guess[i] = 1
        productionCon = {'type' : 'eq', 'fun' : self.production_e, 'args' : (r,w)}
        constraint = [productionCon]
        solution = minimize(self.profit, guess, args = (p,r,w), method = 'SLSQP', constraints = constraint)
        
        print solution

#Second Part: Conusmer 
class consumer(object):
    def __init__ (self,pr,A,K,L,omega,alpha,theta):
        #2 Producers in economy 
        #Producer 1 produces good 1, producer 2 producers good 2
        #pr :producers
        #K: capital
        #L:labor
        #A: techonology
        
        self.pr = pr
        self.A = A
        self.K = K
        self.L = L
        self.alpha = alpha
        self.theta = theta
        self.omega = omega
        production_array = np.array(self.pr)
        self.pa = production_array
        
    def utility(self):
        
    
    
        
        
         
       
    
            
        
           