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


#The first part:producer 
class Producer(object):
    def __init__ (self,pr,A,omega,alpha,theta):
        #2 Producers in economy 
        #Producer 1 produces good 1, producer 2 producers good 2
        #pr :producers
        #L:labor
        #A: techonology
        
        self.pr = pr
        self.A = A
        self.alpha = alpha
        self.theta = theta
        self.omega = omega
        production_array = np.array(self.pr) 
        self.pa = production_array
        
        #The length of input list = 2, the first two item are labor used by firm 1 & 2 respectively, the third and fourth are goods1 and goods 2 produced 
        input_list_f = np.empty(2 + self.pa)
        for i in range(0,2 + self.pa):
            input[i] = 1
        
            
#Production function
    def production_ne(self,input):
        #Production function for firm1 and firm2 (no reaction when no externality)
        for i in range(0,self.pa):
            self.pa[i] = self.A * (input[i]**self.alpha) 
            
            return np.sum(self.pa)
    
    
    def production_e(self,w,input):
        #Production function for firm1
        self.pa[0] = self.A * (input[0]**self.alpha) 
        #Define beta
        beta = self.omega * (1/abs(1-self.pa[0])) + (1 - self.omega)
        #Production function for firm2
        self.pa[1] = (self.A * (input[1]**self.theta)) * beta
        
        return self.pa - self.pa
    
    def profit(self,p,w,input):
        goods_earnings = np.empty(2)
        for i in range(0,2):
            goods_earnings[i] = self.pa[i]*p
            
        factors_earnings = np.empty(2)
        for j in range(0,2):
            factors_earnings[j] = w*input[j]
            
        profit_array = np.empty(2)
        for p in range(0,2):
            profit_array[p] = goods_earnings[p] - factors_earnings[p] 
    
            return (-1) * profit_array
      
           #Maximized under no externality
    def maximize_ne(self,p,w):
        
        guess = np.empty(4)
        for i in range(0,4):
            guess[i] = 1
        productionCon = {'type' : 'eq', 'fun' : self.production_ne, 'args' : (self.pa)}
        constraint = [productionCon]
        solution = minimize(self.profit, guess, args = (self.pa), method = 'SLSQP', constraints = constraint)
        
        print "Total Labor Demanded_ne:" + sum(solution.x[0] + solution.x[1])
        
        print "G1 quantity supplied_ne:" + solution.x[2]
        
        print "G2 quantity supplied_ne:" + solution.x[3]
    
    #Maximized with externality 
    def maximize_e(self,p,w,input):
        
        guess = np.empty(4)
        for i in range(0,4):
            guess[i] = 1
        productionCon = {'type' : 'eq', 'fun' : self.production_e, 'args' : (input)}
        constraint = [productionCon]
        solution = minimize(self.profit, guess, args = (input), method = 'SLSQP', constraints = constraint)
        
        print "Total Labor Demanded_e:" + sum(solution.x[0] + solution.x[1])
        
        print "G1 quantity supplied_e:" + solution.x[2]
        
        print "G2 quantity supplied_e:" + solution.x[3]
   
#The second Part: Conusmer 
class Consumer(object):
    
    def _init_(self,psi,beta,delta):
        # Only one consumer in the economy 
        # ngp: no of types of goods produced by each producer
        # lc: types of labour provided by consumer
        # kc: types of capital provided by consumer
        # p:random input a price list for each good 

        #The array storing the good1 and good2 produced by producers 
        self.psi = psi
        self.beta
        self.delta = delta
        consumption_goods_array = np.array(2)
        self.ca = consumption_goods_array 
        
        #An array storing labor and capital offered, [labor,capital]
        factors_array = np.array(1)
        self.fa = factors_array        
        
        input = np.empty(self.ca + self.fa)
        for i in range(0,self.ca + self.fa):
          input[i] = 1
        
        
    def utility(self,input):
         # utility of the goods consumed in a array form
         for i in range(0,input):
        # B is a random variable from [0,1], append each consumer's utility from both goods.
            goods_utility = sum(input[0:2]) ** self.psi
    
        
        # factor dis-utility when providing k & l to producers for production
         labors_utility = np.empty(1)
         for i in range(0,2):
             labors_utility[i] = self.beta* input[i] ** self.delta
                
        # sum the utility - sum of factor disutility
         total_utility = goods_utility - np.sum(labors_utility)
        
         return total_utility

          
    def constraint(self,input,pi,p,w):
        
        # pi: profits shares from two producers
        # total amount paid for goods
        goods_paid = []
        for i in range(0,2):
            goods_paid[i] = input[i]*p
        
        total_paid = np.sum(goods_paid)
        
        # total amount earned from factors
        
        #An array storing each conumser's labor unit 
        factors_earned = input[2] * w 
    
        #Budget constraint
        budget = pi + factors_earned - total_paid
        
        return budget  
    
       
       
    def MaxUtility(self,p,w,guess):
        #the maximize utility
        guess = np.empty(4)
        for i in range(0,4):
            guess[i] = 1
        budgetCon = {'type' : 'eq', 'fun' : self.constraint, 'args' : (input)}
        constraint = [budgetCon]
        solution = minimize(self.utility, guess, args = (input), method = 'SLSQP', constraints = constraint)
        
        print solution
        #Return the quantity demanded on good 1 
        print "G1 quantity demanded: " + solution.x[0] 
        
        print "G2 quantity demanded: " + solution.x[1]
        
        print "Labor_supplied:" + solution.x[2]
  
                     
 
#The third part:

class economy(self):
    
#Decide the condition under which the externality exists
#Price
#Find the optimal tax(or subsidy) under the contraint: total goods produced = goods consumed, total factors supplied = total factors demanded
   def __init__(self)
        pass
    
    def excess_demand(self):
        
    
    
   def social_welfare(self):
       #Total social welfare under no externality 
       #Concern: whether adding profit and utility directly suffers the doulbe-counting?
       total_social_welfare = Producer.profit
       
       return (-1)*total_social_welfare

#constraint two market equal(quantity comes from )
   def constraint_ne(self):

       
       
       
       
       
       
        
        
    
    
    
       
    
        
    
            
    

    
    
    
  