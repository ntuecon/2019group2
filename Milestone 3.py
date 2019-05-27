#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 09:11:19 2019

@author: JKong
"""

import numpy as np
from scipy.optimize import minimize



#The first part:producer 
class Producer(object):
    def __init__ (self,pr,cr,dict):
        """ We will be concern with only 2 producers and 1 consumers when building this code"""
        #pr: number of producers (where each producer produces only 1 type of good)
        #cr: number of consumers (where each consumer provide only 1 type of factor)
        #dict is the collection of all the parameters
        #A: techonology
        #alpha, theta and omega is the production parameters
        
        self.pr = pr
        self.cr = cr
        para_dict = {"A": 0.5, "alpha": 0.6, "theta": 0.7, "omega": 0.5 }
        self.dict = para_dict
        """input_list will be an array of all the values we are interested in.
        (SG1,SG2,DG1,DG2,SF1,DF1) where S: Supplied, D: Demanded, G: Goods, F: Factors"""
        input_list = []
        for i in range(0, 2 * (self.pr + self.cr)):
            input_list[i] += 1
        
        #To separate the input_list array
        Qty_supplied = input_list[0: self.pr * self.cr]
        self.QS = Qty_supplied
        Qty_demanded = input_list[self.pr * self.cr: self.pr * self.cr + self.pr * self.cr]
        self.QD = Qty_demanded
        factor_supplied = input_list[self.pr * self.cr + self.pr * self.cr: self.pr * self.cr + self.pr * self.cr + self.cr]
        self.fs = factor_supplied
        factor_demanded = input_list[self.pr * self.cr + self.pr * self.cr + self.cr + self.cr]
        self.fd = factor_demanded
        
#Production function
    #production without externality
    def production_ne(self,QS):
        #Production function for firm1 and firm2 (no reaction when no externality)
        for i in range(0,self.pr):
            for x in range(0,self.cr):
                self.QS[i] = self.dict["A"] * (self.fs[x]**self.dict["alpha"]) 
            return self.QS
    
    #production with externality
    def production_e(self,w,QS):
        #Production function for firm1
        self.QS[0] = self.dict["A"] * (self.fs[0]**self.dict["alpha"]) 
        #Define beta, which will be the effect of externality on the 2nd producer
        beta = self.dict["omega"] * (1/abs(1-self.QS[0])) + (1 - self.dict["omega"])
        #Production function for firm2
        self.QS[1] = (self.dict["A"] * (self.fs[0]**self.dict["theta"])) * beta
        return self.QS
    
    #Different prices for the goods, to introduce the price variables
    def price(self):
        price = np.array(self.pr)
        for i in range(0,self.pr):
            price[i] = 1
        return price
        
    def profit(self, price, w, QS):
        #calculate the earnings from sale of goods
        goods_earnings = np.empty(self.pr)
        for i in range(0,self.pr):
            goods_earnings[i] = self.QS[i] * price[i]
        #calculate the earnings from providing factor    
        factors_earnings = np.empty(self.pr)
        for j in range(0,self.pr):
            factors_earnings[j] = self.QS[j] * w
            
        profit_array = np.empty(self.pr)
        for p in range(0,self.pr):
            profit_array[p] = goods_earnings[p] - factors_earnings[p] 
            return (-1) * profit_array
      
    #Maximized under no externality
    def maximize_ne(self,price,w):
        
        guess = np.empty(len(self.QS) + len(price) + 1)
        for i in range(0,len(self.QS) + len(price) + 1):
            guess[i] = 1
        productionCon = {'type' : 'eq', 'fun' : self.production_ne, 'args' : (self.QS)}
        constraint = [productionCon]
        solution = minimize(self.profit, guess, args = (self.pa), method = 'SLSQP', constraints = constraint)
        
        print "G1 quantity supplied_ne: " + solution.x[0]
        
        print "G2 quantity supplied_ne: " + solution.x[1]
        
        print "Price of G1: " + solution.x[2]
        
        print "Price of G2: " + solution.x[3]
        
        print "Wage of labour: " + solution.x[4]
        
        return solution.x
    
    #Maximized with externality 
    def maximize_e(self,p,w,input):
        
        guess = np.empty(len(self.QS) + len(self.price) + 1)
        for i in range(0,len(self.QS) + len(self.price) + 1):
            guess[i] = 1
        productionCon = {'type' : 'eq', 'fun' : self.production_e, 'args' : (self.QS, self.price, w)}
        constraint = [productionCon]
        solution = minimize(self.profit, guess, args = (input), method = 'SLSQP', constraints = constraint)
        
        print "G1 quantity supplied_e: " + solution.x[0]
        
        print "G2 quantity supplied_e: " + solution.x[1]
        
        print "Price of G1: " + solution.x[2]
        
        print "Price of G2: " + solution.x[3]
        
        print "Wage of labour: " + solution.x[4]
        
        return solution.x
   
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
        
        return solution.x
  
                     
 

class economy(self):
    
#Decide the condition under which the externality exists
#Price
#Find the optimal tax(or subsidy) under the contraint: total goods produced = goods consumed, total factors supplied = total factors demanded
   def __init__(self):
        pass
    
   def excess_demand_ne_g1(self,w,p,Producer,Consumer):
        w = Symbol('w')
        Producer.price = Symbol('p')
        solve(Consumer.MaxUtility[0] - Producer.maximize_ne[0],w,Producer.price)

   def excess_demand_ne_g2(self,w,p,Producer,Consumer):
       w = Symbol('w')
       Producer.price = Symbol('p')
       solve(Consumer.MaxUtility[1] - Producer.maximize_ne[1],w,Producer.price)     

        
   def excess_demand_onlabor_ne(self,w,r,Producer,Consumer):
        w = Symbol('w')
        Producer.price = Symbol('p')
        solve(Consumer.MaxUtility[2] - Producer.maximize_ne[2],w,Producer.price)
    
   def excess_demand_e_g1(self,w,r,Producer, Consumer):
        w = Symbol('w')
        Producer.price = Symbol('p')
        solve(Consumer.MaxUtility[0] - Producer.maximize_e[0],w,Producer.price)
        
   def excess_demand_e_g2(self,w,r,Producer, Consumer):
        w = Symbol('w')
        Producer.price = Symbol('p')
        solve(Consumer.MaxUtility[1] - Producer.maximize_e[1],w,Producer.price)


   def excess_demand_onlabor_e(self,w,r,Producer,Consumer):
        w = Symbol('w')
        Producer.price = Symbol('p')
        solve(Consumer.MaxUtility[2] - Producer.maximize_ne[2],w,Producer.price)
        
class government(self):
    
    def __init__(self):
        pass
    
    def optimal_tax_rate(self):