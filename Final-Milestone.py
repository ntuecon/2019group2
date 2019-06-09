#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  9 16:59:21 2019

@author:user
"""

import numpy as np

from scipy.optimize import minimize



#The first part:producer 
class Producer(object):
    def __init__ (self,pr,cr):
        """ We will be concern with only 2 producers and 1 consumers when building this code"""
        #pr: number of producers (where each producer produces only 1 type of good)
        #cr: number of consumers (where each consumer provide only 1 type of factor)
        #dict is the collection of all the parameters
        #A: techonology
        #alpha, theta and omega is the production parameters
        
        self.pr = pr
        self.cr = cr
        para_dict = {'A': 0.5, 'alpha': 0.6, 'theta':0.5}
        self.dict = para_dict
        """input_list will be an array two values we are interested in.
        (quantity_supplied,labor_demanded)"""
        general_list = np.empty(self.pr + self.pr * self.cr)
        for i in range(0, self.pr + self.pr * self.cr):
            general_list[i] = 1
        self.gl = general_list
        
#Production function

    def prodftn(self):
        #Production function for firm1
        externality_array = np.empty(self.pr)
        externality_array[0] = 1
        externality_array[1] = self.dict['theta']
        goods_supplyarray = np.empty(self.pr)
        goods_supplyarray[0] = 1 
        for i in range(0,self.pr):
            goods_supplyarray[i] = (self.dict['A'] * (self.gl[i + self.pr]**self.dict['alpha']))* (1/(goods_supplyarray[0]* externality_array[i]))
        return goods_supplyarray        
                         
    def profit(self, p, w, t):
        #P for price 1/ price 2
        #W for wage
        #t for tax
        #calculate the earnings from sale of goods
        goods_earnings = np.empty(self.pr)
        for i in range(0,self.pr):
            goods_earnings[i] = self.prodftn[i] * p[i]
        #calculate the earnings from providing factor    
        factors_earnings = np.empty(self.pr)
        for j in range(0,self.pr):
            factors_earnings[j] = self.gl[j + self.pr] * w
        #calculate profit  
        #Sum the total profit, t stands for tax 
        profit = np.sum(goods_earnings) - np.sum(factors_earnings) - t * (self.prodftn[0])
        z = (-1)
        return z * profit 
      
    #Maximized under no externality
    def maximize(self,p,w,t):         
        guess = np.empty(self.pr)
        for i in range(0,self.pr):
            guess[i] = 1
        
        constraint = np.empty(self.pr)
        for i in range(0,self.pr):
            constraint[i] = self.prodftn[i] - self.gl[i]
            
        productionCon = {'type' : 'eq', 'fun' : constraint , 'args' : (p, w, t)}
        solution = minimize(self.profit, guess, args = (p,w,t), method = 'SLSQP', constraints = productionCon)
        
        return solution.x

        print "G1 quantity supplied: " + solution.x[0]
        
        print "G2 quantity supplied  " + solution.x[1]
        
        print "Firm1 labor demanded " + solution.x[2]
        
        print "Firm2 labor demanded " + solution.x[3]
        
            
#The second Part: Conusmer 
class Consumer(object):
    
    def _init_(self,pr,cr):
        #pr: number of producers (where each producer produces only 1 type of good)
        #cr: number of consumers (where each consumer provide only 1 type of factor)
        #dict is the collection of all the parameters
        #psi, beta and delta is the production parameters
        
        self.pr = pr
        self.cr = cr
        self.dict = {'psi': 1.2, 'beta': 0.3, 'delta': 0.7}
        consumption_goods_array = np.array(2)
        self.ca = consumption_goods_array 
        
        #An array storing labor and capital offered, [labor,capital]
        factors_array = np.array(1)
        self.fa = factors_array    
        
        general_list = np.empty(self.pr + self.pr + self.pr * self.cr)
        for i in range(0, self.pr + self.pr * self.cr):
            general_list[i] += 1
        self.gl = general_list
        
        
        
    def utility(self,general_list):
         # utility of the goods consumed in a array form
         goods_utility = 0
         for i in range(0,self.pr):
        # B is a random variable from [0,1], append each consumer's utility from both goods.
            goods_utility += self.gl[i] ** self.psi
    
        
        # factor dis-utility when providing k & l to producers for production
        
         labors_utility = 0
         for i in range(0 + self.pr, self.pr +self.pr):
             labors_utility[i] = self.beta* self.gl[i] ** self.delta
                
        # sum the utility - sum of factor disutility
         total_utility = goods_utility - labors_utility
        
         return (-1)* total_utility

          
    def constraint(self, general_list, pi,p,w,t):
        
        # pi: profits shares from two producers
        # total amount paid for goods
        goods_paid = 0
        for i in range(0,self.pr):
            goods_paid += self.gl[i] * p[i]
        
        # total amount earned from factors
        
        #An array storing each conumser's labor unit 
        factors_earned = 0
        for j in range(self.pr, self.pr + self.pr):
            factors_earned += self.gl[j] * w 
    
        #Budget constraint
        z = (-1)
        a = Producer
        pi = (a.profit) * z
        budget = pi + factors_earned - goods_paid  + t * (a.prodftn[0])
        
        return budget  
          
 
    def MaxUtility(self,p,w,t):
        #the maximize utility
        x0 = 0
        budgetCon = {'type' : 'eq', 'fun' : self.constraint, 'args' : (p,w,t)}
        solution = minimize(self.utility, x0, args = (p,w,t), method = 'SLSQP', constraints = budgetCon)
        
        return solution.x
    
        #Return the quantity demanded on good 1 
        print "G1 quantity demanded: " + solution.x[0] 
        
        print "G2 quantity demanded: " + solution.x[1]
        
        print "Labor_supplied:" + solution.x[2]
        
        
        
class Economy(object):
       
    def __init__(self,pr,cr):
        self.pr = pr
        self.cr = cr
    
    def objective(self):  
        a = Producer
        b = Consumer
        ex_post_array_pro = np.empty(self.pr + self.pr * self.cr)
        for i in range(0, self.pr + self.pr * self.cr):
            ex_post_array_pro[i] = a.maximize[i]
    
        ex_post_array_con = np.empty(self.pr + self.cr)
        for j in range(0, self.pr + self.cr):
            ex_post_array_con[i] = b.MaxUtility[i]
        
        excess_demand_array = np.empty(self.pr + self.cr)
        
        excess_demand_goods1 =  (ex_post_array_pro[0] - ex_post_array_con[0])
        excess_demand_goods2 =  (ex_post_array_pro[1] - ex_post_array_con[1])
        excess_demand_labor  =  (ex_post_array_pro[2] + ex_post_array_pro[3]) - ex_post_array_con[2]
    
        excess_demand_array[0] = excess_demand_goods1 ** 2
        excess_demand_array[1] = excess_demand_goods2 ** 2
        excess_demand_array[2] = excess_demand_labor ** 2
        
        return excess_demand_array
    
    def equilibrium(self):
        guess = np.empty(self.pr + self.cr)
        for i in range(0, self.pr + self.cr):
            guess[i] = 100
        
        sol = minimize(self.objective, guess, args = (), method ='Nelder-Mead' )
        
        return sol

#simulation of the production externality

print "Welcome to the world of production externality"
#Ask the user for the number the number of consumers and producers

#Ask the number of producer
pr = int(raw_input("Number of producers: "))
#Ask the number of consumer
cr = int(raw_input("Number of consumers: "))

test = Economy(pr,cr)
print test.objective()

    
    
