nc# -*- coding: utf-8 -*-
"""
Created on Wed May  1 17:46:25 2019

@author:
"""

#Introduce the package 
import numpy as np
import random
from scipy.optimize import minimize
from sympy import diff

class Economy(object):

    def _init_(self,nc,npro,ngp,lc,kc,T,A,lp,kp,externality):
        # nc: no of consumers in an economy
        # npro: no of producers in an economy
        # ngp: no of types of goods produced by each producer
        # lc: unit of labor each consumer has
        # kc: unit of capital each consumer has 
        # T: technology for a producer
        # A: alpha, the production parameters for a producer
        # lp: unit of labor used by producer
        # kp: unit of capital used by producer
        # externality: whether the externality exists or not 
        # p:random input a price list for each good 
        
        self.nc = nc
        self.npro = npro
        self.ngp = ngp
        self.lc = lc
        self.kc = kc
        self.T = T
        self.A = A
        self.lp = lp
        self.kp = kp
        self.externality = externality 
        goods_length = (self.npro * self.ngp)
        self.gl = goods_length
        
#The first part: Consumer
        
    def utility(self,input):
        
         goods = np.array(input[0:self.gl])
         factors = np.array(input[self.gl:(self.gl + self.lc * self.nc + self.kc* self.nc)])
        
        # utility of the goods consumed in a array form
        
         goods_utility = np.empty(self.nc * self.gl)
        
         for i in goods:
        # B will be a random variable from [0,1], append each consumer's utility from a certain good.
            goods_utility[i] = goods[i] ** random.randint(0,1)
    
        
        # factor dis-utility when providing k & l to producers for production
         factors_utility = np.empty(self.nc * self.fc * self.nc)
         for i in factors:
             factors_utility[i] = random.randint(0,1) * factors[i] ** random.randint(0,1) 
                
        # sum the utility - sum of factor disutility
         total_utility = np.sum(goods_utility) - np.sum(factors_utility)
        
         return total_utility
    
    def invertedutility(self,input):
        # inverts the utility function for minimisation
        return (-1) * self.utility(input)
    
    def price(self):
        #generate the price array for each good
        for x in range(self.ngp):
            price = np.array(input[random.randint(1,1000)])
        return price 
            
    def constraint(self,input, pi, r):
        # pr: profit earned from shares, r: price of factors
        
        # separate the goods from the factors
        goods = np.array(input[0:self.gl])
        
        # total amount paid for goods
        goods_paid = np.empty(self.nc * self.gl)
        for i in goods:
            goods_paid[i] = self.price[i] * goods[i]
            
        total_paid = np.sum(goods_paid)
        
        # total amount earned from factors
        
        #An array storing each conumser's labor unit 
        factor_labor_earned = np.empty(self.nc * self.lc )
        
        for l in factor_labor_earned:
            factor_labor_earned[l] += w * factor_labor_earned[l]
        
        #An array storing each conusmer's capital unit 
        factor_capital_earned = np.empty(self.nc * self.kc)
        
        for k in factor_capital_earned:
            factor_capital_earned[k] += r * factor_capital_earned[l]
        
        total_earned = np.sum(factor_capital_earned + factor_labor_earned)
        
        budget = pi + total_earned - total_paid
        
        return budget 
    
    def MaxUtility(self,pi,r,guess):
        # for maximisation of utility function
        # please input the maximisation properly, i am not sure how it is supposed to be
        
        budgetCon = {'type' : 'eq', 'fun' : self.constraint, 'args' : (pi,r,w)}
        constraint = [budgetCon]
        solution = minimize(self.invertedUtility, guess, args = (pi,r,w), method = 'SLSQP', constraints = constraint)
        
        return solution.x
     
    def consumption_quantity(self):
        
        
     

# The second part:Producer
        

     def production(self,input):
          factors_labor = np.array(input[self.lp * self.npro])
          factors_capital = np.array(input[self.kp * self.npro])
        
        # giving a matrix of (G1,G2) where G1 is number of G1 produced
          production = np.empty(self.gl)
        # where each producer uses 2 factors to produce, simulate k & l
        
        
          if self.externality == False:
        # Production of goods under no externality in an economy
              
           for i in range(0,self.gl):
             production[i] = self.T * (factors_labor[i] ** self.A) * (factors_capital[i] ** (1-self.A))
        # The profit array under no externality
        
           for i in production:
            total_profit_noe=0
            total_profit_noe += production[i] * price[i] - (w * factors_labor[i]) - (r * factors_capital[i])
 
            return total_profit_noe

            def maximizedprofit_noe(self,w,r):
                factors_labor = np.array(input[self.lp * self.npro])
                factors_capital = np.array(input[self.kp * self.npro])
        # creating the maximized profit array
                maximized_profit_noe = np.empty(self.npro)
        #Fill this array this each firm's maximized profit via F.O.C
                for j in self.production_noexternality:
        # y = MPL*l + MPK *K
                    maximized_profit_noe[i] =((diff(total_profit_noe[i],w))*factors_labor[i] + (diff(total_profit_noe[i],r)*factors_capital[i] )
        #The maximized profit array 
                   
        
        #Return the quantity of output under maximization for each producer
        
             def total_output_noe(self,w,r):
                 total_output_noe = np.empty(self.npro)
                 for j in self.maximizedprofit_noe:
                      total_output_noe[j] = (maximized_profit_noe[j] + (w * factors_labor[j]) + (r * factors_capital[j]))/self.price[j]
                      return np.sum(total_output_noe)
          
          else:
        #Production of goods under externality in an economy
     
            for j in range(0,self.gl):
                production[j] = self.T*((0.5)*factors_labor[j]**self.A)*(factors_capital[j]**(1-self.A))
        #Profit array under externlaity 
        
            for j in self.production_externality:
                total_profit_ex=0
                total_profit_ex += self.production_externality[j]*price[j]-(w*factors_labor[j])-(r*factors_capital[j])
                
                return total_profit_noe
        
    
                def maximizedprofit_externality(self,w,r):
                     factors_labor = np.array(input[self.lp * self.npro])
                     factors_capital = np.array(input[self.kp * self.npro])
        # creating the maximized profit array
                     maximized_profit_e = np.empty(self.npro)
        #Fill this array this each firm's maximized profit via F.O.C
                     for j in self.production_externality:
        # y = MPL*l + MPK *K
                         maximized_profit_e[i] =((diff(total_profit_e[i],w))*factors_labor[i] + (diff(total_profit_e[i],r)*factors_capital[i]))
        #The maximized profit array 
                          return maximized_profit_e
        
        #Return the total output under maximized profit
                 def total_output_e(self,w,r):
                     total_output_e = np.empty(self.npro)
            
                      for j in self.maximizedprofit_externality:
                          total_output_e[j] = (maximized_profit_e[j] + (w * factors_labor[j]) + (r * factors_capital[j]))/self.price[j]
                          return np.sum(total_output_e)
                      
                 def total_factors_used_e(self,w,r):
                     total_factors_labor_e = np.empty(self.npro)
                      
                      for i in self.maximizedprofit_externality:
                          total_factors_labor_e[i] = (maximized_profit_e[i] + (r * factors_capital[i]))
                     

class Equlibrium(Economy):
    def _init_(self,nc,npro,ngp,fc,T,A,fp,externality):
        # nc: no of consumers in an economy
        # npro: no of producers in an economy
        # ngp: no of goods produced by each producer
        # fc: no of factors each consumer has
        # T: technology for a producer
        # A: alpha, the production parameters for a producer
        # fp: no of factors use by producer
        # externality: whether the externality exists or not 
        # p:random input a price list for each good 
        
        self.nc = nc
        self.npro = npro
        self.ngp = ngp
        self.fc = fc
        self.T = T
        self.A = A
        self.fp = 2
        self.externality = externality 
        goods_length = (self.npro * self.ngp)
        self.gl = goods_length
        
#What mentioned below are only the contour of the idea to reach the equilibrium(excess demand =0), not the precise notataion and code
    
    def excess_demand(w,r,goods_consumed,total_output):
        return goods_consumed - total_output
        price = optimize.fsolve(lambda x:excess_demand(w,r,goods_consumed,total_output))
        print price
    
    def excess_demand_onlabor(w,r,labor_supplied,labor_demanded):
        return labor_demanded - labor_supplied
        wage = optimize.fsolve(lambda x:excess_demand_onlabor(w,r,labor_supplied,labor_demanded))
        print wage
    
    def excess_demand_oncapital(w,r,capital_supplied,capital_demanded):
        return capital_demanded - capital_supplied
        rental = optimize.fsolve(lambda x:excess_demand_oncapital(w,r,capital_supplied,capital_demanded))
        print rental
     
         
    
             

        