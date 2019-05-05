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
        # ngp: no of goods produced by each producer
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

#The second part:Producer
  if self.externality == False:
                           
        # Production of goods under no externality in an economy
      def production_noexternality(self, input):
          factors_labor = np.array(input[self.lp * self.npro])
          factors_capital = np.array(input[self.kp * self.npro])
        
        # giving a matrix of (G1,G2) where G1 is number of G1 produced
         production_noe = np.empty(self.gl)
        # where each producer uses 2 factors to produce, simulate k & l
         for i in range(0,self.gl):
             production_noe[i] = self.T * (factors_labor[i] ** self.A) * (factors_capital[i] ** (1-self.A))
        # The profit array
        
         for i in production_noe:
            total_profit_noe=0
            total_profit_noe += production_noe[i] * price[i] - (w * factors_labor[i]) - (r * factors_capital[i])
 
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
            return maximized_profit_noe
            
   else:
      def production_externality(self,input):
          factors_labor = np.array(input[self.lp * self.npro])
          factors_capital = np.array(input[self.kp * self.npro])
        
        #Production of goods under externality in an economy
        production_e = np.empty(self.gl)
        
        #Production function with externality
        for j in range(0,self.gl):
            production_e[j] = self.T*((0.5)*factors_labor[j]**self.A)*(factors_capital[j]**(1-self.A))
        #Profit array under externlaity 
        
        total_profit_ex=0
        
        for j in self.production_externality:
            total_profit_e += self.production_externality[j]*price[j]-(w*factors_labor[j])-(r*factors_capital[j])
        
    
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
        
     def objective(self,w,r):
         
         
    
             

        