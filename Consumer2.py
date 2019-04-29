#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 14:18:52 2019

@author: JKong
"""

import numpy as np
import random
from scipy.optimize import minimize

class consumer(object):
    
    def _init(self,nc,np,ngp,fc):
        # nc: no of consumers in an economy
        # np: no of producers in an economy
        # ngp: no of goods produced by each producer
        # fc: no of factors each consumer has
        
        self.nc = nc
        self.np = np
        self.ngp = ngp
        self.fc = fc
        goods_length = self.ngp * self.np
        self.gl = goods_length
        
        
    def utility(self,input):
        
        # separate the goods from the factors
        goods = np.array(input[0:self.gl])
        factors = np.array(input[self.gl:(self.gl + self.fc * self.nc)])
        # utility of the goods consumed in a array form
        goods_utility = np.empty(self.nc * self.gl)
        
        for i in goods:
        # B will be a random variable from [0,1], append each consumer's utility from a certain good.
            goods_utility[i] = goods[i] ** random.randint(0,1)
    
        
        # factor disutility when providing k & l to producers for production
        factors_utility = np.empty(self.nc * self.fc * self.nc)
        for i in factors:
            factors_utility[i] = random.randint(0,1) * factors[i] ** random.randint(0,1) 
                
        # sum the utility - sum of factor disutility
        total_utility = np.sum(goods_utility) - np.sum(factors_utility)
        
        return total_utility
    
    def invertedutility(self,input):
        # inverts the utility function for minimisation
        return (-1) * self.utility(input)
        
        
    def constraint(self,input, pi, p, r):
        # pr: profit earned from shares, p: price of goods, r: price of factors
        
        # separate the goods from the factors
        goods = np.array(input[0:self.gl])
        factors = np.array(input[self.gl:(self.gl + self.f)])
        
        # total amount paid for goods
        goods_paid = np.empty(self.nc * self.gl)
        for i in goods:
            goods_paid[i] = p * goods[i]
            
        total_paid = np.sum(goods_paid)
        
        # total amount earned from factors
        factors_earned = np.empty(self.nc * self.f)
        for i in factors:
            factors_earned[i] = r * factors[i]
        total_earned = np.sum(factors_earned)
        
        budget = pi + total_earned - total_paid
        
        return budget
    
    def MaxUtility(self,pi,p,r,guess):
        # for maximisation of utility function
        # please input the maximisation properly, i am not sure how it is supposed to be
        
        budgetCon = {'type' : 'eq', 'fun' : self.constraint, 'args' : (pi,p,r,)}
        constraint = [budgetCon]
        solution = minimize(self.invertedUtility, guess, args = (pi,p,r), method = 'SLSQP', constraints = constraint)
        return solution.x


        
        
        
            
                
        
        