



class Economics_Agents(object):

    import random
    
    #Setting the fixed variables
    w = raw_input("Enter the wage: ")
    px = raw_input("Enter the price of X: ")
    py = raw_input("Enter the price of Y: ")
    
#First step: production side#
    
    psi = float(random.uniform(0, 1))
    print "The psi is " + str(psi)
    
    def __init__(self, L1, L2, X1, X2, Y1, Y2):
        self.L1 = L1 #Labour of first consumer
        self.L2 = L2 #Labour of second consumer
        self.X1 = X1 #quantity of X consumed by consumer 1
        self.X2 = X2 #quantity of X consumed by consumer 2
        self.Y1 = Y1 #quantity of Y consumed by consumer 1
        self.Y2 = Y2 #quantity of Y consumed by consumer 2
    
    
    def profit_ftnX(self,L1, px, w):
        profitX = px * self.L1 ** self.psi - w * self.L1
        return profitX
        print "The profit for X is " + str(profitX)
    
    def profit_ftnY(self, L2, py, w):
        profitY = py * self.L2 ** self.psi - w * self.L2
        return profitY
        print "The profit for Y is " + str(profitY)

    from sympy import diff
    
    Pi_1_Star = diff(profit_ftnX(L1, px, w), L1) #invalid due to unknown L1
    Pi_2_Star = diff(profit_ftnY(L2), L2) #invalid due to unknown L2
    
    labor_1 = ((Pi_1_Star + w)/(psi*px)) ** (1/(psi-1))
    print "The labor input for 1 when maximized profit is" + str(labor_1)
    
    labor_2 =  ((Pi_2_Star + w)/(psi*py)) ** (1/(psi - 1))
    print "The labor input for 2 when maximized profit is" + str(labor_2)


#what listed above are production side#

#Second step : consumption side#

    def income_1(self):
         income1 =(the_optimal_labor_input1(L1)) * self.w #invalid due to unknown L1
         return income1
         print "The budget for 1 is " + str(income1)

    def income_2(self):
         income2 =(the_optimal_labor_input(L2)) * self.w #invalid due to unknown L2
         return income2
         print "The budget for 2 is" + str(income2)

    
    import numpy as np
    from scipy.optimize import minimize

    alpha = float(random.uniform(0, 1))

    def objective1(self, alpha,  X1,  Y1):
        
        return ((float(1.5)*X1)**alpha) * (Y1**(1-alpha))

    def constraint1(self, px, X1, py, Y1, income1):
        
        return (self.px * self.X1 + self.py * self.Y1) - income1

    optimal_utility1= minimize(objective1, method = 'SLSQP',constraints = constraint1)

    print optimal_utility1

    def objective2(self, alpha, X2, Y2):
        
        return ((self.X2) ** alpha) * ((float(1.5) * (self.Y2)) ** (1-alpha))

    def constraint2(self, px, X2, py, Y2, income2):
        
        return (px * (self.X2) + py * (self.Y2))-income2
    
    optimal_utility2= minimize(objective2, method = 'SLSQP',constraints = constraint2)

    print optimal_utility2
