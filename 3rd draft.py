import random
    
#Setting the fixed variables
w = raw_input("Enter the wage: ")
px = raw_input("Enter the price of X: ")
py = raw_input("Enter the price of Y: ")

#1st Step: production side#

psi = float(random.uniform(0, 1))

print "The psi is " + str(psi)

def maximized_profit1(Px,L1,psi,w):
    profit1 = Px*((L1)**psi)-w*L1

    from sympy import diff

    FOC1= diff(profit1, L1)

    return FOC1
 
def maximized_profit2(Py,L2,psi,w):
    profit2 = Py*((L2)**psi)-w*L2

    FOC2= diff(profit2, L2)

    return FOC2

def maximized_labor_input1(FOC1,w,psi,px):
    labor_input1 = ((FOC1+w)/(psi*px)) **(1/(psi-1))

    return labor_input1

def maximized_labor_input2(FOC2,w,psi,py):
    labor_input2 = ((FOC2+w)/(psi*py)) **(1/(psi-1))

    return labor_input2

#Each agent's total income#

def income_1(labor_input1,w):
    budget1= (labor_input1)*w

    return budget1

def income_2(labor_input2,w):
    budget2= (labor_input2)*w

    return budget2

import numpy as np
from scipy.optimize import minimize

alpha = float( random.uniform(0,1))
    
#Maximized agent1's utilities given budget constraints#

def objective1(alpha,x1,y1):
    utility1 = ((float(1.5)*x1)**alpha)*((y1)**(1-alpha))

    return utility1

def constraint1(budget1,px,x1,py,y1):
    budget_constraint1 = (px*x1 + py*y1) - budget1

    return budget_constraint1

    def optimal_utilities_1(utility1,budget_constraint1):
        maximized_utilities_1 =  minimize(utility1, method = 'SLSQP', constraints = budget_constraint1)
        
    return maximized_utilities_1

#Maximized agent2's utilities given budget constraints#

def objective2(alpha,x2,y2):
    utility2 = ((x2)**alpha)*((float(1.5)*y2)**(1-alpha))
        
    return utility2
    
def constraint2(budget2,px,x2,py,y2):
    budget_constraint2 = (px*x2 + py*y2)- budget2

    return budget_constraint2

def optimal_utilities_2(utility2,budget_constraint2):
    maximized_utilities_2 = minimize(utility2, method = 'SLSQP' , constraints = budget_constraint2)
        
    return  maximized_utilities_2

    
    
    
