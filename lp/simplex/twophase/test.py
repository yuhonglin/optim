from scipy.optimize import linprog

c = [1, 2, 1, 0, 0]

A = [[1,1,0,1,0], [0,-1,1,0,1]]

b = [-1,3]

lb = -2

ub = 3

print linprog(c, None, None, A, b, (lb,ub))
