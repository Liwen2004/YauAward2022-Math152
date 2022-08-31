import cvxpy as cp

print(cp.__version__)

N = 6        
T = 7   

d_n_t =(
           (    68, 	68, 	68, 	68, 	68, 	68, 	68),
           (    284,	279,	273,	268,	263,	257,	252),
           (    468,	462,	456,	449,	442,	436,	429),
           (    28, 	28, 	27, 	26, 	26, 	25, 	24),
           (    896,	884, 	872,	860,	847,	834,	821),
           (    2988,	3008,	3027,	3044,	3060,	3073,	3086)
           )
T_n_n = (    
           (0,2,3,1,1,2),
           (2,0,1,1,2,2),
           (3,1,0,2,4,5),
           (1,1,2,0,1,2),
           (1,2,4,1,0,1),
           (2,2,5,2,1,0)            
        )    
        
dist_rate = 700    

D_n_n = [\
           [T_n_n[i][j]*dist_rate for j in range(N)] for i in range(N)\
        ] 


Y_n =   (    # the initial inventory level in state n at time period t=0
           6575,
           4912,
           1076,
           9690,
           7585,
           16669           
        )     


gama_n = (     
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5
    )

tao_n =(    # the maximum fraction of its ventilator supply that each state is willing to share
        0.10,
        0.10,
        0.10,
        0.10,
        0.10,
        0.10        
        )
   
b_n = [ (1-gama_n[i])*Y_n[i] for i in range(N)]      
   
namda = 0.01
        
y_n_t = [    
           [0 for i in range(T)]
               for j in range(N)                     
        ]
        

x=cp.Variable( N*(N-1)*T,    
              integer=True)     
               #nonneg=True)     

print("Now,  set the computation method of the objective function...")

Z_n_n_t=[]     
W_n_t=[]       

y1_n_t = [
           [0 for i in range(T)]
               for j in range(N)                     
         ]

def objective_function(x1):
    global run_cycle
    global N  
    global T  
    global d_n_t  
    global T_n_n  
    global dist_rate
    global D_n_n        
    global Y_n
    global gama_n
    global tao_n   
    global b_n 
    global namda   
    global y_n_t
    global max_rate
    global Z_n_n_t
    global W_n_t
    global y1_n_t 
  
    print("N="+str(N)+" T="+str(T))
    Z_n_n_t = [  [
                    [0 for k in range(T)] \
                        for j in range(N) \
                 ]\
                            for i in range(N)\
              ]
    cursor_in_x = 0          
    for i in range(N):
        for j in range(N):
            for t in range(T):
                if(i != j):
                    Z_n_n_t[i][j][t] = x1[cursor_in_x]
                    cursor_in_x = cursor_in_x + 1
                    #print("hi")
            
    assert ( cursor_in_x == N*(N-1)*T ) 
   
    print("------- 1 ---------------")
    for i in range(N):
        y_n_t[i][0] = b_n[i]     
        y_n_t[i][0] = y_n_t[i][0] - sum( [Z_n_n_t[i][j][0] for j in range(N) ] )  
        y1_n_t[i][0] = y_n_t[i][0]    
        y_n_t[i][0] +=  sum(  [  (Z_n_n_t[k][i][  0 - T_n_n[k][i] ])  \
                                     if (0- T_n_n[k][i])>=0 \
                                         else 0\
                                             for k in range(N) ]\
                                 )

    print("------- 2 ---------------")
    for n in range(N):
        for t in range(1,T):     
            y_n_t[n][t] = y_n_t[n][t-1] - \
                          sum( [Z_n_n_t[n][k][t] for k in range(N) ] )  
            y1_n_t[i][t] = y_n_t[i][t]      
            y_n_t[n][t] +=  sum(\
                                 [   Z_n_n_t[k][n][ t-T_n_n[k][n] ]\
                                     if (t- T_n_n[k][n])>=0 \
                                         else 0\
                                             for k in range(N) ]\
                                 )
  
    print("------- 3 ---------------")
    W_n_t = [
               [ (d_n_t[n][t] - y_n_t[n][t])**2
                     for t in range(T) ]   \
                           for n in range(N)\
            ]

    print("------- 4 ---------------")
    ret = sum( [ sum( W_n_t[n] ) for n in range(N) ] )       
    ret += namda * sum( [ (D_n_n[i][j] * Z_n_n_t[i][j][t]) for t in range(T) for j in range(N) for i in range(N) ] )   
    print("------- return la ---------------")
    return ret
    

result = objective_function(x)
constrains = [x >= 0]
dynamic_constrains=[ (y1_n_t[i][t]>=0) for t in range(T) for i in range(N)]
dynamic_constrains.extend(  [  (y_n_t[i][t] >=  ( (1-tao_n[i]) * b_n[i] ))  \
                                         for t in range(T) for i in range(N) ]   );

constrains.extend(dynamic_constrains)

print("Now, start to run optimization, please wait...")
print(cp.installed_solvers())
prob = cp.Problem(cp.Minimize(result),constrains)
result = prob.solve(solver=cp.CPLEX)

print("solve() returns:",result)
print("status:", prob.status)     
print("optimal value", prob.value) 
print("optimal var:",x.value);

print("#-------------------------------------------------------------------------------------------")
print("#  everyday Demand and actual                                                               ")
print("#-------------------------------------------------------------------------------------------")
for t in range(T):
    print("--------------- day %d,[demand: actual supply] ---------------------" % t)
    for i in range(N):
        print("[%6d,%6d]" % (d_n_t[i][t],int(y_n_t[i][t].value)),end=" ")
    print("")
    

print("#-------------------------------------------------------------------------------------------")
print("#  TRANSFERRING MATRIX                                                                      ")
print("#-------------------------------------------------------------------------------------------")
print("optimized matrix for shift from n to n':")
for t in range(T):
    print("--------------- day %d ---------------------" % t)
    for i in range(N):
        for j in range(N):
            if( i!=j):
                print("%8d" % int(Z_n_n_t[i][j][t].value),end = "  ")
            else:
                print("%8d" % Z_n_n_t[i][j][t],end = "  ")
        print("")
    print("")
print("")
print("")
        
