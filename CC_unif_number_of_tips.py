import math
import numpy
import networkx as nx
import time
import copy
import scipy.stats
from scipy.misc import logsumexp
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D, get_test_data
from matplotlib import cm
import matplotlib

N = 50

def binomial(n,k):
    """Compute n factorial by an additive method."""
    if k > n-k:
        k = n-k  # Use symmetry of Pascal's triangle
    thediag = [i+1 for i in range(k+1)]
    for i in range(n-k-1):
        for j in range(1,k+1):
            thediag[j] += thediag[j-1]
    return thediag[k]

C = [[1.0 for i in range(2*N+1)] for j in range(2*N+1)]
for i in range(len(C)):
    for j in range(len(C[i])):
        if i < j:
            C[i][j] = 0
for i in range(len(C)):
    for j in range(len(C[i])):
        if C[i][j] == 1:
            C[i][j] = binomial(i,j)

shots = [[[0 for i in range(2*N+1)] for j in range(N+1)] for k in range(N+1)]

for totalbins in range(N+1):
    for number_of_shots in range(2*N+1):
        for emptybins in range(N-1,-1,-1):
            if totalbins == 0:
                shots[emptybins][totalbins][number_of_shots] = 0
            elif number_of_shots == 0:
                if emptybins == totalbins:
                    shots[emptybins][totalbins][number_of_shots] = 1.0
                else:
                    shots[emptybins][totalbins][number_of_shots] = 0
            elif number_of_shots == 1:
                if totalbins-emptybins == 1:
                    shots[emptybins][totalbins][number_of_shots] = 1
                else:
                    shots[emptybins][totalbins][number_of_shots] = 0
            elif emptybins+number_of_shots>=totalbins:
                if emptybins < totalbins:
                    shots[emptybins][totalbins][number_of_shots] = (1-float(emptybins)/float(totalbins))*shots[emptybins][totalbins][number_of_shots-1] + float(1+emptybins)/float(totalbins)*shots[emptybins+1][totalbins][number_of_shots-1]


def Pn(totalballs):
    P = math.exp(-lambdh)*lambdh**totalballs/math.factorial(totalballs)
    return P

average_tips = []
x_axis = []
for count in range(1,N/2):
    x_axis.append(count)
    lambdh = float(count)
    P_next = []
    for N_before in range(N):
            P_next_aux = [0 for i in range(N)]
            for i in range(N):
                for j in range(i+1):
                    P_next_aux[i] = P_next_aux[i]+Pn(i-j)*shots[j][N_before][2*(i-j)]
            P_next.append(P_next_aux)
    P_next[0][0] = 1


    option = 1

    
    ##
    ## FIND THE STATIONARY PROPABILITY OF THE MARKOV CHAIN
    ##
    
    if option == 1:
        M = numpy.matrix(P_next)
        accum_P = [[sum(P_next[i][:j+1]) for j in range(len(P_next[i]))] for i in range(len(P_next))]
        #print M
        M500 = numpy.linalg.matrix_power(M, 500)
        P500 = M500[2].tolist()
        P500 = P500[0]
        #print M500
        average_tips.append(sum([i*P500[i] for i in range(len(P500))]))
        #print [sum(P_next[i]) for i in range(len(P_next))]


        
    ##
    ## PLOT ORBITS
    ##
        
    if option == 2:
        v = []
        v_average = [math.sqrt(2)*lambdh for i in range (100)]
        plt.plot(v_average)
        for k in range(10):
            N0 = random.randrange(1, 4*lambdh)
            v_aux = []
            for i in range(100):
                p = random.random()
                j = 0
                while accum_P[N0][j] < p and j < len(accum_P[N0]):
                    j = j + 1
                N0 = j
                v_aux.append(N0)
            plt.plot(v_aux)
            v.append(v_aux)
        p = 0.5
        v_aux = []
        for i in range(len(accum_P)):
            j = 0
            while accum_P[i][j] < p and j < len(accum_P[N0]):
                j = j + 1
            v_aux.append(j-i)
        print v_aux
        plt.show()



plt.plot(average_tips, x_axis)

##
## LEAST SQUARE METHOD TO FIT A LINE INTO THE GRAPH: EXPECTED NUMBER OF TIPS VS LAMBDA*H
##

a = float(len(x_axis))
b = float(sum(x_axis))
c = float(sum(average_tips))
d = float(sum([x_axis[i]*average_tips[i] for i in range(len(x_axis))]))
e = float(sum([x_axis[i]*x_axis[i] for i in range(len(x_axis))]))
f = float(sum([average_tips[i]*average_tips[i] for i in range(len(average_tips))]))

slope = (a*d-b*c)/(a*e-b*b)
a0 = (d-slope*e)/b

print a0, slope

#plt.plot([a0 + x_axis[i]*slope for i in range(len(x_axis))], x_axis)

plt.show()

