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
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

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
    

    ##
    ## FIND AND PLOT THE STATIONARY PROBABILITY OF THE MARKOV CHAIN
    ##
    
    option = 1
    option2 = 3
    
    if option == 1:
        n_tips = [i for i in range(N)]
        M = numpy.matrix(P_next)
        M500 = numpy.linalg.matrix_power(M, 500)
        P500 = M500[2].tolist()
        P500 = P500[0]
        average_tips.append(sum([i*P500[i] for i in range(len(P500))]))
        if option2 == 0:
            plt.xlabel(r'Number of tips - $n$')
            plt.ylabel('Stationary Distribution')
            if count == 1:
                print lambdh
                plt.plot(n_tips,P500, label=r'$\lambda h = 1$')
                plt.legend()
            elif count == int(N/2)-1:
                plt.plot(n_tips,P500, label=r'$\lambda h = 24$')
                plt.legend()
                print lambdh
            elif count == int((int(N/2)-1)/2):
                plt.plot(n_tips,P500, label=r'$\lambda h = 12$')
                plt.legend()
                print lambdh
            else:
                plt.plot(n_tips,P500, '0.5')

        if option2 == 1:
            for i in range(len(P500)):
                n_tips[i] = (n_tips[i]-0.72364808551)/float(1.24786770351*lambdh)
            plt.plot(n_tips,P500)
            #     plt.plot([1.24786770351,1.24786770351],[0,1],linestyle='dashed')
            plt.xlabel(r'Reescaled Number of Tips - $(n-0.723)/(1.247\lambda h)$')
            plt.ylabel('Stationary Distribution')
            if count == 1:
                print lambdh
                plt.plot(n_tips,P500, label=r'$\lambda h = 1$')
                plt.legend()
            elif count == int(N/2)-1:
                plt.plot(n_tips,P500, label=r'$\lambda h = 24$')
                plt.legend()
                print lambdh
            elif count == int((int(N/2)-1)/2):
                plt.plot(n_tips,P500, label=r'$\lambda h = 12$', )
                plt.legend()
                print lambdh
            else:
                plt.plot(n_tips,P500, '0.5')



    ##
    ## PLOT ORBITS
    ##
    
    if option == 2:
        accum_P = [[sum(P_next[i][:j+1]) for j in range(len(P_next[i]))] for i in range(len(P_next))]
        j = 0
        for k in range(10):
            N0 = random.randrange(1, min(4*lambdh,N))
            v_aux = []
            for i in range(100):
                p = random.random()
                j = 0
                while accum_P[N0][j] < p and j < len(accum_P[N0])-1:
                    j = j + 1
                N0 = j
                v_aux.append(N0)
            plt.plot(v_aux)
        plt.show()

if option2 == 3:
    plt.plot(x_axis,average_tips)
    plt.xlabel(r'$\lambda h$')
    plt.ylabel('Average number of tips')
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
plt.show()
