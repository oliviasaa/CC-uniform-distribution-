import math
import numpy
import networkx as nx
import time
import copy
import scipy.stats
from scipy.misc import logsumexp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D, get_test_data
from matplotlib import cm
import matplotlib

def binomial(n,k):
    """Compute n factorial by an additive method."""
    if k > n-k:
        k = n-k  # Use symmetry of Pascal's triangle
    thediag = [i+1 for i in range(k+1)]
    for i in range(n-k-1):
        for j in range(1,k+1):
            thediag[j] += thediag[j-1]
    return thediag[k]

for nn in range(1,25):
    lambdh = float(nn)
    N = 3*nn

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
            for emptybins in range(N,-1,-1):
                if totalbins == 0:
                    if emptybins == 0:
                        shots[emptybins][totalbins][number_of_shots] = 1
                elif number_of_shots == 0:
                    if emptybins == totalbins:
                        shots[emptybins][totalbins][number_of_shots] = 1.0
                    else:
                        shots[emptybins][totalbins][number_of_shots] = 0
                elif emptybins+number_of_shots>=totalbins:
                    if emptybins < totalbins:
                        shots[emptybins][totalbins][number_of_shots] = (1-float(emptybins)/float(totalbins))*shots[emptybins][totalbins][number_of_shots-1] + float(1+emptybins)/float(totalbins)*shots[emptybins+1][totalbins][number_of_shots-1]
    shots[0][0][0] = 1

    def tra(nA,nAB,nB,NA,NB):
        if NA+NB == 0:
            r = 0
        else:
            p1 = (float(NA)/float(NA+NB))**2
            p2 = (float(NB)/float(NA+NB))**2
            p3 = 1-p1-p2
            r = binomial(nA+nAB+nB,nA)*binomial(nB+nAB,nB)*(p1)**(nA)*(p2)**(nB)*(p3)**(nAB)
        return r

    def check_sum(aux):
        for i in range(len(aux)):
            for j in range(len(aux[i])):
                sum1 = 0
                for k in range(len(aux[i][j])):
                    sum1 = sum1 + sum(aux[i][j][k])
                if sum1>1.000000001 or sum1<0.9:
                    print "checksum:", "lambda*h:",lambdh,  "i=", i,"j=",  j,"sum=",  sum1
                sum2 = 0
                for k in range(N):
                    for l in range(N):
                        aux[i][j][k][l] = aux[i][j][k][l]/sum1
                sum2 = sum2 + sum(aux[i][j][k])

    def multiply(matrix):
        list_squared = []
        N = len(matrix[0])
        aux = [[[[0 for m in range(N)] for j in range(N)] for k in range(N)] for l in range(N)]
        for i in range(len(matrix)):
            for j in range(N):
                for k in range(N):
                    for l in range(N):
                        for m in range(N):
                            for n in range(N):
                                aux[i][j][k][l] = aux[i][j][k][l] + matrix[i][j][m][n]*matrix[m][n][k][l]
        return aux
    
    def Pn(totalballs):
        P = math.exp(-lambdh)*lambdh**totalballs/math.factorial(totalballs)
        return P

    P_next = []
    for N_before in range(N):
        P_next_aux = [0 for i in range(N)]
        for i in range(N):
            for j in range(i+1):
                P_next_aux[i] = P_next_aux[i]+Pn(i-j)*shots[j][N_before][2*(i-j)]
        P_next.append(P_next_aux)
    P_next[0][0] = 1
    M = numpy.matrix(P_next)
    accum_P = [[sum(P_next[i][:j+1]) for j in range(len(P_next[i]))] for i in range(len(P_next))]
    M500 = numpy.linalg.matrix_power(M, 500)
    P500 = M500[1].tolist()
    P500 = list(P500[0])
    a = sum(P500)
    for i in range(len(P500)):
        P500[i] = P500[i]/a

    matrix = [[[[0 for i in range(N)] for count in range(N)] for j in range(N)] for k in range(N)]
    for Na in range(1,N):
        print Na
        for Nb in range(N):
            P_next = [[0 for i in range(N)] for j in range(N)]
            summ = 0
            for i in range(N):
                for j in range(N):
                    for k in range(N+1): #nA
                        for l in range(N-k+1): #nAB
                            for m in range(N-k-l+1): #nB
                                if (i-k-l)>=0 and j-m>=0 and 2*k+l<2*N+1 and 2*m+l<2*N+1:
                                    P_next[i][j] = P_next[i][j]+Pn(k+l+m)*tra(k,l,m,Na,Nb)*shots[i-k-l][Na][2*k+l]*shots[j-m][Nb][2*m+l]
                    matrix[Na][Nb][i][j] = P_next[i][j]
                    summ  = summ + matrix[Na][Nb][i][j]

    for Nb in range(1,N):
        summ = 0
        P_next = [[0 for i in range(N)] for j in range(N)]
        for i in range(N):
            for j in range(N):
                for k in range(N+1): #nA
                    for l in range(N-k+1): #nAB
                        for m in range(N-k-l+1): #nB
                            if (i-k-l)==0 and j-m>=0 and 2*k+l<2*N+1 and 2*m+l<2*N+1:
                                P_next[i][j] = P_next[i][j]+Pn(k+l+m)*tra(k,l,m,1,Nb)*shots[j-m][Nb][2*m+l]
                matrix[0][Nb][i][j] = P_next[i][j]
                summ  = summ + matrix[0][Nb][i][j]

    matrix[0][0][0][0] = 1
    check_sum(matrix)



    def min_P_going_to_NB_0(matrix):
        print "P go to NB=0", lambdh
        min = 1
        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                sum1 = 0
                for k in range(len(matrix[i][j])):
                    sum1 = sum1 + matrix[i][j][k][0]
                if sum1 < min:
                    min = sum1
        return min


    matrix2 =  multiply(matrix)
    check_sum(matrix2)
    matrix4 =  multiply(matrix2)
    check_sum(matrix4)
    matrix8 =  multiply(matrix4)
    check_sum(matrix8)
    matrix16 =  multiply(matrix8)
    check_sum(matrix16)
    matrix32 =  multiply(matrix16)
    check_sum(matrix32)
    matrix64 =  multiply(matrix32)
    check_sum(matrix64)
    matrix128 =  multiply(matrix64)
    check_sum(matrix128)
    matrix256 =  multiply(matrix128)
    check_sum(matrix256)
    matrix512 =  multiply(matrix256)
    check_sum(matrix512)
    matrix1024 =  multiply(matrix512)
    check_sum(matrix1024)

    print "1 second", min_P_going_to_NB_0(matrix)
    print "2 second", min_P_going_to_NB_0(matrix2)
    print "4 second", min_P_going_to_NB_0(matrix4)
    print "8 second", min_P_going_to_NB_0(matrix8)
    print "16 second", min_P_going_to_NB_0(matrix16)
    print "32 second", min_P_going_to_NB_0(matrix32)
    print "64 second", min_P_going_to_NB_0(matrix64)
    print "128 second", min_P_going_to_NB_0(matrix128)
    print "256 second", min_P_going_to_NB_0(matrix256)
    print "512 second", min_P_going_to_NB_0(matrix512)
    print "1024 second", min_P_going_to_NB_0(matrix1024)

