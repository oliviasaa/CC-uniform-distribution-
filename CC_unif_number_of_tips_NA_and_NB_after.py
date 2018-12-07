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
    

for count1 in range(1,25):
    lambdh = float(count1)
    N = int(3*count1)



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
                    if emptybins == 0 :
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

    def Pn(totalballs):
        P = math.exp(-lambdh)*lambdh**totalballs/math.factorial(totalballs)
        return P

    def tra(nA,nB,NA,NB):
        if NA+NB == 0:
            r = 0
        else:
            p1 = (float(NA**2)/float((NA+NB)**2-2*NA*NB))
            p2 = (float(NB**2)/float((NA+NB)**2-2*NA*NB))
            r = binomial(nA+nB,nA)*(p1)**(nA)*(p2)**(nB)
        return r


    P_next = []
    for N_before in range(N):
        P_next_aux = [0 for i in range(N)]
        for i in range(N):
            for j in range(i+1):
                P_next_aux[i] = P_next_aux[i]+Pn(i-j)*shots[j][N_before][2*(i-j)]
        P_next.append(P_next_aux)
    P_next[0][0] = 1
    M = numpy.matrix(P_next)
    M500 = numpy.linalg.matrix_power(M, 500)
    P500 = M500[2].tolist()
    P500 = list(P500[0])

    a = sum(P500)
    for i in range(len(P500)):
        P500[i] = P500[i]/a

    matrix = [[[[0 for i in range(N)] for count in range(N)] for j in range(N)] for k in range(N)]
    for NA in range(1,N):
        print NA
        for NB in range(1,N):
            P_next = [[0 for i in range(N)] for j in range(N)]
            summ = 0
            for i in range(N):   #NA
                for j in range(N):    #NB
                    for l in range(N+1): #nA
                        for m in range(N-l+1): #nB
                            if i>=l and j>=m and 2*l<2*N+1 and 2*m<2*N+1:
                                P_next[i][j] = P_next[i][j]+Pn(l+m)*tra(l,m,NA,NB)*shots[i-l][NA][2*l]*shots[j-m][NB][2*m]
                    matrix[NA][NB][i][j] = P_next[i][j]


    for NA in range(1,N):
        P_next = [[0 for i in range(N)] for j in range(N)]
        summ = 0
        for i in range(N):   #NA
            for j in range(N):    #NB
                for l in range(N+1): #nA
                    for m in range(N-l+1): #nB
                        if i>=l and j==m and 2*l<2*N+1 and 2*m<2*N+1:
                            P_next[i][j] = P_next[i][j]+Pn(l+m)*tra(l,m,NA,1)*shots[i-l][NA][2*l]
                matrix[NA][0][i][j] = P_next[i][j]


    for NB in range(1,N):
        P_next = [[0 for i in range(N)] for j in range(N)]
        summ = 0
        for i in range(N):   #NA
            for j in range(N):    #NB
                for l in range(N+1): #nA
                    for m in range(N-l+1): #nB
                        if i==l and j>=m and 2*l<2*N+1 and 2*m<2*N+1:
                            P_next[i][j] = P_next[i][j]+Pn(l+m)*tra(l,m,1,NB)*shots[j-m][NB][2*m]
                matrix[0][NB][i][j] = P_next[i][j]

    matrix[0][0][0][0] = 1



    def multiply(list1,N):
        aux2 = [[[[0 for m in range(N)] for j in range(N)] for k in range(N)] for l in range(N)]
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    for l in range(N):
                        for m in range(N):
                            for n in range(N):
                                aux2[i][j][k][l] = aux2[i][j][k][l]+ list1[i][j][m][n]*list1[m][n][k][l]
        return aux2



    def normalize(list1,N):
        for i in range(N):
            for j in range(N):
                summ = 0
                for l in range(N):
                    for m in range(N):
                        summ = summ + list1[i][j][l][m]
                for l in range(N):
                    for m in range(N):
                        list1[i][j][l][m] = list1[i][j][l][m]/summ
        return list1



    def check_matrix(list1,N):
        summ = 1.0
        for i in range(N):
            for j in range(N):
    #            print sum([sum(list1[i][j][l])  for l in range(N)])
                summ = summ * sum([sum(list1[i][j][l])  for l in range(N)])
        if summ-1.0<-0.00000001 or summ-1.0>0.00000001:
            return 0
        else:
            return 1






    matrix = normalize(matrix,N)
    check_matrix(matrix,N)

    matrix2 = multiply(matrix,N)
    matrix2 =  normalize(matrix2,N)
    check_matrix(matrix2,N)

    matrix4 = multiply(matrix2,N)
    matrix4 =  normalize(matrix4,N)
    check_matrix(matrix4,N)

    matrix8 = multiply(matrix4,N)
    matrix8 =  normalize(matrix8,N)
    check_matrix(matrix8,N)

    matrix16 = multiply(matrix8,N)
    matrix16 =  normalize(matrix16,N)
    check_matrix(matrix16,N)

    matrix32 = multiply(matrix16,N)
    matrix32 =  normalize(matrix32,N)
    check_matrix(matrix32,N)

    matrix64 = multiply(matrix32,N)
    matrix64 =  normalize(matrix64,N)
    check_matrix(matrix64,N)

    matrix128 = multiply(matrix64,N)
    matrix128 =  normalize(matrix128,N)
    check_matrix(matrix128,N)

    matrix256 = multiply(matrix128,N)
    matrix256 =  normalize(matrix256,N)
    check_matrix(matrix256,N)

    matrix512 = multiply(matrix256,N)
    matrix512 =  normalize(matrix512,N)
    check_matrix(matrix512,N)

    matrix1024 = multiply(matrix512,N)
    matrix1024 = normalize(matrix1024,N)
    check_matrix(matrix1024,N)

    matrix2048 = multiply(matrix1024,N)
    matrix2048 = normalize(matrix2048,N)
    check_matrix(matrix2048,N)

    matrix4096 = multiply(matrix2048,N)
    matrix4096 = normalize(matrix4096,N)
    check_matrix(matrix4096,N)


    stationary = matrix1024[1][0]


    def P_unconfirm(CC, list1, N, stationary):
        max = 0
        if CC<0.5:
            return -1
        else:
            for i in range(N):
                for j in range(N):
                    sum1 = 0
                    if (i+j)>0:
                        if float(i)/float(i+j) > CC:
                            for k in range(N):
                                for l in range(N):
                                    if (k+l)>0:
                                        if float(l)/float(k+l) > CC:
                                            sum1 = sum1 + list1[i][j][k][l]
                    if sum1 > max:
                        max = sum1
            if max==0:
                return -1
            else:
                return max



    for i in range(20):
        print "CC=", 0.05+0.05*i, "t=1", P_unconfirm(0.05+0.05*i, matrix, N, stationary)

    for i in range(20):
        print "CC=", 0.05+0.05*i, "t=2", P_unconfirm(0.05+0.05*i, matrix2, N, stationary)

    for i in range(20):
        print "CC=", 0.05+0.05*i, "t=4", P_unconfirm(0.05+0.05*i, matrix4, N, stationary)

    for i in range(20):
        print "CC=", 0.05+0.05*i, "t=8", P_unconfirm(0.05+0.05*i, matrix8, N, stationary)

    for i in range(20):
        print "CC=", 0.05+0.05*i, "t=16", P_unconfirm(0.05+0.05*i, matrix16, N, stationary)

    for i in range(20):
        print "CC=", 0.05+0.05*i, "t=32", P_unconfirm(0.05+0.05*i, matrix32, N, stationary)

    for i in range(20):
        print "CC=", 0.05+0.05*i, "t=64", P_unconfirm(0.05+0.05*i, matrix64, N, stationary)

    for i in range(20):
        print "CC=", 0.05+0.05*i, "t=128", P_unconfirm(0.05+0.05*i, matrix128, N, stationary)

    for i in range(20):
        print "CC=", 0.05+0.05*i, "t=256", P_unconfirm(0.05+0.05*i, matrix256, N, stationary)

    for i in range(20):
        print "CC=", 0.05+0.05*i, "t=512", P_unconfirm(0.05+0.05*i, matrix512, N, stationary)

    for i in range(20):
        print "CC=", 0.05+0.05*i, "t=1024", P_unconfirm(0.05+0.05*i, matrix1024, N, stationary)




