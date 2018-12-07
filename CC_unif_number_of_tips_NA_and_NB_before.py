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


for count1 in range (1,25):
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

    def tra(nA,nB,nN,nAN,nBN,NA,NB,NN):
        if NA+NB+NN == 0:
            r = 0
        else:
            p1 = (float(NA**2)/float((NA+NB+NN)**2-2*NA*NB))
            p2 = (float(NB**2)/float((NA+NB+NN)**2-2*NA*NB))
            p3 = (float(NN**2)/float((NA+NB+NN)**2-2*NA*NB))
            p4 = (float(2*NA*NN)/float((NA+NB+NN)**2-2*NA*NB))
            p5 = (float(2*NB*NN)/float((NA+NB+NN)**2-2*NA*NB))
            r = binomial(nA+nB+nN+nAN+nBN,nA)*binomial(nB+nN+nAN+nBN,nB)*binomial(nN+nAN+nBN,nN)*binomial(nAN+nBN,nAN)*(p1)**(nA)*(p2)**(nB)*(p3)**(nN)*(p4)**(nAN)*(p5)**(nBN)
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

    matrix = [[[[[[0 for i in range(N)] for count in range(N)] for j in range(N)] for k in range(N)] for l in range(N)] for m in range(N)]
    for NA in range(1,N):
        print NA
        for NB in range(1,N):
            for NN in range(0,N):
                P_next = [[[0 for i in range(N)] for j in range(N)]  for k in range(N)]
                summ = 0
                for i in range(N):   #NA
                    for j in range(N):    #NB
                        for k in range(N):   #NN
                            for l in range(N+1): #nA
                                for m in range(N-l+1): #nB
                                    for n in range(N-l-m+1): #nN
                                        for o in range(N-l-m-n+1): #nAN
                                            for p in range(N-l-m-n-o+1): #nBN
                                                if i>=l+o and j>=m+p and k>=n and 2*l+o<2*N+1 and 2*m+p<2*N+1 and 2*n+o+p<2*N+1:
                                                    P_next[i][j][k] = P_next[i][j][k]+Pn(l+m+n+o+p)*tra(l,m,n,o,p,NA,NB,NN)*shots[i-l-o][NA][2*l+o]*shots[j-m-p][NB][2*m+p]*shots[k-n][NN][2*n+o+p]
                            matrix[NA][NB][NN][i][j][k] = P_next[i][j][k]


    for NA in range(1,N):
        for NN in range(N):
            P_next = [[[0 for i in range(N)] for j in range(N)]  for k in range(N)]
            summ = 0
            for i in range(N):   #NA
                for j in range(N):    #NB
                    for k in range(N):   #NN
                        for l in range(N+1): #nA
                            for m in range(N-l+1): #nB
                                for n in range(N-l-m+1): #nN
                                    for o in range(N-l-m-n+1): #nAN
                                        for p in range(N-l-m-n-o+1): #nBN
                                            if i>=l+o and j==m+p and k>=n and 2*l+o<2*N+1 and 2*m+p<2*N+1 and 2*n+o+p<2*N+1:
                                                P_next[i][j][k] = P_next[i][j][k]+Pn(l+m+n+o+p)*tra(l,m,n,o,p,NA,1,NN)*shots[i-l-o][NA][2*l+o]*shots[k-n][NN][2*n+o+p]
                        matrix[NA][0][NN][i][j][k] = P_next[i][j][k]


    for NB in range(1,N):
        for NN in range(N):
            P_next = [[[0 for i in range(N)] for j in range(N)]  for k in range(N)]
            summ = 0
            for i in range(N):   #NA
                for j in range(N):    #NB
                    for k in range(N):   #NN
                        for l in range(N+1): #nA
                            for m in range(N-l+1): #nB
                                for n in range(N-l-m+1): #nN
                                    for o in range(N-l-m-n+1): #nAN
                                        for p in range(N-l-m-n-o+1): #nBN
                                            if i==l+o and j>=m+p and k>=n and 2*l+o<2*N+1 and 2*m+p<2*N+1 and 2*n+o+p<2*N+1:
                                                P_next[i][j][k] = P_next[i][j][k]+Pn(l+m+n+o+p)*tra(l,m,n,o,p,1,NB,NN)*shots[j-m-p][NB][2*m+p]*shots[k-n][NN][2*n+o+p]
                        matrix[0][NB][NN][i][j][k] = P_next[i][j][k]


    for NN in range(N):
        P_next = [[[0 for i in range(N)] for j in range(N)]  for k in range(N)]
        summ = 0
        for i in range(N):   #NA
            for j in range(N):    #NB
                for k in range(N):   #NN
                    for l in range(N+1): #nA
                        for m in range(N-l+1): #nB
                            for n in range(N-l-m+1): #nN
                                for o in range(N-l-m-n+1): #nAN
                                    for p in range(N-l-m-n-o+1): #nBN
                                        if i==l+o and j==m+p and k>=n and 2*l+o<2*N+1 and 2*m+p<2*N+1 and 2*n+o+p<2*N+1:
                                            P_next[i][j][k] = P_next[i][j][k]+Pn(l+m+n+o+p)*tra(l,m,n,o,p,1,1,NN)*shots[k-n][NN][2*n+o+p]
                    matrix[0][0][NN][i][j][k] = P_next[i][j][k]



    def multiply(list1,N):
        aux2 = [[[[[[0 for m in range(N)] for j in range(N)] for k in range(N)] for l in range(N)] for m in range(N)] for n in range(N)]
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    for l in range(N):
                        for m in range(N):
                            for n in range(N):
                                for o in range(N):
                                    for p in range(N):
                                        for q in range(N):
                                            aux2[i][j][k][l][m][n] = aux2[i][j][k][l][m][n] + list1[i][j][k][o][p][q]*list1[o][p][q][l][m][n]
        return aux2



    def normalize(list1,N):
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    summ = 0
                    for l in range(N):
                        for m in range(N):
                            for n in range(N):
                                summ = summ + list1[i][j][k][l][m][n]
                    for l in range(N):
                        for m in range(N):
                            for n in range(N):
                                list1[i][j][k][l][m][n] = list1[i][j][k][l][m][n]/summ
        return list1



    def check_matrix(list1,N):
        summ = 1.0
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    sum1 = 0
                    for l in range(N):
                        for m in range(N):
                            for n in range(N):
                                sum1 = sum1 + list1[i][j][k][l][m][n]
                    summ = summ * sum1
        if summ-1.0<-0.00000001 or summ-1.0>0.00000001:
            return 0
        else:
            return 1



    
    def min_P_going_to_NN_0(matrix):
        print "P go to NB=0", lambdh
        min = 1
        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                for k in range(len(matrix[i][j])):
                    sum1 = 0
                    for l in range(len(matrix[i][j][k])):
                        for m in range(len(matrix[i][j][k][l])):
                            sum1 = sum1 + matrix[i][j][k][l][m][0]
                    if sum1 < min:
                        min = sum1
        return min

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


    print "1 second", min_P_going_to_NN_0(matrix)
    print "2 second", min_P_going_to_NN_0(matrix2)
    print "4 second", min_P_going_to_NN_0(matrix4)
    print "8 second", min_P_going_to_NN_0(matrix8)
    print "16 second", min_P_going_to_NN_0(matrix16)
    print "32 second", min_P_going_to_NN_0(matrix32)
    print "64 second", min_P_going_to_NN_0(matrix64)
    print "128 second", min_P_going_to_NN_0(matrix128)
    print "256 second", min_P_going_to_NN_0(matrix256)
    print "512 second", min_P_going_to_NN_0(matrix512)
    print "1024 second", min_P_going_to_NN_0(matrix1024)

