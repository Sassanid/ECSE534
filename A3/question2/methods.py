import math
from scipy import random
import csv
# Function to check the number of columns of a matrix
def numColumnCheck (A):
    numOfColumuns = 0
    try:
        numOfColumuns = len(A[0])
        return A
    except TypeError:
        B = [[0] for a in range(len(A))]
        for i in range(0, len(A)):
            B[i][0] = A[i]
        return B

# Function to multiply a scalar and a matrix
def scalarmultiplier(a, A):
    A = numColumnCheck(A)
    B = [[0 for i in range(len(A[0]))]for k in range(len(A))]
    for i in range(len(A)):
        for j in range(len(A[0])):
            B[i][j] = a*A[i][j]
    return B

# Function to multiply two matrices
def multiplyMatrix (A, B):
    A = numColumnCheck(A)
    B = numColumnCheck(B)
    if len(A[0]) == len(B):
        C = [[0 for i in range(len(B[0]))]for k in range(len(A))]
        for i in range(len(A)):
            for j in range(len(B[0])):
                for k in range(len(A[0])):
                    C[i][j] += A[i][k]*B[k][j]
        return C
    else:
        print('cannot multiply this two matrices, incorrect dimensions')


# Function to transpose a matrix
def transposeMatrix (A):
   numOfRows = len(A)
   numOfColumns = len(A[0])
   C = [[0 for i in range(numOfRows)]for k in range(numOfColumns)]
   for i in range(numOfRows):
       for j in range(numOfColumns):
           C[j][i] = A[i][j]
   return C

# Function to create a symmetric matrix
def symmetricMatrix(size, n):
    A = [[0 for i in range(size)] for k in range(size)]
    # assign the lower part of A a value
    for i in range(len(A)):
        for j in range(0, i + 1):
            A[i][j] = n * random.random() - n
    B = transposeMatrix (A)
    C = multiplyMatrix (A, B)
    return C


# Function to use the choleski decomposition to find L and y
def choleski(A, b, halfBandwidth=None):
    A = numColumnCheck(A)
    b = numColumnCheck(b)
    if len(b[0])!= 1:
        print('invalid b input')
        return
    try:
        numOfColumuns = len(A[0])
    except TypeError:
        print('A only has one column')
        return
    if len(A) != len(A[0]):
        print('A is not a nxn matrix')
        return
    size = len(A)
    for j in range (size):
        if A[j][j] < 0:
            print("the matrix A is not positive definite")
            return
        A[j][j] = math.sqrt(A[j][j])
        b[j][0] = b[j][0]/A[j][j]
        for i in range (j+1, size):
            if halfBandwidth and i >= j + halfBandwidth:
                break
            A[i][j] = A[i][j]/A[j][j]
            b[i][0] = b[i][0]-A[i][j]*b[j][0]
            for k in range (j+1, i+1):
                if halfBandwidth and k >= j + halfBandwidth:
                    break
                A[i][k] = A[i][k]-A[i][j]*A[k][j]
    return [b,A]

# Function to find the solution through backward elimination, notice here L should be a lower matrix
def backwardElim(y, L):
    y = numColumnCheck(y)
    L = numColumnCheck(L)
    x = [0 for a in range(len(y))]
    for i in range(len(L)-1, -1, -1):
        for j in range(len(L)-1, i, -1):
            y[i][0] = y[i][0] - L[j][i]*x[j]
        x[i] = y[i][0] / L[i][i]
    return x

def matrixAddOrSub(A, B, option):
    A = numColumnCheck(A)
    B = numColumnCheck(B)
    if len(A)!= len(B) or len(A[0])!= len(B[0]):
        print('cannot add or subtract two matrices with different sizes!')
        return
    C = [[0 for a in range(len(A[0]))] for b in range(len(A))]
    if option == 'add':
        for i in range(0, len(A)):
            for j in range(0, len(A[0])):
                C[i][j] = A[i][j] + B[i][j]
    elif option == 'sub':
        for i in range(0, len(A)):
            for j in range(0, len(A[0])):
                C[i][j] = A[i][j] - B[i][j]
    return C


def getCircuit(r):
    with open('test_circuit.csv')as circuitData:
        reader = csv.reader(circuitData)
        for n in reader:
            if (n[0].startswith('#')):
                cirNumber = int(n[0].replace('#', ''))
                if cirNumber == r:
                    A_pre = n[1].split(';')
                    J_pre = n[2].split(';')
                    R_pre = n[3].split(';')
                    E_pre = n[4].split(';')
                    A = [0 for i in range(len(A_pre))]
                    for i in range(len(A)):
                        rowA_pre = A_pre[i].split(',')
                        rowA = []
                        for j in range(len(rowA_pre)):
                            rowA.append(int(rowA_pre[j]))
                        A[i] = rowA
                    J, E = [], []
                    y = [[0 for a in range(len(R_pre))] for b in range(len(R_pre))]
                    for i in range(len(J_pre)):
                        J.append(int(J_pre[i]))
                        E.append(int(E_pre[i]))
                        y[i][i] = 1/int(R_pre[i])
                    return [A, J, y, E]

def solveCircuitProblem(A,J, y, E, halfBandwidth=None):
    A = numColumnCheck(A)
    J = numColumnCheck(J)
    y = numColumnCheck(y)
    E = numColumnCheck(E)
    A_final = multiplyMatrix(A, multiplyMatrix (y, transposeMatrix(A)))
    b_final = multiplyMatrix(A,  matrixAddOrSub(J, multiplyMatrix(y, E), 'sub'))
    choleskiOutput = choleski(A_final, b_final, halfBandwidth)
    voltage = backwardElim(choleskiOutput[0], choleskiOutput[1])
    return voltage

