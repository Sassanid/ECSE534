from methods import *
from scipy import random

# allow the user to specify the matrix size
matrixSize = int(input("Please enter the size of the testing matrix: "))
n = int(input("Please enter the multiplier of the random function(since random only generates value from [0, 1)): "))
A = symmetricMatrix(matrixSize, n)
x = [n*random.random() for i in range(matrixSize)]
b = multiplyMatrix(A, x)
print('the x we generated is:')
for line in x:
    print(line)
print('A = ')
for line in A:
    print(line)
print('b =')
for line in b:
    print(line)
choleskiOutput = choleski(A, b)
solution = backwardElim(choleskiOutput[0], choleskiOutput[1])
print('solution = ')
for line in solution:
    print(line)
