from potentialSolver import *
from methods import *
import math

def generateAandb(mesh, numNode, innerPotential, outerPotential):
    A = [[-4 if a == b else 0 for a in range(numNode)] for b in range(numNode)]
    b = [0 for a in range(numNode)]
    k = 0
    for i in range(0, len(mesh) - 1):
        for j in range(0, len(mesh[0]) - 1):
            if j > 1 and mesh[i][j] == 0 and mesh[i][j - 1] == innerPotential:
                if i == 0:
                    A[k][k + 1] = 1
                    A[k][k + 2] = 2
                    b[k] = -innerPotential
                elif i == 1:
                    A[k][k + 1] = A[k][k - 2] = A[k][k + 5] = 1
                    b[k] = -innerPotential
                k += 1
            elif j + 2 == len(mesh[0]):
                if i == 0:
                    A[k][k - 1] = 1
                    A[k][k + 2] = 2
                    b[k] = -outerPotential
                elif i == 1:
                    A[k][k - 1] = A[k][k + 5] = A[k][k - 2] = 1
                    b[k] == 0
                elif i == len(mesh) - 2:
                    A[k][k - 1] = A[k][k - 5] = 1
                    b[k] = -outerPotential * 2
                else:
                    A[k][k - 1] = A[k][k + 5] = A[k][k - 5] = 1
                k += 1
            elif j == 0 and i > 1:
                if mesh[i - 1][j] == innerPotential:
                    A[k][k + 1] = 2
                    A[k][k + 5] = 1
                    b[k] = -innerPotential
                elif i + 2 == len(mesh):
                    A[k][k + 1] = 2
                    A[k][k - 5] = 1
                    b[k] = -outerPotential
                else:
                    A[k][k + 1] = 2
                    A[k][k + 5] = A[k][k - 5] = 1
                    b[k] = 0
                k += 1
            elif i == 2 and mesh[i - 1][j] == innerPotential:
                A[k][k - 1] = A[k][k + 1] = A[k][k + 5] = 1
                b[k] = -innerPotential
                k += 1
            elif i + 2 == len(mesh):
                A[k][k - 1] = A[k][k + 1] = A[k][k - 5] = 1
                b[k] = -outerPotential
                k += 1
            elif 1 < i and 1 <= j:
                A[k][k - 1] = A[k][k + 1] = A[k][k - 5] = A[k][k + 5] = 1
                b[k] = 0
                k += 1
    return A, b

def conjugateGradient(A, b, numNode):
    x = numColumnCheck([0 for a in range(numNode)])
    r = matrixAddOrSub(b, multiplyMatrix(A, x), 'sub')
    p = [0 for a in range(len(r))]
    for i in range(0, len(r)):
        p[i] = r[i]
    r = numColumnCheck(r)
    p = numColumnCheck(p)
    infNorm_ini = 0
    twoNorm_ini = 0
    print(r)
    for l in range(numNode):
        if abs(r[l][0]) > infNorm_ini:
            infNorm_ini = abs(r[l][0])
        twoNorm_ini += r[l][0] ** 2
    twoNorm_ini = math.sqrt(twoNorm_ini)
    print("iteration number: 0" + "  infinity norm: " + str(infNorm_ini) + "  2-norm: " + str(twoNorm_ini))
    for k in range(numNode):
        alpha = multiplyMatrix(transposeMatrix(p), r)[0][0]/(multiplyMatrix(transposeMatrix(p), multiplyMatrix(A, p))[0][0])
        x = matrixAddOrSub(x, scalarmultiplier(alpha, p), 'add')
        r = matrixAddOrSub(b, multiplyMatrix(A, x), 'sub')
        beta = -((multiplyMatrix(transposeMatrix(p), multiplyMatrix(A, r))[0][0])/(multiplyMatrix(transposeMatrix(p), multiplyMatrix(A, p))[0][0]))
        p = matrixAddOrSub(r, scalarmultiplier(beta, p), 'add')
        # finding the norms
        infNorm = 0
        twoNorm = 0
        for j in range(numNode):
            if abs(r[j][0]) > infNorm:
                infNorm = abs(r[j][0])
            twoNorm += r[j][0] ** 2
        twoNorm = math.sqrt(twoNorm)
        print("iteration number: " + str(k + 1) + "  infinity norm: " + str(infNorm) + "  2-norm: " + str(twoNorm))
    return x

h = 0.02
innerPotential = 110
outerPotential = 0
potentials = potentialMesh(h, 0)
mesh = potentials.mesh
print(mesh)
numNode = 19
(A, b) = generateAandb(mesh, numNode, innerPotential, outerPotential)
print(b)
for n in range(0, numNode):
    print(A[n])
choleskiTest = choleski(A, b)
Afinal = multiplyMatrix (transposeMatrix(A), A)
bfinal = multiplyMatrix(transposeMatrix(A), b)
conjugateSolution = conjugateGradient(Afinal, bfinal, numNode)
choleskiOutput = choleski(Afinal, bfinal)
choleskiSolution = backwardElim(choleskiOutput[0], choleskiOutput[1])
print("Choleski result of the potential equals " + str(choleskiSolution[11]) + " V" )
print("Conjugate result of the potential equals " + str(conjugateSolution[11][0]) + " V" )



