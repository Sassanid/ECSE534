from methods import *
import numpy as np
import math
E = 0.2
R = 512
Vt = 25 * 10**(-3)
Isa = 0.8 * 10**(-6)
Isb = 1.1 * 10**(-6)
k = 0
counter = 0


def newtonRap (V1, V2):
    f1 = Isa * (math.exp((V1 - V2) / Vt) - 1) - Isb * (math.exp(V2 / Vt) - 1)
    f2 = V1- E + R * Isb * (math.exp(V2 / Vt) - 1)
    f1V1Partial = Isa/Vt * math.exp((V1 - V2) / Vt)
    f1V2Partial = -Isa/Vt * math.exp((V1 - V2) / Vt) - Isb/Vt * math.exp(V2 / Vt)
    f2V1Partial = 1
    f2V2Partial = -Isb/Vt * math.exp(V2 / Vt)
    V = [V1, V2]
    f = [f1, f2]
    # F is the Jacobian Matrix
    F = [[f1V1Partial, f1V2Partial], [f2V1Partial, f2V2Partial]]
    invF = np.linalg.inv(F)
    V = matrixAddOrSub(scalarmultiplier(-1, multiplyMatrix(invF, f)), V,'add')
    f1_final = Isa * (math.exp((V[0][0] - V[1][0]) / Vt) - 1) - Isb * (math.exp(V[1][0] / Vt) - 1)
    f2_final = V[0][0] - E + R * Isb * (math.exp(V[1][0] / Vt) - 1)
    return f1_final, f2_final, V

# initial guess

V1 = 0
V2 = 0
f1_initial = Isa * (math.exp((V1 - V2) / Vt) - 1) - Isb * (math.exp(V2 / Vt) - 1)
f2_initial = V1 - E + R * Isb * (math.exp(V2 / Vt) - 1)
print("number of iteration: " + str(counter))
print("f1 = " + str(abs(f1_initial)) + "; " + "f2 = " + str(abs(f2_initial)))
print("V1 = " + str(V1) + "V" + "; " + "V2 = " + str(V2) + "V")
while abs((V1 - E + R * Isb * (math.exp(V2 / Vt) - 1))/(0 - E + R * Isb * (math.exp(0 / Vt) - 1))) >= 10**(-5):
    counter = counter + 1
    (f1, f2, V) = newtonRap(V1, V2)
    V1 = V[0][0]
    V2 = V[1][0]
    print("number of iteration: " + str(counter))
    print("f1 = " + str(abs(f1)) + "; " + "f2 = " + str(abs(f2)))
    print("V1 = " + str(V1) + "V" + "; " + "V2 = " + str(V2) + "V")
