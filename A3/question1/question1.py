from sympy import *
import numpy
from interpolate import *
from M19 import *
import matplotlib.pyplot as plt


B = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9]
H = [0.0, 14.7, 36.5, 71.7, 121.4, 197.4, 256.2, 348.7, 540.6, 1062.8, 2318.0, 4781.9, 8687.4, 13924.3, 22650.2]
Bb = [0.0, 1.3, 1.4, 1.7, 1.8, 1.9]
Hb = [0.0, 540.6, 1062.8, 8687.4, 13924.3, 22650.2]

b_a = numpy.arange(0, 1.2, 0.01)
h_a = [0 for i in range(len(b_a))]

b_b = numpy.arange(0, 2.0, 0.02)
h_b = [0 for i in range(len(b_b))]

b_c = numpy.arange(0, 1.9, 0.01)
h_c = [0 for i in range(len(b_c))]

variable = Symbol('B')
problem_a = Lagrange(variable, B[:6], H[:6])
problem_b = Lagrange(variable, Bb, Hb)
problem_c = cubicHermite(variable, Bb, Hb)

print("Interpolate the first 6 points with full-domain Lagrange polynomials result is: \n" + "H = " + str(problem_a) + "\n")
print("Interpolate the selecting points with full-domain Lagrange polynomials result is: \n" + "H = " + str(problem_b) + "\n")

#plotting block
for i in range(len(b_a)):
    h_a[i] = problem_a.subs({variable: b_a[i]})
pa = plt.figure(1)
plot(b_a, B[:6], h_a, H[:6])

for i in range(len(b_b)):
    h_b[i] = problem_b.subs({variable: b_b[i]})
pb = plt.figure(2)
plot(b_b, Bb, h_b, Hb)

for i in range(len(b_c)):
    for j in range(5):
        if Bb[j] <= b_c[i] < Bb[j + 1]:
            h_c[i] = problem_c[j].subs({variable: b_c[i]})
            break
pb = plt.figure(3)
plot(b_c, Bb, h_c, Hb)

plt.show()

(phi, counter) = NewtonRaphson(B, H)
print("Total iteration number: " + str(counter))
print("The flux is: " + str(phi) + " Wb")

(phi2, counter2) = successivesub(B, H)
print("Total iteration number: " + str(counter2))
print("The flux is: " + str(phi2) + " Wb")

