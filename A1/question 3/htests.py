from potentialSolver import *
h = 0.02
potentials = potentialMesh(h, 0.00001)

for i in range(6):
    print('h equals ' + str(h) + 'm')
    potentials = potentialMesh(h, 0.00001)
    # divide w by 10 as the SOR parameter

    final_mesh = potentials.potentials_SOR(1.3)
    print('value of w = ' + str(1.3))
    print('value of potential at the point (x,y)= (0.06, 0.04) is: ' + str(final_mesh[int(0.06 / h)][int(0.04 / h)]))
    h = h/2
