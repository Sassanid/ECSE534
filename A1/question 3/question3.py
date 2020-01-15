from potentialSolver import *
method = input("S for 'SOR' and j for 'Jacobi': ")
if method == 'S' or method == 's':
    nonUniSpacing = int(input("0 for uniform spacing, 1 for non uniform: "))
    if nonUniSpacing == 0:
        h = float(input("non-spacing constant h: "))
        for w in range(10, 20, 1):
            potentials = potentialMesh(h, 0.00001)
            # divide w by 10 as the SOR parameter

            final_mesh = potentials.potentials_SOR(w/10)
            print('value of w = '+ str(w/10))
            print('value of potential at the point (x,y)= (0.06, 0.04) is: ' + str(final_mesh[int(0.06/h)][int(0.04/h)]))
    else:
        # set h as 0.01 to allow the potentialSolver calculate the number of rows and columns
        width_index = [0, 0.02, 0.04, 0.055, 0.059, 0.06, 0.065, 0.07, 0.08, 0.09, 0.1]
        height_index = [0, 0.01, 0.02, 0.03, 0.035, 0.04, 0.045, 0.055, 0.06, 0.08, 0.1]
        # width_index = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
        # height_index = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
        potentials = potentialMesh(0.01, 0.00001, width_index, height_index)
        final_mesh = potentials.potentials_SOR(1.2)
        print('value of potential at the point (x,y)= (0.06, 0.04) is: '
              + str(final_mesh[potentials.x_interest][potentials.y_interest]))
elif method == 'J' or method == 'j':
    h = 0.02
    for i in range(6):
        print('h equals ' + str(h) + 'm')
        potentials = potentialMesh(h, 0.00001)
        # divide w by 10 as the SOR parameter

        final_mesh = potentials.potentials_jacobi()
        print('value of potential at the point (x,y)= (0.06, 0.04) is: ' + str(final_mesh[int(0.06 / h)][int(0.04 / h)]))
        h = h / 2


