class potentialMesh:
    def __init__(self, h, residual_limit, width_index = None, height_index = None):
        # define the size of the symmetry plane
        self.h = h
        self.residual_limit = residual_limit
        self.outerLength = 0.1
        self.innerHeight = 0.02
        self.innerWidth = 0.04
        self.outerPotential = 0
        self.innerPotential = 110
        self.numColumns = int(self.outerLength/h + 1)
        self.numRows = int(self.outerLength/h + 1)
        self.mesh = None
        self.x_interest = None
        self.y_interest = None
        self.x_inner = None
        self.y_inner = None
        self.width_index = width_index
        self.height_index = height_index
        if width_index and height_index:
            self.x_interest = width_index.index(0.06)
            self.y_interest = height_index.index(0.04)
            for i in range(len(width_index)):
                if width_index[i] > self.innerWidth:
                    self.x_inner = i
                    break
            for j in range(len(width_index)):
                if height_index[j] > self.innerHeight:
                    self.y_inner = j
                    break
            self.mesh = [[self.innerPotential if x < self.x_inner and y < self.y_inner
                else self.outerPotential if x == self.numColumns - 1 and y == self.numRows - 1 else 0.0 for x
                in range(self.numColumns)] for y in range(self.numRows)]
        else:
            self.mesh = [[self.innerPotential if x <= self.innerWidth / self.h and y <= self.innerHeight / self.h
                else self.outerPotential if x == self.numColumns - 1 and y == self.numRows - 1 else 0.0 for x
                in range(self.numColumns)] for y in range(self.numRows)]

    def SOR(self, w):
        # we should keep in mind that the most 'outer' node has V = 0
        for y in range(self.numRows - 1):
            for x in range(int(self.numColumns - 1)):
                if x == 0 and y > int(self.innerHeight/ self.h):
                    # we can assume this formula since when x = 0 d(potential)/dx = 0, we have mesh[x - 1][y] = mesh[x+1}[y]
                    self.mesh[y][x] = (1 - w) * self.mesh[y][x] + (w / 4) * (2*self.mesh[y][x + 1] + self.mesh[y - 1][x] + self.mesh[y + 1][x])
                elif y == 0 and x > int(self.innerWidth/ self.h):
                    # same argument as above
                    self.mesh[y][x] = (1 - w) * self.mesh[y][x] + (w / 4) * (
                    self.mesh[y][x - 1] + self.mesh[y][x + 1] + 2*self.mesh[y + 1][x])
                elif x > int(self.innerWidth/ self.h) or y > int(self.innerHeight/ self.h):
                    self.mesh[y][x] = (1 - w) * self.mesh[y][x] + (w / 4) * (
                    self.mesh[y][x - 1] + self.mesh[y][x + 1] + self.mesh[y - 1][x] + self.mesh[y + 1][x])
        return self.mesh

        # Equation that computes the residue

    def residual(self):
        res = 0
        finalRes = 0
        if self.width_index and self.height_index:
            for y in range(1, self.numRows - 1):
                for x in range(1, self.numColumns - 1):
                    if x == 0 and y >= self.y_inner:
                        a2 = self.width_index[x + 1] - self.width_index[x]
                        a1 = a2
                        b1 = self.height_index[y] - self.height_index[y - 1]
                        b2 = self.height_index[y + 1] - self.height_index[y]
                        # we can assume this formula since when x = 0 d(potential)/dx = 0, we have mesh[x - 1][y] = mesh[x+1}[y]
                        res = (1 / (a1 * a2) + 1 / (b1 * b2)) * self.mesh[y][x] - (
                                2 * self.mesh[y][x + 1] / (a2 * (a1 + a2))
                                + self.mesh[y - 1][x] / (b1 * (b1 + b2)) + self.mesh[y + 1][x] / (b2 * (b1 + b2)))
                        print(res)
                    elif y == 0 and x >= self.x_inner:
                        # same argument as above
                        a1 = self.width_index[x] - self.width_index[x - 1]
                        a2 = self.width_index[x + 1] - self.height_index[x]
                        b2 = self.height_index[y + 1] - self.height_index[y]
                        b1 = b2
                        res = (1 / (a1 * a2) + 1 / (b1 * b2)) * self.mesh[y][x] - (
                                self.mesh[y][x - 1] / (a1 * (a1 + a2)) + self.mesh[y][x + 1] / (a2 * (a1 + a2))
                                + 2 * self.mesh[y + 1][x] / (b2 * (b1 + b2)))
                        print(res)
                    elif x >= self.x_inner or y >= self.y_inner:
                        a1 = self.width_index[x] - self.width_index[x - 1]
                        a2 = self.width_index[x + 1] - self.width_index[x]
                        b1 = self.height_index[y] - self.height_index[y - 1]
                        b2 = self.height_index[y + 1] - self.height_index[y]
                        res = (1 / (a1 * a2) + 1 / (b1 * b2)) * self.mesh[y][x] - (
                                self.mesh[y][x - 1] / (a1 * (a1 + a2)) + self.mesh[y][x + 1] / (a2 * (a1 + a2))
                                + self.mesh[y - 1][x] / (b1 * (b1 + b2)) + self.mesh[y + 1][x] / (b2 * (b1 + b2)))
                    res = abs(res)
                    if res > finalRes:
                        # Updates variable with the biggest residue amongst the free point
                        finalRes = res
        else:
            for y in range(1, self.numRows - 1):
                for x in range(1, self.numColumns - 1):
                    if x == 0 and y > int(self.innerHeight/ self.h):
                        # we can assume this formula since when x = 0 d(potential)/dx = 0, we have mesh[x - 1][y] = mesh[x+1}[y]
                        res = 2 * self.mesh[y][x + 1] + self.mesh[y - 1][x] + self.mesh[y + 1][x] - 4 * self.mesh[y][x]
                    elif y == 0 and x > int(self.innerWidth/ self.h):
                        res = self.mesh[y][x - 1] + self.mesh[y][x + 1] + 2*self.mesh[y + 1][x] - 4 * self.mesh[y][x]
                    elif x > int(self.innerWidth/ self.h) or y > int(self.innerHeight/ self.h):
                        res = self.mesh[y][x - 1] + self.mesh[y][x + 1] + self.mesh[y - 1][x] + self.mesh[y + 1][x] - 4 * self.mesh[y][x]
                    res = abs(res)
                    if res > finalRes:
                        # Updates variable with the biggest residue amongst the free point
                        finalRes = res
        return finalRes

    # The Equation that calculates Jacobian
    def jacobi(self):
        # we should keep in mind that the most 'outer' node has V = 0
        for y in range(self.numRows - 1):
            for x in range(int(self.numColumns - 1)):
                if x == 0 and y > int(self.innerHeight / self.h):
                    # we can assume this formula since when x = 0 d(potential)/dx = 0, we have mesh[x - 1][y] = mesh[x+1}[y]
                    self.mesh[y][x] = 1/4* (2 * self.mesh[y][x + 1] + self.mesh[y - 1][x] + self.mesh[y + 1][x])
                elif y == 0 and x > int(self.innerWidth / self.h):
                    # same argument as above
                    self.mesh[y][x] = 1/4 * (self.mesh[y][x - 1] + self.mesh[y][x + 1] + 2 * self.mesh[y + 1][x])
                elif x > int(self.innerWidth / self.h) or y > int(self.innerHeight / self.h):
                    self.mesh[y][x] = 1/4 * (self.mesh[y][x - 1] + self.mesh[y][x + 1] + self.mesh[y - 1][x] + self.mesh[y + 1][x])
        return self.mesh

    def potentials_SOR(self, w):
        iteration = 0
        if self.width_index and self.height_index:
            self.SOR_non_uniform(w)
            while self.residual() >= self.residual_limit:
                self.SOR_non_uniform(w)
                iteration = iteration + 1
            print('total iteration is: ' + str(iteration))
        else:
            self.SOR(w)
            while self.residual() >= self.residual_limit:
                self.SOR(w)
                iteration = iteration + 1
            print('total iteration is: ' + str(iteration))
        return self.mesh

    def potentials_jacobi(self):
        iteration = 0
        self.jacobi()
        while self.residual() >= self.residual_limit:
            self.jacobi()
            iteration = iteration + 1
        print('total iteration is: ' + str(iteration))
        return self.mesh


    def SOR_non_uniform(self, w):
        # we should keep in mind that the most 'outer' node has V = 0
        # since we should have equal
        for y in range(self.numRows - 1):
            for x in range(int(self.numColumns - 1)):
                if x == 0 and y >= self.y_inner:
                    a2 = self.width_index[x + 1] - self.width_index[x]
                    a1 = a2
                    b1 = self.height_index[y] - self.height_index[y - 1]
                    b2 = self.height_index[y + 1] - self.height_index[y]
                    # we can assume this formula since when x = 0 d(potential)/dx = 0, we have mesh[x - 1][y] = mesh[x+1}[y]
                    self.mesh[y][x] = (1 - w) * self.mesh[y][x] + w * (2 * self.mesh[y][x + 1] / (a2 * (a1 + a2))
                        + self.mesh[y - 1][x] / (b1 * (b1 + b2)) + self.mesh[y + 1][x] / (b2 * (b1 + b2)))/(1 / (a1 * a2) + 1 / (b1 * b2))
                elif y == 0 and x >= self.x_inner:
                    # same argument as above
                    a1 = self.width_index[x] - self.width_index[x - 1]
                    a2 = self.width_index[x + 1] - self.width_index[x]
                    b2 = self.height_index[y + 1] - self.height_index[y]
                    b1 = b2
                    self.mesh[y][x] = (1 - w) * self.mesh[y][x] + w * (self.mesh[y][x - 1] / (a1 * (a1 + a2)) + self.mesh[y][x + 1] / (a2 * (a1 + a2))
                        + 2 * self.mesh[y + 1][x] / (b2 * (b1 + b2)))/(1 / (a1 * a2) + 1 / (b1 * b2))
                elif x >= self.x_inner or y >= self.y_inner:
                    a1 = self.width_index[x] - self.width_index[x - 1]
                    a2 = self.width_index[x + 1] - self.width_index[x]
                    b1 = self.height_index[y] - self.height_index[y - 1]
                    b2 = self.height_index[y + 1] - self.height_index[y]
                    self.mesh[y][x] = (1 - w) * self.mesh[y][x] + w * (self.mesh[y][x - 1]/(a1*(a1 + a2)) + self.mesh[y][x + 1]/(a2*(a1 + a2))
                        + self.mesh[y - 1][x]/(b1*(b1 + b2)) + self.mesh[y + 1][x]/(b2*(b1 + b2)))/(1 / (a1 * a2) + 1 / (b1 * b2))
        return self.mesh


