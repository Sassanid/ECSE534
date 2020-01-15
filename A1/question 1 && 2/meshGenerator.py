class mesh:
    def __init__(self, meshSize, testVoltage):
        self.meshSize = meshSize
        self.numNodes = (meshSize + 1)**2
        self.numBranches = 2*meshSize*(meshSize + 1)
        self.testVoltage = testVoltage

# for voltage testing purpose, we add a last branch
# that connects the bottom left corner and the top right corner with a test voltage source
# in series with a 10kOhm resistor
# the final resistance we get should be in order of KOhm
    def yMatrix(self):
        y = [[0 for a in range(self.numBranches + 1)] for b in range(self.numBranches + 1)]
        for i in range(len(y)):
            y[i][i] = 1/10
        return y

    def JMatrix(self):
        J = [0 for a in range(self.numBranches + 1)]
        return J

    def EMatrix(self):
        E = [0 for a in range(self.numBranches)]
        # suppose we add a 10V test voltage at the last branch
        E.append(self.testVoltage)
        return E

    def AMatrix(self):
        # we choose the last node, which is the most bottom right node as the reference
        # otherwise rows of A will sum to 0
        # in this way the A we build will not contain this last row
        A = [[0 for a in range(self.numBranches + 1)] for b in range(self.numNodes)]
        global j
        j = 0  # column index

        for i in range(len(A)):
            # first, construct the circuits as if there is no voltage source connected
            # if the node is not at the most top
            if i > self.meshSize:
                A[i][j - (self.meshSize + 1)] = -1

            # if the node is not at the most bottom
            if i < self.numNodes - 1 - self.meshSize:
                A[i][j + self.meshSize] = 1

            # if the node is not at the most left
            if not i % (self.meshSize + 1) == 0:
                A[i][j - 1] = -1

            # if the node is not at the most right
            if not (i - self.meshSize) % (self.meshSize + 1) == 0:
                A[i][j] = 1
                j = j + 1
            else:
                j = j + (self.meshSize + 1)

        # finally, consider the special nodes where we connect the voltage source
        A[self.meshSize][self.numBranches] = 1
        A[self.numNodes - 1 - self.meshSize][self.numBranches] = -1

        # we can choose one of the node as our reference, here I choose the top-right node
        del A[self.meshSize]
        return A
