import time
from methods import *
from meshGenerator import *
meshSize = int(input("Please enter the mesh size: "))
ifBandwidth = int(input("1 for half-bandwidth computation, 0 for normal: "))
# set a 2 volts test voltage
testVoltage = 2
testCircuit = mesh(meshSize, testVoltage)
A = testCircuit.AMatrix()
J = testCircuit.JMatrix()
y = testCircuit.yMatrix()
E = testCircuit.EMatrix()

B = multiplyMatrix(A, multiplyMatrix(y,transposeMatrix(A)))
halfBandwidth = None
if ifBandwidth == 1:
    halfBandwidth = meshSize + 2

startTime = time.time()  # starts to account the time for solving the circuit
voltages = solveCircuitProblem(A, J, y, E, halfBandwidth)
endTime = time.time()  # Stops the clock

voltagesAcrossTheCircuit = voltages[testCircuit.numNodes - 2 - testCircuit.meshSize]

# use the voltage divider to find the total resistance of the mesh
# from bottom left to top right
resistance = voltagesAcrossTheCircuit/(testVoltage - voltagesAcrossTheCircuit) * 10
print("resistance in kOhm is:" + str(round(resistance, 2)))
print("time taken in s:" + str(endTime - startTime))
