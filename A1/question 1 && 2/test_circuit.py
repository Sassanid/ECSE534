from methods import *

circuitNumber = int(input("Please type in the circuit number: "))
circuit = getCircuit(circuitNumber)
voltage = solveCircuitProblem(circuit[0], circuit[1], circuit[2], circuit[3])
for n in range(len(voltage)):
    # just keep the decimal to one digit
    print('the voltage at node ' + str(n + 1) + ' is ' + str(round(voltage[n], 2)))

