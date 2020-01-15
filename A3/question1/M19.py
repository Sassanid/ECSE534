import math
A = 1/100**2
Lc = 30/100
La = 0.5/100
u0 = 4*math.pi/10**7
N = 800
I = 10
mmf = N*I
k = La/(u0*A)
def getbandh(phi, B, H):
    b = phi/A
    h = 0
    h_derivative = 0
    for i in range(len(B) - 1):
        if B[i] <= b < B[i + 1]:
            h_derivative = (H[i + 1] - H[i])/((B[i + 1] - B[i]) * A)
            h = H[i] + (b - B[i]) * ((H[i + 1] - H[i])/(B[i + 1] - B[i]))
            break
        elif b > B[len(B) - 1]:
            h_derivative = (H[len(B) - 1] - H[len(B) - 2]) / ((B[len(B) - 1] - B[len(B) - 2]) * A)
            h = H[len(B) - 1] + (b - B[len(B) - 1]) * ((H[len(B) - 1] - H[len(B) - 2]) / (B[len(B) - 1] - B[len(B) - 2]))
            break
    return h, h_derivative
def getInverseH(phi, B, H):
    h = (mmf - k*phi)/Lc
    b = 0
    for i in range(len(H) - 1):
        if H[i] <= h < H[i + 1]:
            b = B[i] + (h - H[i]) * ((B[i + 1] - B[i])/(H[i + 1] - H[i]))
            break
        elif h > H[len(B) - 1]:
            b = B[len(B) - 1] + (h - H[len(B) - 1]) * ((B[len(B) - 1] - B[len(B) - 2])/(H[len(B) - 1] - H[len(B) - 2]))
            break
    return b
def fphi(phi, h):
    return phi + (Lc * h - mmf)/k
def fphiderivative(h_derivative):
    return 1 + Lc * h_derivative/k
def NewtonRaphson(B, H):
    phi = 0
    (h0, h_derivative0) = getbandh(0, B, H)
    counter = 0
    h = h0
    h_derivative = h_derivative0
    while abs(fphi(phi, h)/fphi(0, h0)) >= 10**(-6):
        phi = -(fphi(phi, h)/fphiderivative(h_derivative)) + phi
        (h, h_derivative) = getbandh(phi, B, H)
        counter = counter + 1
    return phi, counter

#this method does not converge
def successivesub_initial(B, H):
    phi = 0
    (h0, h_derivative0) = getbandh(0, B, H)
    counter = 0
    h = h0
    while abs(fphi(phi, h)/fphi(0, h0)) >= 10**(-6):
        phi = -fphi(phi, h) + phi
        (h, h_derivative) = getbandh(phi, B, H)
        counter = counter + 1
    return phi, counter

def successivesub(B, H):
    phi = 0
    (h0, h_derivative0) = getbandh(0, B, H)
    counter = 0
    h = h0
    while abs(fphi(phi, h)/fphi(0, h0)) >= 10**(-6):
        phi = A*getInverseH(phi, B, H)
        (h, h_derivative) = getbandh(phi, B, H)
        counter = counter + 1
    return phi, counter

