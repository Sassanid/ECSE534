import matplotlib.pyplot as plt

def Lagrange(variable, xPoints, yPoints):
    l = len(xPoints)
    L = [None for i in range(l)]
    yf = 0
    for i in range(l):
        F1 = 1
        F2 = 1
        for j in range(l):
            if j != i:
                F1 *= (variable - xPoints[j])
                F2 *= (xPoints[i] - xPoints[j])
        L[i] = F1/F2
        yf += yPoints[i] * L[i]
    yf = yf.expand()
    return yf

def cubicHermite(variable, xPoints, yPoints):
    l = len(xPoints)
    yf = [0 for i in range(l - 1)]
    for i in range(l - 1):
        ydiff1 = (yPoints[i + 1] - yPoints[i])/(xPoints[i + 1] - xPoints[i])
        if i == l - 2:
            ydiff2 = yPoints[i + 1] / xPoints[i]
        else:
            ydiff2 = (yPoints[i + 2] - yPoints[i + 1]) / (xPoints[i + 2] - xPoints[i + 1])
        U1 = (1 - 2 * (variable - xPoints[i])/(xPoints[i] - xPoints[i + 1]))*((variable - xPoints[i + 1])/(xPoints[i] - xPoints[i + 1]))**2
        U2 = (1 - 2 * (variable - xPoints[i + 1])/(xPoints[i + 1] - xPoints[i]))*((variable - xPoints[i])/(xPoints[i + 1] - xPoints[i]))**2
        V1 = (variable - xPoints[i])*((variable - xPoints[i + 1])/(xPoints[i] - xPoints[i + 1]))**2
        V2 = (variable - xPoints[i + 1])*((variable - xPoints[i])/(xPoints[i + 1] - xPoints[i]))**2
        yf[i] = (yPoints[i]*U1 + yPoints[i + 1]*U2 + ydiff1*V1 + ydiff2*V2).expand()
        print(yf[i])
    return yf


def plot(b, B, h, H):
    plt.plot(H, B, "r.", label="Table1 Data")
    plt.plot(h, b, "k", label="Interpolated Curve")
    plt.xlabel("H(A/m)")
    plt.ylabel("B(T)")
    plt.legend()


