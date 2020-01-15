import math


def gaussLegendre(start, end, function, N, h = None):
    if h == None:
        h = [(end - start)/N for i in range(N)]
    integral_sum = 0
    xstart = start
    xend = start
    for i in range(0, N):
        xend = xend + h[i]
        xi = (xstart + xend)/2
        integral_sum += h[i] * function(xi)
        xstart = xstart + h[i]
    return integral_sum


start = 0
end = 1
# integral of sin(x)
actual1 = 0.45969769413
print("integral of sin(x), where the actual value is " + str(actual1) + ":\n")
function1 = math.sin
for N in range(1, 21):
    gaussLegendreResult = gaussLegendre(start, end, function1, N)
    print("N = " + str(N) + "  the Gauss-Legendre solution: " + str(gaussLegendreResult) + "  the error is: " + str(abs(gaussLegendreResult - actual1)))
print("\n")
# integral of ln(x)
actual2 = -1
print("integral of ln(x), where the actual value is " + str(actual2) + ":\n")
function2 = math.log
for N in range(10, 210, 10):
    gaussLegendreResult = gaussLegendre(start, end, function2, N)
    print("N = " + str(N) + "  the Gauss-Legendre solution: " + str(gaussLegendreResult) + "  the error is: " + str(abs(gaussLegendreResult - actual2)))
print("\n")
# integral of ln(0.2|sin(x)|)
actual3 = -2.66616
print("integral of ln(0.2|sin(x)|, where the actual value is " + str(actual3) + ":\n")


def function3(x):
    result = math.log(0.2*abs(math.sin(x)))
    return result


for N in range(10, 210, 10):
    gaussLegendreResult = gaussLegendre(start, end, function3, N)
    print("N = " + str(N) + "  the Gauss-Legendre solution: " + str(gaussLegendreResult) + "  the error is: " + str(abs(gaussLegendreResult - actual3)))
print("\n")
# integral of of ln(x) and ln(0.2|sin(x)|) with inequal interval
Ntest = 10
h = [0 for i in range(Ntest)]
smallest_interval = 2*(end - start)/(Ntest*(Ntest + 1))
for i in range(Ntest):
    h[i] = smallest_interval * (i + 1)
print("integral of ln(x), where the actual value is " + str(actual2) + ":\n")
gaussLegendreResult = gaussLegendre(start, end, function2, Ntest, h)
print("N = " + str(Ntest) + "  the Gauss-Legendre solution: " + str(gaussLegendreResult) + "  the error is: " + str(abs(gaussLegendreResult - actual2)) + "\n")
print("integral of ln(0.2|sin(x)|, where the actual value is " + str(actual3) + ":\n")
gaussLegendreResult = gaussLegendre(start, end, function3, Ntest, h)
print("N = " + str(Ntest) + "  the Gauss-Legendre solution: " + str(gaussLegendreResult) + "  the error is: " + str(abs(gaussLegendreResult - actual3)))
