# FiniteDifferenceSixthOrderProgram.py

from sympy import *

def taylorSeriesExpansion(f, x, h):
    factorialList = [0 for k in range(7)]
    factorialList[0] = 1
    for k in range(1, 7):
        factorialList[k] = k*factorialList[k-1]
    listOfDerivs = [0 for k in range(7)]
    listOfDerivs[0] = f
    for k in range(1, 7): listOfDerivs[k] = diff(listOfDerivs[k-1], x)
    powerOfH = [0 for k in range(7)]
    powerOfH[0] = 1
    for k in range(1, 7): powerOfH[k] = powerOfH[k-1]*h
    seriesExpand = 0
    for k in range(7):
        seriesExpand = expand(
            seriesExpand \
            + Rational(1, factorialList[k])*listOfDerivs[k]*powerOfH[k]
        )
    return seriesExpand

x = Symbol("x")
h = Symbol("h")
f = Function("f")(x)
nOrder = 6
nOrderBy2 = nOrder//2
latticeSites = [
    taylorSeriesExpansion(f, x, (k - nOrder//2)*h) \
        for k in range(nOrder+1)
    ]
aCoeffList = [Symbol("a["+str(k)+"]") for k in range(nOrder+1)]
bCoeffList = [Symbol("b["+str(k)+"]") for k in range(nOrder+1)]
cCoeffList = [Symbol("c["+str(k)+"]") for k in range(nOrder+1)]
firstOrderDeriv = 0
secondOrderDeriv = 0
thirdOrderDeriv = 0
for k in range(nOrder+1):
    firstOrderDeriv = expand(
        firstOrderDeriv + aCoeffList[k]*latticeSites[k]
        )
    secondOrderDeriv = expand(
        secondOrderDeriv + bCoeffList[k]*latticeSites[k]
        )
    thirdOrderDeriv = expand(
        thirdOrderDeriv + cCoeffList[k]*latticeSites[k]
    )
listOfEquationsForFirst = [
    firstOrderDeriv.coeff(h, k) for k in range(nOrder+1)
    ]
listOfEquationsForSecond = [
    secondOrderDeriv.coeff(h, k) for k in range(nOrder+1)
    ]
listOfEquationsForThird = [
    thirdOrderDeriv.coeff(h, k) for k in range(nOrder+1)
    ]
listOfEquationsForFirst[1] = listOfEquationsForFirst[1] - diff(f, x)
listOfEquationsForSecond[2] = listOfEquationsForSecond[2] - diff(f, x, x)
listOfEquationsForThird[3] = listOfEquationsForThird[3] - diff(f, x, x, x)
assert(
    solve(listOfEquationsForFirst, aCoeffList) == {
        aCoeffList[0]: -Rational(1, 60),
        aCoeffList[1]: Rational(3, 20),
        aCoeffList[2]: -Rational(3, 4),
        aCoeffList[3]: 0,
        aCoeffList[4]: Rational(3, 4),
        aCoeffList[5]: -Rational(3, 20),
        aCoeffList[6]: Rational(1, 60)
        }
)
assert(
    solve(listOfEquationsForSecond, bCoeffList) == {
        bCoeffList[0]: Rational(1, 90),
        bCoeffList[1]: -Rational(3, 20),
        bCoeffList[2]: Rational(3, 2),
        bCoeffList[3]: -Rational(49, 18),
        bCoeffList[4]: Rational(3, 2),
        bCoeffList[5]: -Rational(3, 20),
        bCoeffList[6]: Rational(1, 90)
        }
)
assert(
    solve(listOfEquationsForThird, cCoeffList) == {
        cCoeffList[0]: Rational(1, 8),
        cCoeffList[1]: -1,
        cCoeffList[2]: Rational(13, 8),
        cCoeffList[3]: 0,
        cCoeffList[4]: -Rational(13, 8),
        cCoeffList[5]: 1,
        cCoeffList[6]: -Rational(1, 8)
    }
)

