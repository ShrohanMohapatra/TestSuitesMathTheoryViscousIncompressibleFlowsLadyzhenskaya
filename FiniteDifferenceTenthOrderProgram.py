# FiniteDifferenceTenthOrderProgram.py

from sympy import *

def taylorSeriesExpansion(f, x, h):
    factorialList = [0 for k in range(11)]
    factorialList[0] = 1
    for k in range(1, 11):
        factorialList[k] = k*factorialList[k-1]
    listOfDerivs = [0 for k in range(11)]
    listOfDerivs[0] = f
    for k in range(1, 11): listOfDerivs[k] = diff(listOfDerivs[k-1], x)
    powerOfH = [0 for k in range(11)]
    powerOfH[0] = 1
    for k in range(1, 11): powerOfH[k] = powerOfH[k-1]*h
    seriesExpand = 0
    for k in range(11):
        seriesExpand = expand(
            seriesExpand \
            + Rational(1, factorialList[k])*listOfDerivs[k]*powerOfH[k]
        )
    return seriesExpand

x = Symbol("x")
h = Symbol("h")
f = Function("f")(x)
nOrder = 10
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
        aCoeffList[0]: -Rational(1, 1260),
        aCoeffList[1]: Rational(5, 504),
        aCoeffList[2]: -Rational(5, 84),
        aCoeffList[3]: Rational(5, 21),
        aCoeffList[4]: -Rational(5, 6),
        aCoeffList[5]: 0,
        aCoeffList[6]: Rational(5, 6),
        aCoeffList[7]: -Rational(5, 21),
        aCoeffList[8]: Rational(5, 84),
        aCoeffList[9]: -Rational(5, 504),
        aCoeffList[10]: Rational(1, 1260)
        }
)

assert(
    solve(listOfEquationsForSecond, bCoeffList) == {
        bCoeffList[0]: Rational(1, 3150),
        bCoeffList[1]: -Rational(5, 1008),
        bCoeffList[2]: Rational(5, 126),
        bCoeffList[3]: -Rational(5, 21),
        bCoeffList[4]: Rational(5, 3),
        bCoeffList[5]: -Rational(5269, 1800),
        bCoeffList[6]: Rational(5, 3),
        bCoeffList[7]: -Rational(5, 21),
        bCoeffList[8]: Rational(5, 126),
        bCoeffList[9]: -Rational(5, 1008),
        bCoeffList[10]: Rational(1, 3150)
        }
)

assert(
    solve(listOfEquationsForThird, cCoeffList) == {
        cCoeffList[0]: Rational(41, 6048),
        cCoeffList[1]: -Rational(1261, 15120),
        cCoeffList[2]: Rational(541, 1120),
        cCoeffList[3]: -Rational(4369, 2520),
        cCoeffList[4]: Rational(1669, 720),
        cCoeffList[5]: 0,
        cCoeffList[6]: -Rational(1669, 720),
        cCoeffList[7]: Rational(4369, 2520),
        cCoeffList[8]: -Rational(541, 1120),
        cCoeffList[9]: Rational(1261, 15120),
        cCoeffList[10]: -Rational(41, 6048)
    }
)

