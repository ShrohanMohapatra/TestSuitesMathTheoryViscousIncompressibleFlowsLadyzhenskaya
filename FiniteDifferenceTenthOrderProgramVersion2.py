# FiniteDifferenceTenthOrderProgramVersion2.py

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

h = Symbol("h")
x = Symbol("x")
f = Function("f")(x)

fXMinus5H = expand(series(f.subs(x, x-5*h).doit(), h, 0, 12)) # taylorSeriesExpansion(f, x, -5*h)
fXMinus4H = expand(series(f.subs(x, x-4*h).doit(), h, 0, 12)) # taylorSeriesExpansion(f, x, -4*h)
fXMinus3H = expand(series(f.subs(x, x-3*h).doit(), h, 0, 12)) # taylorSeriesExpansion(f, x, -3*h)
fXMinus2H = expand(series(f.subs(x, x-2*h).doit(), h, 0, 12)) # taylorSeriesExpansion(f, x, -2*h)
fXMinusH  = expand(series(f.subs(x, x-h).doit(), h, 0, 12)) # taylorSeriesExpansion(f, x, -h)

fXPlusH = expand(series(f.subs(x, x+h).doit(), h, 0, 12)) # taylorSeriesExpansion(f, x, h)
fXPlus2H = expand(series(f.subs(x, x+2*h).doit(), h, 0, 12)) # taylorSeriesExpansion(f, x, 2*h)
fXPlus3H = expand(series(f.subs(x, x+3*h).doit(), h, 0, 12)) # taylorSeriesExpansion(f, x, 3*h)
fXPlus4H = expand(series(f.subs(x, x+4*h).doit(), h, 0, 12)) # taylorSeriesExpansion(f, x, 4*h)
fXPlus5H = expand(series(f.subs(x, x+5*h).doit(), h, 0, 12))  # taylorSeriesExpansion(f, x, 5*h)

firstOrderDerivCheck = expand(
    Rational(1, 1260)/h*(fXPlus5H - fXMinus5H) \
    - Rational(5, 504)/h*(fXPlus4H - fXMinus4H) \
    + Rational(5, 84)/h*(fXPlus3H - fXMinus3H) \
    - Rational(5, 21)/h*(fXPlus2H - fXMinus2H) \
    + Rational(5, 6)/h*(fXPlusH - fXMinusH)
    )

thirdOrderDerivCheck = expand(
    -Rational(41, 6048)/h**3*(fXPlus5H - fXMinus5H) \
    +Rational(1261, 15120)/h**3*(fXPlus4H - fXMinus4H) \
    -Rational(541, 1120)/h**3*(fXPlus3H - fXMinus3H) \
    +Rational(4369, 2520)/h**3*(fXPlus2H - fXMinus2H) \
    -Rational(1669, 720)/h**3*(fXPlusH - fXMinusH)
    )

print("-"*50)
pprint(firstOrderDerivCheck)
print("-"*50)
pprint(thirdOrderDerivCheck)
print("-"*50)
