# SixPointFiniteDifferenceCheck.jl

using LinearAlgebra

function firstOrderDerivSixthOrder(xAxis, numPixel)
    fAxis = zeros(Float64, 1, numPixel)
    for k = 1:1:numPixel
        fAxis[k] = xAxis[k]^4 + 6*xAxis[k]^3 + 16*xAxis[k]^2 + 8*xAxis[k] - 144
    end
    fDerivAxisFiniteDiff = zeros(Float64, 1, numPixel)
    fDerivCheck = zeros(Float64, 1, numPixel)
    h = xAxis[2] - xAxis[1]
    fDerivAxisFiniteDiff[1] = 4*xAxis[1]^3 + 18*xAxis[1]^2 + 32*xAxis[1] + 8
    fDerivAxisFiniteDiff[2] = 4*xAxis[2]^3 + 18*xAxis[2]^2 + 32*xAxis[2] + 8
    fDerivAxisFiniteDiff[3] = 4*xAxis[3]^3 + 18*xAxis[3]^2 + 32*xAxis[3] + 8
    for k = 4:1:numPixel-3
        fDerivAxisFiniteDiff[k] = -1/60/h*fAxis[k-3] + 3/20/h*fAxis[k-2] - 3/4/h*fAxis[k-1] + 3/4/h*fAxis[k+1] - 3/20/h*fAxis[k+2] + 1/60/h*fAxis[k+3] 
    end
    fDerivAxisFiniteDiff[numPixel-2] = 4*xAxis[numPixel-2]^3 + 18*xAxis[numPixel-2]^2 + 32*xAxis[numPixel-2] + 8
    fDerivAxisFiniteDiff[numPixel-1] = 4*xAxis[numPixel-1]^3 + 18*xAxis[numPixel-1]^2 + 32*xAxis[numPixel-1] + 8
    fDerivAxisFiniteDiff[numPixel] = 4*xAxis[numPixel]^3 + 18*xAxis[numPixel]^2 + 32*xAxis[numPixel] + 8
    for k = 1:1:numPixel
        fDerivCheck[k] = 4*xAxis[k]^3 + 18*xAxis[k]^2 + 32*xAxis[k] + 8
    end
    deltaFderivCheck = zeros(Float64, 1, numPixel)
    for k = 1:1:numPixel
        deltaFderivCheck[k] = abs(fDerivCheck[k] - fDerivAxisFiniteDiff[k])
    end
    println(norm(deltaFderivCheck)/numPixel)
end

function thirdOrderDerivSixthOrder(xAxis, numPixel)
    fAxis = zeros(Float64, 1, numPixel)
    for k = 1:1:numPixel
        fAxis[k] = xAxis[k]^4 + 6*xAxis[k]^3 + 16*xAxis[k]^2 + 8*xAxis[k] - 144
    end
    fThirdDerivAxisFiniteDiff = zeros(Float64, 1, numPixel)
    fThirdDerivCheck = zeros(Float64, 1, numPixel)
    h = xAxis[2] - xAxis[1]
    fThirdDerivAxisFiniteDiff[1] = 24*xAxis[1] + 36
    fThirdDerivAxisFiniteDiff[2] = 24*xAxis[2] + 36
    fThirdDerivAxisFiniteDiff[3] = 24*xAxis[3] + 36
    for k = 4:1:numPixel-3
        fThirdDerivAxisFiniteDiff[k] = 1/8/h^3*fAxis[k-3] - 1/h^3*fAxis[k-2] + 13/8/h^3*fAxis[k-1] - 13/8/h^3*fAxis[k+1] + 1/h^3*fAxis[k+2] - 1/8/h^3*fAxis[k+3] 
    end
    fThirdDerivAxisFiniteDiff[numPixel-2] = 24*xAxis[numPixel-2] + 36
    fThirdDerivAxisFiniteDiff[numPixel-1] = 24*xAxis[numPixel-1] + 36
    fThirdDerivAxisFiniteDiff[numPixel] = 24*xAxis[numPixel] + 36
    for k = 1:1:numPixel
        fThirdDerivCheck[k] = 24*xAxis[k] + 36
    end
    deltaFThirdDerivCheck = zeros(Float64, 1, numPixel)
    for k = 1:1:numPixel
        deltaFThirdDerivCheck[k] = abs(fThirdDerivCheck[k] - fThirdDerivAxisFiniteDiff[k])
    end
    println(norm(deltaFThirdDerivCheck)/numPixel)
end

function firstOrderDerivTenthOrder(xAxis, numPixel)
    fAxis = zeros(Float64, 1, numPixel)
    for k = 1:1:numPixel
        fAxis[k] = xAxis[k]^4 + 6*xAxis[k]^3 + 16*xAxis[k]^2 + 8*xAxis[k] - 144
    end
    fDerivAxisFiniteDiff = zeros(Float64, 1, numPixel)
    fDerivCheck = zeros(Float64, 1, numPixel)
    h = xAxis[2] - xAxis[1]
    fDerivAxisFiniteDiff[1] = 4*xAxis[1]^3 + 18*xAxis[1]^2 + 32*xAxis[1] + 8
    fDerivAxisFiniteDiff[2] = 4*xAxis[2]^3 + 18*xAxis[2]^2 + 32*xAxis[2] + 8
    fDerivAxisFiniteDiff[3] = 4*xAxis[3]^3 + 18*xAxis[3]^2 + 32*xAxis[3] + 8
    fDerivAxisFiniteDiff[4] = 4*xAxis[4]^3 + 18*xAxis[4]^2 + 32*xAxis[4] + 8
    fDerivAxisFiniteDiff[5] = 4*xAxis[5]^3 + 18*xAxis[5]^2 + 32*xAxis[5] + 8
    for k = 6:1:numPixel-5
        fDerivAxisFiniteDiff[k] = 1/1260/h*(fAxis[k+5] - fAxis[k-5]) - 5/504/h*(fAxis[k+4] - fAxis[k-4]) + 5/84/h*(fAxis[k+3] - fAxis[k-3]) - 5/21/h*(fAxis[k+2] - fAxis[k-2]) + 5/6/h*(fAxis[k+1] - fAxis[k-1])
    end
    fDerivAxisFiniteDiff[numPixel-4] = 4*xAxis[numPixel-4]^3 + 18*xAxis[numPixel-4]^2 + 32*xAxis[numPixel-4] + 8
    fDerivAxisFiniteDiff[numPixel-3] = 4*xAxis[numPixel-3]^3 + 18*xAxis[numPixel-3]^2 + 32*xAxis[numPixel-3] + 8
    fDerivAxisFiniteDiff[numPixel-2] = 4*xAxis[numPixel-2]^3 + 18*xAxis[numPixel-2]^2 + 32*xAxis[numPixel-2] + 8
    fDerivAxisFiniteDiff[numPixel-1] = 4*xAxis[numPixel-1]^3 + 18*xAxis[numPixel-1]^2 + 32*xAxis[numPixel-1] + 8
    fDerivAxisFiniteDiff[numPixel] = 4*xAxis[numPixel]^3 + 18*xAxis[numPixel]^2 + 32*xAxis[numPixel] + 8
    for k = 1:1:numPixel
        fDerivCheck[k] = 4*xAxis[k]^3 + 18*xAxis[k]^2 + 32*xAxis[k] + 8
    end
    deltaFFirstDerivCheck = zeros(Float64, 1, numPixel)
    for k = 1:1:numPixel
        deltaFFirstDerivCheck[k] = abs(fDerivCheck[k] - fDerivAxisFiniteDiff[k])
    end
    println(norm(deltaFFirstDerivCheck)/numPixel)
end

function thirdOrderDerivTenthOrder(xAxis, numPixel)
    fAxis = zeros(Float64, 1, numPixel)
    for k = 1:1:numPixel
        fAxis[k] = xAxis[k]^4 + 6*xAxis[k]^3 + 16*xAxis[k]^2 + 8*xAxis[k] - 144
    end
    fThirdDerivAxisFiniteDiff = zeros(Float64, 1, numPixel)
    fThirdDerivCheck = zeros(Float64, 1, numPixel)
    h = xAxis[2] - xAxis[1]
    fThirdDerivAxisFiniteDiff[1] = 24*xAxis[1] + 36
    fThirdDerivAxisFiniteDiff[2] = 24*xAxis[2] + 36
    fThirdDerivAxisFiniteDiff[3] = 24*xAxis[3] + 36
    fThirdDerivAxisFiniteDiff[4] = 24*xAxis[4] + 36
    fThirdDerivAxisFiniteDiff[5] = 24*xAxis[5] + 36
    for k = 6:1:numPixel-5
        fThirdDerivAxisFiniteDiff[k] = -41/6048/h^3*(fAxis[k+5] - fAxis[k-5]) + 1261/15120/h^3*(fAxis[k+4] - fAxis[k-4]) - 541/1120/h^3*(fAxis[k+3] - fAxis[k-3]) + 4369/2520/h^3*(fAxis[k+2] - fAxis[k-2]) - 1669/720/h^3*(fAxis[k+1] - fAxis[k-1]) 
    end
    fThirdDerivAxisFiniteDiff[numPixel-4] = 24*xAxis[numPixel-4] + 36
    fThirdDerivAxisFiniteDiff[numPixel-3] = 24*xAxis[numPixel-3] + 36
    fThirdDerivAxisFiniteDiff[numPixel-2] = 24*xAxis[numPixel-2] + 36
    fThirdDerivAxisFiniteDiff[numPixel-1] = 24*xAxis[numPixel-1] + 36
    fThirdDerivAxisFiniteDiff[numPixel] = 24*xAxis[numPixel] + 36
    for k = 1:1:numPixel
        fThirdDerivCheck[k] = 24*xAxis[k] + 36
    end
    deltaFThirdDerivCheck = zeros(Float64, 1, numPixel)
    for k = 1:1:numPixel
        deltaFThirdDerivCheck[k] = abs(fThirdDerivCheck[k] - fThirdDerivAxisFiniteDiff[k])
    end
    println(norm(deltaFThirdDerivCheck)/numPixel)
end

Lmax = 10.0000
numPixel = 1000
xAxis = zeros(Float64, 1, numPixel)
for k = 1:1:numPixel
    xAxis[k] = k*Lmax/numPixel - Lmax/2
end
println("h = ", Lmax/numPixel)
firstOrderDerivSixthOrder(xAxis, numPixel)
thirdOrderDerivSixthOrder(xAxis, numPixel)
firstOrderDerivTenthOrder(xAxis, numPixel)
thirdOrderDerivTenthOrder(xAxis, numPixel)