# TenthOrderNumericalIntegrationExample.jl

using Test

function tenthOrderTrapezoidalIntegrator(xAxis, fAxis, nPixels)
    finalIntegral = 0.0000
    h = (xAxis[nPixels] - xAxis[1])/nPixels
    extractedFaxis = zeros(Float64, 1, 10)
    for k = 1:1:nPixels-1
        for m = 0:1:9
            extractedFaxis[m+1] = ((10 - m)/10*fAxis[k] + m/10*fAxis[k+1]) 
        end
        finalIntegral = finalIntegral + h/89600*(
            2857*(extractedFaxis[1] + extractedFaxis[10])
            + 15741*(extractedFaxis[2] + extractedFaxis[9])
            + 1080*(extractedFaxis[3] + extractedFaxis[8])
            + 19344*(extractedFaxis[4] + extractedFaxis[7])
            + 5778*(extractedFaxis[5] + extractedFaxis[6])
            )
    end
    return finalIntegral
end

nPixels = 40000
xMax = 100
xAxis = zeros(Float64, 1, nPixels)
fAxis = zeros(Float64, 1, nPixels)
for k = 1:1:nPixels
    xAxis[k] = -xMax/2 + xMax/nPixels*(k-1)
    fAxis[k] = exp(-xAxis[k]^2)
end

numericalIntegral = tenthOrderTrapezoidalIntegrator(
                        xAxis, fAxis, nPixels
                        )
convergedIntegral = sqrt(pi)

@test(
    abs(convergedIntegral - numericalIntegral) < 1.e-4
    )