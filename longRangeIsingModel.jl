using LaTeXStrings
using Dates

function norm(arr) # Calculates L1 norm of an array
    arrNorm = sum(abs.(arr))
    arrNorm
end;

function coinflip()
    return 2*(rand()>0.5) - 1
end

function evolveIsing(nCurves, flips, downsizeFactor)

    m = [ [ ] for curve in 1:nCurves]


    for curve in 1:nCurves # number of monte carlo chains to run

        connectionsSet = calculations.generateConnectionSet()
        connectionsArray = calculations.connectionSetToConnectivityArrayConversion(connectionsSet)

        spins = [ coinflip() for _ in 1:nPoints ]
        
        global flipsPlot # to pull it out of the for loop
        flipsPlot = []

        for flip in 1:flips # time
            spins = calculations.metropolis.singleStep(calculations.nPoints, calculations.T, connectionsArray, spins)
            
            if flip%downsizeFactor == 1 # It's not desired to save every single step
                append!(flipsPlot, flip)
                append!(m[curve], sum(spins) / nPoints)
            end
        end

    end
    return flipsPlot, m
end

function msquared(m)
    nCurves = length(m)

    m2Mean = [ 0. for _ in m[1] ]
    m2STD =  [ 0. for _ in m[1] ]


    for flip in 1:length(m[1])
        mean = 0.
        stdSquared = 0.
        for curve in m
            mean += curve[flip] ^ 2
        end
        mean /= nCurves
        
        for (i, curve) in enumerate(m)
            stdSquared += (curve[flip]^2 - mean) ^2
        end
        std = sqrt(stdSquared / nCurves)
        m2Mean[flip] = mean
        m2STD[flip] = std
    end
    m2Mean, m2STD
end
#print(2lkjsadf)

function plotIsing(m, flipsPlot)
    p = plot(grid = false) 
    
    title = L"N = %$(calculations.nPoints), N_l = %$(calculations.nBonds), \sigma = %$(calculations.sigma), T = %$T"
    
    plot!(p, flipsPlot, m, 
        legend=false, xlabel=L"t", ylabel=L"m" )
    plot!(p, title = title)

    
    savefig("figures/isingRuns/T_$(T)_nPoints_$(calculations.nPoints)_nBonds_$(calculations.nBonds)_sigma_$(calculations.sigma).pdf")
end

function plotm2(m2stats, flipsPlot)
    m2Mean, m2STD = m2stats
    plot(xaxis=:log, yaxis=:log, grid=false)

    plot!(flipsPlot, m2Mean, ribbon=m2STD/2, grid=false, legend=false)
    plot!(xlabel=L"t", ylabel=L"\langle m^2 \rangle")
    title = L"N = %$(calculations.nPoints), N_l = %$(calculations.nBonds), \sigma = %$(calculations.sigma), T = %$T"
    plot!(title = title)
    
    savefig("figures/m2Runs/T_$(T)_nPoints_$(calculations.nPoints)_nBonds_$(calculations.nBonds)_sigma_$(calculations.sigma).pdf")
end

function evolveAndPlot(params, nCurves, flips, downsizeFactor)
    nPoints, nBonds, sigma, T = params
        
    calculations.nPoints = nPoints
    calculations.nBonds = nBonds
    calculations.sigma = sigma
    calculations.T = T

    flipsPlot, m = evolveIsing(nCurves, flips, downsizeFactor)
    m2stats= msquared(m)

    plotIsing(m, flipsPlot)
    plotm2(m2stats, flipsPlot)
end