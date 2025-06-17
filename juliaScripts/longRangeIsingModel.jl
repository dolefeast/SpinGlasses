using Dates
using DataFrames
using CSV
using LaTeXStrings

function norm(arr) # Calculates L1 norm of an array
    arrNorm = sum(abs.(arr))
    arrNorm
end

function coinflip()
    return 2*(rand()>0.5) - 1
end

function evolveIsing(params, nCurves=50, flips=10^6, downsizeFactor=1000)

    calculations.T, calculations.nPoints, calculations.nBonds, calculations.sigma = params
    local m = [ [ ] for curve in 1:nCurves]


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
    df = DataFrame(collect(eachrow(reduce(hcat, m))), string.(flipsPlot))
    return df
end

function msquared(m)
    nCurves = length(m)

    return calculations.arrayStatistics(m)
end

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

function evolveAndPlot(params, nCurves=50, flips=100000, downsizeFactor=100)
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

# Need: Ising runs for the parameters in params
nPointsArray = [1024, 4096, 8192, 16384]
coordNumberArray = [2, 4, 8]
sigmaArray = [0.2, 0.8, 1.5]
TArray = [0.1, 0.2]

paramsArray = calculations.IOutils.buildParamsArray(TArray, nPointsArray, coordNumberArray,
sigmaArray)

for params in paramsArray
    m = @timed evolveIsing(params)
    CSV.write("/home/dein/spinGlasses/juliaScripts/results/Ising Runs/$(calculations.IOutils.filename(params...)).csv", m[1])
    write("/home/dein/spinGlasses/juliaScripts/results/Ising Runs/Time Logs/$(calculations.IOutils.filename(params...)).csv", "Date: $(now()), time elapsed = $(m[2])")
    
end