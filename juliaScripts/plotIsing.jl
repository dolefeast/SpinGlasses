using CSV
using DataFrames
using DataFramesMeta
using Glob
using LaTeXStrings
using Printf
using Plots
using calculations

# Define the path
directory = "results/Ising Runs/"
function readIsingRuns(directory)

# Dictionary to hold results
    results = DataFrame()

    # Define regex pattern to extract parameters from filenames
    pattern = r"T_([0-9.]+)_nPoints_(\d+)_nBonds_(\d+)_sigma_([0-9.]+)\.csv"

    # Get list of all .csv files in the directory
    for filepath in Glob.glob("*.csv", directory)
        filename = split(filepath, '/')[end]

        m = match(pattern, filename)
        if isnothing(m)
            @warn "Filename did not match expected format: $filename"
            continue
        end

        # Extract parameters
        T = parse(Float64, m.captures[1])
        nPoints = parse(Int, m.captures[2])
        nBonds = parse(Int, m.captures[3])
        sigma = parse(Float64, m.captures[4])

        # Load CSV into matrix
        df = CSV.File(filepath, header=0) |> DataFrame
        isingrun = Matrix(df)

        # Store in the results list
        push!(results, (
            T = T,
            nPoints = nPoints,
            nBonds = nBonds,
            sigma = sigma,
            isingrun = isingrun
        ))
    end
    results
end
function plotIsingRuns(results)
    for (i, _isingRuns) in enumerate(results[!, :isingrun])
        
        totalPlot = plot(bottom_margin=15Plots.mm, grid=false)# , ylims=(0,0.5))
        statsPlot = plot(bottom_margin=15Plots.mm, grid=false)# , ylims=(0,0.5))

        T = results[i,:T]
        nPoints = results[i,:nPoints]
        nBonds = results[i,:nBonds]
        sigma = results[i,:sigma]

        flipsPlot = _isingRuns[1, :]

        isingRuns = [ [ _isingRuns[i,j]
            for j in 1:size(_isingRuns)[2]]
            for i in 2:size(_isingRuns)[1]]


        m2Stats = calculations.arrayStatistics(isingRuns, x->x^2)

        plot!(totalPlot, (flipsPlot ./ nPoints, isingRuns), legend=false)
        plot!(statsPlot, (flipsPlot ./ nPoints, m2Stats[1]), ribbon=m2Stats[2]/2, legend=false, xaxis=:log, yaxis=:log)
            #legend=L"N_l = %$(nBonds)")
                
                plot!(statsPlot
        , title=L"N = %$(nPoints), σ = %$(sigma), N_l =%$(nBonds), T = %$(T)", xlabel=L"t",
        ylabel=L"\langle m^2 \rangle",
        xticks=10. .^ ( floor(log(flipsPlot[1]/nPoints) +1):floor(1+log10(flipsPlot[end]/nPoints))),
        minorticks=10)
                plot!(totalPlot
        , title=L"N = %$(nPoints), σ = %$(sigma), N_l =%$(nBonds), T = %$(T)", xlabel=L"t", ylabel=L"m")

        savefig(totalPlot, "/home/dein/spinGlasses/figures/Ising Runs/T_$(T)_nPoints_$(nPoints)_nBonds_$(nBonds)_sigma_$(sigma).pdf")
        savefig(statsPlot, "/home/dein/spinGlasses/figures/m2 Runs/T_$(T)_nPoints_$(nPoints)_nBonds_$(nBonds)_sigma_$(sigma).pdf")
    end
end

# Repeats code; should have a way to more comfortably plot the data
function plotIsingRunsSigmaGroupings(results)

    TGroupings = groupby(results, :T)
    for fixedTdf in TGroupings
    nPointsGroupings = groupby(fixedTdf, :nPoints)
    T = fixedTdf[1, :T]
        for fixednPointsdf in nPointsGroupings
            nPoints = fixednPointsdf[1, :nPoints]
            sigmaGroupings = groupby(fixednPointsdf, :sigma)

            for fixedSigmadf in sigmaGroupings
                totalPlot = plot(bottom_margin=15Plots.mm, grid=false , ylims=(0.00001,1))
                statsPlot = plot(bottom_margin=15Plots.mm, grid=false , ylims=(0.00001,1))
            
                sigma = fixedSigmadf[1, :sigma]
                for (i, _isingRuns) in enumerate(sort!(fixedSigmadf, :nBonds)[:, :isingrun])

                    flipsPlot = _isingRuns[1, :]

                    isingRuns = [ [ _isingRuns[i,j]
                        for j in 1:size(_isingRuns)[2]]
                        for i in 2:size(_isingRuns)[1]]

                    m2Stats = calculations.arrayStatistics(isingRuns, x->x^2)

                    plot!(totalPlot, (flipsPlot / nPoints, isingRuns), legend=false)

                    plot!(statsPlot, (flipsPlot / nPoints, m2Stats[1]),ribbon=m2Stats[2]/2, xaxis=:log, yaxis=:log,
                     xticks=10. .^ ( floor(log(flipsPlot[1]/nPoints) +1):floor(1+log10(flipsPlot[end]/nPoints))),label=L"N_l = %$(fixedSigmadf[i, :nBonds])")
                    #legend=L"N_l = %$(nBonds)")
                end
                        
                        plot!(statsPlot
                , title=L"N = %$(nPoints), \sigma = %$(sigma), T = %$(T)", xlabel=L"t", ylabel=L"\langle m^2 \rangle", legend=:topleft, minorticks=10, mayorticks=5)

                #savefig(totalPlot, "/home/dein/spinGlasses/figures/Ising Runs/Sigma groupings/T_$(T)_nPoints_$(nPoints)_sigma_$(sigma).pdf")
                savefig(statsPlot, "/home/dein/spinGlasses/figures/m2 Runs/Sigma groupings/T_$(T)_nPoints_$(nPoints)_sigma_$(sigma).pdf")
            end
        end
    end
end


function searchIndexFromBelow(x, arr)
    """
    Returns the maximum i::Int such that arr[i]<x,
    such that el < x. arr must be an ordered set
    """
    if x > arr[end] || x < arr[1]
        return nothing
    end
    for (i, el) in enumerate(arr)
        if el > x
            return i-1
        end
    end
end