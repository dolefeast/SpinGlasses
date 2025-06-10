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
    df = CSV.File(filepath) |> DataFrame
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

function plotIsingRuns(results)
    for (i, _isingRuns) in enumerate(results[!, :isingrun])
        
        totalPlot = plot(bottom_margin=15Plots.mm, grid=false)# , ylims=(0,0.5))
        statsPlot = plot(bottom_margin=15Plots.mm, grid=false)# , ylims=(0,0.5))

        T = results[i,:T]
        nPoints = results[i,:nPoints]
        nBonds = results[i,:nBonds]
        sigma = results[i,:sigma]
        
        isingRuns = [ [ _isingRuns[i,j] 
            for j in 1:size(_isingRuns)[2]] 
            for i in 1:size(_isingRuns)[1]]


        m2Stats = calculations.arrayStatistics(isingRuns, x->x^2)

        plot!(totalPlot, (1:100:100000, isingRuns), legend=false) 
        plot!(statsPlot, (1:100:100000, m2Stats[1]), ribbon=m2Stats[2]/2, legend=false, xaxis=:log, yaxis=:log) 
            #legend=L"N_l = %$(nBonds)")
                
                plot!(statsPlot
        , title=L"N = %$(nPoints), σ = %$(sigma), N_l =%$(nBonds), T = %$(T)", xlabel=L"t", ylabel=L"\langle m^2 \rangle")
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
                for (i, _isingRuns) in enumerate(fixedSigmadf[:, :isingrun])
                                
                isingRuns = [ [ _isingRuns[i,j] 
                    for j in 1:size(_isingRuns)[2]] 
                    for i in 1:size(_isingRuns)[1]]

                m2Stats = calculations.arrayStatistics(isingRuns, x->x^2)

                plot!(totalPlot, (1:100:100000, isingRuns), legend=false) 
                plot!(statsPlot, (1:100:100000, m2Stats[1]), ribbon=m2Stats[2]/2, xaxis=:log, yaxis=:log, label=L"N_l = %$(fixedSigmadf[i, :nBonds])") 
                    #legend=L"N_l = %$(nBonds)")
                    end
                        
                        plot!(statsPlot
                , title=L"N = %$(nPoints), \sigma = %$(sigma), T = %$(T)", xlabel=L"t", ylabel=L"\langle m^2 \rangle", legend=:topleft)

                #savefig(totalPlot, "/home/dein/spinGlasses/figures/Ising Runs/Sigma groupings/T_$(T)_nPoints_$(nPoints)_sigma_$(sigma).pdf")
                savefig(statsPlot, "/home/dein/spinGlasses/figures/m2 Runs/Sigma groupings/T_$(T)_nPoints_$(nPoints)_sigma_$(sigma).pdf")
            end
        end
    end
end