using CSV
using DataFrames
using DataFramesMeta
using Glob
using LaTeXStrings
using Printf
using Plots
using calculations

# Define the path
directory = "results/Cluster Counting"

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
    clustersizes = Matrix(df)

    # Store in the results list
    push!(results, (
        T = T,
        nPoints = nPoints,
        nBonds = nBonds,
        sigma = sigma,
        clustersizes = clustersizes
    ))
end

nPointsGroupings = groupby(results, :nPoints)
for fixedNPoints in nPointsGroupings
    sigmaGroupings = groupby(fixedNPoints, :sigma)
    
    totalPlot = plot(layout = (length(sigmaGroupings),1), size=(300, 1500), bottom_margin=15Plots.mm, grid=false)# , ylims=(0,0.5))

    statsPlot = plot(layout = (length(sigmaGroupings),1), size=(300, 1500), bottom_margin=15Plots.mm, grid=false)# , ylims=(0,0.5))
    leftZoom = plot(layout = (length(sigmaGroupings),1), size=(300, 1500), bottom_margin=15Plots.mm, grid=false)# , ylims=(0,0.5))
    rightZoom = plot(layout = (length(sigmaGroupings),1), size=(300, 1500), bottom_margin=15Plots.mm, grid=false)# , ylims=(0,0.5))

    nPoints = fixedNPoints[1, :nPoints]
    for (j, fixedsigma) in enumerate(sigmaGroupings)

        sigma = fixedsigma[1, :sigma]
        nBondsGroupings = groupby(sort!(fixedsigma, :nBonds), :nBonds)

        for (i, fixedNBonds) in enumerate(nBondsGroupings)
            nBonds = fixedNBonds[1, :nBonds]
            for _clustersizes in  fixedNBonds[!,:clustersizes]
                clustersizes = [ [ _clustersizes[j,i] 
                    for j in 1:size(_clustersizes)[1]] 
                    for i in 1:size(_clustersizes)[2]]
                if i%4==0
                    plot!(totalPlot, 
                    clustersizes, sp=j, legend=false) #, legend=L"N_l = $(fixedNBonds[1, :nBonds])",
                end
                clustersizesStats = calculations.arrayStatistics(clustersizes)
                norm = sum(clustersizesStats[1])

                plot!(statsPlot, clustersizesStats[1]/norm, ribbon=clustersizesStats[2]/2 / norm, label=L"N_l = %$(nBonds)", sp=j) 
              
                fraction  = 10/nPoints
                leftRegion = 1:Int(floor(nPoints * fraction))
                rightRegion = Int(floor(nPoints * (1 - fraction))):(nPoints - 1)
                plot!(leftZoom, (leftRegion, clustersizesStats[1][leftRegion]/norm), ribbon=clustersizesStats[2][leftRegion]/2 / norm, label=L"N_l = %$(nBonds)", sp=j)
                plot!(rightZoom, (rightRegion, clustersizesStats[1][rightRegion]/norm), ribbon=clustersizesStats[2][rightRegion]/2 / norm, label=L"N_l = %$(nBonds)", sp=j)
                #legend=L"N_l = %$(nBonds)")
            end
        end


        plot!(statsPlot,
 title=L"N = %$(nPoints), \sigma = %$(sigma)", sp=j, xlabel=L"\textrm{Cluster\ Size}")
        plot!(totalPlot,
 title=L"N = %$(nPoints), \sigma = %$(sigma)", sp=j, xlabel=L"\textrm{Cluster\ Size}")
        plot!(leftZoom,
        title=L"N = %$(nPoints), \sigma = %$(sigma)", sp=j, xlabel=L"\textrm{Cluster\ Size}")
        plot!(rightZoom,
        title=L"N = %$(nPoints), \sigma = %$(sigma)", sp=j, xlabel=L"\textrm{Cluster\ Size}")
    end
    savefig(leftZoom, "/home/dein/spinGlasses/figures/ClusterLengthDistributionWithNBonds/zoom-in/Left/nPoints_$(nPoints).pdf")
    savefig(rightZoom, "/home/dein/spinGlasses/figures/ClusterLengthDistributionWithNBonds/zoom-in/Right/nPoints_$(nPoints).pdf")
    savefig(statsPlot, "/home/dein/spinGlasses/figures/ClusterLengthDistributionWithNBonds/nPoints_$(nPoints).pdf")
    savefig(totalPlot, "/home/dein/spinGlasses/figures/ClusterLengthDistributionWithNBonds/Before Averaging/nPoints_$(nPoints).pdf")
end