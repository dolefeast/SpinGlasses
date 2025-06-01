using Pkg

Pkg.activate("/home/sanz/spinGlasses/calculations/")

using calculations
    
function buildParamsArray(nPointsArray, 
    coordNumberArray,
    sigmaArray,
    TArray
    )
    """
    Builds the parameter array to 
    perform calculations in bulk
    """
    params = []
    for nPoints in nPointsArray
        for z in coordNumberArray
            for sigma in sigmaArray
                for T in TArray
                    append!(params, 
                    (
                        nPoints,
                        z * nPoints,
                        sigma,
                        T
                    )
                    )
                end
            end
        end
    end
end
nPointsArray = [2^10, 2^12, 2^14, 2^16]
coordNumberArray = [2, 4, 8]
sigmaArray = [0.2, 0.8, 1.5]
TArray = [0.1, 0.2]

params = buildParamsArray(nPointsArray, coordNumberArray,
sigmaArray, TArray)
print(params)