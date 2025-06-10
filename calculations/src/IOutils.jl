module IOutils

function buildParamsArray(
    TArray,
    nPointsArray,
    coordNumberArray,
    sigmaArray,
        )
    """
    Builds the parameter array to 
    perform calculations in bulk. 
    Calculates the cartesian product of the given arrays.
    """

    params = []
    for nPoints in nPointsArray
        for z in coordNumberArray
            for sigma in sigmaArray
                for T in TArray
                    push!(params,                    
                    ((  
                        T,
                        Int(nPoints),
                        Int(z * nPoints),
                        sigma,
                    ))
                    )
                end
            end
        end
    end
    params
end

function filename(T::Float64, nPoints::Int64, nBonds::Int64, sigma::Float64)
    return "T_$(T)_nPoints_$(nPoints)_nBonds_$(nBonds)_sigma_$(sigma)"
end

end