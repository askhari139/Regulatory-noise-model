using CSV, DataFrames, JLD2

include("lambda_sampler.jl")

function save_lambda(dt, noiseType)
    fl = joinpath("data", join(["lambda_hist_", noiseType, "_dt",string(dt), ".jld2"]))
    if (!ispath(fl))
        println(string(dt)*"_"*noiseType*" not found.")
        return()
    end
    @load fl lambda_inh lambda_act
    df = DataFrame(LambdaType=String[], NoiseLevel=Float64[], LambdaInit=Float64[],
               BinMid=Float64[], Count=Int[])

    for (cache, ltype) in zip([lambda_inh, lambda_act], ["inh", "act"])
        for (noise, inner) in cache
            for (lI, hist) in inner
                mids = (hist.edges[1:end-1] .+ hist.edges[2:end]) ./ 2
                n = length(mids)
                append!(df, DataFrame(LambdaType=fill(ltype,n), NoiseLevel=fill(noise,n),
                                    LambdaInit=fill(lI,n), BinMid=mids, Count=hist.counts))
            end
        end
    end

    CSV.write("data/lambda_hist_flat_"*noiseType*"_dt"*string(dt)*".csv", df)
end

noiseTypes = ["Additive", "Multiplicative", "Fluctuating"]
dt = 1.0
for ns in noiseTypes
    save_lambda(dt, ns)
end
