export AICc, BIC, compute_statistics

function AICc(RSS,k,N)
# Returns the AICc using the residual sum of squares (RSS), the number of parameters (k), and the number of observations (N).
        
    if k < N/40
        AICc = N*log(RSS/N) + 2*k;
    end 
    if k >= N/40
        AICc = N*log(RSS/N) + (2*k*N)/(N-k-1); 
    end
    return AICc
    
end

function BIC(RSS,k,N)
    # Returns the BIC using the residual sum of squares (RSS), the number of parameters (k), and the number of observations (N).
        
    BIC = N*log(RSS/N) + k*log(N);
    return BIC
        
end

function compute_statistics(mp_fits, mpb1_fits, models, conds)
    results = Dict{String, DataFrame}()

    for cond in conds
        fits = cond == "MP" ? mp_fits :
               cond == "MPB1" ? mpb1_fits :
               error("Unknown condition: $cond")

        Stat = zeros(length(models), 5)

        for i in eachindex(models)
            fit = fits[i]
            k_cond = length(fit.param)
            RSS_cond = sum(fit.resid .^ 2)
            Nobs = length(fit.resid)

            Stat[i, 1] = k_cond
            Stat[i, 2] = AICc(RSS_cond, k_cond, Nobs)
            Stat[i, 3] = BIC(RSS_cond, k_cond, Nobs)
            Stat[i, 4] = RSS_cond
        end

        minAICc = minimum(Stat[:, 2])
        minBIC = minimum(Stat[:, 3])

        dAICc = Stat[:, 2] .- minAICc
        dBIC = Stat[:, 3] .- minBIC

        statTable = DataFrame(
            "Model" => models,
            "N" => Int.(Stat[:, 1]),
            "AICc" => Stat[:, 2],
            "ΔAICc" => dAICc,
            "BIC" => Stat[:, 3],
            "ΔBIC" => dBIC,
            "RSS" => Stat[:, 4]
        )

        results[cond] = statTable
    end

    return results
end