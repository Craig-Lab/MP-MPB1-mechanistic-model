module Results

using DataFrames
using LsqFit: stderror

export build_param_df

# build_param_df(dataset_names, fits, param_names, indices)
# param_names: Vector{Symbol} e.g. [:r, :k]
# indices: Vector{Int} matching param positions in fit.param
function build_param_df(dataset_names::Vector{String}, fits::Vector, param_names::Vector{Symbol}, indices::Vector{Int})
    df = DataFrame(Dataset = dataset_names)
    for (name, idx) in zip(param_names, indices)
        vals = [fit.param[idx] for fit in fits]
        errs = [stderror(fit)[idx] for fit in fits]
        df[!, name] = vals
        df[!, Symbol(string(name, "_stderr"))] = errs
    end
    return df
end

end # module
