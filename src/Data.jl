module DataTools

using CSV, DataFrames, Statistics

export load_dataset, load_tumor_datasets, calculate_error_bars, load_immune_datasets, load_immune_datasets_cis

"""
    calculate_error_bars(ydata::AbstractVector, ymin::AbstractVector, ymax::AbstractVector)

Calculate error bars as tuples of (upper_error, lower_error) from observed data and min/max bounds.

# Arguments
- `ydata`: Vector of observed values
- `ymin`: Vector of minimum bounds
- `ymax`: Vector of maximum bounds

# Returns
Tuple of (upper_error, lower_error) where:
- `upper_error = ymax - ydata`
- `lower_error = ydata - ymin`
"""
function calculate_error_bars(ydata::AbstractVector, ymin::AbstractVector, ymax::AbstractVector)
    upper_error = ymax .- ydata
    lower_error = ydata .- ymin
    return (upper_error, lower_error)
end

function load_dataset(filename::AbstractString; data_dir::AbstractString = joinpath(@__DIR__, "..", "..", "Data"))
    path = joinpath(data_dir, filename)
    df = CSV.read(path, DataFrame)
    t = df[:, 1]
    y = df[:, 2]
    u0 = [y[1]]
    
    # Extract error bounds if they exist (columns 3 and 4)
    ymin = ncol(df) >= 3 ? df[:, 3] : nothing
    ymax = ncol(df) >= 4 ? df[:, 4] : nothing
    
    # Calculate error bars if bounds exist
    yerr = (ymin !== nothing && ymax !== nothing) ? calculate_error_bars(y, ymin, ymax) : nothing
    
    return (df = df, t = t, y = y, u0 = u0, ymin = ymin, ymax = ymax, yerr = yerr)
end

function load_tumor_datasets(; data_dir::AbstractString = joinpath(@__DIR__, "..", "..", "Data"))
    return (
        mp = load_dataset("MP_Vehicle_Data.csv", data_dir = data_dir),
        mp_cis = load_dataset("MP_Cisplatin_Data.csv", data_dir = data_dir),
        mp_icb = load_dataset("MP_ICB_Data.csv", data_dir = data_dir),
        mp_cis_icb = load_dataset("MP_Cisplatin_ICB_Data.csv", data_dir = data_dir),
        mpb1 = load_dataset("MPB1_Vehicle_Data.csv", data_dir = data_dir),
        mpb1_cis = load_dataset("MPB1_Cisplatin_Data.csv", data_dir = data_dir),
        mpb1_icb = load_dataset("MPB1_ICB_Data.csv", data_dir = data_dir),
        mpb1_cis_icb = load_dataset("MPB1_Cisplatin_ICB_Data.csv", data_dir = data_dir)
    )
end

function load_immune_datasets(; data_dir::AbstractString = joinpath(@__DIR__, "..", "..", "Data"))
    df_immune_mp = CSV.read(joinpath(data_dir, "MP_Vehicle_Immune_Data.csv"), DataFrame)
    df_immune_mpb1 = CSV.read(joinpath(data_dir, "MPB1_Vehicle_Immune_Data.csv"), DataFrame)

    nk_cells_prop_mp = mean(df_immune_mp[!, "Proportion of NK cells (vehicle)"]) / 100
    nk_cells_prop_mp_sem = (std(df_immune_mp[!, "Proportion of NK cells (vehicle)"]) / sqrt(nrow(df_immune_mp)))/100
    nk_cells_prop_mpb1 = mean(df_immune_mpb1[!, "Proportion of NK cells (vehicle)"]) / 100
    nk_cells_prop_mpb1_sem = (std(df_immune_mpb1[!, "Proportion of NK cells (vehicle)"]) / sqrt(nrow(df_immune_mpb1)))/100

    t_cells_prop_mp = mean(df_immune_mp[!, "Proportion of T cells (vehicle)"]) / 100
    t_cells_prop_mp_sem = (std(df_immune_mp[!, "Proportion of T cells (vehicle)"]) / sqrt(nrow(df_immune_mp)))/100
    t_cells_prop_mpb1 = mean(df_immune_mpb1[!, "Proportion of T cells (vehicle)"]) / 100
    t_cells_prop_mpb1_sem = (std(df_immune_mpb1[!, "Proportion of T cells (vehicle)"]) / sqrt(nrow(df_immune_mpb1)))/100

    return (
        df_immune_mp = df_immune_mp,
        df_immune_mpb1 = df_immune_mpb1,
        nk_cells_prop_mp = nk_cells_prop_mp,
        nk_cells_prop_mp_sem = nk_cells_prop_mp_sem,
        nk_cells_prop_mpb1 = nk_cells_prop_mpb1,
        nk_cells_prop_mpb1_sem = nk_cells_prop_mpb1_sem,
        t_cells_prop_mp = t_cells_prop_mp,
        t_cells_prop_mp_sem = t_cells_prop_mp_sem,
        t_cells_prop_mpb1 = t_cells_prop_mpb1,
        t_cells_prop_mpb1_sem = t_cells_prop_mpb1_sem
    )
end

function load_immune_datasets_cis(; data_dir::AbstractString = joinpath(@__DIR__, "..", "..", "Data"))
    df_immune_cis_mp = CSV.read(joinpath(data_dir, "MP_Cisplatin_Immune_Data.csv"), DataFrame)
    df_immune_cis_mpb1 = CSV.read(joinpath(data_dir, "MPB1_Cisplatin_Immune_Data.csv"), DataFrame)

    nk_cells_prop_cis_mp = mean(df_immune_cis_mp[!, "Proportion of NK cells (cisplatin)"]) / 100
    nk_cells_prop_cis_mp_sem = (std(df_immune_cis_mp[!, "Proportion of NK cells (cisplatin)"]) / sqrt(nrow(df_immune_cis_mp)))/100
    nk_cells_prop_cis_mpb1 = mean(df_immune_cis_mpb1[!, "Proportion of NK cells (cisplatin)"]) / 100
    nk_cells_prop_cis_mpb1_sem = (std(df_immune_cis_mpb1[!, "Proportion of NK cells (cisplatin)"]) / sqrt(nrow(df_immune_cis_mpb1)))/100

    t_cells_prop_cis_mp = mean(df_immune_cis_mp[!, "Proportion of T cells (cisplatin)"]) / 100
    t_cells_prop_cis_mp_sem = (std(df_immune_cis_mp[!, "Proportion of T cells (cisplatin)"]) / sqrt(nrow(df_immune_cis_mp)))/100
    t_cells_prop_cis_mpb1 = mean(df_immune_cis_mpb1[!, "Proportion of T cells (cisplatin)"]) / 100
    t_cells_prop_cis_mpb1_sem = (std(df_immune_cis_mpb1[!, "Proportion of T cells (cisplatin)"]) / sqrt(nrow(df_immune_cis_mpb1)))/100

    return (
        df_immune_cis_mp = df_immune_cis_mp,
        df_immune_cis_mpb1 = df_immune_cis_mpb1,
        nk_cells_prop_cis_mp = nk_cells_prop_cis_mp,
        nk_cells_prop_cis_mp_sem = nk_cells_prop_cis_mp_sem,
        nk_cells_prop_cis_mpb1 = nk_cells_prop_cis_mpb1,
        nk_cells_prop_cis_mpb1_sem = nk_cells_prop_cis_mpb1_sem,
        t_cells_prop_cis_mp = t_cells_prop_cis_mp,
        t_cells_prop_cis_mp_sem = t_cells_prop_cis_mp_sem,
        t_cells_prop_cis_mpb1 = t_cells_prop_cis_mpb1,
        t_cells_prop_cis_mpb1_sem = t_cells_prop_cis_mpb1_sem
    )
end

end
