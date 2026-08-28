using DifferentialEquations
using LsqFit
using CSV, DataFrames
include("../src/Models.jl");
include("../src/Data.jl");
include("../src/Assumptions.jl");
include("../src/Fixed_params.jl");
include("../src/Predictions.jl")
using .Predictions


data = DataTools.load_tumor_datasets(data_dir = joinpath(@__DIR__, "..", "Data"))

tdata_mp = data.mp.t
ydata_mp = data.mp.y
tdata_mpb1 = data.mpb1.t
ydata_mpb1 = data.mpb1.y

# Load immune data and proportions via DataTools
immune = DataTools.load_immune_datasets(data_dir = joinpath(@__DIR__, "..", "Data"))
df_immune_mp = immune.df_immune_mp
df_immune_mpb1 = immune.df_immune_mpb1
nk_cells_prop_mp = immune.nk_cells_prop_mp
nk_cells_prop_mpb1 = immune.nk_cells_prop_mpb1
t_cells_prop_mp = immune.t_cells_prop_mp
t_cells_prop_mpb1 = immune.t_cells_prop_mpb1

# Cell counts for NK and T cells, using the assumptions defined in Assumptions.jl
nk_cell_count_mp = nk_cells_prop_mp * prop_cd45_cells * total_cell_count
nk_cell_count_mpb1 = nk_cells_prop_mpb1 * prop_cd45_cells * total_cell_count
t_cell_count_mp = t_cells_prop_mp * prop_cd8cells_in_t_cells * prop_cd45_cells * total_cell_count
t_cell_count_mpb1 = t_cells_prop_mpb1 * prop_cd8cells_in_t_cells * prop_cd45_cells * total_cell_count

# Combine tumour volume, NK cell count, and T cell count data for fitting
# Scale each component so residuals are on a comparable numerical scale.
tumour_scale_mp = maximum(ydata_mp)
nk_scale_mp = maximum([nk_cell_count_mp, 1.0])
t_scale_mp = maximum([t_cell_count_mp, 1.0])

tumour_scale_mpb1 = maximum(ydata_mpb1)
nk_scale_mpb1 = maximum([nk_cell_count_mpb1, 1.0])
t_scale_mpb1 = maximum([t_cell_count_mpb1, 1.0])

ydata_mp_combined = vcat(ydata_mp, [nk_cell_count_mp, t_cell_count_mp])
ydata_mpb1_combined = vcat(ydata_mpb1, [nk_cell_count_mpb1, t_cell_count_mpb1])
ydata_mp_combined_scaled = vcat(ydata_mp ./ tumour_scale_mp, [nk_cell_count_mp / nk_scale_mp, t_cell_count_mp / t_scale_mp])
ydata_mpb1_combined_scaled = vcat(ydata_mpb1 ./ tumour_scale_mpb1, [nk_cell_count_mpb1 / nk_scale_mpb1, t_cell_count_mpb1 / t_scale_mpb1])

# Definition of the model fixed parameters
fixed_params_mp = fixed_params_immune_mp
fixed_params_mpb1 = fixed_params_immune_mpb1


# Load fitted parameters from timeshift fitting (will be initial estimates for fitting the immune model)
time_shift_result = CSV.read(joinpath(@__DIR__, "..", "Fitted_params_results/fitted_parameters_normalised_timeshift.csv"), DataFrame)

# Case where we initialise the tumour at treatment initiation 
# Initial conditions for immune cells and tumour growth model (vehicle)
# Where T, N, E = u
# with T in mm3, N and E in cells
N0_treatment_initiation_mp = time_shift_result[time_shift_result.Dataset .== "MP", :nk_cells_at_treatment_initiation][1]
N0_treatment_initiation_mpb1 = time_shift_result[time_shift_result.Dataset .== "MPB1", :nk_cells_at_treatment_initiation][1]
E0_treatment_initiation_mp = time_shift_result[time_shift_result.Dataset .== "MP", :t_cells_at_shift_at_treatment_initiation][1]
E0_treatment_initiation_mpb1 = time_shift_result[time_shift_result.Dataset .== "MPB1", :t_cells_at_shift_at_treatment_initiation][1]


# Initial estimates for fitted params
r_0_mp = time_shift_result[time_shift_result.Dataset .== "MP", :r][1]
r_0_mpb1 = time_shift_result[time_shift_result.Dataset .== "MPB1", :r][1]
pN_0_mp = time_shift_result[time_shift_result.Dataset .== "MP", :pN][1]
pN_0_mpb1 = time_shift_result[time_shift_result.Dataset .== "MPB1", :pN][1]
q_0_mp = time_shift_result[time_shift_result.Dataset .== "MP", :q][1]
q_0_mpb1 = time_shift_result[time_shift_result.Dataset .== "MPB1", :q][1]

# Bounds
r_min = 0.0
r_max = 5.0
pN_min = 0.00001
pN_max = 1.0
q_min = 0.00001
q_max = 1.0

# Time points for immune cell data (experimental t=10, shifted to model time)
t_immune_mp = 10.0
t_immune_mpb1 = 10.0
u0_immune_mp = [ydata_mp[1], N0_treatment_initiation_mp, E0_treatment_initiation_mp]
u0_immune_mpb1 = [ydata_mpb1[1], N0_treatment_initiation_mpb1, E0_treatment_initiation_mpb1]

tumour_immune_model_mp_with_immune = (t, p) -> Predictions.predict_tumour_immune(t, t_immune_mp, p, u0_immune_mp, fixed_params_mp)
tumour_immune_model_mpb1_with_immune = (t, p) -> Predictions.predict_tumour_immune(t, t_immune_mpb1, p, u0_immune_mpb1, fixed_params_mpb1)

# Fit the immune model to combined tumour and immune cell data for MP and MPB1
# The model predictions are returned on the same scaled basis as the data.
tumour_immune_model_mp_with_immune = (t, p) -> begin
    pred = Predictions.predict_tumour_immune(t, t_immune_mp, p, u0_immune_mp, fixed_params_mp)
    return vcat(pred[1:length(ydata_mp)] ./ tumour_scale_mp,
                pred[end - 1] / nk_scale_mp,
                pred[end] / t_scale_mp)
end

tumour_immune_model_mpb1_with_immune = (t, p) -> begin
    pred = Predictions.predict_tumour_immune(t, t_immune_mpb1, p, u0_immune_mpb1, fixed_params_mpb1)
    return vcat(pred[1:length(ydata_mpb1)] ./ tumour_scale_mpb1,
                pred[end - 1] / nk_scale_mpb1,
                pred[end] / t_scale_mpb1)
end

fit_immune_mp = curve_fit(tumour_immune_model_mp_with_immune, tdata_mp, ydata_mp_combined_scaled, [r_0_mp, pN_0_mp, q_0_mp], 
                           lower=[r_min, pN_min, q_min], upper=[r_max, pN_max, q_max])
fit_immune_mpb1 = curve_fit(tumour_immune_model_mpb1_with_immune, tdata_mpb1, ydata_mpb1_combined_scaled, [r_0_mpb1, pN_0_mpb1, q_0_mpb1], 
                             lower=[r_min, pN_min, q_min], upper=[r_max, pN_max, q_max])


# Create DataFrame for parameters and errors
param_df_immune = DataFrame(
    Dataset = ["MP", "MPB1"],
    r = [fit_immune_mp.param[1], fit_immune_mpb1.param[1]],
    pN = [fit_immune_mp.param[2], fit_immune_mpb1.param[2]],
    q = [fit_immune_mp.param[3], fit_immune_mpb1.param[3]],
    r_stderr = [stderror(fit_immune_mp)[1], stderror(fit_immune_mpb1)[1]],
    pN_stderr = [stderror(fit_immune_mp)[2], stderror(fit_immune_mpb1)[2]],
    q_stderr = [stderror(fit_immune_mp)[3], stderror(fit_immune_mpb1)[3]],
    prop_cd45_cells = fill(prop_cd45_cells, length(["MP", "MPB1"])),
    prop_cd8cells_in_t_cells = fill(prop_cd8cells_in_t_cells, length(["MP", "MPB1"]))
)


# Write to CSV file
CSV.write("./Fitted_params_results/fitted_parameters_normalised_immune.csv", param_df_immune)


println("Immune model fitting complete. Parameters saved to fitted_parameters_normalised_immune.csv")