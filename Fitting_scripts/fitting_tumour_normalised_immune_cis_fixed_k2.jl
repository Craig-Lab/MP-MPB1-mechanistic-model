using DifferentialEquations
using LsqFit
using CSV, DataFrames
include("../src/Models.jl");
include("../src/Data.jl");
include("../src/Assumptions.jl");
include("../src/Fixed_params.jl");
include("../src/Predictions.jl");
using .Predictions

data = DataTools.load_tumor_datasets(data_dir = joinpath(@__DIR__, "..", "Data"));

tdata_mp_cis = data.mp_cis.t;
ydata_mp_cis = data.mp_cis.y;
tdata_mpb1_cis = data.mpb1_cis.t;
ydata_mpb1_cis = data.mpb1_cis.y;

# Fixed bolus dosing setup
BW = 0.025;
dose_amt = 3.0 * BW;
initial_day_dose = 0.01;  # Start dosing at 0.01 days to avoid dosing at t=0 which is causing issues
dose_starts = [initial_day_dose, initial_day_dose + 7.0];  # Doses 1 week apart

# Load immune data and proportions via DataTools
immune_cis = DataTools.load_immune_datasets_cis(data_dir = joinpath(@__DIR__, "..", "Data"));
df_immune_cis_mp = immune_cis.df_immune_cis_mp;
df_immune_cis_mpb1 = immune_cis.df_immune_cis_mpb1;
nk_cells_prop_cis_mp = immune_cis.nk_cells_prop_cis_mp;
nk_cells_prop_cis_mpb1 = immune_cis.nk_cells_prop_cis_mpb1;
t_cells_prop_cis_mp = immune_cis.t_cells_prop_cis_mp;
t_cells_prop_cis_mpb1 = immune_cis.t_cells_prop_cis_mpb1;

# Cell counts for NK and T cells, using the assumptions defined in Assumptions.jl
nk_cell_count_mp = nk_cells_prop_cis_mp * prop_cd45_cells * total_cell_count;
nk_cell_count_mpb1 = nk_cells_prop_cis_mpb1 * prop_cd45_cells * total_cell_count;
t_cell_count_mp = t_cells_prop_cis_mp * prop_cd8cells_in_t_cells * prop_cd45_cells * total_cell_count;
t_cell_count_mpb1 = t_cells_prop_cis_mpb1 * prop_cd8cells_in_t_cells * prop_cd45_cells * total_cell_count;

# Scale each component so residuals are on a comparable numerical scale.
tumour_scale_mp = maximum(ydata_mp_cis)
nk_scale_mp = maximum([nk_cell_count_mp, 1.0])
t_scale_mp = maximum([t_cell_count_mp, 1.0])

tumour_scale_mpb1 = maximum(ydata_mpb1_cis)
nk_scale_mpb1 = maximum([nk_cell_count_mpb1, 1.0])
t_scale_mpb1 = maximum([t_cell_count_mpb1, 1.0])

ydata_mp_combined = vcat(ydata_mp_cis, [nk_cell_count_mp, t_cell_count_mp])
ydata_mpb1_combined = vcat(ydata_mpb1_cis, [nk_cell_count_mpb1, t_cell_count_mpb1])
ydata_mp_combined_scaled = vcat(ydata_mp_cis ./ tumour_scale_mp, [nk_cell_count_mp / nk_scale_mp, t_cell_count_mp / t_scale_mp])
ydata_mpb1_combined_scaled = vcat(ydata_mpb1_cis ./ tumour_scale_mpb1, [nk_cell_count_mpb1 / nk_scale_mpb1, t_cell_count_mpb1 / t_scale_mpb1])


# Definition of the model fixed parameters
fixed_params_mp = fixed_params_immune_cis_mp;
fixed_params_mpb1 = fixed_params_immune_cis_mpb1;

# Load fitted parameters from timeshift fitting (will be initial estimates for fitting the immune model)
time_shift_result = CSV.read(joinpath(@__DIR__, "..", "Fitted_params_results/fitted_parameters_normalised_timeshift.csv"), DataFrame)

# Case where we initialise the tumour at treatment initiation 
# Initial conditions for immune cells and tumour growth model (vehicle)
# Where T, N, E = u
# with T in mm3, N and E in cells
N0_treatment_initiation_mp = time_shift_result[time_shift_result.Dataset .== "MP", :nk_cells_at_treatment_initiation][1];
N0_treatment_initiation_mpb1 = time_shift_result[time_shift_result.Dataset .== "MPB1", :nk_cells_at_treatment_initiation][1];
E0_treatment_initiation_mp = time_shift_result[time_shift_result.Dataset .== "MP", :t_cells_at_shift_at_treatment_initiation][1];
E0_treatment_initiation_mpb1 = time_shift_result[time_shift_result.Dataset .== "MPB1", :t_cells_at_shift_at_treatment_initiation][1];

# Time at which the immune cells were measured
t_immune_mp = 10.0;
t_immune_mpb1 = 10.0;
# Initial conditions for immune cells and tumour growth model, from vehicle data
# Where C1, C2, T1, T2, T3, T4, N, E = u
u0_immune_cis_mp = [0, 0, ydata_mp_cis[1], 0, 0, 0, N0_treatment_initiation_mp, E0_treatment_initiation_mp];
u0_immune_cis_mpb1 = [0, 0, ydata_mpb1_cis[1], 0, 0, 0, N0_treatment_initiation_mpb1, E0_treatment_initiation_mpb1];

# Initial estimates for fitted params
k2N_0_mp = 1.0;  # (L/mg*day)
k2N_0_mpb1 = k2N_0_mp; # (L/mg*day)
k2E_0_mp = 1.0;  # (L/mg*day)
k2E_0_mpb1 = k2E_0_mp; # (L/mg*day)

# Bounds
k2N_min = 0.0;
k2N_max = 500.0;
k2E_min = 0.0;
k2E_max = 500.0;

# The model predictions are returned on the same scaled basis as the data.
tumour_immune_model_mp_cis = (t, p) -> begin
    pred = Predictions.predict_tumour_immune_cis_pkpd(t, t_immune_mp, p, u0_immune_cis_mp, fixed_params_mp, dose_amt, dose_starts)
    return vcat(pred[1:length(ydata_mp_cis)] ./ tumour_scale_mp,
                pred[end - 1] / nk_scale_mp,
                pred[end] / t_scale_mp)
end

tumour_immune_model_mpb1_cis = (t, p) -> begin
    pred = Predictions.predict_tumour_immune_cis_pkpd(t, t_immune_mpb1, p, u0_immune_cis_mpb1, fixed_params_mpb1, dose_amt, dose_starts)
    return vcat(pred[1:length(ydata_mpb1_cis)] ./ tumour_scale_mpb1,
                pred[end - 1] / nk_scale_mpb1,
                pred[end] / t_scale_mpb1)
end

# Fit the immune model to combined tumour and immune cell data for MP and MPB1
fit_immune_mp = curve_fit(tumour_immune_model_mp_cis, tdata_mp_cis, ydata_mp_combined_scaled, [k2N_0_mp, k2E_0_mp], 
                           lower=[k2N_min, k2E_min], upper=[k2N_max, k2E_max])
fit_immune_mpb1 = curve_fit(tumour_immune_model_mpb1_cis, tdata_mpb1_cis, ydata_mpb1_combined_scaled, [k2N_0_mpb1, k2E_0_mpb1], 
                             lower=[k2N_min, k2E_min], upper=[k2N_max, k2E_max])


# Create DataFrame for parameters and errors
param_df_immune = DataFrame(
    Dataset = ["MP", "MPB1"],
    k2N = [fit_immune_mp.param[1], fit_immune_mpb1.param[1]],
    k2E = [fit_immune_mp.param[2], fit_immune_mpb1.param[2]],
    k2N_stderr = [stderror(fit_immune_mp)[1], stderror(fit_immune_mpb1)[1]],
    k2E_stderr = [stderror(fit_immune_mp)[2], stderror(fit_immune_mpb1)[2]],
    prop_cd45_cells = fill(prop_cd45_cells, length(["MP", "MPB1"])),
    prop_cd8cells_in_t_cells = fill(prop_cd8cells_in_t_cells, length(["MP", "MPB1"]))
);

# Write to CSV file
CSV.write("./Fitted_params_results/fitted_parameters_normalised_immune_cisplatin.csv", param_df_immune)

println("Immune model fitting complete. Parameters saved to fitted_parameters_normalised_immune_cisplatin.csv")