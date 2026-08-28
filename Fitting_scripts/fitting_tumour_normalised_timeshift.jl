using DifferentialEquations
using LsqFit
using CSV, DataFrames
include("../src/Models.jl")
include("../src/Data.jl")
include("../src/Assumptions.jl")
include("../src/Fixed_params.jl")


data = DataTools.load_tumor_datasets()

tdata_mp = data.mp.t
ydata_mp = data.mp.y
tdata_mpb1 = data.mpb1.t
ydata_mpb1 = data.mpb1.y
u0_mp = [tumour_initial_size]
u0_mpb1 = [tumour_initial_size]


# Defining the parameters for the logistic fit
logistic_params = CSV.read(joinpath(@__DIR__, "..", "Fitted_params_results", "fitted_parameters_logistic.csv"), DataFrame)
r_0_mp_logistic_fit = logistic_params[logistic_params.Dataset .== "MP", :r][1]
r_0_mpb1_logistic_fit = logistic_params[logistic_params.Dataset .== "MPB1", :r][1]
k_0_mp_logistic_fit = logistic_params[logistic_params.Dataset .== "MP", :k][1]
k_0_mpb1_logistic_fit = logistic_params[logistic_params.Dataset .== "MPB1", :k][1]

# Parameters for the logistic model fits (r, k)
p_log_mp = (r_0_mp_logistic_fit, k_0_mp_logistic_fit)
p_log_mpb1 = (r_0_mpb1_logistic_fit, k_0_mpb1_logistic_fit)

# Finding at what time the tumour reaches 106 mm3 for MP and 62 mm3 for MPB1 using the logistic model with the fitted parameters
prob_mp = ODEProblem(tumour_logistic!, u0_mp, (0.0, 300.0), p_log_mp)
prob_mpb1 = ODEProblem(tumour_logistic!, u0_mpb1, (0.0, 300.0), p_log_mpb1)
T0_mp = 106.0  # mm3
T0_mpb1 = 62.0  # mm3
cb_mp = ContinuousCallback((u, t, integrator) -> u[1] - T0_mp, terminate!; rootfind = true, save_positions = (true, false))
cb_mpb1 = ContinuousCallback((u, t, integrator) -> u[1] - T0_mpb1, terminate!; rootfind = true, save_positions = (true, false))
sol_mp = solve(prob_mp, Vern7(), callback = cb_mp, abstol = 1e-8, reltol = 1e-6)
sol_mpb1 = solve(prob_mpb1, Vern7(), callback = cb_mpb1, abstol = 1e-8, reltol = 1e-6)

# Debug output to check if the tumour reaches the target sizes within 300 days
# if sol_mp.t[end] < 300.0 && sol_mp.u[end][1] >= T0_mp - 1e-8
#     println("MP tumour reached 106.0 mm^3 at t = ", sol_mp.t[end], " days")
# else
#     println("MP tumour did not reach 106.0 mm^3 within 300 days.")
#     println("  max T reached = ", sol_mp.u[end][1])
# end

# if sol_mpb1.t[end] < 300.0 && sol_mpb1.u[end][1] >= T0_mpb1 - 1e-8
#     println("MPB1 tumour reached 62.0 mm^3 at t = ", sol_mpb1.t[end], " days")
# else
#     println("MPB1 tumour did not reach 62.0 mm^3 within 300 days.")
#     println("  max T reached = ", sol_mpb1.u[end][1])
# end


# Shift time data so that the first data point (experimental t=0) is treated as being at t=43 days in the model
# This means shifting all data times forward by the time it takes to reach T0
time_to_T0_mp = sol_mp.t[end]
time_to_T0_mpb1 = sol_mpb1.t[end]

tdata_mp_shifted = tdata_mp .+ time_to_T0_mp
tdata_mpb1_shifted = tdata_mpb1 .+ time_to_T0_mpb1

println("\nTime data shifted:")
println("MP: First data point (experimental t=0) now treated as model t = ", time_to_T0_mp, " days")
println("MPB1: First data point (experimental t=0) now treated as model t = ", time_to_T0_mpb1, " days")


# Initial fit of the tumour + immune cells model to estimate initial values for NK and T cells at treatment initiation time
# With parameters r, p and q to be fitted, and the rest of the parameters fixed based on assumptions and literature values (see Assumptions.jl and Fixed_params.jl)

# Load immune data and proportions via DataTools
immune = DataTools.load_immune_datasets()
df_immune_mp = immune.df_immune_mp
df_immune_mpb1 = immune.df_immune_mpb1
nk_cells_prop_mp = immune.nk_cells_prop_mp
nk_cells_prop_mpb1 = immune.nk_cells_prop_mpb1
t_cells_prop_mp = immune.t_cells_prop_mp
t_cells_prop_mpb1 = immune.t_cells_prop_mpb1

# #
# # Normally in Assumptions, here the code for looking at impact of cell count proportions
# prop_cd45_cells = 0.6  # assuming 20, 40 or 60% infiltration of CD45+ cells
# prop_cd8cells_in_t_cells = 0.4  # assuming 20, 30 or 40% of T cells are CD8+
# #

# Cell counts for NK and T cells, using the assumptions defined in Assumptions.jl
nk_cell_count_mp = nk_cells_prop_mp * prop_cd45_cells * total_cell_count
nk_cell_count_mpb1 = nk_cells_prop_mpb1 * prop_cd45_cells * total_cell_count
t_cell_count_mp = t_cells_prop_mp * prop_cd8cells_in_t_cells * prop_cd45_cells * total_cell_count
t_cell_count_mpb1 = t_cells_prop_mpb1 * prop_cd8cells_in_t_cells * prop_cd45_cells * total_cell_count

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

# Time points for immune cell data (experimental t=10, shifted to model time)
t_immune_mp = 10.0 + time_to_T0_mp
t_immune_mpb1 = 10.0 + time_to_T0_mpb1

# Initial estimates for fitted params
r_0_mp = 1.7      # (1/days), informed on prior fitting runs of tumour-immune cells model
# should be 0.305 from prior fitting of logistic growth to MP vehicle data
r_0_mpb1 = 1.8    # (1/days), informed on prior fitting runs of tumour-immune cells model
# should be 0.423 from prior fitting of logistic growth to MPB1 vehicle data
pN_0 = 0.001      # (1/mm3*day), informed on prior fitting runs of tumour-immune cells model
# was p = 0.1 (1/mm3*day) adapted from de Pillis et al. 2005
q_0 = 0.001       # (1/mm3*day),
# was q = 3.42e-4 (1/mm3*day), adapted from de Pillis et al. 2005

# Bounds
r_min = 0.0
r_max = 5.0
pN_min = 0.00001
pN_max = 1.0
q_min = 0.00001
q_max = 1.0


# Case where we initialise the tumour close to 0
# Initial conditions for immune cells and tumour growth model (vehicle)
# Where T, N, E = u
# with T in mm3, N and E in cells
u0_immune_mp = [tumour_initial_size, (fixed_params_mp.N0 / fixed_params_mp.dN), 0.0]
u0_immune_mpb1 = [tumour_initial_size, (fixed_params_mpb1.N0 / fixed_params_mpb1.dN), 0.0]

include("../src/Predictions.jl")
using .Predictions
#tumour_immune_model_mp_with_immune = (t, p) -> Predictions.predict_tumour_immune(t, t_immune_mp, p, u0_immune_mp, fixed_params_mp)
#tumour_immune_model_mpb1_with_immune = (t, p) -> Predictions.predict_tumour_immune(t, t_immune_mpb1, p, u0_immune_mpb1, fixed_params_mpb1)

# Fit the immune model to combined tumour and immune cell data for MP and MPB1
# fit_immune_mp = curve_fit(tumour_immune_model_mp_with_immune, tdata_mp_shifted, ydata_mp_combined, [r_0_mp, pN_0, q_0], 
#                            lower=[r_min, pN_min, q_min], upper=[r_max, pN_max, q_max])
# fit_immune_mpb1 = curve_fit(tumour_immune_model_mpb1_with_immune, tdata_mpb1_shifted, ydata_mpb1_combined, [r_0_mpb1, pN_0, q_0], 
#                              lower=[r_min, pN_min, q_min], upper=[r_max, pN_max, q_max])

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

fit_immune_mp = curve_fit(tumour_immune_model_mp_with_immune, tdata_mp_shifted, ydata_mp_combined_scaled, [r_0_mp, pN_0, q_0], 
                           lower=[r_min, pN_min, q_min], upper=[r_max, pN_max, q_max])
fit_immune_mpb1 = curve_fit(tumour_immune_model_mpb1_with_immune, tdata_mpb1_shifted, ydata_mpb1_combined_scaled, [r_0_mpb1, pN_0, q_0], 
                             lower=[r_min, pN_min, q_min], upper=[r_max, pN_max, q_max])




# Finding the NK cell and T cell counts at treatment initiation time (model time) based on fitted parameters
t_fine_mp = range(0, maximum(tdata_mp_shifted), length=200)
t_fine_mpb1 = range(0, maximum(tdata_mpb1_shifted), length=200)


predictions_mp = Predictions.predict_tumour_immune_states(t_fine_mp, [fit_immune_mp.param[1], fit_immune_mp.param[2], fit_immune_mp.param[3]], u0_immune_mp, fixed_params_mp)
predictions_mpb1 = Predictions.predict_tumour_immune_states(t_fine_mpb1, [fit_immune_mpb1.param[1], fit_immune_mpb1.param[2], fit_immune_mpb1.param[3]], u0_immune_mpb1, fixed_params_mpb1)

# Extract T, N, E for each dataset
yfit_T_mp = predictions_mp[1, :]
yfit_N_mp = predictions_mp[2, :]
yfit_E_mp = predictions_mp[3, :]

yfit_T_mpb1 = predictions_mpb1[1, :]
yfit_N_mpb1 = predictions_mpb1[2, :]
yfit_E_mpb1 = predictions_mpb1[3, :]

# Compute cell counts at the shift time (first data point in T)
shift_time_mp = time_to_T0_mp
idx_shift_mp = findmin(abs.(t_fine_mp .- shift_time_mp))[2]
shift_N_mp = yfit_N_mp[idx_shift_mp]
shift_E_mp = yfit_E_mp[idx_shift_mp]

shift_time_mpb1 = time_to_T0_mpb1
idx_shift_mpb1 = findmin(abs.(t_fine_mpb1 .- shift_time_mpb1))[2]
shift_N_mpb1 = yfit_N_mpb1[idx_shift_mpb1]
shift_E_mpb1 = yfit_E_mpb1[idx_shift_mpb1]

print("\nEstimated NK cell count at treatment initiation time (model time) for MP: ", shift_N_mp)
print("\nEstimated NK cell count at treatment initiation time (model time) for MPB1: ", shift_N_mpb1)
print("\nEstimated T cell count at treatment initiation time (model time) for MP: ", shift_E_mp)
print("\nEstimated T cell count at treatment initiation time (model time) for MPB1: ", shift_E_mpb1)

# Create DataFrame for parameters and errors
param_df_immune = DataFrame(
    Dataset = ["MP", "MPB1"],
    r = [fit_immune_mp.param[1], fit_immune_mpb1.param[1]],
    pN = [fit_immune_mp.param[2], fit_immune_mpb1.param[2]],
    q = [fit_immune_mp.param[3], fit_immune_mpb1.param[3]],
    r_stderr = [stderror(fit_immune_mp)[1], stderror(fit_immune_mpb1)[1]],
    pN_stderr = [stderror(fit_immune_mp)[2], stderror(fit_immune_mpb1)[2]],
    q_stderr = [stderror(fit_immune_mp)[3], stderror(fit_immune_mpb1)[3]],
    timeshift = [time_to_T0_mp, time_to_T0_mpb1],
    nk_cells_at_treatment_initiation = [shift_N_mp, shift_N_mpb1],
    t_cells_at_shift_at_treatment_initiation = [shift_E_mp, shift_E_mpb1],
    prop_cd45_cells = fill(prop_cd45_cells, length(["MP", "MPB1"])),
    prop_cd8cells_in_t_cells = fill(prop_cd8cells_in_t_cells, length(["MP", "MPB1"]))
)

param_cell_count = DataFrame(
    Dataset = ["MP", "MPB1"],
    nk_cells_at_treatment_initiation = [shift_N_mp, shift_N_mpb1],
    t_cells_at_shift_at_treatment_initiation = [shift_E_mp, shift_E_mpb1],
    prop_cd45_cells = fill(prop_cd45_cells, length(["MP", "MPB1"])),
    prop_cd8cells_in_t_cells = fill(prop_cd8cells_in_t_cells, length(["MP", "MPB1"]))
)

# Write to CSV file (for fitting)
CSV.write("../Fitted_params_results/fitted_parameters_normalised_timeshift.csv", param_df_immune)

# Write to file (for looking at hypothesis on cell counts proportion)
#CSV.write("../Fitted_params_results/si_immune_cell_count_hypothesis/si_fitted_parameters_normalised_timeshift_06_04.csv", param_cell_count)

println("Immune model initial fitting from tumour initiation complete. Parameters saved to fitted_parameters_normalised_timeshift.csv")