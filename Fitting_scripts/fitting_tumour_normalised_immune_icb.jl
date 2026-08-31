using DifferentialEquations
using LsqFit
using CSV, DataFrames
include("../src/Models.jl")
include("../src/Data.jl")
include("../src/Assumptions.jl")
include("../src/Fixed_params.jl")

# Load tumour volume data
data = DataTools.load_tumor_datasets()
tdata_mp_icb = data.mp_icb.t
ydata_mp_icb = data.mp_icb.y
tdata_mpb1_icb = data.mpb1_icb.t
ydata_mpb1_icb = data.mpb1_icb.y

# No immune data

# Definition of the model fixed parameters
fixed_params_mp = fixed_params_mp_icb
fixed_params_mpb1 = fixed_params_mpb1_icb


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

u0_immune_icb_mp = [ydata_mp_icb[1], N0_treatment_initiation_mp, E0_treatment_initiation_mp]
u0_immune_icb_mpb1 = [ydata_mpb1_icb[1], N0_treatment_initiation_mpb1, E0_treatment_initiation_mpb1]

# Initial estimates for fitted params
q_0_mp = 2.00e-3    # (1/mm3*days), from prior fitting of tumour-immune cells model to vehicle
q_0_mpb1 = 7.74e-4  # (1/mm3*days), from prior fitting of tumour-immune cells model to vehicle

# Bounds
q_min = 0.00001
q_max = 1.0

include("../src/Predictions.jl")
using .Predictions
tumour_immune_model_mp_with_immune_icb = (t, p) -> begin
    states = Predictions.predict_tumour_immune_icb_q_states(t, p, u0_immune_icb_mp, fixed_params_mp_icb)
    return states[1, :]  # Only return the tumour volume predictions, as fitting is done to tumour volume data only
end
tumour_immune_model_mpb1_with_immune_icb = (t, p) -> begin
    states = Predictions.predict_tumour_immune_icb_q_states(t, p, u0_immune_icb_mpb1, fixed_params_mpb1_icb)
    return states[1, :]  # Only return the tumour volume predictions, as fitting is done to tumour volume data only
end


fit_immune_mp = curve_fit(tumour_immune_model_mp_with_immune_icb, tdata_mp_icb, ydata_mp_icb, [q_0_mp], 
                           lower=[q_min], upper=[q_max])
fit_immune_mpb1 = curve_fit(tumour_immune_model_mpb1_with_immune_icb, tdata_mpb1_icb, ydata_mpb1_icb, [q_0_mpb1], 
                             lower=[q_min], upper=[q_max])


# Create DataFrame for parameters and errors
param_df_immune_icb = DataFrame(
    Dataset = ["MP", "MPB1"],
    q = [fit_immune_mp.param[1], fit_immune_mpb1.param[1]],
    q_stderr = [stderror(fit_immune_mp)[1], stderror(fit_immune_mpb1)[1]],
)

# Write to CSV file
CSV.write("../Fitted_params_results/fitted_parameters_normalised_immune_icb.csv", param_df_immune_icb)

println("Immune model fitting complete. Parameters saved to fitted_parameters_normalised_immune_icb.csv")