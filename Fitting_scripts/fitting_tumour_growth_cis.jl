using DifferentialEquations
using LsqFit
using CSV, DataFrames
include("../src/Models.jl")
include("../src/Data.jl")
include("../src/Fixed_params.jl")
include("../src/Predictions.jl")
using .Predictions

# Importing the data and defining time, tumor volumes, and initial conditions for each dataset
data = DataTools.load_tumor_datasets(data_dir = joinpath(@__DIR__, "..", "Data"))

tdata_mp_cis = data.mp_cis.t
ydata_mp_cis = data.mp_cis.y
tdata_mpb1_cis = data.mpb1_cis.t
ydata_mpb1_cis = data.mpb1_cis.y

# Fixed bolus dosing setup
BW = 0.025
dose_amt = 3.0 * BW
initial_day_dose = 0.01  # Start dosing at 0.01 days to avoid dosing at t=0 which is causing issues
dose_starts = [initial_day_dose, initial_day_dose + 7.0]  # Doses 1 week apart

# Initial conditions for PK/PD model
# With the variables: C1, C2, T1, T2, T3, T4
# Initial condition is 0 for all except T1
u0_pkpd_mp_cis = [0.0, 0.0, ydata_mp_cis[1], 0.0, 0.0, 0.0]
u0_pkpd_mpb1_cis = [0.0, 0.0, ydata_mpb1_cis[1], 0.0, 0.0, 0.0]

# Use fixed parameters defined in src/Fixed_params.jl
fixed_params_mp = fixed_params_cis_mp
fixed_params_mpb1 = fixed_params_cis_mpb1

# PK/PD prediction wrappers using the centralized Predictions helper
tumour_cis_pkpd_mp = (t, p) -> Predictions.predict_cis_pkpd(t, p, u0_pkpd_mp_cis, fixed_params_mp, dose_amt, dose_starts)
tumour_cis_pkpd_mpb1 = (t, p) -> Predictions.predict_cis_pkpd(t, p, u0_pkpd_mpb1_cis, fixed_params_mpb1, dose_amt, dose_starts)

# Initial estimate
k2_0 = 15.0

# Bounds
k2min = 0.0
k2max = 150.0

# Fit PK/PD cisplatin bolus model
fit_cis_pkpd_mp = curve_fit(
    tumour_cis_pkpd_mp,
    tdata_mp_cis,
    ydata_mp_cis,
    [k2_0],
    lower = [k2min],
    upper = [k2max]
)

fit_cis_pkpd_mpb1 = curve_fit(
    tumour_cis_pkpd_mpb1,
    tdata_mpb1_cis,
    ydata_mpb1_cis,
    [k2_0],
    lower = [k2min],
    upper = [k2max]
)

param_df_pkpd_cis = DataFrame(
    Dataset = ["MP Cisplatin", "MPB1 Cisplatin"],
    k2 = [
        fit_cis_pkpd_mp.param[1],
        fit_cis_pkpd_mpb1.param[1]
    ],
    k2_stderr = [
        stderror(fit_cis_pkpd_mp)[1],
        stderror(fit_cis_pkpd_mpb1)[1]
    ],
)

CSV.write(
    "./Fitted_params_results/fitted_parameters_pkpd_cisplatin_bolus.csv",
    param_df_pkpd_cis
)