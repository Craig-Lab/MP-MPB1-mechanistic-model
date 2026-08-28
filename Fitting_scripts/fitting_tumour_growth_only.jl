using DifferentialEquations
using LsqFit
using CSV, DataFrames
include("../src/Models.jl")
include("../src/Data.jl")
include("../src/Predictions.jl")
using .Predictions
include("../src/Results.jl")
using .Results

# Importing the data and defining time, tumor volumes, and initial conditions for each dataset
data = DataTools.load_tumor_datasets(data_dir = joinpath(@__DIR__, "..", "Data"))

tdata_mp = data.mp.t
ydata_mp = data.mp.y
tdata_mp_cis = data.mp_cis.t
ydata_mp_cis = data.mp_cis.y
tdata_mp_icb = data.mp_icb.t
ydata_mp_icb = data.mp_icb.y
tdata_mp_cis_icb = data.mp_cis_icb.t
ydata_mp_cis_icb = data.mp_cis_icb.y

tdata_mpb1 = data.mpb1.t
ydata_mpb1 = data.mpb1.y
tdata_mpb1_cis = data.mpb1_cis.t
ydata_mpb1_cis = data.mpb1_cis.y
tdata_mpb1_icb = data.mpb1_icb.t
ydata_mpb1_icb = data.mpb1_icb.y
tdata_mpb1_cis_icb = data.mpb1_cis_icb.t
ydata_mpb1_cis_icb = data.mpb1_cis_icb.y

u0_mp = data.mp.u0
u0_mpb1 = data.mpb1.u0
u0_mp_cis = data.mp_cis.u0
u0_mpb1_cis = data.mpb1_cis.u0
u0_mp_icb = data.mp_icb.u0
u0_mpb1_icb = data.mpb1_icb.u0
u0_mp_cis_icb = data.mp_cis_icb.u0
u0_mpb1_cis_icb = data.mpb1_cis_icb.u0


tumour_model_exponential_mp = (t, p) -> Predictions.predict_exponential(t, p, u0_mp)
tumour_model_exponential_mpb1 = (t, p) -> Predictions.predict_exponential(t, p, u0_mpb1)

tumour_model_logistic_mp = (t, p) -> Predictions.predict_logistic(t, p, u0_mp)
tumour_model_logistic_mp_cis = (t, p) -> Predictions.predict_logistic(t, p, u0_mp_cis)
tumour_model_logistic_mpb1 = (t, p) -> Predictions.predict_logistic(t, p, u0_mpb1)
tumour_model_logistic_mpb1_cis = (t, p) -> Predictions.predict_logistic(t, p, u0_mpb1_cis)
tumour_model_logistic_mp_icb = (t, p) -> Predictions.predict_logistic(t, p, u0_mp_icb)
tumour_model_logistic_mpb1_icb = (t, p) -> Predictions.predict_logistic(t, p, u0_mpb1_icb)
tumour_model_logistic_mp_cis_icb = (t, p) -> Predictions.predict_logistic(t, p, u0_mp_cis_icb)
tumour_model_logistic_mpb1_cis_icb = (t, p) -> Predictions.predict_logistic(t, p, u0_mpb1_cis_icb)

tumour_model_gompertz_mp = (t, p) -> Predictions.predict_gompertz(t, p, u0_mp)
tumour_model_gompertz_mpb1 = (t, p) -> Predictions.predict_gompertz(t, p, u0_mpb1)


# Initial estimates for [r, k]
r0 = 0.514  # (1/days), Tumor growth rate from de Pillis
k0 = 2000.0  # (mm^3), Tumor carrying capacity (could be more fine-tuned)

rmin = 0.0
rmax = 1.0
kmin = 0.0
kmax = 20000.0

# Fitting the models to vehicle data for MP and MPB1
fit_exponential_mp = curve_fit(tumour_model_exponential_mp, tdata_mp, ydata_mp, [r0], lower=[rmin], upper=[rmax])
fit_exponential_mpb1 = curve_fit(tumour_model_exponential_mpb1, tdata_mpb1, ydata_mpb1, [r0], lower=[rmin], upper=[rmax])
fit_logistic_mp = curve_fit(tumour_model_logistic_mp, tdata_mp, ydata_mp, [r0, k0], lower=[rmin, kmin], upper=[rmax, kmax])
fit_logistic_mpb1 = curve_fit(tumour_model_logistic_mpb1, tdata_mpb1, ydata_mpb1, [r0, k0], lower=[rmin, kmin], upper=[rmax, kmax])
fit_gompertz_mp = curve_fit(tumour_model_gompertz_mp, tdata_mp, ydata_mp, [r0, k0], lower=[rmin, kmin], upper=[rmax, kmax])
fit_gompertz_mpb1 = curve_fit(tumour_model_gompertz_mpb1, tdata_mpb1, ydata_mpb1, [r0, k0], lower=[rmin, kmin], upper=[rmax, kmax])

# Fitting logistic model to cisplatin-treated, ICB, and Cis+ICB data for MP and MPB1
fit_logistic_mp_cis = curve_fit(tumour_model_logistic_mp_cis, tdata_mp_cis, ydata_mp_cis, [r0, k0], lower=[rmin, kmin], upper=[rmax, kmax])
fit_logistic_mpb1_cis = curve_fit(tumour_model_logistic_mpb1_cis, tdata_mpb1_cis, ydata_mpb1_cis, [r0, k0], lower=[rmin, kmin], upper=[rmax, kmax])
fit_logistic_mp_icb = curve_fit(tumour_model_logistic_mp_icb, tdata_mp_icb, ydata_mp_icb, [r0, k0], lower=[rmin, kmin], upper=[rmax, kmax])
fit_logistic_mpb1_icb = curve_fit(tumour_model_logistic_mpb1_icb, tdata_mpb1_icb, ydata_mpb1_icb, [r0, k0], lower=[rmin, kmin], upper=[rmax, kmax])
fit_logistic_mp_cis_icb = curve_fit(tumour_model_logistic_mp_cis_icb, tdata_mp_cis_icb, ydata_mp_cis_icb, [r0, k0], lower=[rmin, kmin], upper=[rmax, kmax])
fit_logistic_mpb1_cis_icb = curve_fit(tumour_model_logistic_mpb1_cis_icb, tdata_mpb1_cis_icb, ydata_mpb1_cis_icb, [r0, k0], lower=[rmin, kmin], upper=[rmax, kmax])


# Saving results in DataFrames and writing to CSV files
datasets = ["MP", "MPB1"]
param_df_exponential = Results.build_param_df(datasets, [fit_exponential_mp, fit_exponential_mpb1], [:r], [1])
param_df_logistic = Results.build_param_df(datasets, [fit_logistic_mp, fit_logistic_mpb1], [:r, :k], [1, 2])
param_df_gompertz = Results.build_param_df(datasets, [fit_gompertz_mp, fit_gompertz_mpb1], [:r, :k], [1, 2])

param_df_logistic_cis = Results.build_param_df(["MP Cisplatin", "MPB1 Cisplatin"], [fit_logistic_mp_cis, fit_logistic_mpb1_cis], [:r, :k], [1, 2])
param_df_logistic_icb = Results.build_param_df(["MP ICB", "MPB1 ICB"], [fit_logistic_mp_icb, fit_logistic_mpb1_icb], [:r, :k], [1, 2])
param_df_logistic_cis_icb = Results.build_param_df(["MP Cisplatin + ICB", "MPB1 Cisplatin + ICB"], [fit_logistic_mp_cis_icb, fit_logistic_mpb1_cis_icb], [:r, :k], [1, 2])

# Write to CSV files
CSV.write("./Fitted_params_results/fitted_parameters_exponential.csv", param_df_exponential)
CSV.write("./Fitted_params_results/fitted_parameters_logistic.csv", param_df_logistic)
CSV.write("./Fitted_params_results/fitted_parameters_gompertz.csv", param_df_gompertz)
CSV.write("./Fitted_params_results/fitted_parameters_logistic_cisplatin.csv", param_df_logistic_cis)
CSV.write("./Fitted_params_results/fitted_parameters_logistic_icb.csv", param_df_logistic_icb)
CSV.write("./Fitted_params_results/fitted_parameters_logistic_cis_icb.csv", param_df_logistic_cis_icb)