using DifferentialEquations
using LsqFit, StatsBase, Distributions
using CSV, DataFrames, JLD2, PlotlyJS

include("../src/craig_lab_template.jl");
include("FigureFunctions.jl");

include("../src/Models.jl");
include("../src/Data.jl");
include("../src/Predictions.jl");
include("../src/Results.jl");
include("../src/Assumptions.jl");
include("../src/Fixed_params.jl");

using .Predictions
using .Results

vehicle_colour = "#c099ba";
cis_colour = "#974855";
icb_colour = "#b8dbe5";
cis_icb_colour =  "#1b325b";

random_color1 = "#d17466";
random_color2 = "#a781a2";
random_color3 = "#9fc2cc";

## -------------------------------------

# Importing the data and defining time, tumor volumes, and initial conditions for each dataset
data = DataTools.load_tumor_datasets(data_dir = joinpath(@__DIR__, "..", "Data"));

tdata_mp = data.mp.t;
ydata_mp = data.mp.y;
yerr_mp = data.mp.yerr;
tdata_mp_cis = data.mp_cis.t;
ydata_mp_cis = data.mp_cis.y;
yerr_mp_cis = data.mp_cis.yerr;
tdata_mp_icb = data.mp_icb.t;
ydata_mp_icb = data.mp_icb.y;
yerr_mp_icb = data.mp_icb.yerr;
tdata_mp_cis_icb = data.mp_cis_icb.t;
ydata_mp_cis_icb = data.mp_cis_icb.y;
yerr_mp_cis_icb = data.mp_cis_icb.yerr;

tdata_mpb1 = data.mpb1.t;
ydata_mpb1 = data.mpb1.y;
yerr_mpb1 = data.mpb1.yerr;
tdata_mpb1_cis = data.mpb1_cis.t;
ydata_mpb1_cis = data.mpb1_cis.y;
yerr_mpb1_cis = data.mpb1_cis.yerr;
tdata_mpb1_icb = data.mpb1_icb.t;
ydata_mpb1_icb = data.mpb1_icb.y;
yerr_mpb1_icb = data.mpb1_icb.yerr;
tdata_mpb1_cis_icb = data.mpb1_cis_icb.t;
ydata_mpb1_cis_icb = data.mpb1_cis_icb.y;
yerr_mpb1_cis_icb = data.mpb1_cis_icb.yerr;

u0_mp   = data.mp.u0;
u0_mpb1 = data.mpb1.u0;
u0_mp_cis   = data.mp_cis.u0;
u0_mpb1_cis = data.mpb1_cis.u0;
u0_mp_icb   = data.mp_icb.u0;
u0_mpb1_icb = data.mpb1_icb.u0;
u0_mp_cis_icb   = data.mp_cis_icb.u0;
u0_mpb1_cis_icb = data.mpb1_cis_icb.u0;

## -------------------------------------
# Tumour growth only

# Define tumor growth models
logistic_mp           = (t, p) -> Predictions.predict_logistic(t, p, u0_mp);
logistic_mp_cis       = (t, p) -> Predictions.predict_logistic(t, p, u0_mp_cis);
logistic_mpb1         = (t, p) -> Predictions.predict_logistic(t, p, u0_mpb1);
logistic_mpb1_cis     = (t, p) -> Predictions.predict_logistic(t, p, u0_mpb1_cis);
logistic_mp_icb       = (t, p) -> Predictions.predict_logistic(t, p, u0_mp_icb);
logistic_mpb1_icb     = (t, p) -> Predictions.predict_logistic(t, p, u0_mpb1_icb);
logistic_mp_cis_icb   = (t, p) -> Predictions.predict_logistic(t, p, u0_mp_cis_icb);
logistic_mpb1_cis_icb = (t, p) -> Predictions.predict_logistic(t, p, u0_mpb1_cis_icb);

exponential_mp   = (t, p) -> Predictions.predict_exponential(t, p, u0_mp);
exponential_mpb1 = (t, p) -> Predictions.predict_exponential(t, p, u0_mpb1);
gompertz_mp      = (t, p) -> Predictions.predict_gompertz(t, p, u0_mp);
gompertz_mpb1    = (t, p) -> Predictions.predict_gompertz(t, p, u0_mpb1);

# Load fitted parameter results
param_exp = CSV.read(joinpath(@__DIR__, "..", "Fitted_params_results/fitted_parameters_exponential.csv"), DataFrame);
param_log = CSV.read(joinpath(@__DIR__, "..", "Fitted_params_results/fitted_parameters_logistic.csv"), DataFrame);
param_gom = CSV.read(joinpath(@__DIR__, "..", "Fitted_params_results/fitted_parameters_gompertz.csv"), DataFrame);
param_log_cis = CSV.read(joinpath(@__DIR__, "..", "Fitted_params_results/fitted_parameters_logistic_cisplatin.csv"), DataFrame);
param_log_icb = CSV.read(joinpath(@__DIR__, "..", "Fitted_params_results/fitted_parameters_logistic_icb.csv"), DataFrame);
param_log_cis_icb = CSV.read(joinpath(@__DIR__, "..", "Fitted_params_results/fitted_parameters_logistic_cis_icb.csv"), DataFrame);

# Extract parameters (tumour growth)
r_exp_mp = param_exp[param_exp.Dataset .== "MP", :r][1];
r_exp_mpb1 = param_exp[param_exp.Dataset .== "MPB1", :r][1];

r_log_mp = param_log[param_log.Dataset .== "MP", :r][1];
k_log_mp = param_log[param_log.Dataset .== "MP", :k][1];
r_log_mpb1 = param_log[param_log.Dataset .== "MPB1", :r][1];
k_log_mpb1 = param_log[param_log.Dataset .== "MPB1", :k][1];

r_gom_mp = param_gom[param_gom.Dataset .== "MP", :r][1];
k_gom_mp = param_gom[param_gom.Dataset .== "MP", :k][1];
r_gom_mpb1 = param_gom[param_gom.Dataset .== "MPB1", :r][1];
k_gom_mpb1 = param_gom[param_gom.Dataset .== "MPB1", :k][1];

r_log_mp_cis = param_log_cis[param_log_cis.Dataset .== "MP Cisplatin", :r][1];
k_log_mp_cis = param_log_cis[param_log_cis.Dataset .== "MP Cisplatin", :k][1];
r_log_mpb1_cis = param_log_cis[param_log_cis.Dataset .== "MPB1 Cisplatin", :r][1];
k_log_mpb1_cis = param_log_cis[param_log_cis.Dataset .== "MPB1 Cisplatin", :k][1];

r_log_mp_icb = param_log_icb[param_log_icb.Dataset .== "MP ICB", :r][1];
k_log_mp_icb = param_log_icb[param_log_icb.Dataset .== "MP ICB", :k][1];
r_log_mpb1_icb = param_log_icb[param_log_icb.Dataset .== "MPB1 ICB", :r][1];
k_log_mpb1_icb = param_log_icb[param_log_icb.Dataset .== "MPB1 ICB", :k][1];

r_log_mp_cis_icb = param_log_cis_icb[param_log_cis_icb.Dataset .== "MP Cisplatin + ICB", :r][1];
k_log_mp_cis_icb = param_log_cis_icb[param_log_cis_icb.Dataset .== "MP Cisplatin + ICB", :k][1];
r_log_mpb1_cis_icb = param_log_cis_icb[param_log_cis_icb.Dataset .== "MPB1 Cisplatin + ICB", :r][1];
k_log_mpb1_cis_icb = param_log_cis_icb[param_log_cis_icb.Dataset .== "MPB1 Cisplatin + ICB", :k][1];


# Generate predictions on a finer time grid for smooth curves
t_fine_mp   = range(tdata_mp[1], tdata_mp[end], length=200);
t_fine_mpb1 = range(tdata_mpb1[1], tdata_mpb1[end], length=200);

yfit_exp_mp   = exponential_mp(t_fine_mp, [r_exp_mp]);
yfit_exp_mpb1 = exponential_mpb1(t_fine_mpb1, [r_exp_mpb1]);

yfit_gom_mp   = gompertz_mp(t_fine_mp, [r_gom_mp, k_gom_mp]);
yfit_gom_mpb1 = gompertz_mpb1(t_fine_mpb1, [r_gom_mpb1, k_gom_mpb1]);

yfit_log_mp   = logistic_mp(t_fine_mp, [r_log_mp, k_log_mp]);
yfit_log_mpb1 = logistic_mpb1(t_fine_mpb1, [r_log_mpb1, k_log_mpb1]);

yfit_log_mp_cis   = logistic_mp_cis(t_fine_mp, [r_log_mp_cis, k_log_mp_cis]);
yfit_log_mpb1_cis = logistic_mpb1_cis(t_fine_mpb1, [r_log_mpb1_cis, k_log_mpb1_cis]);

yfit_log_mp_icb   = logistic_mp_icb(t_fine_mp, [r_log_mp_icb, k_log_mp_icb]);
yfit_log_mpb1_icb = logistic_mpb1_icb(t_fine_mpb1, [r_log_mpb1_icb, k_log_mpb1_icb]);

yfit_log_mp_cis_icb   = logistic_mp_cis_icb(t_fine_mp, [r_log_mp_cis_icb, k_log_mp_cis_icb]);
yfit_log_mpb1_cis_icb = logistic_mpb1_cis_icb(t_fine_mpb1, [r_log_mpb1_cis_icb, k_log_mpb1_cis_icb]);

# ====================================
# Figure S1
# ====================================

plotall_growth_only([tdata_mp, ydata_mp], [t_fine_mp, yfit_exp_mp, yfit_log_mp, yfit_gom_mp], yerr_mp, "All Models - MP Vehicle", [random_color1, random_color2, random_color3]);
plotall_growth_only([tdata_mpb1, ydata_mpb1], [t_fine_mpb1, yfit_exp_mpb1, yfit_log_mpb1, yfit_gom_mpb1], yerr_mpb1, "All Models - MPB1 Vehicle", [random_color1, random_color2, random_color3]; seelegend = false);

## -------------------------------------
# Logistic growth with treatment 

# ====================================
# Figure 2
# ====================================

fig_color = [vehicle_colour, cis_colour, icb_colour, cis_icb_colour];

fig_mp_data = [tdata_mp, ydata_mp, ydata_mp_cis, ydata_mp_icb, ydata_mp_cis_icb];
fig_mp_sol  = [t_fine_mp, yfit_log_mp, yfit_log_mp_cis, yfit_log_mp_icb, yfit_log_mp_cis_icb];
fig_mp_err  = [yerr_mp, yerr_mp_cis, yerr_mp_icb, yerr_mp_cis_icb];

fig_mpb1_data = [tdata_mpb1, ydata_mpb1, ydata_mpb1_cis, ydata_mpb1_icb, ydata_mpb1_cis_icb];
fig_mpb1_sol  = [t_fine_mp, yfit_log_mpb1, yfit_log_mpb1_cis, yfit_log_mpb1_icb, yfit_log_mpb1_cis_icb];
fig_mpb1_err  = [yerr_mpb1, yerr_mpb1_cis, yerr_mpb1_icb, yerr_mpb1_cis_icb];

plot_log_treatments(fig_mp_data, fig_mp_sol, fig_mp_err, "Logistic Model Across Treatment - MP", fig_color);
plot_log_treatments(fig_mpb1_data, fig_mpb1_sol, fig_mpb1_err, "Logistic Model Across Treatment - MPB1", fig_color; seelegend = false);

## -------------------------------------
# Logistic tumour growth and cisplatin PKPD

# Import of params and functions
param_pkpd_cis = CSV.read(joinpath(@__DIR__, "..", "Fitted_params_results/fitted_parameters_pkpd_cisplatin_bolus.csv"), DataFrame);
k2_mp_cis = param_pkpd_cis[param_pkpd_cis.Dataset .== "MP Cisplatin", :k2][1];
k2_mpb1_cis = param_pkpd_cis[param_pkpd_cis.Dataset .== "MPB1 Cisplatin", :k2][1];

# Use fixed parameters defined in src/Fixed_params.jl
fixed_params_mp = fixed_params_cis_mp;
fixed_params_mpb1 = fixed_params_cis_mpb1;

# Fixed bolus dosing setup
BW = 0.025;
dose_amt = 3.0 * BW;
initial_day_dose = 0.01;
dose_starts = [initial_day_dose, initial_day_dose + 7.0];

# Initial conditions for PK/PD model
u0_pkpd_mp_cis = [0.0, 0.0, ydata_mp_cis[1], 0.0, 0.0, 0.0];
u0_pkpd_mpb1_cis = [0.0, 0.0, ydata_mpb1_cis[1], 0.0, 0.0, 0.0];

yfit_cis_pkpd_mp = Predictions.predict_cis_pkpd(t_fine_mp, [k2_mp_cis], u0_pkpd_mp_cis, fixed_params_mp, dose_amt, dose_starts);
yfit_cis_pkpd_mpb1 = Predictions.predict_cis_pkpd(t_fine_mpb1, [k2_mpb1_cis], u0_pkpd_mpb1_cis, fixed_params_mpb1, dose_amt, dose_starts);

# ====================================
# Figure 3
# ====================================

fig_color = [vehicle_colour, cis_colour];

fig_mp_data = [tdata_mp, ydata_mp, ydata_mp_cis];
fig_mp_sol  = [t_fine_mp, yfit_log_mp, yfit_cis_pkpd_mp];
fig_mp_err  = [yerr_mp, yerr_mp_cis];

fig_mpb1_data = [tdata_mpb1, ydata_mpb1, ydata_mpb1_cis];
fig_mpb1_sol  = [t_fine_mp, yfit_log_mpb1, yfit_cis_pkpd_mpb1];
fig_mpb1_err  = [yerr_mpb1, yerr_mpb1_cis];

plot_pkpd_cis(fig_mp_data, fig_mp_sol, fig_mp_err, "Logistic vs Cisplatin PKPD - MP", fig_color);
plot_pkpd_cis(fig_mpb1_data, fig_mpb1_sol, fig_mpb1_err, "Logistic vs Cisplatin PKPD - MPB1", fig_color; seelegend = false);

## -------------------------------------
# Timeshift for vehicle tumour + immune cells initial conditions

# Load immune data and proportions via DataTools
immune = DataTools.load_immune_datasets(data_dir = joinpath(@__DIR__, "..", "Data"));

df_immune_mp = immune.df_immune_mp;
df_immune_mpb1 = immune.df_immune_mpb1;
nk_cells_prop_mp = immune.nk_cells_prop_mp;
nk_cells_prop_mpb1 = immune.nk_cells_prop_mpb1;
t_cells_prop_mp = immune.t_cells_prop_mp;
t_cells_prop_mpb1 = immune.t_cells_prop_mpb1;

# Error bars for immune cells
nk_cells_prop_mp_sem = immune.nk_cells_prop_mp_sem;
nk_cells_prop_mpb1_sem = immune.nk_cells_prop_mpb1_sem;
t_cells_prop_mp_sem = immune.t_cells_prop_mp_sem;
t_cells_prop_mpb1_sem = immune.t_cells_prop_mpb1_sem;

# Cell counts for NK and T cells
nk_cell_count_mp = nk_cells_prop_mp * prop_cd45_cells * total_cell_count;
nk_cell_count_mpb1 = nk_cells_prop_mpb1 * prop_cd45_cells * total_cell_count;
t_cell_count_mp = t_cells_prop_mp * prop_cd8cells_in_t_cells * prop_cd45_cells * total_cell_count;
t_cell_count_mpb1 = t_cells_prop_mpb1 * prop_cd8cells_in_t_cells * prop_cd45_cells * total_cell_count;

nk_cell_count_mp_sem = nk_cells_prop_mp_sem * prop_cd45_cells * total_cell_count;
nk_cell_count_mpb1_sem = nk_cells_prop_mpb1_sem * prop_cd45_cells * total_cell_count;
t_cell_count_mp_sem = t_cells_prop_mp_sem * prop_cd8cells_in_t_cells * prop_cd45_cells * total_cell_count;
t_cell_count_mpb1_sem = t_cells_prop_mpb1_sem * prop_cd8cells_in_t_cells * prop_cd45_cells * total_cell_count;

# Load fitted parameters
param_immuneTS = CSV.read(joinpath(@__DIR__, "..", "Fitted_params_results", "fitted_parameters_normalised_timeshift.csv"), DataFrame);

# Extract immune model parameters for MP and MPB1
r_immune_mp = param_immuneTS[param_immuneTS.Dataset .== "MP", :r][1];
pN_immune_mp = param_immuneTS[param_immuneTS.Dataset .== "MP", :pN][1];
q_immune_mp = param_immuneTS[param_immuneTS.Dataset .== "MP", :q][1];

r_immune_mpb1 = param_immuneTS[param_immuneTS.Dataset .== "MPB1", :r][1];
pN_immune_mpb1 = param_immuneTS[param_immuneTS.Dataset .== "MPB1", :pN][1];
q_immune_mpb1 = param_immuneTS[param_immuneTS.Dataset .== "MPB1", :q][1];

fixed_params_mp = fixed_params_immune_mp;
fixed_params_mpb1 = fixed_params_immune_mpb1;

# Initial conditions for immune cells and tumour growth model (vehicle)
u0_immune_mp = [tumour_initial_size, (fixed_params_mp.N0 / fixed_params_mp.dN), 0.0];
u0_immune_mpb1 = [tumour_initial_size, (fixed_params_mpb1.N0 / fixed_params_mpb1.dN), 0.0];

time_to_T0_mp = param_immuneTS[param_immuneTS.Dataset .== "MP", :timeshift][1];
time_to_T0_mpb1 = param_immuneTS[param_immuneTS.Dataset .== "MPB1", :timeshift][1];

tdata_mp_shifted = tdata_mp .+ time_to_T0_mp;
tdata_mpb1_shifted = tdata_mpb1 .+ time_to_T0_mpb1;

t_fine_mp = range(0, maximum(tdata_mp_shifted), length=200);
t_fine_mpb1 = range(0, maximum(tdata_mpb1_shifted), length=200);

predictions_mp = Predictions.predict_tumour_immune_states(t_fine_mp, [r_immune_mp, pN_immune_mp, q_immune_mp], u0_immune_mp, fixed_params_mp);
predictions_mpb1 = Predictions.predict_tumour_immune_states(t_fine_mpb1, [r_immune_mpb1, pN_immune_mpb1, q_immune_mpb1], u0_immune_mpb1, fixed_params_mpb1);

# Extract T, N, E for each dataset
yfit_T_mp = predictions_mp[1, :];
yfit_N_mp = predictions_mp[2, :];
yfit_E_mp = predictions_mp[3, :];

yfit_T_mpb1 = predictions_mpb1[1, :];
yfit_N_mpb1 = predictions_mpb1[2, :];
yfit_E_mpb1 = predictions_mpb1[3, :];

# ====================================
# Figure S2
# ====================================

mp_data = [tdata_mp_shifted, ydata_mp, t_cell_count_mp, nk_cell_count_mp];
mp_error = [yerr_mp, t_cell_count_mp_sem, nk_cell_count_mp_sem];
mp_sol = [t_fine_mp, yfit_T_mp, yfit_E_mp, yfit_N_mp];

plot_immune(mp_data, mp_sol, mp_error, "MP (Time shift)", vehicle_colour; seelegend = true, t2treatment = time_to_T0_mp);

mpb1_data = [tdata_mpb1_shifted, ydata_mpb1, t_cell_count_mpb1, nk_cell_count_mpb1];
mpb1_error = [yerr_mpb1, t_cell_count_mpb1_sem, nk_cell_count_mpb1_sem];
mpb1_sol = [t_fine_mpb1, yfit_T_mpb1, yfit_E_mpb1, yfit_N_mpb1];

plot_immune(mpb1_data, mpb1_sol, mpb1_error, "MPB1 (Time shift)", vehicle_colour; seelegend = false, t2treatment = time_to_T0_mpb1);

# ====================================
# Figure 4
# ====================================

# Load fitted parameters
param_immune = CSV.read(joinpath(@__DIR__, "..", "Fitted_params_results", "fitted_parameters_normalised_immune.csv"), DataFrame);

# Extract immune model parameters for MP and MPB1
r_immune_mp = param_immune[param_immune.Dataset .== "MP", :r][1];
pN_immune_mp = param_immune[param_immune.Dataset .== "MP", :pN][1];
q_immune_mp = param_immune[param_immune.Dataset .== "MP", :q][1];

r_immune_mpb1 = param_immune[param_immune.Dataset .== "MPB1", :r][1];
pN_immune_mpb1 = param_immune[param_immune.Dataset .== "MPB1", :pN][1];
q_immune_mpb1 = param_immune[param_immune.Dataset .== "MPB1", :q][1];

# Load fitted parameters from timeshift fitting (will be initial estimates for fitting the immune model)
time_shift_result = CSV.read(joinpath(@__DIR__, "..", "Fitted_params_results", "fitted_parameters_normalised_timeshift.csv"), DataFrame);
N0_treatment_initiation_mp = time_shift_result[time_shift_result.Dataset .== "MP", :nk_cells_at_treatment_initiation][1];
N0_treatment_initiation_mpb1 = time_shift_result[time_shift_result.Dataset .== "MPB1", :nk_cells_at_treatment_initiation][1];
E0_treatment_initiation_mp = time_shift_result[time_shift_result.Dataset .== "MP", :t_cells_at_shift_at_treatment_initiation][1];
E0_treatment_initiation_mpb1 = time_shift_result[time_shift_result.Dataset .== "MPB1", :t_cells_at_shift_at_treatment_initiation][1];

t_immune_mp = 10.0;
t_immune_mpb1 = 10.0;
u0_immune_mp = [ydata_mp[1], N0_treatment_initiation_mp, E0_treatment_initiation_mp];
u0_immune_mpb1 = [ydata_mpb1[1], N0_treatment_initiation_mpb1, E0_treatment_initiation_mpb1];

t_fine_mp = range(0, maximum(tdata_mp), length=200);
t_fine_mpb1 = range(0, maximum(tdata_mpb1), length=200);
predictions_mp = Predictions.predict_tumour_immune_states(t_fine_mp, [r_immune_mp, pN_immune_mp, q_immune_mp], u0_immune_mp, fixed_params_mp);
predictions_mpb1 = Predictions.predict_tumour_immune_states(t_fine_mpb1, [r_immune_mpb1, pN_immune_mpb1, q_immune_mpb1], u0_immune_mpb1, fixed_params_mpb1);

# Extract T, N, E for each dataset
yfit_T_mp = predictions_mp[1, :];
yfit_N_mp = predictions_mp[2, :];
yfit_E_mp = predictions_mp[3, :];

yfit_T_mpb1 = predictions_mpb1[1, :];
yfit_N_mpb1 = predictions_mpb1[2, :];
yfit_E_mpb1 = predictions_mpb1[3, :];

# plot figures
mp_data = [tdata_mp, ydata_mp, t_cell_count_mp, nk_cell_count_mp];
mp_error = [yerr_mp, t_cell_count_mp_sem, nk_cell_count_mp_sem];
mp_sol = [t_fine_mp, yfit_T_mp, yfit_E_mp, yfit_N_mp];

plot_immune(mp_data, mp_sol, mp_error, "MP (Vehicle)", vehicle_colour; seelegend = true);

mpb1_data = [tdata_mpb1, ydata_mpb1, t_cell_count_mpb1, nk_cell_count_mpb1];
mpb1_error = [yerr_mpb1, t_cell_count_mpb1_sem, nk_cell_count_mpb1_sem];
mpb1_sol = [t_fine_mpb1, yfit_T_mpb1, yfit_E_mpb1, yfit_N_mpb1];

plot_immune(mpb1_data, mpb1_sol, mpb1_error, "MPB1 (Vehicle)", vehicle_colour; seelegend = false);

## -------------------------------------
# Tumour + immune cells with cisplatin PKPD model 

# Load immune data and proportions
immune_cis = DataTools.load_immune_datasets_cis(data_dir = joinpath(@__DIR__, "..", "Data"));

df_immune_cis_mp = immune_cis.df_immune_cis_mp;
df_immune_cis_mpb1 = immune_cis.df_immune_cis_mpb1;
nk_cells_prop_cis_mp = immune_cis.nk_cells_prop_cis_mp;
nk_cells_prop_cis_mpb1 = immune_cis.nk_cells_prop_cis_mpb1;
t_cells_prop_cis_mp = immune_cis.t_cells_prop_cis_mp;
t_cells_prop_cis_mpb1 = immune_cis.t_cells_prop_cis_mpb1;

# Error bars for immune cells
nk_cells_prop_cis_mp_sem = immune_cis.nk_cells_prop_cis_mp_sem;
nk_cells_prop_cis_mpb1_sem = immune_cis.nk_cells_prop_cis_mpb1_sem;
t_cells_prop_cis_mp_sem = immune_cis.t_cells_prop_cis_mp_sem;
t_cells_prop_cis_mpb1_sem = immune_cis.t_cells_prop_cis_mpb1_sem;

# Fixed bolus dosing setup
BW = 0.025;
dose_amt = 3.0 * BW;
initial_day_dose = 0.01;
dose_starts = [initial_day_dose, initial_day_dose + 7.0];  # Doses 1 week apart

# Load fitted parameters
param_immune_cis = CSV.read(joinpath(@__DIR__, "..", "Fitted_params_results", "fitted_parameters_normalised_immune_cisplatin.csv"), DataFrame);

# Cell counts for NK and T cells, using the assumptions defined in Assumptions.jl
nk_cell_count_mp_cis = nk_cells_prop_cis_mp * prop_cd45_cells * total_cell_count;
nk_cell_count_mpb1_cis = nk_cells_prop_cis_mpb1 * prop_cd45_cells * total_cell_count;
t_cell_count_mp_cis = t_cells_prop_cis_mp * prop_cd8cells_in_t_cells * prop_cd45_cells * total_cell_count;
t_cell_count_mpb1_cis = t_cells_prop_cis_mpb1 * prop_cd8cells_in_t_cells * prop_cd45_cells * total_cell_count;

# Error bars for immune cells
nk_cell_count_mp_cis_sem = nk_cells_prop_cis_mp_sem * prop_cd45_cells * total_cell_count;
nk_cell_count_mpb1_cis_sem = nk_cells_prop_cis_mpb1_sem * prop_cd45_cells * total_cell_count;
t_cell_count_mp_cis_sem = t_cells_prop_cis_mp_sem * prop_cd8cells_in_t_cells * prop_cd45_cells * total_cell_count;
t_cell_count_mpb1_cis_sem = t_cells_prop_cis_mpb1_sem * prop_cd8cells_in_t_cells * prop_cd45_cells * total_cell_count;

# Extract immune model parameters for MP and MPB1
#k2_immune_mp = param_immune_cis[param_immune_cis.Dataset .== "MP", :k2][1];
k2_immune_mp = 8.07
k2N_immune_mp = param_immune_cis[param_immune_cis.Dataset .== "MP", :k2N][1];
k2E_immune_mp = param_immune_cis[param_immune_cis.Dataset .== "MP", :k2E][1];
#k2_immune_mpb1 = param_immune_cis[param_immune_cis.Dataset .== "MPB1", :k2][1];
k2_immune_mpb1 = 63.7
k2N_immune_mpb1 = param_immune_cis[param_immune_cis.Dataset .== "MPB1", :k2N][1];
k2E_immune_mpb1 = param_immune_cis[param_immune_cis.Dataset .== "MPB1", :k2E][1];

# Load fitted parameters from timeshift fitting (will be initial estimates for fitting the immune model)
time_shift_result = CSV.read(joinpath(@__DIR__, "..", "Fitted_params_results", "fitted_parameters_normalised_timeshift.csv"), DataFrame);
N0_treatment_initiation_mp = time_shift_result[time_shift_result.Dataset .== "MP", :nk_cells_at_treatment_initiation][1];
N0_treatment_initiation_mpb1 = time_shift_result[time_shift_result.Dataset .== "MPB1", :nk_cells_at_treatment_initiation][1];
E0_treatment_initiation_mp = time_shift_result[time_shift_result.Dataset .== "MP", :t_cells_at_shift_at_treatment_initiation][1];
E0_treatment_initiation_mpb1 = time_shift_result[time_shift_result.Dataset .== "MPB1", :t_cells_at_shift_at_treatment_initiation][1];

# Definition of the model fixed parameters
fixed_params_mp = fixed_params_immune_cis_mp;
fixed_params_mpb1 = fixed_params_immune_cis_mpb1;

t_fine_mp = range(0, maximum(tdata_mp_cis), length=200);
t_fine_mpb1 = range(0, maximum(tdata_mpb1_cis), length=200);
t_immune_mp = 10.0;
t_immune_mpb1 = 10.0;
u0_immune_cis_mp = [0, 0, ydata_mp_cis[1], 0, 0, 0, N0_treatment_initiation_mp, E0_treatment_initiation_mp];
u0_immune_cis_mpb1 = [0, 0, ydata_mpb1_cis[1], 0, 0, 0, N0_treatment_initiation_mpb1, E0_treatment_initiation_mpb1];

predictions_cisPKPD_mp = Predictions.predict_tumour_immune_cis_pkpd_states(t_fine_mp, [k2_immune_mp, k2N_immune_mp, k2E_immune_mp], u0_immune_cis_mp, fixed_params_mp, dose_amt, dose_starts);
predictions_cisPKPD_mpb1 = Predictions.predict_tumour_immune_cis_pkpd_states(t_fine_mpb1, [k2_immune_mpb1, k2N_immune_mpb1, k2E_immune_mpb1], u0_immune_cis_mpb1, fixed_params_mpb1, dose_amt, dose_starts);

# Extract T, N, E for each dataset
yfit_T_cisPKPD_mp = predictions_cisPKPD_mp[1, :];
yfit_N_cisPKPD_mp = predictions_cisPKPD_mp[2, :];
yfit_E_cisPKPD_mp = predictions_cisPKPD_mp[3, :];
yfit_T_cisPKPD_mpb1 = predictions_cisPKPD_mpb1[1, :];
yfit_N_cisPKPD_mpb1 = predictions_cisPKPD_mpb1[2, :];
yfit_E_cisPKPD_mpb1 = predictions_cisPKPD_mpb1[3, :];

# Compute shared y-limits per row so the two columns use the same y-axis range
t_cell_count_mp_cisplatin = t_cells_prop_cis_mp * prop_cd8cells_in_t_cells * prop_cd45_cells * total_cell_count;
t_cell_count_mpb1_cisplatin = t_cells_prop_cis_mpb1 * prop_cd8cells_in_t_cells * prop_cd45_cells * total_cell_count;

# -------------------------------------

# Fixing k2 - Load fitted parameters
param_immune_cis_fixed_k2 = CSV.read(joinpath(@__DIR__, "..", "Fitted_params_results", "fitted_parameters_normalised_immune_cisplatin.csv"), DataFrame);

# Extract immune model parameters for MP and MPB1
k2_immune_mp_fixed_k2 = 8.07;
k2N_immune_mp_fixed_k2 = param_immune_cis_fixed_k2[param_immune_cis_fixed_k2.Dataset .== "MP", :k2N][1];
k2E_immune_mp_fixed_k2 = param_immune_cis_fixed_k2[param_immune_cis_fixed_k2.Dataset .== "MP", :k2E][1];
k2_immune_mpb1_fixed_k2 = 63.7;
k2N_immune_mpb1_fixed_k2 = param_immune_cis_fixed_k2[param_immune_cis_fixed_k2.Dataset .== "MPB1", :k2N][1];
k2E_immune_mpb1_fixed_k2 = param_immune_cis_fixed_k2[param_immune_cis_fixed_k2.Dataset .== "MPB1", :k2E][1];

predictions_cis_mp_fixed_k2 = Predictions.predict_tumour_immune_cis_pkpd_states(t_fine_mp, [k2_immune_mp_fixed_k2, k2N_immune_mp_fixed_k2, k2E_immune_mp_fixed_k2], u0_immune_cis_mp, fixed_params_mp, dose_amt, dose_starts);
predictions_cis_mpb1_fixed_k2 = Predictions.predict_tumour_immune_cis_pkpd_states(t_fine_mpb1, [k2_immune_mpb1_fixed_k2, k2N_immune_mpb1_fixed_k2, k2E_immune_mpb1_fixed_k2], u0_immune_cis_mpb1, fixed_params_mpb1, dose_amt, dose_starts);

# Extract T, N, E for each dataset
yfit_T_cis_mp_fixed_k2 = predictions_cis_mp_fixed_k2[1, :];
yfit_N_cis_mp_fixed_k2 = predictions_cis_mp_fixed_k2[2, :];
yfit_E_cis_mp_fixed_k2 = predictions_cis_mp_fixed_k2[3, :];
yfit_T_cis_mpb1_fixed_k2 = predictions_cis_mpb1_fixed_k2[1, :];
yfit_N_cis_mpb1_fixed_k2 = predictions_cis_mpb1_fixed_k2[2, :];
yfit_E_cis_mpb1_fixed_k2 = predictions_cis_mpb1_fixed_k2[3, :];

# -------------------------------------

# No refitting, no direct effect of cisplatin on immune cells (recruitment = death term)
k2_immune_mp_fixed_k2 = 8.07;
k2N_immune_mp_fixed_k2 = k2_immune_mp_fixed_k2;
k2E_immune_mp_fixed_k2 = k2_immune_mp_fixed_k2;
k2_immune_mpb1_fixed_k2 = 63.7;
k2N_immune_mpb1_fixed_k2 = k2_immune_mpb1_fixed_k2;
k2E_immune_mpb1_fixed_k2 = k2_immune_mpb1_fixed_k2;

predictions_cis_mp_norefit = Predictions.predict_tumour_immune_cis_pkpd_states(t_fine_mp, [k2_immune_mp_fixed_k2, k2N_immune_mp_fixed_k2, k2E_immune_mp_fixed_k2], u0_immune_cis_mp, fixed_params_mp, dose_amt, dose_starts);
predictions_cis_mpb1_norefit = Predictions.predict_tumour_immune_cis_pkpd_states(t_fine_mpb1, [k2_immune_mpb1_fixed_k2, k2N_immune_mpb1_fixed_k2, k2E_immune_mpb1_fixed_k2], u0_immune_cis_mpb1, fixed_params_mpb1, dose_amt, dose_starts);

# Extract T, N, E for each dataset
yfit_T_cis_mp_norefit = predictions_cis_mp_norefit[1, :];
yfit_N_cis_mp_norefit = predictions_cis_mp_norefit[2, :];
yfit_E_cis_mp_norefit = predictions_cis_mp_norefit[3, :];
yfit_T_cis_mpb1_norefit = predictions_cis_mpb1_norefit[1, :];
yfit_N_cis_mpb1_norefit = predictions_cis_mpb1_norefit[2, :];
yfit_E_cis_mpb1_norefit = predictions_cis_mpb1_norefit[3, :];

# ====================================
# Figure 5
# ====================================

mp_data = [tdata_mp, 
             ydata_mp, ydata_mp_cis,
             t_cell_count_mp, t_cell_count_mp_cis,
             nk_cell_count_mp, nk_cell_count_mp_cis];
mp_error = [yerr_mp, yerr_mp_cis,
              t_cell_count_mp_sem, t_cell_count_mp_cis_sem, 
              nk_cell_count_mp_sem, nk_cell_count_mp_cis_sem];
mp_sol = [t_fine_mp, 
            yfit_T_mp, yfit_T_cis_mp_fixed_k2, yfit_T_cis_mp_norefit,
            yfit_E_mp, yfit_E_cis_mp_fixed_k2, yfit_E_cis_mp_norefit,
            yfit_N_mp, yfit_N_cis_mp_fixed_k2, yfit_N_cis_mp_norefit];

plot_immune_cis(mp_data, mp_sol, mp_error, "MP (Cisplatin)", [vehicle_colour, cis_colour]; seelegend = true);

mpb1_data = [tdata_mpb1, 
             ydata_mpb1, ydata_mpb1_cis,
             t_cell_count_mpb1, t_cell_count_mpb1_cis,
             nk_cell_count_mpb1, nk_cell_count_mpb1_cis];
mpb1_error = [yerr_mpb1, yerr_mpb1_cis,
              t_cell_count_mpb1_sem, t_cell_count_mpb1_cis_sem, 
              nk_cell_count_mpb1_sem, nk_cell_count_mpb1_cis_sem];
mpb1_sol = [t_fine_mpb1, 
            yfit_T_mpb1, yfit_T_cis_mpb1_fixed_k2, yfit_T_cis_mpb1_norefit,
            yfit_E_mpb1, yfit_E_cis_mpb1_fixed_k2, yfit_E_cis_mpb1_norefit,
            yfit_N_mpb1, yfit_N_cis_mpb1_fixed_k2, yfit_N_cis_mpb1_norefit];

plot_immune_cis(mpb1_data, mpb1_sol, mpb1_error, "MPB1 (Cisplatin)", [vehicle_colour, cis_colour]; seelegend = false);

## -------------------------------------
# Tumour + immune cells with ICB

# Definition of the model fixed parameters
fixed_params_mp_icb = fixed_params_mp_icb;
fixed_params_mpb1_icb = fixed_params_mpb1_icb;

# Load fitted parameters from timeshift fitting (will be initial estimates for fitting the immune model)
time_shift_result = CSV.read(joinpath(@__DIR__, "..", "Fitted_params_results", "fitted_parameters_normalised_timeshift.csv"), DataFrame);

N0_treatment_initiation_mp = time_shift_result[time_shift_result.Dataset .== "MP", :nk_cells_at_treatment_initiation][1];
N0_treatment_initiation_mpb1 = time_shift_result[time_shift_result.Dataset .== "MPB1", :nk_cells_at_treatment_initiation][1];
E0_treatment_initiation_mp = time_shift_result[time_shift_result.Dataset .== "MP", :t_cells_at_shift_at_treatment_initiation][1];
E0_treatment_initiation_mpb1 = time_shift_result[time_shift_result.Dataset .== "MPB1", :t_cells_at_shift_at_treatment_initiation][1];

u0_immune_icb_mp = [ydata_mp_icb[1], N0_treatment_initiation_mp, E0_treatment_initiation_mp];
u0_immune_icb_mpb1 = [ydata_mpb1_icb[1], N0_treatment_initiation_mpb1, E0_treatment_initiation_mpb1];

# Load fitted parameters
param_immune_icb = CSV.read(joinpath(@__DIR__, "..", "Fitted_params_results/fitted_parameters_normalised_immune_icb.csv"), DataFrame);

# Extract immune model parameters for MP and MPB1
q_immune_mp = param_immune_icb[param_immune_icb.Dataset .== "MP", :q][1];
q_immune_mpb1 = param_immune_icb[param_immune_icb.Dataset .== "MPB1", :q][1];

# SHIFTED TIME SCALE
t_fine_mp = range(0, maximum(tdata_mp_icb), length=200);
t_fine_mpb1 = range(0, maximum(tdata_mpb1_icb), length=200);
predictions_mp_icb = Predictions.predict_tumour_immune_icb_q_states(t_fine_mp, [q_immune_mp], u0_immune_icb_mp, fixed_params_mp_icb);
predictions_mpb1_icb = Predictions.predict_tumour_immune_icb_q_states(t_fine_mpb1, [q_immune_mpb1], u0_immune_icb_mpb1, fixed_params_mpb1_icb);

# Extract T, N, E for each dataset
yfit_T_icb_mp = predictions_mp_icb[1, :];
yfit_N_icb_mp = predictions_mp_icb[2, :];
yfit_E_icb_mp = predictions_mp_icb[3, :];

yfit_T_icb_mpb1 = predictions_mpb1_icb[1, :];
yfit_N_icb_mpb1 = predictions_mpb1_icb[2, :];
yfit_E_icb_mpb1 = predictions_mpb1_icb[3, :];

# ====================================
# Figure 6
# ====================================

mp_data = [tdata_mp, 
             ydata_mp, ydata_mp_icb,
             t_cell_count_mp,
             nk_cell_count_mp];
mp_error = [yerr_mp, yerr_mp_icb,
              t_cell_count_mp_sem, 
              nk_cell_count_mp_sem];
mp_sol = [t_fine_mp, 
            yfit_T_mp, yfit_T_icb_mp,
            yfit_E_mp, yfit_E_icb_mp,
            yfit_N_mp, yfit_N_icb_mp];

plot_immune_icb(mp_data, mp_sol, mp_error, "MP (ICB)", [vehicle_colour, icb_colour]; seelegend = true);

mpb1_data = [tdata_mpb1, 
             ydata_mpb1, ydata_mpb1_icb,
             t_cell_count_mpb1,
             nk_cell_count_mpb1];
mpb1_error = [yerr_mpb1, yerr_mpb1_icb,
              t_cell_count_mpb1_sem, 
              nk_cell_count_mpb1_sem];
mpb1_sol = [t_fine_mpb1, 
            yfit_T_mpb1, yfit_T_icb_mpb1,
            yfit_E_mpb1, yfit_E_icb_mpb1,
            yfit_N_mpb1, yfit_N_icb_mpb1];

plot_immune_icb(mpb1_data, mpb1_sol, mpb1_error, "MPB1 (ICB)", [vehicle_colour, icb_colour]; seelegend = false);

## -------------------------------------
# Tumour + immune cells with all treatment models (Cis, ICB & Cis + ICB)

# Initial conditions for immune cells:
time_shift_result = CSV.read(joinpath(@__DIR__, "..", "Fitted_params_results", "fitted_parameters_normalised_timeshift.csv"), DataFrame);
N0_treatment_initiation_mp = time_shift_result[time_shift_result.Dataset .== "MP", :nk_cells_at_treatment_initiation][1];
N0_treatment_initiation_mpb1 = time_shift_result[time_shift_result.Dataset .== "MPB1", :nk_cells_at_treatment_initiation][1];
E0_treatment_initiation_mp = time_shift_result[time_shift_result.Dataset .== "MP", :t_cells_at_shift_at_treatment_initiation][1];
E0_treatment_initiation_mpb1 = time_shift_result[time_shift_result.Dataset .== "MPB1", :t_cells_at_shift_at_treatment_initiation][1];

# Load fitted parameters
param_immune_cis = CSV.read(joinpath(@__DIR__, "..", "Fitted_params_results", "fitted_parameters_normalised_immune_cisplatin.csv"), DataFrame);

# Extract immune model parameters for MP and MPB1
#k2_immune_mp = param_immune_cis[param_immune_cis.Dataset .== "MP", :k2][1];
k2_immune_mp = 8.07;
k2N_immune_mp = param_immune_cis[param_immune_cis.Dataset .== "MP", :k2N][1];
k2E_immune_mp = param_immune_cis[param_immune_cis.Dataset .== "MP", :k2E][1];
#k2_immune_mpb1 = param_immune_cis[param_immune_cis.Dataset .== "MPB1", :k2][1];
k2_immune_mpb1 = 63.7;
k2N_immune_mpb1 = param_immune_cis[param_immune_cis.Dataset .== "MPB1", :k2N][1];
k2E_immune_mpb1 = param_immune_cis[param_immune_cis.Dataset .== "MPB1", :k2E][1];

# Load fitted parameters
param_immune_icb = CSV.read(joinpath(@__DIR__, "..", "Fitted_params_results/fitted_parameters_normalised_immune_icb.csv"), DataFrame);

# Extract immune model parameters for MP and MPB1
q_immune_mp = param_immune_icb[param_immune_icb.Dataset .== "MP", :q][1];
q_immune_mpb1 = param_immune_icb[param_immune_icb.Dataset .== "MPB1", :q][1];

# Definition of the model fixed parameters
fixed_params_mp_cis_icb = deepcopy(fixed_params_icb_cis_mp);
fixed_params_mpb1_cis_icb = deepcopy(fixed_params_icb_cis_mpb1);

# Fixed bolus dosing setup
BW = 0.025;
dose_amt = 3.0 * BW;
initial_day_dose = 0.01;
dose_starts = [initial_day_dose, initial_day_dose + 7.0];  # Doses 1 week apart

# Initial conditions
u0_immune_cis_icb_mp = [0, 0, ydata_mp_cis_icb[1], 0, 0, 0, N0_treatment_initiation_mp, E0_treatment_initiation_mp];
u0_immune_cis_icb_mpb1 = [0, 0, ydata_mpb1_cis_icb[1], 0, 0, 0, N0_treatment_initiation_mpb1, E0_treatment_initiation_mpb1];

t_fine_mp = range(0, maximum(tdata_mp_icb), length=200);
t_fine_mpb1 = range(0, maximum(tdata_mpb1_icb), length=200);
predictions_cis_icb_mp = Predictions.predict_tumour_immune_cis_pkpd_states(t_fine_mp, [k2_immune_mp, k2N_immune_mp, k2E_immune_mp], u0_immune_cis_icb_mp, fixed_params_mp_cis_icb, dose_amt, dose_starts);
predictions_cis_icb_mpb1 = Predictions.predict_tumour_immune_cis_pkpd_states(t_fine_mpb1, [k2_immune_mpb1, k2N_immune_mpb1, k2E_immune_mpb1], u0_immune_cis_icb_mpb1, fixed_params_mpb1_cis_icb, dose_amt, dose_starts);

# Extract T, N, E for each dataset
yfit_T_cis_icb_mp = predictions_cis_icb_mp[1, :];
yfit_N_cis_icb_mp = predictions_cis_icb_mp[2, :];
yfit_E_cis_icb_mp = predictions_cis_icb_mp[3, :];
yfit_T_cis_icb_mpb1 = predictions_cis_icb_mpb1[1, :];
yfit_N_cis_icb_mpb1 = predictions_cis_icb_mpb1[2, :];
yfit_E_cis_icb_mpb1 = predictions_cis_icb_mpb1[3, :];

# ====================================
# Figure 7
# ====================================

mp_data = [tdata_mp, 
             ydata_mp, ydata_mp_cis, ydata_mp_icb, ydata_mp_cis_icb,
             t_cell_count_mp, t_cell_count_mp_cis,
             nk_cell_count_mp, nk_cell_count_mp_cis];
mp_error = [yerr_mp, yerr_mp_cis, yerr_mp_icb, yerr_mp_cis_icb,
              t_cell_count_mp_sem, t_cell_count_mp_cis_sem, 
              nk_cell_count_mp_sem, nk_cell_count_mp_cis_sem];
mp_sol = [t_fine_mp, 
            yfit_T_mp, yfit_T_cis_mp_fixed_k2, yfit_T_icb_mp, yfit_T_cis_icb_mp,
            yfit_E_mp, yfit_E_cis_mp_fixed_k2, yfit_E_icb_mp, yfit_E_cis_icb_mp,
            yfit_N_mp, yfit_N_cis_mp_fixed_k2, yfit_N_icb_mp, yfit_N_cis_icb_mp];

plot_immune_all(mp_data, mp_sol, mp_error, "MP (All treatments)", [vehicle_colour, cis_colour, icb_colour, cis_icb_colour]; seelegend = true);

mpb1_data = [tdata_mpb1, 
             ydata_mpb1,  ydata_mpb1_cis, ydata_mpb1_icb, ydata_mpb1_cis_icb,
             t_cell_count_mpb1, t_cell_count_mpb1_cis,
             nk_cell_count_mpb1, nk_cell_count_mpb1_cis];
mpb1_error = [yerr_mpb1, yerr_mpb1_cis, yerr_mpb1_icb, yerr_mpb1_cis_icb,
              t_cell_count_mpb1_sem, t_cell_count_mpb1_cis_sem,  
              nk_cell_count_mpb1_sem, nk_cell_count_mpb1_cis_sem];
mpb1_sol = [t_fine_mpb1, 
            yfit_T_mpb1, yfit_T_cis_mpb1_fixed_k2, yfit_T_icb_mpb1, yfit_T_cis_icb_mpb1,
            yfit_E_mpb1, yfit_E_cis_mpb1_fixed_k2, yfit_E_icb_mpb1, yfit_E_cis_icb_mpb1,
            yfit_N_mpb1, yfit_N_cis_mpb1_fixed_k2, yfit_N_icb_mpb1, yfit_N_cis_icb_mpb1];

plot_immune_all(mpb1_data, mpb1_sol, mpb1_error, "MPB1 (All treatments)", [vehicle_colour, cis_colour, icb_colour, cis_icb_colour]; seelegend = false);