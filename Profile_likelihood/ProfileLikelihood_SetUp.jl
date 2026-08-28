using DifferentialEquations
using LsqFit, StatsBase, Distributions
using CSV, DataFrames, JLD2, PlotlyJS

include("../src/craig_lab_template.jl");
include("../src/Models.jl");
include("../src/Data.jl");
include("../src/Predictions.jl");
include("../src/Results.jl");
include("../src/Assumptions.jl");
include("../src/Fixed_params.jl");

using .Predictions
using .Results

function PL(pidx, prange, pbest, lower, upper, model, tdata, ydata)
    # Return the profile likelihood for the parameter with index "pidx"

    N_points = length(prange)
    N_param  = length(pbest)

    PLL_CF = zeros(N_points)
    PLL_p  = zeros(N_points, N_param)

    # indices of free parameters to optimize
    free_idx = setdiff(1:N_param, [pidx])

    for i in eachindex(prange)

        pfix = prange[i]

        # initial guess for free parameters
        p0_free = pbest[free_idx]

        # bounds for free parameters
        lb_free = lower[free_idx]
        ub_free = upper[free_idx]

        # model wrapper for curve_fit which fix parameter at index pfix
        model_fixed = function(t, pfree)

            pfull = similar(pfree, N_param)

            pfull[pidx] = pfix
            pfull[free_idx] .= pfree

            return model(t, pfull)

        end

        # optimize the free parameters
        fit = curve_fit(model_fixed, tdata, ydata, p0_free, lower = lb_free, upper = ub_free)

        # reconstruct full parameter vector
        pfull = similar(pbest)

        pfull[pidx] = pfix
        pfull[free_idx] .= coef(fit)

        PLL_CF[i] = sum(abs2, residuals(fit))
        PLL_p[i,:] .= pfull

    end

    return PLL_CF, PLL_p
end

#-------------------------------------

function profile_parameter(model, tdata, ydata;
    p0,
    lower,
    upper,
    pidx,
    pminfold=0.99,
    pmaxfold=1.0,
    # xlabel="Parameter",
    # title="Profile likelihood",
    # filename=nothing,
)

    # Fit model
    fit = curve_fit(model, tdata, ydata, p0; lower=lower, upper=upper)

    pbest = coef(fit)
    CFbest = sum(abs2, residuals(fit))

    n = length(ydata)
    N_param = length(pbest)

    if pbest[pidx] == 0
        pmin = -0.01; #-pminfold*1.0
        pmax = 100; #pmaxfold*1.0
    else   
        pmin = pbest[pidx] - pminfold*abs(pbest[pidx])
        pmax = pbest[pidx] + pmaxfold*abs(pbest[pidx])
    end
    prange = range(pmin, pmax,length=501);

    # Profile likelihood
    CF_PLL, p_PLL = PL(pidx, prange, pbest, lower, upper, model, tdata, ydata)

    m2LL_best = n*log((2π/(n-N_param))*CFbest) + (n-N_param)
    CI = m2LL_best + quantile(Chisq(1),0.95)

    m2LL = n .* log.((2π/(n-N_param)) .* CF_PLL) .+ (n-N_param)

    idx = findall(i -> (m2LL[i] - CI) * (m2LL[i+1] - CI) <= 0, 1:length(m2LL)-1)
    
    if length(idx) == 2 
        CI95_min = prange[idx[1]]; 
        CI95_max = prange[idx[2]];
    elseif length(idx) == 1  
        if prange[idx[1]] < pbest[pidx]
            CI95_min = prange[idx[1]]; 
            CI95_max = Inf;
        elseif prange[idx[1]] > pbest[pidx]
            CI95_max = prange[idx[1]]; 
            CI95_min = -Inf;
        end 
    elseif length(idx) > 2
        CI95_max = Inf;
        CI95_min = -Inf;

        for i in 1:2:length(idx)-1
            left  = idx[i]
            right = idx[i+1]

            mid = div(left + right, 2)

            if prange[left] <= pbest[pidx] <= prange[right] && m2LL[mid] <= CI
                CI95_min = prange[left]
                CI95_max = prange[right]
                break
            end
        end
    else
        CI95_max = Inf;
        CI95_min = -Inf;
    end

    return (
        prange = prange,    
        m2LL = m2LL,
        CI = CI,
        pbest = pbest[pidx], 
        CI95 = (CI95_min, CI95_max)
    )
end

#-------------------------------------

function plot_profile_likelihood(
    models, pminfold, pmaxfold;
    subtitles = nothing,
    filename = nothing,
    template = craig_lab_template,
    seeaxis = true
)

    nModels = length(models)
    nParams = length(models[1].p0)

    # Confidence intervals
    CI95 = DataFrame(
        Model = String[],
        Parameter = String[],
        CI95 = Tuple{Float64,Float64}[]
    )

    # Make sure Figures directory exists
    if filename !== nothing
        mkpath(joinpath(pwd(), "Figures"))
    end

    for (model_idx, m) in enumerate(models)

        for param_idx in 1:nParams

            res = profile_parameter(
                m.model,
                m.t,
                m.y;
                p0 = m.p0,
                lower = m.lower,
                upper = m.upper,
                pidx = param_idx,
                pminfold = pminfold[param_idx],
                pmaxfold = pmaxfold[param_idx]
            )

            # Store confidence interval
            push!(
                CI95,
                (
                    m.name,
                    m.param_names[param_idx],
                    res.CI95
                )
            )

            m2ll_trace = scatter(
                x = collect(res.prange),
                y = res.m2LL,
                mode = "lines",
                name = "-2LL",
                showlegend = false,
                line = attr(width = 6, color = "#835258")
            )

            cutoff_trace = scatter(
                x = collect(res.prange),
                y = fill(res.CI, length(res.prange)),
                mode = "lines",
                name = "95% cutoff",
                showlegend = false,
                line = attr(width = 6, color = "#545B77")
            )

            bestfit_trace = scatter(
                x = [res.pbest, res.pbest],
                y = [minimum(res.m2LL) - 5, maximum(vcat(res.m2LL, [res.CI])) + 5],
                mode = "lines",
                name = "Best fit",
                showlegend = false,
                line = attr(width = 6, color = "#dbdbdb", dash = "dot")
            )
            

            xatt = attr(
                title = m.param_labels[param_idx],
                titlefont = attr(size = 28),
                tickfont = attr(size = 28),
                nticks = 5
            )

            if seeaxis == true

                yatt = attr(
                    title = "-2 Log-Likelihood",
                    titlefont = attr(size = 28),
                    tickfont = attr(size = 28),
                    range = [minimum(res.m2LL) - 5, maximum(vcat(res.m2LL, [res.CI])) + 5],
                    nticks = 5
                )

            else

                yatt = attr(
                    title = "",
                    titlefont = attr(size = 28),
                    tickfont = attr(size = 28),
                    range = [minimum(res.m2LL) - 5, maximum(vcat(res.m2LL, [res.CI])) + 5],
                    nticks = 5
                )

            end

            annotations = Any[]

            if subtitles !== nothing

                subtitle_idx = (model_idx - 1) * nParams + param_idx
                push!(
                    annotations,
                    attr(
                        text = "<b>$(subtitles[subtitle_idx])</b>",
                        x = 0.5,
                        y = 1.02,
                        xref = "paper",
                        yref = "paper",
                        xanchor = "center",
                        yanchor = "bottom",
                        showarrow = false,
                        font = attr(size = 32)
                    )
                )

            end

            layout = Layout(
                template = template,
                annotations = annotations,
                xaxis = xatt,
                yaxis = yatt,
                width = 900,
                height = 600,
                margin = attr(l = 130, r = 20, t = 80, b = 80),
                showlegend = false
            )

            fig = plot([m2ll_trace, cutoff_trace, bestfit_trace], layout)
            display(fig)

            if filename !== nothing

                parameter_name = m.param_names[param_idx]
                plot_filename = joinpath(pwd(), "Profile Likelihood", "$(filename) - $(m.name) - $(parameter_name).svg")
                savefig(fig, plot_filename; width = 900, height = 600)

            end

        end
    end

    return CI95

end

#-------------------------------------

# Importing the data and defining time, tumor volumes, and initial conditions for each dataset
data = DataTools.load_tumor_datasets(data_dir = joinpath(@__DIR__, "..", "Data"));

tdata_mp = data.mp.t;
ydata_mp = data.mp.y;
tdata_mp_cis = data.mp_cis.t;
ydata_mp_cis = data.mp_cis.y;
tdata_mp_icb = data.mp_icb.t;
ydata_mp_icb = data.mp_icb.y;
tdata_mp_cis_icb = data.mp_cis_icb.t;
ydata_mp_cis_icb = data.mp_cis_icb.y;

tdata_mpb1 = data.mpb1.t;
ydata_mpb1 = data.mpb1.y;
tdata_mpb1_cis = data.mpb1_cis.t;
ydata_mpb1_cis = data.mpb1_cis.y;
tdata_mpb1_icb = data.mpb1_icb.t;
ydata_mpb1_icb = data.mpb1_icb.y;
tdata_mpb1_cis_icb = data.mpb1_cis_icb.t;
ydata_mpb1_cis_icb = data.mpb1_cis_icb.y;

u0_mp   = data.mp.u0;
u0_mpb1 = data.mpb1.u0;
u0_mp_cis   = data.mp_cis.u0;
u0_mpb1_cis = data.mpb1_cis.u0;
u0_mp_icb   = data.mp_icb.u0;
u0_mpb1_icb = data.mpb1_icb.u0;
u0_mp_cis_icb   = data.mp_cis_icb.u0;
u0_mpb1_cis_icb = data.mpb1_cis_icb.u0;

#-------------------------------------
# Tumour growth only
#-------------------------------------

# Define tumor growth models
tumour_model_logistic_mp           = (t, p) -> Predictions.predict_logistic(t, p, u0_mp);
tumour_model_logistic_mp_cis       = (t, p) -> Predictions.predict_logistic(t, p, u0_mp_cis);
tumour_model_logistic_mpb1         = (t, p) -> Predictions.predict_logistic(t, p, u0_mpb1);
tumour_model_logistic_mpb1_cis     = (t, p) -> Predictions.predict_logistic(t, p, u0_mpb1_cis);
tumour_model_logistic_mp_icb       = (t, p) -> Predictions.predict_logistic(t, p, u0_mp_icb);
tumour_model_logistic_mpb1_icb     = (t, p) -> Predictions.predict_logistic(t, p, u0_mpb1_icb);
tumour_model_logistic_mp_cis_icb   = (t, p) -> Predictions.predict_logistic(t, p, u0_mp_cis_icb);
tumour_model_logistic_mpb1_cis_icb = (t, p) -> Predictions.predict_logistic(t, p, u0_mpb1_cis_icb);

tumour_model_gompertz_mp   = (t, p) -> Predictions.predict_gompertz(t, p, u0_mp);
tumour_model_gompertz_mpb1 = (t, p) -> Predictions.predict_gompertz(t, p, u0_mpb1);

# Initial estimates for [r, k]
r0 = 0.514;  # (1/days), Tumor growth rate from de Pillis
k0 = 2000.0; # (mm^3), Tumor carrying capacity (could be more fine-tuned)

# Bounds
rmin = 0.0;
rmax = 1.0;
kmin = 0.0;
kmax = 20000.0;

#-------------------------------------
# Tumour growth with immune environment
#-------------------------------------

# Load immune data and proportions via DataTools
immune = DataTools.load_immune_datasets(data_dir = joinpath(@__DIR__, "..", "Data"));

df_immune_mp   = immune.df_immune_mp;
df_immune_mpb1 = immune.df_immune_mpb1;
nk_cells_prop_mp   = immune.nk_cells_prop_mp;
nk_cells_prop_mpb1 = immune.nk_cells_prop_mpb1;
t_cells_prop_mp   = immune.t_cells_prop_mp;
t_cells_prop_mpb1 = immune.t_cells_prop_mpb1;

# Cell counts for NK and T cells, using the assumptions defined in Assumptions.jl
nk_cell_count_mp   = nk_cells_prop_mp * prop_cd45_cells * total_cell_count;
nk_cell_count_mpb1 = nk_cells_prop_mpb1 * prop_cd45_cells * total_cell_count;
t_cell_count_mp    = t_cells_prop_mp * prop_cd8cells_in_t_cells * prop_cd45_cells * total_cell_count;
t_cell_count_mpb1  = t_cells_prop_mpb1 * prop_cd8cells_in_t_cells * prop_cd45_cells * total_cell_count;

# Combine tumour volume, NK cell count, and T cell count data for fitting
ydata_mp_combined   = vcat(ydata_mp, [nk_cell_count_mp, t_cell_count_mp]);
ydata_mpb1_combined = vcat(ydata_mpb1, [nk_cell_count_mpb1, t_cell_count_mpb1]);

# Load fitted parameters from timeshift fitting (will be initial estimates for fitting the immune model)
time_shift = CSV.read(joinpath(@__DIR__, "..", "Fitted_params_results", "fitted_parameters_normalised_timeshift.csv"), DataFrame);

# Case where we initialise the tumour at treatment initiation 
N0_treatment_initiation_mp   = time_shift[time_shift.Dataset .== "MP", :nk_cells_at_treatment_initiation][1];
N0_treatment_initiation_mpb1 = time_shift[time_shift.Dataset .== "MPB1", :nk_cells_at_treatment_initiation][1];
E0_treatment_initiation_mp   = time_shift[time_shift.Dataset .== "MP", :t_cells_at_shift_at_treatment_initiation][1];
E0_treatment_initiation_mpb1 = time_shift[time_shift.Dataset .== "MPB1", :t_cells_at_shift_at_treatment_initiation][1];

# Time at which the immune cells were measured
t_immune_mp   = 10.0;
t_immune_mpb1 = 10.0;

# Initial conditions for immune cells and tumour growth model
u0_immune_mp   = [ydata_mp[1], N0_treatment_initiation_mp, E0_treatment_initiation_mp];
u0_immune_mpb1 = [ydata_mpb1[1], N0_treatment_initiation_mpb1, E0_treatment_initiation_mpb1];

# Define tumor growth models
tumour_immune_model_mp   = (t, p) -> Predictions.predict_tumour_immune(t, t_immune_mp, p, u0_immune_mp, fixed_params_immune_mp);
tumour_immune_model_mpb1 = (t, p) -> Predictions.predict_tumour_immune(t, t_immune_mpb1, p, u0_immune_mpb1, fixed_params_immune_mpb1);

# Initial estimates for fitted params
r0_immune_mp    = time_shift[time_shift.Dataset .== "MP", :r][1];
r0_immune_mpb1  = time_shift[time_shift.Dataset .== "MPB1", :r][1];
pN0_immune_mp   = time_shift[time_shift.Dataset .== "MP", :pN][1];
pN0_immune_mpb1 = time_shift[time_shift.Dataset .== "MPB1", :pN][1];
q0_immune_mp    = time_shift[time_shift.Dataset .== "MP", :q][1];
q0_immune_mpb1  = time_shift[time_shift.Dataset .== "MPB1", :q][1];

# Bounds
rmin_immune  = 0.0;
rmax_immune  = 5.0;
pNmin_immune = 0.00001;
pNmax_immune = 1.0;
qmin_immune  = 0.00001;
qmax_immune  = 1.0;

#-------------------------------------
# Tumour growth with immune environment - NORMALIZED
#-------------------------------------

# Scaling factors
tumour_scale_mp = maximum(ydata_mp);
nk_scale_mp     = maximum([nk_cell_count_mp, 1.0]);
t_scale_mp      = maximum([t_cell_count_mp, 1.0]);

tumour_scale_mpb1 = maximum(ydata_mpb1);
nk_scale_mpb1     = maximum([nk_cell_count_mpb1, 1.0]);
t_scale_mpb1      = maximum([t_cell_count_mpb1, 1.0]);

# Scaled data
ydata_mp_combined_scaled   = vcat(ydata_mp ./ tumour_scale_mp, [nk_cell_count_mp / nk_scale_mp, t_cell_count_mp / t_scale_mp]);
ydata_mpb1_combined_scaled = vcat(ydata_mpb1 ./ tumour_scale_mpb1, [nk_cell_count_mpb1 / nk_scale_mpb1, t_cell_count_mpb1 / t_scale_mpb1]);

time_shift_scaled = CSV.read(joinpath(@__DIR__, "..", "Fitted_params_results", "fitted_parameters_normalised_timeshift.csv"), DataFrame);

# Case where we initialise the tumour at treatment initiation 
N0_treatment_initiation_mp_scaled   = time_shift_scaled[time_shift_scaled.Dataset .== "MP", :nk_cells_at_treatment_initiation][1];
N0_treatment_initiation_mpb1_scaled = time_shift_scaled[time_shift_scaled.Dataset .== "MPB1", :nk_cells_at_treatment_initiation][1];
E0_treatment_initiation_mp_scaled   = time_shift_scaled[time_shift_scaled.Dataset .== "MP", :t_cells_at_shift_at_treatment_initiation][1];
E0_treatment_initiation_mpb1_scaled = time_shift_scaled[time_shift_scaled.Dataset .== "MPB1", :t_cells_at_shift_at_treatment_initiation][1];

# Time at which the immune cells were measured
t_immune_mp   = 10.0;
t_immune_mpb1 = 10.0;

# Initial conditions for immune cells and tumour growth model
u0_immune_mp_scaled   = [ydata_mp[1], N0_treatment_initiation_mp_scaled, E0_treatment_initiation_mp_scaled];
u0_immune_mpb1_scaled = [ydata_mpb1[1], N0_treatment_initiation_mpb1_scaled, E0_treatment_initiation_mpb1_scaled];

# Define models
tumour_immune_model_mp_scaled = (t, p) -> begin
    pred = Predictions.predict_tumour_immune(t, t_immune_mp, p, u0_immune_mp_scaled, fixed_params_immune_mp)
    return vcat(pred[1:length(ydata_mp)] ./ tumour_scale_mp,
                pred[end - 1] / nk_scale_mp,
                pred[end] / t_scale_mp)
end;

tumour_immune_model_mpb1_scaled = (t, p) -> begin
    pred = Predictions.predict_tumour_immune(t, t_immune_mpb1, p, u0_immune_mpb1_scaled, fixed_params_immune_mpb1)
    return vcat(pred[1:length(ydata_mpb1)] ./ tumour_scale_mpb1,
                pred[end - 1] / nk_scale_mpb1,
                pred[end] / t_scale_mpb1)
end;

# Initial estimates for fitted params
r0_immune_mp_scaled    = time_shift_scaled[time_shift_scaled.Dataset .== "MP", :r][1];
r0_immune_mpb1_scaled  = time_shift_scaled[time_shift_scaled.Dataset .== "MPB1", :r][1];
pN0_immune_mp_scaled   = time_shift_scaled[time_shift_scaled.Dataset .== "MP", :pN][1];
pN0_immune_mpb1_scaled = time_shift_scaled[time_shift_scaled.Dataset .== "MPB1", :pN][1];
q0_immune_mp_scaled    = time_shift_scaled[time_shift_scaled.Dataset .== "MP", :q][1];
q0_immune_mpb1_scaled  = time_shift_scaled[time_shift_scaled.Dataset .== "MPB1", :q][1];

# Bounds
rmin_immune  = 0.0;
rmax_immune  = 5.0;
pNmin_immune = 0.00001;
pNmax_immune = 1.0;
qmin_immune  = 0.00001;
qmax_immune  = 1.0;

fit_immune_mp_scaled   = curve_fit(tumour_immune_model_mp_scaled, tdata_mp, ydata_mp_combined_scaled, [r0_immune_mp_scaled, pN0_immune_mp_scaled, q0_immune_mp_scaled], 
                           lower=[rmin_immune, pNmin_immune, qmin_immune], upper=[rmax_immune, pNmax_immune, qmax_immune])
fit_immune_mpb1_scaled = curve_fit(tumour_immune_model_mpb1_scaled, tdata_mpb1, ydata_mpb1_combined_scaled, [r0_immune_mpb1_scaled, pN0_immune_mpb1_scaled, q0_immune_mpb1_scaled], 
                             lower=[rmin_immune, pNmin_immune, qmin_immune], upper=[rmax_immune, pNmax_immune, qmax_immune])

#-------------------------------------
# Cisplatin PKPD with immune environment
#-------------------------------------

# Load immune data and proportions via DataTools
immune_cis = DataTools.load_immune_datasets_cis(data_dir = joinpath(@__DIR__, "..", "Data"));

df_immune_cis_mp       = immune_cis.df_immune_cis_mp;
df_immune_cis_mpb1     = immune_cis.df_immune_cis_mpb1;
nk_cells_prop_cis_mp   = immune_cis.nk_cells_prop_cis_mp;
nk_cells_prop_cis_mpb1 = immune_cis.nk_cells_prop_cis_mpb1;
t_cells_prop_cis_mp    = immune_cis.t_cells_prop_cis_mp;
t_cells_prop_cis_mpb1  = immune_cis.t_cells_prop_cis_mpb1;

# Cell counts for NK and T cells, using the assumptions defined in Assumptions.jl
nk_cell_count_cis_mp   = nk_cells_prop_cis_mp * prop_cd45_cells * total_cell_count;
nk_cell_count_cis_mpb1 = nk_cells_prop_cis_mpb1 * prop_cd45_cells * total_cell_count;
t_cell_count_cis_mp    = t_cells_prop_cis_mp * prop_cd8cells_in_t_cells * prop_cd45_cells * total_cell_count;
t_cell_count_cis_mpb1  = t_cells_prop_cis_mpb1 * prop_cd8cells_in_t_cells * prop_cd45_cells * total_cell_count;

ydata_mp_cis_combined   = vcat(ydata_mp_cis, [nk_cell_count_cis_mp, t_cell_count_cis_mp]);
ydata_mpb1_cis_combined = vcat(ydata_mpb1_cis, [nk_cell_count_cis_mpb1, t_cell_count_cis_mpb1]);

# Fixed bolus dosing setup
BW = 0.025;
dose_amt = 3.0 * BW;
initial_day_dose = 0.01;  # Start dosing at 0.01 days to avoid dosing at t=0 which is causing issues
dose_starts = [initial_day_dose, initial_day_dose + 7.0];  # Doses 1 week apart

# Initial conditions for immune cells and tumour growth model
u0_immune_cis_mp   = [0, 0, ydata_mp_cis[1], 0, 0, 0, N0_treatment_initiation_mp, E0_treatment_initiation_mp];
u0_immune_cis_mpb1 = [0, 0, ydata_mpb1_cis[1], 0, 0, 0, N0_treatment_initiation_mpb1, E0_treatment_initiation_mpb1];

# Define models
tumour_immune_model_mp_cis   = (t, p) -> Predictions.predict_tumour_immune_cis_pkpd(t, t_immune_mp, p, u0_immune_cis_mp, fixed_params_immune_cis_mp, dose_amt, dose_starts);
tumour_immune_model_mpb1_cis = (t, p) -> Predictions.predict_tumour_immune_cis_pkpd(t, t_immune_mpb1, p, u0_immune_cis_mpb1, fixed_params_immune_cis_mpb1, dose_amt, dose_starts);

# Initial estimates for fitted params
k2N0_mp   = 1.0;      # (L/mg*day)
k2N0_mpb1 = k2N0_mp;  # (L/mg*day)
k2E0_mp   = 1.0;      # (L/mg*day)
k2E0_mpb1 = k2E0_mp;  # (L/mg*day)

# Bounds
k2Nmin = 0.0;
k2Nmax = 500.0;
k2Emin = 0.0;
k2Emax = 500.0;

#-------------------------------------
# Cisplatin PKPD with immune environment - NORMALIZED
#-------------------------------------

# Scaling factors
tumour_scale_mp_cis = maximum(ydata_mp_cis);
nk_scale_mp_cis = maximum([nk_cell_count_cis_mp, 1.0]);
t_scale_mp_cis  = maximum([t_cell_count_cis_mp, 1.0]);

tumour_scale_mpb1_cis = maximum(ydata_mpb1_cis);
nk_scale_mpb1_cis = maximum([nk_cell_count_cis_mpb1, 1.0]);
t_scale_mpb1_cis  = maximum([t_cell_count_cis_mpb1, 1.0]);

# Scaled dataydata_mp_combined_scaled = vcat(ydata_mp_cis ./ tumour_scale_mp, [nk_cell_count_mp / nk_scale_mp, t_cell_count_mp / t_scale_mp])
ydata_mp_cis_combined_scaled   = vcat(ydata_mp_cis ./ tumour_scale_mp_cis, [nk_cell_count_cis_mp / nk_scale_mp_cis, t_cell_count_cis_mp / t_scale_mp_cis]);
ydata_mpb1_cis_combined_scaled = vcat(ydata_mpb1_cis ./ tumour_scale_mpb1_cis, [nk_cell_count_cis_mpb1 / nk_scale_mpb1_cis, t_cell_count_cis_mpb1 / t_scale_mpb1_cis]);

# Initial conditions for immune cells and tumour growth model
u0_immune_cis_mp_scaled   = [0, 0, ydata_mp_cis[1], 0, 0, 0, N0_treatment_initiation_mp_scaled, E0_treatment_initiation_mp_scaled];
u0_immune_cis_mpb1_scaled = [0, 0, ydata_mpb1_cis[1], 0, 0, 0, N0_treatment_initiation_mpb1_scaled, E0_treatment_initiation_mpb1_scaled];

# Define models
tumour_immune_model_mp_cis_scaled = (t, p) -> begin
    pred = Predictions.predict_tumour_immune_cis_pkpd(t, t_immune_mp, p, u0_immune_cis_mp_scaled, fixed_params_immune_cis_mp, dose_amt, dose_starts)
    return vcat(pred[1:length(ydata_mp_cis)] ./ tumour_scale_mp_cis,
                pred[end - 1] / nk_scale_mp_cis,
                pred[end] / t_scale_mp_cis)
end

tumour_immune_model_mpb1_cis_scaled = (t, p) -> begin
    pred = Predictions.predict_tumour_immune_cis_pkpd(t, t_immune_mpb1, p, u0_immune_cis_mpb1_scaled, fixed_params_immune_cis_mpb1, dose_amt, dose_starts)
    return vcat(pred[1:length(ydata_mpb1_cis)] ./ tumour_scale_mpb1_cis,
                pred[end - 1] / nk_scale_mpb1_cis,
                pred[end] / t_scale_mpb1_cis)
end

fit_immune_mp_cis_scaled   = curve_fit(tumour_immune_model_mp_cis_scaled, tdata_mp_cis, ydata_mp_cis_combined_scaled, [k2N0_mp, k2E0_mp], 
                           lower=[k2Nmin, k2Emin], upper=[k2Nmax, k2Emax]);
fit_immune_mpb1_cis_scaled = curve_fit(tumour_immune_model_mpb1_cis_scaled, tdata_mpb1_cis, ydata_mpb1_cis_combined_scaled, [k2N0_mpb1, k2E0_mpb1], 
                             lower=[k2Nmin, k2Emin], upper=[k2Nmax, k2Emax]);
