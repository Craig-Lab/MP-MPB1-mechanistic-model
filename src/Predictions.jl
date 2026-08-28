module Predictions

include("Models.jl")
using DifferentialEquations
using DiffEqCallbacks

export predict_exponential, predict_logistic, predict_gompertz
export predict_cis_pkpd, predict_tumour_immune_cis_pkpd
export predict_tumour_immune, predict_tumour_immune_states, predict_tumour_immune_cis_pkpd_states, predict_tumour_immune_icb_q_states

# predict_exponential(t, p, u0) -> Vector predictions at times t
function predict_exponential(t, p, u0)
    prob = ODEProblem(tumour_exponential!, u0, (t[1], t[end]), p)
    sol = solve(prob, Tsit5(); saveat=t)
    return sol[1, :]
end


function predict_logistic(t, p, u0)
    r, k = p
    prob = ODEProblem(tumour_logistic!, u0, (t[1], t[end]), (r, k))
    sol = solve(prob, Tsit5(); saveat=t)
    return sol[1, :]
end

function predict_gompertz(t, p, u0)
    r, k = p
    prob = ODEProblem(tumour_gompertz!, u0, (t[1], t[end]), (r, k))
    sol = solve(prob, Tsit5(); saveat=t)
    return sol[1, :]
end

function predict_cis_pkpd(t, p_fit, u0, p_fixed, dose_amt, dose_starts)
    k2 = p_fit[1]

    p = (;
        p_fixed.V1,
        p_fixed.V2,
        p_fixed.CL,
        p_fixed.Q,
        p_fixed.r,
        p_fixed.k,
        p_fixed.k1,
        k2,
    )

    function bolus!(integrator)
        integrator.u[1] += dose_amt / p_fixed.V1
    end

    cb = PresetTimeCallback(dose_starts, bolus!)

    prob = ODEProblem(tumour_cis_pkpd!, u0, (t[1], t[end]), p)

    # Do not use saveat=t here because callback times can add extra saved points.
    sol = solve(prob, Tsit5(); callback=cb)

    # Sample exactly at the observed data times so output length matches ydata.
    return [sum(sol(ti)[3:6]) for ti in t]
end


function predict_tumour_immune(t_tumour, t_immune, p, u0, fixed_p)
    r, pN, q = p
    p_full = (; fixed_p..., r = r, pN = pN, q = q)
    # Combine all times needed and sort
    t_all = vcat(t_tumour, [t_immune])
    t_unique = sort(unique(t_all))
    t_max = maximum(t_unique)
    prob = ODEProblem(tumour_immune!, u0, (0.0, t_max), p_full)
    sol = solve(prob, Tsit5(); saveat=t_unique)
    # Extract T values at tumour observation times
    T_vals = sol(t_tumour)[1, :]
    # Extract N and E values at immune cell observation time
    N_val = sol(t_immune)[2]
    E_val = sol(t_immune)[3]
    # Return combined vector [T_obs_1, T_obs_2, ..., N_obs, E_obs]
    return vcat(T_vals, [N_val, E_val])
end

function predict_tumour_immune_states(t, p, u0, fixed_p)
    r, pN, q = p
    p_full = (; fixed_p..., r = r, pN = pN, q = q)
    prob = ODEProblem(tumour_immune!, u0, (0.0, maximum(t)), p_full)
    sol = solve(prob, Tsit5(); saveat=t)
    return sol[1:3, :]
end

function predict_tumour_immune_cis_pkpd(t_tumour, t_immune, p, u0, fixed_p, dose_amt, dose_starts)
    # k2, k2N, k2E = p
    # p_full = (; fixed_p..., k2 = k2, k2N = k2N, k2E = k2E, infusions = Tuple{Float64, Float64, Float64}[])

    k2N, k2E = p
    p_full = (; fixed_p..., k2N = k2N, k2E = k2E, infusions = Tuple{Float64, Float64, Float64}[])

    function bolus!(integrator)
        integrator.u[1] += dose_amt / p_full.V1
    end

    cb = PresetTimeCallback(dose_starts, bolus!)

    # Combine all times needed and sort
    t_all = vcat(t_tumour, [t_immune])
    t_unique = sort(unique(t_all))
    t_max = maximum(t_unique)

    prob = ODEProblem(tumour_immune_cis_pkpd!, u0, (0.0, t_max), p_full)
    sol = solve(prob, Tsit5(); callback=cb, saveat=t_unique)

    # Extract T values at tumour observation times
    T_vals = [sum(sol(ti)[3:6]) for ti in t_tumour]
    # Extract N and E values at immune cell observation time
    N_val = sol(t_immune)[7]
    E_val = sol(t_immune)[8]
    # Return combined vector [T_obs_1, T_obs_2, ..., N_obs, E_obs]
    return vcat(T_vals, [N_val, E_val])
end

function predict_tumour_immune_cis_pkpd_states(t, params, u0, fixed_p, dose_amt, dose_starts)
    k2, k2N, k2E = params
    p_full = (; fixed_p..., k2 = k2, k2N = k2N, k2E = k2E, infusions = Tuple{Float64, Float64, Float64}[])
    
    # Define bolus callback
    function bolus!(integrator)
        integrator.u[1] += dose_amt / p_full.V1
    end
    cb = PresetTimeCallback(dose_starts, bolus!)
    
    # Always solve from t=0 to the maximum time in the data
    t_max = maximum(t)
    prob = ODEProblem(tumour_immune_cis_pkpd!, u0, (0.0, t_max), p_full)
    sol = solve(prob, Tsit5(); callback=cb, saveat=t)
    
    # Extract tumor (sum of T1, T2, T3, T4), NK cells (N), and effector T cells (E)
    T_total = [sum(sol(ti)[3:6]) for ti in t]  # Sum T1+T2+T3+T4 at each time
    N_cells = [sol(ti)[7] for ti in t]  # NK cells at each time
    E_cells = [sol(ti)[8] for ti in t]  # Effector T cells at each time
    return vcat(T_total', N_cells', E_cells')  # Return as 3 x length(t) matrix
end


function predict_tumour_immune_icb_q_states(t, params, u0, fixed_p)
    q = params[1]
    p_full = (; fixed_p..., q = q)
    # Always solve from t=0 to the maximum time in the data
    t_max = maximum(t)
    prob = ODEProblem(tumour_immune!, u0, (0.0, t_max), p_full)
    sol = solve(prob, Tsit5(); saveat=t)
    return sol[1:3, :]  # Return T, N, E at the requested times
end

end # module
