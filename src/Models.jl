#  tumour growth only
function tumour_exponential!(du, u, p, t)
    T = u[1]
    r = p[1]
    du[1] = r * T
end

function tumour_logistic!(du, u, p, t)
    T = u[1]
    r, k = p
    du[1] = r * T * (1 - T / k)
end

function tumour_gompertz!(du, u, p, t)
    T = u[1]
    r, k = p
    du[1] = r * T * log(k / T)
end

function tumour_cis_pkpd!(du, u, p, t)
    C1, C2, T1, T2, T3, T4 = u


    V1 = p.V1
    V2 = p.V2
    CL = p.CL
    Q  = p.Q
    r  = p.r
    k  = p.k
    k1 = p.k1
    k2 = p.k2
    k12 = Q / V1
    k21 = Q / V2
    k10 = CL / V1

    du[1] = - (k12 + k10) * C1 + k21 * C2
    du[2] = k12 * C1 - k21 * C2
    du[3] = r * T1 * (1 - T1 / k) - k2 * T1 * C1
    du[4] = k2 * T1 * C1 - k1 * T2
    du[5] = k1 * T2 - k1 * T3
    du[6] = k1 * T3 - k1 * T4
end

function tumour_immune!(du, u, p, t)
    T, N, E = u
    r   = p.r
    k   = p.k
    d  = p.d
    c   = p.c
    dN  = p.dN
    αN  = p.αN
    hN  = p.hN
    pN  = p.pN
    dE  = p.dE
    αE  = p.αE
    hE  = p.hE
    q   = p.q
    hT  = p.hT
    N0  = p.N0

    du[1] = r * T * (1 - T/k) - c*N*T - d*(E/(hT+E))*T
    du[2] = N0 - dN * N + αN * (T / (hN + T)) * N - pN * N * T
    du[3] = -dE * E + αE * (T / (hE + T)) * E - q * E * T + (1.1*pN) * N * T

end

function tumour_immune_cis_pkpd!(du, u, p, t)
    C1, C2, T1, T2, T3, T4, N, E = u
    V1 = p.V1
    V2 = p.V2
    CL = p.CL
    Q  = p.Q
    r   = p.r
    k   = p.k
    k1 = p.k1
    k2 = p.k2
    k12 = Q / V1
    k21 = Q / V2
    k10 = CL / V1

    # Are named r_C,N and r_C,E in manuscript
    k2N = p.k2N
    k2E = p.k2E

    d  = p.d
    c   = p.c
    dN  = p.dN
    αN  = p.αN
    hN  = p.hN
    pN  = p.pN
    dE  = p.dE
    αE  = p.αE
    hE  = p.hE
    q   = p.q
    hT  = p.hT
    N0  = p.N0

    du[1] = - (k12 + k10)*C1 + k21*C2
    du[2] = k12*C1 - k21 * C2

    du[3] = r*T1*(1 - T1/k) - c*N*T1 - d*(E/(hT+E))*T1 - k2*T1*C1
    du[4] = k2*T1*C1 - k1*T2
    du[5] = k1*T2 - k1*T3
    du[6] = k1*T3 - k1*T4

    # If only T1 has impact on immune cells dynamics (we consider the rest is slowly dying)
    # and we model cell recruitment differences instead of killing rate differences
    # where k2N and k2E represent the impact of cisplatin on the recruitment of NK and CD8 cells, respectively
    du[7] = N0 - dN*N + αN*(T1 / (hN + T1))*N - pN*N*T1 + (k2N - k2)*N*C1
    du[8] = -dE*E + αE*(T1 / (hE + T1))*E - q*E*T1 + (1.1*pN)*N*T1 + (k2E - k2)*E*C1

end

#ICB-style tumour + effector + PD1/PDL1
function tumour_icb!(du, u, p, t)
    T, N, E, I, PD1U, PD1B, PDL1 = u
    V1 = p.V1
    r   = p.r
    k   = p.k
    d  = p.d
    c   = p.c
    dN  = p.dN
    αN  = p.αN
    hN  = p.hN
    pN  = p.pN
    dE  = p.dE
    αE  = p.αE
    hE  = p.hE
    q   = p.q
    hT  = p.hT
    N0  = p.N0

    muI = p.muI
    muP = p.muP
    dI = p.dI
    ρE  = p.ρE
    ρT  = p.ρT
    kYQ = p.kYQ

    du[1] = r * T * (1 - T/k) - c*N*T - d*(E/(hT+E))*T
    du[2] = N0 - dN * N + αN * (T / (hN + T)) * N - pN * N * T
    du[3] = -dE * E + αE * (T / (hE + T)) * E * (1 / (1 + (PD1U * PDL1) / kYQ)) - q * E * T + (1.1*pN) * N * T
    du[4] = - muI * PD1U - dI * I
    du[5] = ρE * (αE * (T / (hE + T)) * E) * (1 / (1 + (PD1U * PDL1) / kYQ)) - (PD1U / (PD1U + PD1B)) * ρE * dE * E - muP * PD1U * I
    du[6] = muP * PD1U * I - (PD1B / (PD1U + PD1B)) * ρE * dE * E
    du[7] = ρT * (r * T * (1 - T/k) - c*N*T - d*(E/(hT+E))*T)
end
