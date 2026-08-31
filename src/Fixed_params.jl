
# Logistic tumour growth and PKPD model of cisplatin
fixed_params_cis_mp = (;
    V1 = 0.0786, # L, from mice and rat data, Fukushima et al. 2016
    V2 = 0.119,  # L, from mice and rat data, Fukushima et al. 2016
    CL = 5.56,   # L/day, from mice and rat data, Fukushima et al. 2016
    Q  = 1.17,   # L/day, from mice and rat data, Fukushima et al. 2016
    r  = 0.305,  # (1/days), from prior fitting of logistic growth to MP vehicle data
    k  = 1592.8, # (mm^3), from prior fitting of logistic growth to MP vehicle data
    k1 = 2.16   # (1/days), from mice and rat data, Fukushima et al. 2016 (before was k1 = 1.64,   # (1/days), from rat data, Cassia Pigatto et al.)
)

fixed_params_cis_mpb1 = (;
    V1 = 0.0786, # L, from mice and rat data, Fukushima et al. 2016
    V2 = 0.119,  # L, from mice and rat data, Fukushima et al. 2016
    CL = 5.56,   # L/day, from mice and rat data, Fukushima et al. 2016
    Q  = 1.17,   # L/day, from mice and rat data, Fukushima et al. 2016
    r  = 0.423,  # (1/days), from prior fitting of logistic growth to MPB1 vehicle data
    k  = 2137.0, # (mm^3), from prior fitting of logistic growth to MPB1 vehicle data
    k1 = 2.16   # (1/days), from mice and rat data, Fukushima et al. 2016 (before was k1 = 1.64,   # (1/days), from rat data, Cassia Pigatto et al.)
)


# Tumour growth and immune cells model, for fitting from tumour initiation
fixed_params_immune_mp = (;
    k  = 1592.8,    # (mm^3), from prior fitting
    αN = 2.5e-2,    # (1/days), de Pillis et al. 2005
    αE = 3.75e-2,   # (1/days), de Pillis et al. 2005
    c = 3.50e-6,    # (1/cell*day), de Pillis et al. 2005
    dN = 4.12e-2,   # (1/days), de Pillis et al. 2005
    dE = 2.0e-2,    # (1/days), de Pillis et al. 2005
    d = 1.43,       # (1/days), de Pillis et al. 2005
    hN = 4.49e-3,   # (mm^3), adapted from de Pillis et al. 2005
    hE = 4.49e-3,   # (mm^3), adapted from de Pillis et al. 2005
    hT = 1.18e6,    # (cells), adapted from Mongeon and Craig 2025
    N0 = 1.3e4,     # (cells/day), de Pillis et al. 2005
    #uNE = 0.11,     # (1/mm3*day), adapted from de Pillis et al. 2005
)

fixed_params_immune_mpb1 = (;
    k  = 2137.0,    # (mm^3), from prior fitting
    αN = 2.5e-2,    # (1/days), de Pillis et al. 2005
    αE = 3.75e-2,   # (1/days), de Pillis et al. 2005
    c = 3.50e-6,    # (1/cell*day), de Pillis et al. 2005
    dN = 4.12e-2,   # (1/days), de Pillis et al. 2005
    dE = 2.0e-2,    # (1/days), de Pillis et al. 2005
    d = 1.43,       # (1/days), de Pillis et al. 2005
    hN = 4.49e-3,   # (mm^3), adapted from de Pillis et al. 2005
    hE = 4.49e-3,   # (mm^3), adapted from de Pillis et al. 2005
    hT = 1.18e6,    # (cells), adapted from Mongeon and Craig 2025
    N0 = 1.3e4,     # (cells/day), de Pillis et al. 2005
    #uNE = 0.11,     # (1/mm3*day), adapted from de Pillis et al. 2005
)

fixed_params_immune_cis_mp = (;
    V1 = 0.0786,    # L, from mice and rat data, Fukushima et al. 2016
    V2 = 0.119,     # L, from mice and rat data, Fukushima et al. 2016
    CL = 5.56,      # L/day, from mice and rat data, Fukushima et al. 2016
    Q  = 1.17,      # L/day, from mice and rat data, Fukushima et al. 2016
    k  = 1592.8,    # (mm^3), from prior fitting
    k1 = 2.16,       # (1/days), from mice and rat data, Fukushima et al. 2016
    k2 = 8.07,      # (L/mg*day), from prior fitting of logistic cisplatin PK/PD model
    αN = 2.5e-2,    # (1/days), de Pillis et al. 2005
    αE = 3.75e-2,   # (1/days), de Pillis et al. 2005
    c = 3.50e-6,    # (1/cell*day), de Pillis et al. 2005
    dN = 4.12e-2,   # (1/days), de Pillis et al. 2005
    dE = 2.0e-2,    # (1/days), de Pillis et al. 2005
    d = 1.43,       # (1/days), de Pillis et al. 2005
    hN = 4.49e-3,   # (mm^3), adapted from de Pillis et al. 2005
    hE = 4.49e-3,   # (mm^3), adapted from de Pillis et al. 2005
    hT = 1.18e6,    # (cells), adapted from Mongeon and Craig 2025
    N0 = 1.3e4,     # (cells/day), de Pillis et al. 2005
    r  = 1.57,      # (1/days), from prior fitting of tumour-immune cells model to vehicle
    pN = 7.49e-4,   # (1/mm3*days), from prior fitting of tumour-immune cells model to vehicle
    q = 2.00e-3,    # (1/mm3*days), from prior fitting of tumour-immune cells model to vehicle
)

fixed_params_immune_cis_mpb1 = (;
    V1 = 0.0786,    # L, from mice and rat data, Fukushima et al. 2016
    V2 = 0.119,     # L, from mice and rat data, Fukushima et al. 2016
    CL = 5.56,      # L/day, from mice and rat data, Fukushima et al. 2016
    Q  = 1.17,      # L/day, from mice and rat data, Fukushima et al. 2016
    k  = 2137.0,    # (mm^3), from prior fitting
    k1 = 2.16,       # (1/days), from mice and rat data, Fukushima et al. 2016 
    k2 = 63.7,      # (L/mg*day), from prior fitting of logistic cisplatin PK/PD model
    αN = 2.5e-2,    # (1/days), de Pillis et al. 2005
    αE = 3.75e-2,   # (1/days), de Pillis et al. 2005
    c = 3.50e-6,    # (1/cell*day), de Pillis et al. 2005
    dN = 4.12e-2,   # (1/days), de Pillis et al. 2005
    dE = 2.0e-2,    # (1/days), de Pillis et al. 2005
    d = 1.43,       # (1/days), de Pillis et al. 2005
    hN = 4.49e-3,   # (mm^3), adapted from de Pillis et al. 2005
    hE = 4.49e-3,   # (mm^3), adapted from de Pillis et al. 2005
    hT = 1.18e6,    # (cells), adapted from Mongeon and Craig 2025
    N0 = 1.3e4,     # (cells/day), de Pillis et al. 2005
    r  = 2.03,      # (1/days), from prior fitting of tumour-immune cells model to vehicle
    pN = 2.62e-4,   # (1/mm3*days), from prior fitting of tumour-immune cells model to vehicle
    q = 7.74e-4,    # (1/mm3*days), from prior fitting of tumour-immune cells model to vehicle
)

fixed_params_mp_icb = (;
    k  = 1592.8,    # (mm^3), from prior fitting
    αN = 2.5e-2,    # (1/days), de Pillis et al. 2005
    αE = 3.75e-2,   # (1/days), de Pillis et al. 2005
    c = 3.50e-6,    # (1/cell*day), de Pillis et al. 2005
    dN = 4.12e-2,   # (1/days), de Pillis et al. 2005
    dE = 2.0e-2,    # (1/days), de Pillis et al. 2005
    d = 1.43,       # (1/days), de Pillis et al. 2005
    hN = 4.49e-3,   # (mm^3), adapted from de Pillis et al. 2005
    hE = 4.49e-3,   # (mm^3), adapted from de Pillis et al. 2005
    hT = 1.18e6,    # (cells), adapted from Mongeon and Craig 2025
    N0 = 1.3e4,     # (cells/day), de Pillis et al. 2005
    #uNE = 0.11,     # (1/mm3*day), adapted from de Pillis et al. 2005
    r  = 1.57,      # (1/days), from prior fitting of tumour-immune cells model to vehicle
    pN = 7.49e-4,   # (1/mm3*days), from prior fitting of tumour-immune cells model to vehicle
)

fixed_params_mpb1_icb = (;
    k  = 2137.0,    # (mm^3), from prior fitting
    αN = 2.5e-2,    # (1/days), de Pillis et al. 2005
    αE = 3.75e-2,   # (1/days), de Pillis et al. 2005
    c = 3.50e-6,    # (1/cell*day), de Pillis et al. 2005
    dN = 4.12e-2,   # (1/days), de Pillis et al. 2005
    dE = 2.0e-2,    # (1/days), de Pillis et al. 2005
    d = 1.43,       # (1/days), de Pillis et al. 2005
    hN = 4.49e-3,   # (mm^3), adapted from de Pillis et al. 2005
    hE = 4.49e-3,   # (mm^3), adapted from de Pillis et al. 2005
    hT = 1.18e6,    # (cells), adapted from Mongeon and Craig 2025
    N0 = 1.3e4,     # (cells/day), de Pillis et al. 2005
    #uNE = 0.11,     # (1/mm3*day), adapted from de Pillis et al. 2005
    r  = 2.03,      # (1/days), from prior fitting of tumour-immune cells model to vehicle
    pN = 2.62e-4,   # (1/mm3*days), from prior fitting of tumour-immune cells model to vehicle
)

fixed_params_icb_cis_mp = (;
    V1 = 0.0786,    # L, from mice and rat data, Fukushima et al. 2016
    V2 = 0.119,     # L, from mice and rat data, Fukushima et al. 2016
    CL = 5.56,      # L/day, from mice and rat data, Fukushima et al. 2016
    Q  = 1.17,      # L/day, from mice and rat data, Fukushima et al. 2016
    k1 = 2.16,       # (1/days), from mice and rat data, Fukushima et al. 2016 
    k2 = 8.07,      # (L/mg*day), from prior fitting of logistic cisplatin PK/PD model
    k  = 1592.8,    # (mm^3), from prior fitting
    αN = 2.5e-2,    # (1/days), de Pillis et al. 2005
    αE = 3.75e-2,   # (1/days), de Pillis et al. 2005
    c = 3.50e-6,    # (1/cell*day), de Pillis et al. 2005
    dN = 4.12e-2,   # (1/days), de Pillis et al. 2005
    dE = 2.0e-2,    # (1/days), de Pillis et al. 2005
    d = 1.43,       # (1/days), de Pillis et al. 2005
    hN = 4.49e-3,   # (mm^3), adapted from de Pillis et al. 2005
    hE = 4.49e-3,   # (mm^3), adapted from de Pillis et al. 2005
    hT = 1.18e6,    # (cells), adapted from Mongeon and Craig 2025
    N0 = 1.3e4,     # (cells/day), de Pillis et al. 2005
    #uNE = 0.11,     # (1/mm3*day), adapted from de Pillis et al. 2005
    r  = 1.57,      # (1/days), from prior fitting of tumour-immune cells model to vehicle
    pN = 7.49e-4,   # (1/mm3*days), from prior fitting of tumour-immune cells model to vehicle
    q = 2.14e-3,    # (1/mm3*days), from prior fitting of tumour-immune cells model to ICB
)

fixed_params_icb_cis_mpb1 = (;
    V1 = 0.0786,    # L, from mice and rat data, Fukushima et al. 2016
    V2 = 0.119,     # L, from mice and rat data, Fukushima et al. 2016
    CL = 5.56,      # L/day, from mice and rat data, Fukushima et al. 2016
    Q  = 1.17,      # L/day, from mice and rat data, Fukushima et al. 2016
    k1 = 2.16,       # (1/days), from mice and rat data, Fukushima et al. 2016 
    k2 = 63.7,      # (L/mg*day), from prior fitting of logistic cisplatin PK/PD model
    k  = 2137.0,    # (mm^3), from prior fitting
    αN = 2.5e-2,    # (1/days), de Pillis et al. 2005
    αE = 3.75e-2,   # (1/days), de Pillis et al. 2005
    c = 3.50e-6,    # (1/cell*day), de Pillis et al. 2005
    dN = 4.12e-2,   # (1/days), de Pillis et al. 2005
    dE = 2.0e-2,    # (1/days), de Pillis et al. 2005
    d = 1.43,       # (1/days), de Pillis et al. 2005
    hN = 4.49e-3,   # (mm^3), adapted from de Pillis et al. 2005
    hE = 4.49e-3,   # (mm^3), adapted from de Pillis et al. 2005
    hT = 1.18e6,    # (cells), adapted from Mongeon and Craig 2025
    N0 = 1.3e4,     # (cells/day), de Pillis et al. 2005
    #uNE = 0.11,     # (1/mm3*day), adapted from de Pillis et al. 2005
    r  = 2.03,      # (1/days), from prior fitting of tumour-immune cells model to vehicle
    pN = 2.62e-4,   # (1/mm3*days), from prior fitting of tumour-immune cells model to vehicle
    q = 1.41e-4,    # (1/mm3*days), from prior fitting of tumour-immune cells model to ICB
)