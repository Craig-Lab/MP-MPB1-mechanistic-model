
include("ProfileLikelihood_SetUp.jl");


## Commands for likelihood profiles and 95% confidence intervals

#-------------------------------------
# Tumour growth only
#-------------------------------------

# MP
pminfold_MP = [0.35, 0.5];
pmaxfold_MP = [0.35, 1.75];

models_MP = [
(name="Logistic MP", model=tumour_model_logistic_mp, t=tdata_mp, y=ydata_mp, p0=[r0,k0], lower=[rmin,kmin], upper=[rmax,kmax], param_names=["r","K"], param_labels=["Growth rate (day⁻¹)","Carrying capacity (cells)"]),
(name="Logistic MP + CIS", model=tumour_model_logistic_mp_cis, t=tdata_mp_cis, y=ydata_mp_cis, p0=[r0,k0], lower=[rmin,kmin], upper=[rmax,kmax], param_names=["r","K"], param_labels=["Growth rate (day⁻¹)","Carrying capacity (cells)"]),
(name="Logistic MP + ICB", model=tumour_model_logistic_mp_icb, t=tdata_mp_icb, y=ydata_mp_icb, p0=[r0,k0], lower=[rmin,kmin], upper=[rmax,kmax], param_names=["r","K"], param_labels=["Growth rate (day⁻¹)","Carrying capacity (cells)"]),
(name="Logistic MP + CIS + ICB", model=tumour_model_logistic_mp_cis_icb, t=tdata_mp_cis_icb, y=ydata_mp_cis_icb, p0=[r0,k0], lower=[rmin,kmin], upper=[rmax,kmax], param_names=["r","K"], param_labels=["Growth rate (day⁻¹)","Carrying capacity (cells)"])
];


subtitle_MP = ["No treatment", "No treatment", "Cisplatin", "Cisplatin", "ICB", "ICB", "Cisplatin + ICB", "Cisplatin + ICB"];
filename_MP = "PLL MP Logistic";

CI95_MP = plot_profile_likelihood(
    models_MP, pminfold_MP, pmaxfold_MP;
    subtitles=subtitle_MP,
    filename=filename_MP,
    template=craig_lab_template
);

CSV.write("Confidence Intervals/95% CI - MP growth.csv", CI95_MP);

# MPB1
pminfold_MPB1 = [0.35, 0.5];
pmaxfold_MPB1 = [0.35, 1.75];

models_MPB1 = [
(name="Logistic MPB1", model=tumour_model_logistic_mpb1, t=tdata_mpb1, y=ydata_mpb1, p0=[r0,k0], lower=[rmin,kmin], upper=[rmax,kmax], param_names=["r","K"], param_labels=["Growth rate (day⁻¹)","Carrying capacity (cells)"]),
(name="Logistic MPB1 + CIS", model=tumour_model_logistic_mpb1_cis, t=tdata_mpb1_cis, y=ydata_mpb1_cis, p0=[r0,k0], lower=[rmin,kmin], upper=[rmax,kmax], param_names=["r","K"], param_labels=["Growth rate (day⁻¹)","Carrying capacity (cells)"]),
(name="Logistic MPB1 + ICB", model=tumour_model_logistic_mpb1_icb, t=tdata_mpb1_icb, y=ydata_mpb1_icb, p0=[r0,k0], lower=[rmin,kmin], upper=[rmax,kmax], param_names=["r","K"], param_labels=["Growth rate (day⁻¹)","Carrying capacity (cells)"]),
(name="Logistic MPB1 + CIS + ICB", model=tumour_model_logistic_mpb1_cis_icb, t=tdata_mpb1_cis_icb, y=ydata_mpb1_cis_icb, p0=[r0,k0], lower=[rmin,kmin], upper=[rmax,kmax], param_names=["r","K"], param_labels=["Growth rate (day⁻¹)","Carrying capacity (cells)"])
];

subtitle_MPB1 = ["No treatment", "No treatment", "Cisplatin", "Cisplatin", "ICB", "ICB", "Cisplatin + ICB", "Cisplatin + ICB"];
filename_MPB1 = "PLL MPB1 Logistic";

CI95_MPB1 = plot_profile_likelihood(
    models_MPB1, pminfold_MPB1, pmaxfold_MPB1;
    subtitles=subtitle_MPB1,
    filename=filename_MPB1,
    template=craig_lab_template
);

CSV.write("Confidence Intervals/95% CI - MPB1 growth.csv", CI95_MPB1);

#-------------------------------------
# Tumour growth with immune environment
#-------------------------------------

pminfold_immune = [0.2, 0.4, 0.4];
pmaxfold_immune = [0.2, 1.0, 1.0];

models_immune = [
(
    name = "Immune MP", 
    model = tumour_immune_model_mp_scaled, 
    t = tdata_mp, y = ydata_mp_combined_scaled, 
    p0 = [r0_immune_mp,pN0_immune_mp,q0_immune_mp], 
    lower = [rmin_immune,pNmin_immune,qmin_immune], 
    upper = [rmax_immune,pNmax_immune,qmax_immune],
    param_names=["r","p","q"], 
    param_labels=["Growth rate (day⁻¹)","NK inactivation rate (mm<sup>-3</sup>day⁻¹)","CD8+ T inactivation rate (mm<sup>-3</sup>day⁻¹)"]
),
(
    name = "Immune MPB1", 
    model = tumour_immune_model_mpb1_scaled, 
    t = tdata_mpb1, y = ydata_mpb1_combined_scaled, 
    p0 = [r0_immune_mpb1,pN0_immune_mpb1,q0_immune_mpb1], 
    lower = [rmin_immune,pNmin_immune,qmin_immune], 
    upper = [rmax_immune,pNmax_immune,qmax_immune],
    param_names=["r","p","q"], 
    param_labels=["Growth rate (day⁻¹)","Rate of NK cell inactivation (mm<sup>-3</sup>day⁻¹)","Rate of CD8+T cells inactivation (mm<sup>-3</sup>day⁻¹)"]
)
];

subtitle_immune = ["Growth rate (r)", "NK inactivation rate (p)", "CD8+ T inactivation rate (q)", "Growth rate (r)", "NK inactivation rate (p)", "CD8+ T inactivation rate (q)"];
filename_immune = "PLL Immune";

CI95_immune = plot_profile_likelihood(
    models_immune, pminfold_immune, pmaxfold_immune;
    subtitles=subtitle_immune,
    filename=filename_immune,
    template=craig_lab_template
);

CSV.write("Confidence Intervals/95% CI - Immune extended.csv", CI95_immune);

#-------------------------------------
# Cisplatin PKPD
#-------------------------------------

pminfold_PKPD = [0.175, 0.25];
pmaxfold_PKPD = [0.175, 0.175];

models_PKPD = [
(
    name = "Cis PKPD MP", 
    model = tumour_immune_model_mp_cis_scaled, 
    t = tdata_mp_cis, y = ydata_mp_cis_combined_scaled, 
    p0 = [k2N0_mp,k2E0_mp], 
    lower = [k2Nmin,k2Emin], 
    upper = [k2Nmax,k2Emax],
    param_names=["r1","r2"], 
    param_labels=["NK recruitement rate (L⋅mg⁻¹day⁻¹)","CD8+ T recruitement rate (L⋅mg⁻¹day⁻¹)"]
),
(
    name = "Cis PKPD MPB1", 
    model = tumour_immune_model_mpb1_cis_scaled, 
    t = tdata_mpb1_cis, y = ydata_mpb1_cis_combined_scaled, 
    p0 = [k2N0_mpb1,k2E0_mpb1], 
    lower = [k2Nmin,k2Emin], 
    upper = [k2Nmax,k2Emax],
    param_names=["r1","r2"], 
    param_labels=["NK recruitement rate (L⋅mg⁻¹day⁻¹)","CD8+ T recruitement rate (L⋅mg⁻¹day⁻¹)"]
)
];

subtitle_PKPD = ["NK recruitement rate", "CD8+ T recruitement rate", "NK recruitement rate", "CD8+ T recruitement rate"];
filename_PKPD = "PLL PKPD";

CI95_PKPD = plot_profile_likelihood(
    models_PKPD, pminfold_PKPD, pmaxfold_PKPD;
    subtitles=subtitle_PKPD,
    filename=filename_PKPD,
    template=craig_lab_template
);

CSV.write("Confidence Intervals/95% CI - PKPD Cisplatin.csv", CI95_PKPD);

println("All likelihood profiles generated")
