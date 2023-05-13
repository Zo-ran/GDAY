a0rhizo <- 0.05                         # minimum allocation to rhizodeposition [0.0-0.1]
a1rhizo <- 0.6                          # slope of allocation to rhizodeposition [0.2-1]
actncmax <- 0.333333                    # Active pool (=1/3) N:C ratio of new SOM - maximum [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology.
actncmin <- 0.066667                    # Active pool (=1/15) N:C of new SOM - when Nmin=Nmin0 [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology.
actpcmax <- 0.333333                    # Active pool (=1/30) P:C ratio of new SOM - maximum [units: gP/gC]. Based on forest version of CENTURY (Parton et al. 1993)
actpcmin <- 0.066667                    # Active pool (=1/80) P:C of new SOM - when Pmin=Pmin0 [units: gP/gC]. Based on forest version of CENTURY (Parton et al. 1993)
ageold <- 10000.0                       # Plant age when max leaf N C ratio is lowest
ageyoung <- 0.0                         # Plant age when max leaf N C ratio is highest
albedo <- 0.2
alpha_c4 <- 0.06                        # quantium efficiency for C4 plants has no Ci and temp dependancy, so if a fixed constant.
alpha_j <- 0.308                        # initial slope of rate of electron transport, used in calculation of quantum yield. Value calculated by Belinda
b_root <- -999.9
b_topsoil <- -999.9
bdecay <- 0.1                           # branch and large root turnover rate (1/yr)
biochemical_p_constant <- 150.0         # Michaelis-Menton constant for biochemical P mineralisation [g N (g P)-1]; Wang et al., 2007, GB1018
branch0 <- 5.61                         # constant in branch-stem allometry (trees)
branch1 <- 0.346                        # exponent in branch-stem allometry
bretrans <- 0.0                         # branch n retranslocation fraction
c_alloc_bmax <- 0.1                     # allocation to branches at branch n_crit and p_crit. If using allometric model this is the max alloc to branches
c_alloc_bmin <- 0.1                     # allocation to branches at zero branch n/c and p/c. If using allometric model this is the min alloc to branches
c_alloc_cmax <- 0.0                     # allocation to coarse roots at n_crit and p_crit. If using allometric model this is the max alloc to coarse roots
c_alloc_fmax <- 0.45                    # allocation to leaves at leaf n_crit and p_crit. If using allometric model this is the max alloc to leaves
c_alloc_fmin <- 0.05                    # allocation to leaves at zero leaf n/c and p/c. If using allometric model this is the min alloc to leaves
c_alloc_rmax <- 0.45                    # allocation to roots at root n_crit and p_crit. If using allometric model this is the max alloc to fine roots
c_alloc_rmin <- 0.05                    # allocation to roots at zero root n/c and p/c. If using allometric model this is the min alloc to fine roots
cfracts <- 0.5                          # carbon fraction of dry biomass
crdecay <- 0.00                         # coarse roots turnover rate (1/yr)
cretrans <- 0.0                         # coarse root n retranslocation fraction
croot0 <- 0.34                          # constant in coarse_root-stem allometry (trees)
croot1 <- 0.84                          # exponent in coarse_root-stem allometry
crit_n_cost_of_p <- 15.0                # Critical value of N cost of root P uptake above which phosphatase production starts [g N (g P)-1]; Wang et al., 2007, GB1018
ctheta_root <- 0.525                    # Fitted parameter based on Landsberg and Waring
ctheta_topsoil <- 0.65                  # Fitted parameter based on Landsberg and Waring
cue <- 0.5                              # carbon use efficiency, or the ratio of NPP to GPP
d0 <- 0.0
d0x <- 0.35                             # Length scale for exponential decline of Umax(z)
d1 <- 0.0
delsj <- 644.4338                       # Deactivation energy for electron transport (J mol-1 k-1)
density <- 800.0                        # sapwood density kg DM m-3 (trees)
direct_frac <- 0.5                      # direct beam fraction of incident radiation - this is only used with the BEWDY model
displace_ratio <- 0.75                  # Value for coniferous forest (0.78) from Jarvis et al 1976, taken from Jones 1992 pg 67. More standard assumption is 2/3
disturbance_doy <- 1
dz0v_dh <- 0.05                         # Rate of change of vegetation roughness length for momentum with height. Value from Jarvis? for conifer 0.075
eac <- 79430.0                          # Activation energy for carboxylation [J mol-1]
eag <- 37830.0                          # Activation energy at CO2 compensation point [J mol-1]
eaj <- 43790.0                          # Activation energy for electron transport (J mol-1)
eao <- 36380.0                          # Activation energy for oxygenation [J mol-1]
eav <- 51560.0                          # Activation energy for Rubisco (J mol-1)
edj <- 200000.0                         # Deactivation energy for electron transport (J mol-1)
faecescn <- 25.0                        # Faeces C:N ratio
faecesn <- 0.0                          # Faeces N content
faecescp <- 300.0                       # Faeces C:P ratio
faecesp <- 0.0                          # Faeces P content
pfdecay <- 0.6                          # foliage turnover rate (1/yr)
fdecaydry <- 0.6                        # Foliage turnover rate - dry soil (1/yr)
fhw <- 0.8                              # n:c ratio of stemwood - immobile pool and new ring
finesoil <- 0.2                         # clay+silt fraction
fmleaf <- 0.0
fmroot <- 0.0
fmfaeces <- 0.0
fix_lai <- -999.9                       # value to fix LAI to, control fixed_lai flag must be set
fracfaeces <- 0.3                       # Fractn of grazd C that ends up in faeces (0..1)
fracteaten <- 0.5                       # Fractn of leaf production eaten by grazers
fractosoil <- 0.85                      # Fractn of grazed N recycled to soil:faeces+urine
fractosoilp <- 0.85                     # Fractn of grazed P recycled to soil:faeces+urine
fractup_soil <- 0.5                     # fraction of uptake from top soil layer
fretrans <- 0.5                         # foliage n retranslocation fraction - 46-57% in young E. globulus trees - see Corbeels et al 2005 ecological modelling 187, pg 463. Roughly 50% from review Aerts '96
fretransp <- 0.5                        # foliage p retranslocation fraction - 39.5-69 in Southern US FACE site - Finzi et al. 2001 Ecology
g1 <- 3.8667                            # stomatal conductance parameter: Slope of reln btw gs and assimilation (fitted by species/pft).
gamstar25 <- 42.75                      # Base rate of CO2 compensation point at 25 deg C [umol mol-1]
growth_efficiency <- 0.7                # growth efficiency (yg) - used only in Bewdy
height0 <- 5.0                          # Height when leaf:sap area ratio = leafsap0 (trees)
height1 <- 25.0                         # Height when leaf:sap area ratio = leafsap1 (trees)
heighto <- 4.826                        # constant in avg tree height (m) - stem (t C/ha) reln
htpower <- 0.35                         # Exponent in avg tree height (m) - stem (t C/ha) reln
intercep_frac <- 0.15                   # Maximum intercepted fraction, values in Oishi et al 2008, AFM, 148, 1719-1732 ~13.9% +/- 4.1, so going to assume 15% following Landsberg and Sands 2011, pg. 193.
jmax <- -999.9                          # maximum rate of electron transport (umol m-2 s-1)
jmaxna <- 49.930                        # slope of the reln btween jmax and leaf N content, units = (umol [gN]-1 s-1) # And for Vcmax-N slopes (vcmaxna) see Table 8.2 in CLM4_tech_note, Oleson et al. 2010.
jmaxnb <- 0.0                           # intercept of jmax vs n, units = (umol [gN]-1 s-1) # And for Vcmax-N slopes (vcmaxna) see Table 8.2 in CLM4_tech_note, Oleson et al. 2010.
jmaxpa <- 933.90                        # slope of the reln btween jmax and leaf P content, units = (umol [gP]-1 s-1)
jmaxpb <- 0.0                           # intercept of jmax vs p, units = (umol [gP]-1 s-1)
jv_intercept <- -999.9                  # Jmax to Vcmax intercept
jv_slope <- -999.9                      # Jmax to Vcmax slope
kc25 <- 404.9                           # Base rate for carboxylation by Rubisco at 25degC [mmol mol-1]
kdec1 <- 3.965571                       # surface structural decay rate (1/yr)
kdec2 <- 14.61                          # surface metabolic decay rate (1/yr)
kdec3 <- 4.904786                       # soil structural decay rate (1/yr)
kdec4 <- 18.262499                      # soil metabolic decay rate(1/yr)
kdec5 <- 7.305                          # active pool decay rate (1/yr)
kdec6 <- 0.198279                       # slow pool decay rate (1/yr)
kdec7 <- 0.006783                       # passive pool decay rate (1/yr)
kext <- 0.5                             # extinction coefficient for light
kn <- 0.3                               # extinction coefficient of nitrogen in the canopy, assumed to be 0.3 by default which comes half from Belinda's head and is supported by fig 10 in Lloyd et al. Biogeosciences, 7, 1833–1859, 2010
kp <- 0.3                               # extinction coefficient of phosphorus in the canopy
ko25 <- 278400.0                        # Base rate for oxygenation by Rubisco at 25degC [umol mol-1]. Note value in Bernacchie 2001 is in mmol!!
kq10 <- 0.08                            # exponential coefficient for Rm vs T
kr <- 0.5                               # N uptake coefficent (0.05 kg C m-2 to 0.5 tonnes/ha) see Silvia's PhD, Dewar and McM, 96.
krp <- 0.00001                          # P uptake coefficent
ks <- 0.5                               # an empirical constant [t P ha-1] - sorption capacity increased with the age of the substrate
lai_closed <- 0.5                       # Leaf angle distribution: 0 = spherical leaf angle distribution; 1 = horizontal leaves; -1 = vertical leaves
latitude <- -33.61                      # LAI of closed canopy (max cover fraction is reached (m2 (leaf) m-2 (ground) ~ 2.5)
leaf_width <- 0.001                     # latitude (degrees, negative for south)
leafsap0 <- 4000.0                      # longitude (degrees, negative for west)
leafsap1 <- 2700.0                      # leaf area  to sapwood cross sectional area ratio when Height = Height0 (mm^2/mm^2)
ligfaeces <- 0.25                       # leaf to sap area ratio when Height = Height1 (mm^2/mm^2)
ligroot <- 0.22                         # Faeces lignin as fraction of biomass
ligshoot <- 0.18                        # lignin-to-biomass ratio in root litter; Values from White et al. = 0.22  - Value in Smith et al. 2013 = 0.16, note subtly difference in eqn C9.
liteffnc <- 0.0
max_intercep_lai <- 3.0                 # canopy LAI at which interception is maximised.
max_p_biochemical <- 0.001              # max rate of biochemical P mineralisation [g P m-2 y-1]; Wang et al., 2007, GB101
measurement_temp <- 25.0                # temperature Vcmax/Jmax are measured at, typical 25.0 (celsius)
ncbnew <- 0.003                         # N alloc param: new branch N C at critical leaf N C
ncbnewz <- 0.003                        # N alloc param: new branch N C at zero leaf N C
nccnew <- 0.003                         # N alloc param: new coarse root N C at critical leaf N C
nccnewz <- 0.003                        # N alloc param: new coarse root N C at zero leaf N C
ncmaxfold <- 0.04                       # max N:C ratio of foliage in old stand, if the same as young=no effect
ncmaxfyoung <- 0.04                     # max N:C ratio of foliage in young stand, if the same as old=no effect
ncmaxr <- 0.03                          # max N:C ratio of roots
ncrfac <- 0.8                           # N:C of fine root prodn / N:C of leaf prodn
ncwimm <- 0.003                         # N alloc param: Immobile stem N C at critical leaf N C
ncwimmz <- 0.003                        # N alloc param: Immobile stem N C at zero leaf N C
ncwnew <- 0.003                         # N alloc param: New stem ring N:C at critical leaf N:C (mob)
ncwnewz <- 0.003                        # N alloc param: New stem ring N:C at zero leaf N:C (mobile)
nf_crit <- 0.015                        # leaf N:C below which N availability limits productivity
nf_min <- 0.005                         # leaf N:C minimum N concentration which allows productivity
nmax <- 0.24
nmin <- 0.95
nmin0 <- 0.0                            # (bewdy) minimum leaf n for +ve p/s (g/m2)
nmincrit <- 2.0                         # mineral N pool corresponding to Actnc0,etc (g/m2)
ntheta_root <- 5.5                      # Critical mineral N pool at max soil N:C (g/m2) (Parton et al 1993, McMurtrie et al 2001).
ntheta_topsoil <- 8.0                   # Fitted parameter based on Landsberg and Waring
nuptakez <- 0.0                         # Fitted parameter based on Landsberg and Waring
oi <- 210000.0                          # constant N uptake per year (1/yr)
p_atm_deposition <- 0.000086            # intercellular concentration of O2 [umol mol-1]
p_rate_par_weather <- 0.0001
passivesoilnz <- 1.0
passivesoilpz <- 1.0
passivesoilz <- 1.0
passncmax <- 0.142857                   # Passive pool (=1/7) N:C ratio of new SOM - maximum [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology.
passncmin <- 0.1                        # Passive pool (=1/10) N:C of new SOM - when Nmin=Nmin0 [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology.
passpcmax <- 0.05                       # Passive pool (=1/20) P:C ratio of new SOM - maximum [units: gP/gC]
passpcmin <- 0.005                      # Passive pool (=1/200) P:C of new SOM - when Pmin=Pmin0 [units: gP/gC]
pcbnew <- 0.0003                        # P alloc param: new branch P C at critical leaf P C
pcbnewz <- 0.0003                       # P alloc param: new branch P C at zero leaf P C
pccnew <- 0.0003                        # P alloc param: new coarse root P C at critical leaf P C
pccnewz <- 0.0003                       # P alloc param: new coarse root P C at zero leaf P C
pcmaxfold <- 0.002                      # max P:C ratio of foliage in old stand, if the same as young=no effect
pcmaxfyoung <- 0.002                    # max P:C ratio of foliage in young stand, if the same as old=no effect
pcmaxr <- 0.0006                        # max P:C ratio of roots
pcrfac <- 0.8                           # P:C of fine root prodp / P:C c of leaf prodp
pcwimm <- 0.00014                       # P alloc param: Immobile stem P C at critical leaf P C
pcwimmz <- 0.00014                      # P alloc param: Immobile stem P C at zero leaf P C
pcwnew <- 0.00014                       # P alloc param: New stem ring P:C at critical leaf P:C (mob)
pcwnewz <- 0.00014                      # P alloc param: New stem ring P:C at zero leaf P:C (mobile)
pf_crit <- 0.002                        # leaf P:C below which P availability limits productivity
pf_min <- 0.0002                        # leaf P:C minimum P concentration which allows productivity
phmax <- 7.6                            # max pH for determining effect on solubility of secondary P
phmin <- 5.0                            # min pH for determining effect on solubility of secondary P
phtextmin <- 0.000008                   # the solubility of secondary P corresponding to min pH (/yr)
phtextmax <- 0.00015                    # the solubility of secondary P corresponding to max pH (/yr)
phtextslope <- 0.00004                  # slope controlling effect of sand on secondary P flow to mineral P
p_lab_avail <- 0.0                      # Fraction of labile P available for plant uptake
pmax <- 0.002
pmin <- 0.01                            # (bewdy) minimum leaf p for +ve p/s (g/m2)
pmin0 <- 0.0                            # mineral P pool corresponding to Actpc0,etc (g/m2)
pmincrit <- 2.0                         # Critical mineral P pool at max soil P:C (g/m2)
prateloss <- 0.05                       # Rate of P loss from mineral P pool (/yr), Ref Wang et al., 2007, GB1018
prateuptake <- 3.6                      # Rate of P uptake from mineral P pool (/yr), guess value
prescribed_leaf_NC <- 0.03              # If the N-Cycle is switched off this needs to be set, e.g. 0.03
prescribed_leaf_PC <- 0.00249           # If the P-Cycle is switched off this needs to be set, e.g. 0.00249
previous_ncd <- 35.0                    # In the first year we don't have last years met_data, so I have precalculated the average of all the november-jan chilling values
psecmnp <- 0.000022                     # controls the flow from secondary to mineral P, used when text_effect_p set to 0
psi_sat_root <- -999.9                  # MPa
psi_sat_topsoil <- -999.9               # MPa
psie_topsoil <- NULL                    # Soil water potential at saturation (m)
psie_root <- NULL                       # Soil water potential at saturation (m)
puptakez <- 0.0255                      # constant P uptake per year (1/yr)
qs <- 1.0                               # exponent in water stress modifier, =1.0 JULES type representation, the smaller the values the more curved the depletion.
r0 <- 0.1325                            # root C at half-maximum N uptake (kg C/m3)
rate_ssorb_occ <- 0.000012              # Rate constant of the transfer of P from strongly sorbed pool to occluded pool, m-1 Yang et al. 2014, Biogeosciences
rate_sorb_ssorb <- 0.048                # Rate constant of the transfer of P from sorbed pool to strongly sorbed pool, m-1 Yang et al. 2014, Biogeosciences
rateloss <- 0.05                        # Rate of N loss from mineral N pool (/yr)
rateuptake <- 1.8                       # Rate of N uptake from mineral N pool (/yr) from here? http://face.ornl.gov/Finzi-PNAS.pdf Seems to correspond to very low NPP values
prdecay <- 0.6                          # root turnover rate (1/yr)
rdecaydry <- 0.6                        # root turnover rate - dry soil (1/yr)
resp_coeff <- 0.2                       # Respiration rate: from LPJ ENF, EBF, C3G = 1.2, Trop EBF, C4G = 0.2
retransmob <- 0.0                       # Fraction stem mobile N retranscd (/yr)
rfmult <- 1.0
root_exu_CUE <- -999.9
rooting_depth <- 2500.0                 # Rooting depth (mm)
rootsoil_type <- "sandy_clay_loam"
soil_order <- "andisol"                 # soil order
rretrans <- 0.0                         # root n retranslocation fraction
sapturnover <- 0.1                      # Sapwood turnover rate: conversion of sapwood to heartwood (1/yr)
psla <- 5.1                             # specific leaf area (m2 one-sided/kg DW)
slamax <- 5.1                           # (if equal slazero=no effect) specific leaf area new fol at max leaf N/C (m2 one-sided/kg DW)
slazero <- 5.1                          # (if equal slamax=no effect) specific leaf area new fol at zero leaf N/C (m2 one-sided/kg DW)
slowncmax <- 0.066666                   # Slow pool (=1/15) N:C ratio of new SOM - maximum [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology.
slowncmin <- 0.025                      # Slow pool (=1/40) N:C of new SOM - when Nmin=Nmin0" [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology.
slowpcmax <- 0.011111                   # Slow pool (=1/90) P:C ratio of new SOM - maximum [units: gP/gC].
slowpcmin <- 0.005                      # Slow pool (=1/200) P:C of new SOM - when Pmin=Pmin0" [units: gP/gC].
smax <- 700.0                           # Maximum amount of sorbed P in the soil [tt P ha-1]
store_transfer_len <- -999.9
structcn <- 150.0                       # C:N ratio of structural bit of litter input
structrat <- 0.0                        # structural input n:c as fraction of metab
structcp <- 5500.0                      # C:P ratio of structural bit of litter input, Ref Attiwill 1980, Aus. J. Bot. 28, 199-222 Table 9 sum of branch, stem, sap and heartwood;
structratp <- 0.0                       # structural input p:c as fraction of metab
soilph <- 4.5                           # soil pH value
sorpmx <- 5.0                           # maximum P sorption potential for a soil
sorpaf <- 1.0                           # slope term which controls the fraction of mineral P that is labile
targ_sens <- 0.5                        # sensitivity of allocation (leaf/branch) to track the target, higher values = less responsive.
theta <- 0.7                            # curvature of photosynthetic light response curve
theta_fc_root <- -999.9
theta_fc_topsoil <- -999.9
theta_sp_root <- -999.9
theta_sp_topsoil <- -999.9
theta_wp_root <- -999.9
theta_wp_topsoil <- -999.9
topsoil_depth <- 450.0                  # Topsoil depth (mm)
topsoil_type <- "loamy_sand"
pvcmax <- -999.9                        # maximum rate of carboxylation (umol m-2 s-1)
vcmaxna <- 27.707                       # slope of the reln btween vcmax and leaf N content, units = (umol [gN]-1 s-1) # And for Vcmax-N slopes (vcmaxna) see Table 8.2 in CLM4_tech_note, Oleson et al. 2010.
vcmaxnb <- 0.0                          # intercept of vcmax vs n, units = (umol [gN]-1 s-1) # And for Vcmax-N slopes (vcmaxna) see Table 8.2 in CLM4_tech_note, Oleson et al. 2010.
vcmaxpa <- 516.83
vcmaxpb <- 0.0
watdecaydry <- 0.0                      # water fractn for dry litterfall rates
watdecaywet <- 0.1                      # water fractn for wet litterfall rates
wcapac_root <- 300.0                    # Max plant avail soil water -root zone, i.e. total (mm) (smc_sat-smc_wilt) * root_depth (750mm) = [mm (water) / m (soil depth)]
wcapac_topsoil <- 67.5                  # Max plant avail soil water -top soil (mm)
wdecay <- 0.1                           # wood turnover rate (1/yr)
wetloss <- 0.5                          # Daily rainfall lost per lai (mm/day)
wretrans <- 0.7                         # mobile wood N retranslocation fraction
z0h_z0m <- 1.0                          # Assume z0m = z0h, probably a big assumption [as z0h often < z0m.], see comment in code!! But 0.1 might be a better assumption
growing_seas_len <- 0
lad <- 0.0                              # spherical leaf angle distribution
decayrate <- rep(0.0, 7)
leaf_abs <- 0.5                         # absorptance of solar radiation (0-1), typically 0.4-0.6
layer_thickness <- 0.1                  # Soil layer thickness (m)
n_layers <- 20                          # Number of soil layers
root_k <- 100.0                         # mass of roots for reaching 50% maximum depth (g m-2)
root_radius <- 0.0005                   # (m)
root_density <- 0.5e6                   # g biomass m-3
max_depth <- 2.0                        # (m)
root_resist <- 20                       # Evergreen value: fine root hydraulic resistivity (MPa s g mmol-1 H2O)
min_lwp <- -2.0                         # minimum leaf water potential (MPa)
potA <- NULL
potB <- NULL
cond1 <- NULL
cond2 <- NULL
cond3 <- NULL
porosity <- NULL
field_capacity <- NULL
wetting <- 10                           # number of wetting layers
burn_specific_yr <- -999.9
hurricane_doy <- -999.9
hurricane_yr <- -999.9
intercep_lai <- 3.0
longitude <- 35.9
prime_y <- 0.0
prime_z <- 0.0
theta_sat_topsoil <- -999.9
theta_sat_root <- -999.9
