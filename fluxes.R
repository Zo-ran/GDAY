# C fluxes
gpp_gCm2 <- 0.0
gpp_am <- 0.0
gpp_pm <- 0.0
npp_gCm2 <- 0.0
gpp <- 0.0
npp <- 0.0
nep <- 0.0
auto_resp <- 0.0
hetero_resp <- 0.0
retrans <- 0.0
retransp <- 0.0
apar <- 0.0

# N fluxes
nuptake <- 0.0
nloss <- 0.0
npassive <- 0.0              # n passive -> active
ngross <- 0.0                # N gross mineralisation
nimmob <- 0.0                # N immobilisation in SOM
nlittrelease <- 0.0          # N rel litter <- struct + metab
activelossf <- 0.0           # frac of active C -> CO2
nmineralisation <- 0.0

# P fluxes
puptake <- 0.0
ploss <- 0.0
ppassive <- 0.0              # p passive -> active
pgross <- 0.0                # P gross mineralisation
pimmob <- 0.0                # P immobilisation in SOM
plittrelease <- 0.0          # P rel litter <- struct + metab
pmineralisation <- 0.0

# water fluxes
wue <- 0.0
et <- 0.0
soil_evap <- 0.0
transpiration <- 0.0
interception <- 0.0
throughfall <- 0.0
canopy_evap <- 0.0
runoff <- 0.0
gs_mol_m2_sec <- 0.0
ga_mol_m2_sec <- 0.0
omega <- 0.0
day_ppt <- 0.0
day_wbal <- 0.0
total_soil_resist <- 0.0

# daily C production
cpleaf <- 0.0
cproot <- 0.0
cpcroot <- 0.0
cpbranch <- 0.0
cpstem <- 0.0

# daily N production
npleaf <- 0.0
nproot <- 0.0
npcroot <- 0.0
npbranch <- 0.0
npstemimm <- 0.0
npstemmob <- 0.0

# daily P production
ppleaf <- 0.0
pproot <- 0.0
ppcroot <- 0.0
ppbranch <- 0.0
ppstemimm <- 0.0
ppstemmob <- 0.0

# dying stuff
deadleaves <- 0.0   # Leaf litter C production (t/ha/yr)
deadroots <- 0.0    # Root litter C production (t/ha/yr)
deadcroots <- 0.0   # Coarse root litter C production (t/ha/yr)
deadbranch <- 0.0   # Branch litter C production (t/ha/yr)
deadstems <- 0.0    # Stem litter C production (t/ha/yr)
deadleafn <- 0.0    # Leaf litter N production (t/ha/yr)
deadrootn <- 0.0    # Root litter N production (t/ha/yr)
deadcrootn <- 0.0   # Coarse root litter N production (t/ha/yr)
deadbranchn <- 0.0  # Branch litter N production (t/ha/yr)
deadstemn <- 0.0    # Stem litter N production (t/ha/yr)
deadleafp <- 0.0    # Leaf litter P production (t/ha/yr)
deadrootp <- 0.0    # Root litter P production (t/ha/yr)
deadcrootp <- 0.0   # Coarse root litter P production (t/ha/yr)
deadbranchp <- 0.0  # Branch litter P production (t/ha/yr)
deadstemp <- 0.0    # Stem litter P production (t/ha/yr)
deadsapwood <- 0.0

# grazing stuff
ceaten <- 0.0       # C consumed by grazers (t C/ha/y)
neaten <- 0.0       # N consumed by grazers (t N/ha/y)
peaten <- 0.0       # P consumed by grazers (t P/ha/y)
faecesc <- 0.0      # Flux determined by faeces C:N
nurine <- 0.0       # Rate of N input to soil in urine (t/ha/y)
purine <- 0.0       # Rate of P input to soil in urine (t/ha/y)

leafretransn <- 0.0
leafretransp <- 0.0

# C N & P Surface litter
surf_struct_litter <- 0.0
surf_metab_litter <- 0.0
n_surf_struct_litter <- 0.0
n_surf_metab_litter <- 0.0
p_surf_struct_litter <- 0.0
p_surf_metab_litter <- 0.0

# C N & P Root Litter
soil_struct_litter <- 0.0
soil_metab_litter <- 0.0
n_soil_struct_litter <- 0.0
n_soil_metab_litter <- 0.0
p_soil_struct_litter <- 0.0
p_soil_metab_litter <- 0.0

# C N & P litter fluxes to slow pool
surf_struct_to_slow <- 0.0
soil_struct_to_slow <- 0.0
n_surf_struct_to_slow <- 0.0
n_soil_struct_to_slow <- 0.0
p_surf_struct_to_slow <- 0.0
p_soil_struct_to_slow <- 0.0

# C N & P litter fluxes to active pool
surf_struct_to_active <- 0.0
soil_struct_to_active <- 0.0
n_surf_struct_to_active <- 0.0
n_soil_struct_to_active <- 0.0
p_surf_struct_to_active <- 0.0
p_soil_struct_to_active <- 0.0

# Metabolic fluxes to active pool
surf_metab_to_active <- 0.0
soil_metab_to_active <- 0.0
n_surf_metab_to_active <- 0.0
n_soil_metab_to_active <- 0.0
p_surf_metab_to_active <- 0.0
p_soil_metab_to_active <- 0.0

# fluxes out of active pool
active_to_slow <- 0.0
active_to_passive <- 0.0
n_active_to_slow <- 0.0
n_active_to_passive <- 0.0
p_active_to_slow <- 0.0
p_active_to_passive <- 0.0

# fluxes out of slow pool
slow_to_active <- 0.0
slow_to_passive <- 0.0
n_slow_to_active <- 0.0
n_slow_to_passive <- 0.0
p_slow_to_active <- 0.0
p_slow_to_passive <- 0.0
p_slow_biochemical <- 0.0

# C N & P fluxes from passive to active pool
passive_to_active <- 0.0
n_passive_to_active <- 0.0
p_passive_to_active <- 0.0

# C source fluxes from the active, slow and passive pools
c_into_active <- 0.0
c_into_slow <- 0.0
c_into_passive <- 0.0

# inorganic P flux exchanges
p_lab_in <- 0.0
p_lab_out <- 0.0
p_sorb_in <- 0.0
p_sorb_out <- 0.0
p_min_to_ssorb <- 0.0
p_ssorb_to_min <- 0.0
p_ssorb_to_occ <- 0.0
p_par_to_min <- 0.0

# CO2 flows to the air
# C flows to the air
co2_to_air <- rep(0.0, 7)

# C allocated fracs
alleaf <- 0.0
alroot <- 0.0
alcroot <- 0.0
albranch <- 0.0
alstem <- 0.0

# Misc stuff
cica_avg <- 0.0 # used in water balance, only when running mate model

rabove <- 0.0
tfac_soil_decomp <- 0.0
co2_rel_from_surf_struct_litter <- 0.0
co2_rel_from_soil_struct_litter <- 0.0
co2_rel_from_surf_metab_litter <- 0.0
co2_rel_from_soil_metab_litter <- 0.0
co2_rel_from_active_pool <- 0.0
co2_rel_from_slow_pool <- 0.0
co2_rel_from_passive_pool <- 0.0

# Hydraulics stuff
soil_conduct <- NULL
swp <- NULL
soilR <- NULL
fraction_uptake <- NULL
ppt_gain <- NULL
water_loss <- NULL
water_gain <- NULL
est_evap <- NULL

# met data
# sub-daily
rain <- 0
wind <- 0
press <- 0
vpd <- 0
tair <- 0
sw_rad <- 0
par <- 0
Ca <- 0
ndep <- 0
nfix <- 0       # N inputs from biological fixation (t/ha/timestep (d/30min))
pdep <- 0
tsoil <- 0

# dailys
tair_am <- 0
tair_pm <- 0
sw_rad_am <- 0
sw_rad_pm <- 0
vpd_am <- 0
vpd_pm <- 0
wind_am <- 0
wind_pm <- 0
Tk_am <- 0
Tk_pm <- 0

root_exc <- 0
root_exn <- 0


#canopy_wk
cos_zenith <- 0
elevation <- 0
diffuse_frac <- 0
kb <- 0
cwdirect_frac <- 0
N0 <- 0
ileaf <- 0
Cs <- 0
dleaf <- 0
an_leaf <- vector("numeric", 2)      # leaf net photosynthesis (umol m-2 s-1)
rd_leaf <- vector("numeric", 2)      # leaf respiration in the light (umol m-2 s-1)
gsc_leaf <- vector("numeric", 2)      # leaf stomatal conductance to CO2 (mol m-2 s-1)
apar_leaf <- vector("numeric", 2)    # leaf abs photosyn. active rad. (umol m-2 s-1)
trans_leaf <- vector("numeric", 2)   # leaf transpiration (mol m-2 s-1)
rnet_leaf <- vector("numeric", 2)    # leaf net radiation (W m-2)
lai_leaf <- vector("numeric", 2)     # sunlit and shaded leaf area (m2 m-2)
omega_leaf <- vector("numeric", 2)   # leaf decoupling coefficient (-)
tleaf <- vector("numeric", 2)        # leaf temperature (deg C)
lwp_leaf <- vector("numeric", 2)     # leaf water potential (MPa)
fwsoil_leaf <- vector("numeric", 2)  # Effective beta
cscalar <- vector("numeric", 2)      # scale from single leaf to canopy
ts_Cs <- 0           # Temporary variable to store Cs
ts_vcmax <- 0        # Temporary variable to store vcmax
ts_km <- 0           # Temporary variable to store km
ts_gamma_star <- 0   # Temporary variable to store gamma_star
ts_rd <- 0           # Temporary variable to store rd
ts_Vj <- 0           # Temporary variable to store Vj
tleaf_new <- 0       # new leaf temperature (deg C)
an_canopy <- 0       # canopy net photosynthesis (umol m-2 s-1)
rd_canopy <- 0       # canopy respiration in the light (umol m-2 s-1)
gsc_canopy <- 0      # canopy stomatal conductance to CO2 (mol m-2 s-1)
apar_canopy <- 0     # canopy abs photosyn. active rad. (umol m-2 s-1)
omega_canopy <- 0    # canopy decoupling coefficient (-)
trans_canopy <- 0    # canopy transpiration (mm 30min-1)
rnet_canopy <- 0     # canopy net radiation (W m-2)
lwp_canopy <- 0      # Leaf water potential for the canopy(MPa)

rexc_cue <- 0
co2_released_exud <- 0
ninflow <- 0
factive <- 0
rtslow <- 0
p_atm_dep <- 0