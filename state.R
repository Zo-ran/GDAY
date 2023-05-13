activesoil <- 0.001                     # active C som pool (t/ha)
activesoiln <- 0.00004                  # active N som pool (t/ha)
activesoilp <- 0.000002                 # active P som pool (t/ha)
age <- 100.0                            # Current stand age (years)
avg_albranch <- 0.0                     # Average branch growing season allocation fractions
avg_alcroot <- 0.0                      # Average coarse root growing season allocation fractions
avg_alleaf <- 0.0                       # Average leaf growing season allocation fractions
avg_alroot <- 0.0                       # Average fine root growing season allocation fractions
avg_alstem <- 0.0                       # Average stem growing season allocation fractions
branch <- 0.001                         # branch c (t/ha)
branchn <- 0.00004                      # branch n (t/ha)
branchp <- 0.000002                     # branch p (t/ha)
canht <- 0.1                            # canopy height (m)
croot <- 0.0                            # coarse root c (t/ha)
crootn <- 0.0                           # coarse root n (t/ha)
crootp <- 0.0                           # coarse root p (t/ha)
cstore <- 0.0                           # C store for deciduous model (t/ha) annual ?
inorgn <- 0.0000                        # Inorganic soil N pool - dynamic (t/ha)
inorgp <- 0.0205                        # Inorganic soil P pool - dynamic (t/ha)
inorgavlp <- 0.096                      # Inorganic soil P pool - available mineral P = lab + sorbed (t/ha)
inorglabp <- 0.0000                     # Inorganic soil P pool - labile P (t/ha)
inorgsorbp <- 0.0                       # Inorganic soil P pool - sorbed P (t/ha)
inorgssorbp <- 0.0                      # Inorganic soil P pool - strongly sorbed P (t/ha)
inorgoccp <- 0.0                        # Inorganic soil P pool - occluded P (t/ha)
inorgparp <- 0.054                      # Inorganic soil P pool - parent P (t/ha)
lai <- NULL                             # leaf area index m2 (leaf) m-2 (ground)
remaining_days <- vector("numeric", 366)
leaf_out_days <- vector("numeric", 366)
growing_days <- vector("numeric", 366)
metabsoil <- 0.0                        # metabolic soil c (t/ha)
metabsoiln <- 0.0                       # metabolic soil n (t/ha)
metabsoilp <- 0.0                       # metabolic soil p (t/ha)
metabsurf <- 0.0                        # metabolic surface c (t/ha)
metabsurfn <- 0.0                       # metabolic surface n (t/ha)
metabsurfp <- 0.0                       # metabolic surface p (t/ha)
nstore <- 0.0                           # N store for deciduous model (t/ha)
pstore <- 0.0                           # P store for deciduous model (t/ha)
passivesoil <- 0.001                    # passive C som pool (t/ha)
passivesoiln <- 0.0004                  # passive N som pool (t/ha)
passivesoilp <- 0.000002                # passive P som pool (t/ha)
pawater_root <- 240.0                   # plant available water - root zone (mm)
pawater_topsoil <- 50.0                 # plant available water - top soil(mm)
prev_sma <- 1.0
root <- 0.001                           # root c (t/ha)
root_depth <- -9999.9                   # rooting depth, Dmax (m)
rootn <- 0.00004                        # root n (t/ha)
rootp <- 0.000002                       # root p (t/ha)
sapwood <- 0.001
shoot <- 0.001                          # shoot c (t/ha)
shootn <- 0.00004                       # shoot n (t/ha)
shootp <- 0.000002                      # shoot p (t/ha)
ssla <- 4.4                             # specific leaf area
slowsoil <- 0.001                       # slow C som pool (t/ha)
slowsoiln <- 0.00004                    # slow N som pool (t/ha)
slowsoilp <- 0.000002                   # slow P som pool (t/ha)
stem <- 0.001
stemn <- 0.00004                        # Stem N (t/ha) = stemnimm + stemnmob
stemnimm <- 0.00004
stemnmob <- 0.0
stemp <- 0.000002                       # Stem P (t/ha) = stempimm + stempmob
stempimm <- 0.000002
stempmob <- 0.0
structsoil <- 0.001                     # soil structural c (t/ha)
structsoiln <- 0.00004                  # soil structural n (t/ha)
structsoilp <- 0.000002                 # soil structural p (t/ha)
structsurf <- 0.001                     # surface structural c (t/ha)
structsurfn <- 0.00004                  # surface structural n (t/ha)
structsurfp <- 0.0000024                # surface structural p (t/ha)
shootnc <- NULL                         # shoot nc ratio
rootnc <- NULL                          # root pn ratio
shootpc <- NULL                         # shoot pc ratio
rootpc <- NULL                          # root pc ratio
anpp <- NULL                            # aboveground NPP
litterc <- NULL                         # litter carbon
littern <- NULL                         # litter nitrogen
litterp <- NULL                         # litter phosphorus
littercbg <- NULL                       # litter C belowground
littercag <- NULL                       # litter C aboveground
litternag <- NULL                       # litter N aboveground
litternbg <- NULL                       # litter N belowground
litterpag <- NULL                       # litter P aboveground
litterpbg <- NULL                       # litter P belowground
plantc <- 0.0                           # plant C
plantn <- 0.0                           # plant N
plantp <- 0.0                           # plant P
totalc <- 0.0                           # total C
totaln <- 0.0                           # total N
totalp <- 0.0                           # total P
soilc <- 0.0                            # Soil C
soiln <- 0.0                            # soil N
soilp <- 0.0                            # soil P
canopy_store <- 0.0
svcmax <- 0.0
twq <- 0.0                              # temperature of warmest quarter
wtfac_root <- 1.0
wtfac_topsoil <- 1.0
thickness <- NULL
root_mass <- NULL
root_length <- NULL
layer_depth <- NULL
water_frac <- NULL
initial_water <- 0.0
dry_thick <- 0.1                        # Thickness of dry soil layer above water table (m)
rooted_layers <- 0
predawn_swp <- 0.0                      # MPa
midday_lwp <- 0.0                       # MPa
max_lai <- -999.9
max_shoot <- -999.9
metabcnmax <- 5500.0
metabcnmin <- 5500.0
metabcpmax <- 5500.0
metabcpmin <- 5500.0
fipar <- 0
delta_sw_store <- 0
weighted_swp <- 0.0
wetting_bot <- NULL
wetting_top <- NULL

