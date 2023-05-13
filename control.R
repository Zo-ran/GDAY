git_hash <- "Err"
ofp_hdr <- NULL
cfg_fname <- "*NOTSET*"
met_fname <- "*NOTSET*"
out_fname <- "*NOTSET*"
out_subdaily_fname <- "*NOTSET*"
out_fname_hdr <- "*NOTSET*"
out_param_fname <- "*NOTSET*"
ifp <- NULL
ofp <- NULL
ofp_sd <- NULL
ofp_hdr <- NULL
day_idx <- 1
alloc_model <- 1                        # C allocation scheme: FIXED, GRASSES, ALLOMETRIC
assim_model <- MATE                     # Photosynthesis model: BEWDY (not coded :p) or MATE
calc_sw_params <- TRUE                  # false=user supplies field capacity and wilting point, true=calculate them based on cosby et al.
deciduous_model <- FALSE                # evergreen_model=False, deciduous_model=True
fixed_stem_nc <- TRUE                   # False=vary stem N:C with foliage, True=fixed stem N:C
fixed_stem_pc <- TRUE                   # False=vary stem P:C with foliage, True=fixed stem P:C
fixed_lai <- FALSE                      # Fix LAI
fixleafnc <- FALSE                      # fixed leaf N C ?
fixleafpc <- FALSE                      # fixed leaf P C ?
grazing <- FALSE                        # Is foliage grazed? 0=No, 1=daily, 2=annual and then set disturbance_doy=doy
gs_model <- MEDLYN                      # Stomatal conductance model, currently only this one is implemented
model_optroot <- FALSE                  # Ross's optimal root model...not sure if this works yet...0=off, 1=on
modeljm <- 1                            # modeljm=0, Jmax and Vcmax parameters are read in, modeljm=1, parameters are calculated from leaf N & P content, modeljm=2, Vcmax is calculated from leaf N & P content but Jmax is related to Vcmax
ncycle <- TRUE                          # Nitrogen cycle on or off?
pcycle <- FALSE                         # Phosphorus cycle on or off?
triose_p <- TRUE                        # Triose phosphates limitation on photosynthesis on or off?
tpu_removed <- TRUE                     # deducting TPU limitation effect
nuptake_model <- 1                      # 0=constant uptake, 1=func of N inorgn, 2=depends on rate of soil N availability
puptake_model <- 1                      # 0=constant uptake, 1=func of P inorgp, 2=depends on rate of soil P availability
output_ascii <- TRUE                    # If this is false you get a binary file as an output.
passiveconst <- FALSE                   # hold passive pool at passivesoil
print_options <- END                    # DAILY=every timestep, END=end of run
ps_pathway <- C3                        # Photosynthetic pathway, c3/c4
respiration_model <- FIXED              # Plant respiration ... Fixed, TEMPERATURE or BIOMASS
strfloat <- 0                           # Structural pool input N:C varies=1, fixed=0
strpfloat <- 0                          # Structural pool input P:C varies=1, fixed=0
sw_stress_model <- 1                    # JULES type linear stress func, or Landsberg and Waring non-linear func
text_effect_p <- 1                      # soil texture effect on strongly sorbed P flow to mineral P = 1 use texture effect = 0 use pre-defined constant
use_eff_nc <- 0                         # use constant leaf n:c for  metfrac s
water_stress <- TRUE                    # water stress modifier turned on=TRUE (default)...ability to turn off to test things without drought stress = FALSE
water_balance <- 0                      # Water calculations: 0=simple 2 layered bucket 1=SPA-style hydraulics
aci_relationship <- WALKER              # Controlling for the relationship between jmax and leaf N/P, and vcmax and leaf N/P, based on either Walker 2014 global synthesis or Ellsworth 2015 EucFACE met_data
spin_up <- TRUE                         # Spin up to a steady state? If False it just runs the model
num_years <- 0                          # Total number of years simulated
num_days <- 0                           # Number of days in a year: 365/366
total_num_days <- 0                     # Total number of days
PRINT_GIT <- FALSE                      # print the git hash to the cmd line and exit? Called from cmd line parsar
sub_daily <- FALSE                      # Run at daily or 30 minute timestep
num_hlf_hrs <- 48
pdebug <- FALSE
adjust_rtslow <- FALSE
disturbance <- FALSE
exudation <- FALSE
frost <- FALSE
hurricane <- FALSE
sma_obj <- NULL
