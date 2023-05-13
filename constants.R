EPSILON <- 1E-08
STRING_LENGTH <- 2000

# Stomatal conductanct models
MEDLYN <- 0

# Photosynthesis models
MATE <- 0
BEWDY <- 1

# A Ci relationship
WALKER <- 0
ELLSWORTH <- 1
BASELINE <- 2

# Respiration models
FIXED <- 0
VARY <- 1

# Allocation models
FIXED <- 0
GRASSES <- 1
ALLOMETRIC <- 2

# Ps pathway
C3 <- 0
C4 <- 1

# output time step, where end = the final state
SUBDAILY <- 0
DAILY <- 1
END <- 2

# Texture identifiers
SILT <- 1
SAND <- 2
CLAY <- 3

# water balance identifiers
BUCKET <- 0
HYDRAULICS <- 1

MOL_2_MMOL <- 1000.0
MMOL_2_MOL <- 1E-03
NDAYS_IN_YR <- 365.25
DEG_TO_KELVIN <- 273.15
RGAS <- 8.314
SECS_IN_HOUR <- 3600.0
UMOL_TO_MOL <- 1E-6
MOL_C_TO_GRAMS_C <- 12.0
G_AS_TONNES <- 1E-6
M2_AS_HA <- 1E-4
KG_AS_G <- 1E+3
WATT_HR_TO_MJ <- 0.0036
MJ_TO_WATT_HR <- 1.0 / 0.0036
MM_TO_M <- 0.001
M_TO_MM <- 1000.0
GRAMS_C_TO_MOL_C <- 1.0 / 12.0
UMOL_TO_MOL <- 1E-6
MOL_TO_UMOL <- 1E6
GRAM_C_2_TONNES_HA <- G_AS_TONNES / M2_AS_HA
UMOL_2_GRAMS_C <- UMOL_TO_MOL * MOL_C_TO_GRAMS_C
TONNES_HA_2_KG_M2 <- 0.1
TONNES_HA_2_G_M2 <- 100.0
YRS_IN_DAYS <- 1.0 / 365.25
DAYS_IN_YRS <- 365.25
G_M2_2_TONNES_HA <- 0.01
KG_M2_2_TONNES_HA <- 10.0
KPA_2_MPA <- 0.001
KG_AS_TONNES <- 1E-3
G_TO_KG <- 0.001
TONNES_AS_KG <- 1.0 / KG_AS_TONNES
CP <- 1010.0              # specific heat of dry air (j kg-1 k-1)
MASS_AIR <- 29.0E-3       # mol mass air (kg mol-1)
H2OLV0 <- 2.501E6         # latent heat H2O (J kg-1)
H2OMW <- 18E-3            # mol mass H20 (kg mol-1)
DHEAT <- 21.5E-6          # molecular diffusivity for heat
GBVGBH <- 1.075           # Ratio of Gbw:Gbh
GSVGSC <- 1.57            # Ratio of Gsw:Gsc
GBHGBC <- 1.32            # Ratio of Gbh:Gbc
SIGMA <- 5.67e-8          # Steffan Boltzman constant (W/m2/K4)
LEAF_EMISSIVITY <- 0.95   # Emissivity of thermal radiation by leaf
KPA_2_PA <- 1000
KPA_2_MPA <- 0.001
METER_OF_HEAD_TO_MPA <- 9.81 * KPA_2_MPA # Height (m) x gravity (m/s2) = pressure (kPa)
PA_2_KPA <- 0.001
CM_2_M <- 0.01

# Solar radiaiton 1 W m-2 ~ 2.3 umol m-2 s-1 PAR Landsberg and Sands, Cp2, pg 20. (1.0 / 2.3)
SW_2_PAR <- 2.3
PAR_2_SW <- 1.0 / SW_2_PAR
J_TO_MJ <- 1.0E-6
MJ_TO_J <- 1.0 / J_TO_MJ
J_2_UMOL <- 4.57               # Conversion from J to umol quanta
UMOL_2_JOL <- 1.0 / J_2_UMOL   # Conversion from umol quanta to J
SEC_2_HLFHR <- 1800
MOLE_WATER_2_G_WATER <- 18.02 # oxygen = 16g/mol, hydrogren = 1.01 g/mol
SUNLIT <- 1
SHADED <- 2
NUM_LEAVES <- 2

MAXSTP <- 10000
TINY <- 1.0e-30
SAFETY <- 0.9
PGROW <- -0.2
PSHRNK <- -0.25
ERRCON <- 1.89e-4
