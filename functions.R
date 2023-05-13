allocate_numerical_libs_stuff <- function () {
    xp <<- vector("numeric", kmax)
    yp <<- matrix(nrow = N, ncol = kmax)
    yscal <<- vector("numeric", N)
    y <<- vector("numeric", N)
    dydx <<- vector("numeric", N)
    ystart <<- vector("numeric", N)
    ak2 <<- vector("numeric", N)
    ak3 <<- vector("numeric", N)
    ak4 <<- vector("numeric", N)
    ak5 <<- vector("numeric", N)
    ak6 <<- vector("numeric", N)
    ytemp <<- vector("numeric", N)
    yerr <<- vector("numeric", N)
}

initialise_roots <- function () {
    # Set up all the rooting arrays for use with the hydraulics assumptions

    thickness <<- vector("numeric", n_layers)
    # root mass is g biomass, i.e. ~twice the C content
    root_mass <<- vector("numeric", n_layers)
    root_length <<- vector("numeric", n_layers)
    layer_depth <<- vector("numeric", n_layers)
    thick <- 0.1
    layer_depth[1] <<- thick
    thickness[1] <<- thick
    for (i in 2:n_layers) {
        thick <- thick + layer_thickness
        layer_depth[i] <<- thick
        thickness[i] <<- layer_thickness
        # made up initalisation, following SPA, get replaced second timestep
        root_mass[i] <<- 0.1
        root_length[i] <<- 0.1
        rooted_layers <<- n_layers
    }
}

write_output_subdaily_header <- function() {
    # Write 30 min fluxes headers to an output CSV file. This is very basic
    # for now...
    #create data frame with 0 rows and 3 columns
    df <- data.frame(matrix(ncol = 9, nrow = 0))
    #provide column names
    colnames(df) <- c('year', 'doy', 'hod', 'an_canopy', 'rd_canopy', 'gsc_canopy', 'apar_canopy', 'trans_canopy', 'tleaf')

    return (df)
}

write_output_header <- function() {
    df <- data.frame(matrix(ncol = 107, nrow = 0))
    #provide column names
    colnames(df) <- c(
        'year', 'doy', 'wtfac_root', 'wtfac_topsoil', 'pawater_root', 'pawater_topsoil',
        'shoot', 'lai', 'branch', 'stem', 'root', 'croot',
        'shootn', 'branchn', 'stemn', 'rootn', 'crootn',
        'shootp', 'branchp', 'stemp', 'rootp', 'crootp',
        'cstore', 'nstore', 'pstore',
        'soilc', 'soiln', 'soilp', 'inorgn',
        'inorgp', 'inorgavlp', 'inorglabp', 'inorgsorbp', 'inorgssorbp', 'inorgoccp', 'inorgparp',
        'litterc', 'littercag', 'littercbg', 'litternag', 'litternbg',
        'litterpag', 'litterpbg',
        'activesoil', 'slowsoil', 'passivesoil',
        'activesoiln', 'slowsoiln', 'passivesoiln', 'activesoilp', 'slowsoilp', 'passivesoilp',
        'et', 'transpiration', 'soil_evap', 'canopy_evap', 'runoff',
        'gs_mol_m2_sec', 'ga_mol_m2_sec',
        'deadleaves', 'deadbranch', 'deadstems', 'deadroots', 'deadcroots',
        'deadleafn', 'deadbranchn', 'deadstemn', 'deadrootn', 'deadcrootn',
        'deadleafp', 'deadbranchp', 'deadstemp', 'deadrootp', 'deadcrootp',
        'nep', 'gpp', 'npp', 'hetero_resp', 'auto_resp', 'apar',
        'cpleaf', 'cpbranch', 'cpstem', 'cproot', 'cpcroot',
        'npleaf', 'npbranch', 'npstemimm', 'npstemmob', 'nproot', 'npcroot',
        'ppleaf', 'ppbranch', 'ppstemimm', 'ppstemmob', 'pproot', 'ppcroot',
        'nuptake', 'ngross', 'nmineralisation', 'nloss',
        'puptake', 'pgross', 'pmineralisation', 'ploss',
        'leafretransn', 'leafretransp'
    )

    return (df)
}

setup_hydraulics_arrays <- function () {
    potA <<- vector("numeric", n_layers)
    potB <<- vector("numeric", n_layers)
    cond1 <<- vector("numeric", n_layers)
    cond2 <<- vector("numeric", n_layers)
    cond3 <<- vector("numeric", n_layers)
    porosity <<- vector("numeric", n_layers)
    field_capacity <<- vector("numeric", n_layers)
    soil_conduct <<- vector("numeric", n_layers)
    swp <<- vector("numeric", n_layers)
    soilR <<- vector("numeric", n_layers)
    fraction_uptake <<- vector("numeric", n_layers)
    ppt_gain <<- vector("numeric", n_layers)
    water_loss <<- vector("numeric", n_layers)
    water_gain <<- vector("numeric", n_layers)
    water_frac <<- vector("numeric", n_layers)
    wetting_bot <<- vector("numeric", wetting)
    wetting_top <<- vector("numeric", wetting)
    est_evap <<- vector("numeric", n_layers)
}

fill_up_solar_arrays <- function() {
    hour_idx <<- 1
    ntimesteps <- total_num_days * 48
    cz_store <<- vector("numeric", ntimesteps)
    ele_store <<- vector("numeric", ntimesteps)
    df_store <<- vector("numeric", ntimesteps)
    for (nyr in 1:num_years) {
        year <- ma$year[hour_idx]
        if (is_leap_year(year)) {
            num_days <- 366
        } else {
            num_days <- 365
        }
        for (doy in 1:num_days) {
            for (hod in 1:num_hlf_hrs) {
                # calculate solar geometry
                hod <- hod / 2
                gamma <- day_angle(doy)
                rdec <- calculate_solar_declination(doy, gamma)
                et <- calculate_eqn_of_time(gamma)
                t0 <- calculate_solar_noon(et, longitude)
                h <- calculate_hour_angle(hod, t0)
                rlat <- DEG2RAD(latitude)

                # A13 - De Pury & Farquhar
                sin_beta <- sin(rlat) * sin(rdec) + cos(rlat) * cos(rdec) * cos(h)

                cos_zenith <<- sin_beta # The same thing, going to use throughout
                if (cos_zenith > 1.0)
                    cos_zenith <<- 1.0
                else if (cos_zenith < 0.0)
                    cos_zenith <<- 0.0

                zenith_angle <- RAD2DEG(acos(cos_zenith))
                elevation <<- 90.0 - zenith_angle


                sw_rad <- ma$par[hour_idx] * PAR_2_SW


                #get diffuse frac
                # sine of the elev of the sun above the horizon is the same as cos_zen
                So <- calc_extra_terrestrial_rad(doy, cos_zenith)

                # atmospheric transmisivity
                tau <- estimate_clearness(sw_rad, So)

                cos_zen_sq <- cos_zenith * cos_zenith

                # For zenith angles > 80 degrees, diffuse_frac = 1.0
                if (cos_zenith > 0.17) {

                    # Spitters formula
                    R <- 0.847 - 1.61 * cos_zenith + 1.04 * cos_zen_sq
                    K <- (1.47 - R) / 1.66
                    if (tau <= 0.22) {
                        diffuse_frac <<- 1.0
                    } else if (tau > 0.22 && tau <= 0.35) {
                        diffuse_frac <<- 1.0 - 6.4 * (tau - 0.22) * (tau - 0.22)
                    } else if (tau > 0.35 && tau <= K) {
                        diffuse_frac <<- 1.47 - 1.66 * tau
                    } else {
                        diffuse_frac <<- R
                    }

                } else {
                    diffuse_frac <<- 1.0
                }

                # doubt we need this, should check
                if (diffuse_frac <= 0.0) {
                    diffuse_frac <<- 0.0
                } else if (diffuse_frac >= 1.0) {
                    diffuse_frac <<- 1.0
                }

                direct_frac <<- 1.0 - diffuse_frac


                cz_store[hour_idx] <<- cos_zenith
                ele_store[hour_idx] <<- elevation
                df_store[hour_idx] <<- diffuse_frac
                hour_idx <<- hour_idx + 1
            }
        }
    }
}

calc_carbon_allocation_fracs <- function (npitfac) {
    # Carbon allocation fractions to move photosynthate through the plant.
    #
    # Parameters:
    # -----------
    # npitfac : float
    #     the smallest value of leaf N:C as a fraction of 'Ncmaxfyoung' (max 1.0) &
    #     leaf P:C as a fraction of "Pcmaxfyoung" (max 1.0)
    #
    # Returns:
    # --------
    # alleaf : float
    #     allocation fraction for shoot
    # alroot : float
    #     allocation fraction for fine roots
    # albranch : float
    #     allocation fraction for branches
    # alstem : float
    #     allocation fraction for stem
    #
    # References:
    # -----------
    # Corbeels, M. et al (2005) Ecological Modelling, 187, 449-474.
    # McMurtrie, R. E. et al (2000) Plant and Soil, 224, 135-152.

    # this is obviously arbitary
    min_stem_alloc <- 0.01

    if (alloc_model == FIXED){
        alleaf <<- (c_alloc_fmax + npitfac *
                     (c_alloc_fmax - c_alloc_fmin))

        alroot <<- (c_alloc_rmax + npitfac *
                     (c_alloc_rmax - c_alloc_rmin))

        albranch <<- (c_alloc_bmax + npitfac *
                       (c_alloc_bmax - c_alloc_bmin))

        # allocate remainder to stem
        alstem <<- 1.0 - alleaf - alroot - albranch

        alcroot <<- c_alloc_cmax * alstem
        alstem <<- alstem - alcroot

    } else if (alloc_model == GRASSES) {

        # First figure out root allocation given available water & nutrients hyperbola shape to allocation
        alroot <<- (c_alloc_rmax * c_alloc_rmin /
                     (c_alloc_rmin + (c_alloc_rmax - c_alloc_rmin) * prev_sma))
        alleaf <<- 1.0 - alroot

        # Now adjust root & leaf allocation to maintain balance, accounting
        #    for stress e.g. Sitch et al. 2003, GCB.
        #
        #  leaf-to-root ratio under non-stressed conditons
        # lr_max = 0.8
        #
        #  Calculate adjustment on lr_max, based on current "stress"
        #    calculated from running mean of N and water stress
        # stress = lr_max * prev_sma
        #
        # calculate new allocation fractions based on imbalance in *biomass*
        # mis_match = shoot / (root * stress)
        #
        #
        # if (mis_match > 1.0) {
        #     reduce leaf allocation fraction
        #     adj = alleaf / mis_match
        #     alleaf = MAX(c_alloc_fmin, MIN(c_alloc_fmax, adj))
        #     alroot = 1.0 - alleaf
        # } else {
        #      reduce root allocation
        #     adj = alroot * mis_match
        #     alroot = MAX(c_alloc_rmin, MIN(c_alloc_rmax, adj))
        #     alleaf = 1.0 - alroot
        # }
        alstem <<- 0.0
        albranch <<- 0.0
        alcroot <<- 0.0

    } else if (alloc_model == ALLOMETRIC) {

         # Calculate tree height: allometric reln using the power function
         #   (Causton, 1985)
        canht <<- heighto * (stem ^ htpower)

        # LAI to stem sapwood cross-sectional area (As m-2 m-2)
        #    (dimensionless)
        #    Assume it varies between LS0 and LS1 as a linear function of tree
        #    height (m)
        arg1 <- sapwood * TONNES_AS_KG * M2_AS_HA
        arg2 <- canht * density * cfracts
        sap_cross_sec_area <- arg1 / arg2
        leaf2sap <- lai / sap_cross_sec_area

        # Allocation to leaves dependant on height. Modification of pipe
        #    theory, leaf-to-sapwood ratio is not constant above a certain
        #    height, due to hydraulic constraints (Magnani et al 2000 Deckmyn
        #    et al. 2006).

        if (canht < height0) {
            leaf2sa_target <- leafsap0
        } else if (float_eq(canht, height1)) {
            leaf2sa_target <- leafsap1
        } else if (canht > height1) {
            leaf2sa_target <- leafsap1
        } else {
            arg1 <- leafsap0
            arg2 <- leafsap1 - leafsap0
            arg3 <- canht - height0
            arg4 <- height1 - height0
            leaf2sa_target <- arg1 + (arg2 * arg3 / arg4)
        }
        alleaf <<- alloc_goal_seek(leaf2sap, leaf2sa_target, c_alloc_fmax, targ_sens)

        # Allocation to branch dependent on relationship between the stem
        #    and branch
        target_branch <- branch0 * (stem ^ branch1)
        albranch <<- alloc_goal_seek(branch, target_branch, c_alloc_bmax, targ_sens)

        coarse_root_target <- croot0 * (stem ^ croot1)
        alcroot <<- alloc_goal_seek(croot, coarse_root_target, c_alloc_cmax, targ_sens)

        # figure out root allocation given available water & nutrients
        #    hyperbola shape to allocation, this is adjusted below as we aim
        #    to maintain a functional balance

        alroot <<- (c_alloc_rmax * c_alloc_rmin /
                     (c_alloc_rmin + (c_alloc_rmax - c_alloc_rmin) * prev_sma))

        alstem <<- 1.0 - alroot - albranch - alleaf - alcroot

        # minimum allocation to leaves - without it tree would die, as this
        #    is done annually.
        if (deciduous_model) {
            if (alleaf < 0.05) {
                min_leaf_alloc <- 0.05
                if (alstem > min_leaf_alloc)
                    alstem <<- alstem - min_leaf_alloc
                else
                    alroot <<- alroot - min_leaf_alloc
                alleaf <<- min_leaf_alloc
            }
        }
    } else {
        stop(paste("Unknown C allocation model:", alloc_model))
    }

    # fprintf(stderr, "alleaf %f, alstem %f, alroot %f,  canht %f\n",
    #        alleaf, albranch + alstem, alroot, canht)
    #
    # Total allocation should be one, if not print warning
    total_alloc <- alroot + alleaf + albranch + alstem + alcroot
    if (total_alloc > 1.0 + EPSILON) {
        stop(paste("Allocation fracs > 1:", total_alloc))
    }
}

allocate_stored_cnp <- function() {
    # Allocate stored C,N and P. This is either down as the model is initialised
    # for the first time or at the end of each year.

    # ========================
    # Carbon - fixed fractions
    # ========================
    c_to_alloc_shoot <<- alleaf * cstore
    c_to_alloc_root <<- alroot * cstore
    c_to_alloc_croot <<- alcroot * cstore
    c_to_alloc_branch <<- albranch * cstore
    c_to_alloc_stem <<- alstem * cstore

    # =========================================================
    #  Nitrogen - Fixed ratios N allocation to woody components.
    # =========================================================

    # N flux into new ring (immobile component structrual components)
    n_to_alloc_stemimm <<- cstore * alstem * ncwimm

    # N flux into new ring (mobile component can be retrans for new woody tissue)
    n_to_alloc_stemmob <<- cstore * alstem * (ncwnew - ncwimm)
    n_to_alloc_branch <<- cstore * albranch * ncbnew
    n_to_alloc_croot <<- cstore * alcroot * nccnew

    # Calculate remaining N left to allocate to leaves and roots
    ntot <- max(0.0, (nstore - n_to_alloc_stemimm - n_to_alloc_stemmob - n_to_alloc_branch))

    # allocate remaining N to flexible-ratio pools
    n_to_alloc_shoot <<- (ntot * alleaf / (alleaf + alroot * ncrfac))
    n_to_alloc_root <<- ntot - n_to_alloc_shoot

    # leaf_NC = n_to_alloc_shoot / c_to_alloc_shoot
    # if leaf_NC > 0.04:
    #     n_to_alloc_shoot = c_to_alloc_shoot * 0.14
    #
    # n_to_alloc_root = ntot - n_to_alloc_shoot
    #
    #
    # root_NC = n_to_alloc_root / c_to_alloc_root
    # ncmaxr = 0.04 * ncrfac
    # if root_NC > ncmaxr:
    #     extrar = (n_to_alloc_root -
    #               (c_to_alloc_root * ncmaxr))
    #
    #     inorgn += extrar
    #     n_to_alloc_root -= extrar
    #
    #
    # =========================================================
    # Phosphorus - Fixed ratios P allocation to woody components.
    # =========================================================
    #
    # P flux into new ring (immobile component structrual components)
    p_to_alloc_stemimm <<- cstore * alstem * pcwimm

    # P flux into new ring (mobile component can be retrans for new woody tissue
    p_to_alloc_stemmob <<- cstore * alstem * (pcwnew - pcwimm)
    p_to_alloc_branch <<- cstore * albranch * pcbnew
    p_to_alloc_croot <<- cstore * alcroot * pccnew

    # Calculate remaining P left to allocate to leaves and roots
    ptot <- max(0.0, (pstore - p_to_alloc_stemimm - p_to_alloc_stemmob - p_to_alloc_branch))

    # allocate remaining P to flexible-ratio pools
    p_to_alloc_shoot <<- (ptot * alleaf / (alleaf + alroot * pcrfac))
    p_to_alloc_root <<- ptot - p_to_alloc_shoot
}

SMA_ADD <- function (value) {
    sma_obj$lv <<- sma_obj$lv + 1
     if (sma_obj$lv <= sma_obj$period) {
         sma_obj$values[sma_obj$lv] <<- value
         sma_obj$sum <<- sma_obj$sum + value
         sma_obj$sma <<- sma_obj$sum / sma_obj$lv
     } else {
         if (sma_obj$lv %% sma_obj$period == 0) {
            idx <- sma_obj$period
         } else {
            idx <- sma_obj$lv %% sma_obj$period
         }
        sma_obj$sum <<- sma_obj$sum - sma_obj$values[idx]
        sma_obj$sum <<- sma_obj$sum + value
        sma_obj$sma <<- sma_obj$sum / sma_obj$period
        sma_obj$values[idx] <<- value
     }
}

correct_rate_constants <- function(output) {
    # adjust rate constants for the number of days in years

    if (output) {
        rateuptake <<- rateuptake * NDAYS_IN_YR
        prateuptake <<- prateuptake * NDAYS_IN_YR
        rateloss <<- rateloss * NDAYS_IN_YR
        prateloss <<- prateloss * NDAYS_IN_YR
        retransmob <<- retransmob * NDAYS_IN_YR
        pfdecay <<- pfdecay * NDAYS_IN_YR
        fdecaydry <<- fdecaydry * NDAYS_IN_YR
        crdecay <<- crdecay * NDAYS_IN_YR
        prdecay <<- prdecay * NDAYS_IN_YR
        rdecaydry <<- rdecaydry * NDAYS_IN_YR
        bdecay <<- bdecay * NDAYS_IN_YR
        wdecay <<- wdecay * NDAYS_IN_YR
        sapturnover <<- sapturnover * NDAYS_IN_YR
        kdec1 <<- kdec1 * NDAYS_IN_YR
        kdec2 <<- kdec2 * NDAYS_IN_YR
        kdec3 <<- kdec3 * NDAYS_IN_YR
        kdec4 <<- kdec4 * NDAYS_IN_YR
        kdec5 <<- kdec5 * NDAYS_IN_YR
        kdec6 <<- kdec6 * NDAYS_IN_YR
        kdec7 <<- kdec7 * NDAYS_IN_YR
        nuptakez <<- nuptakez * NDAYS_IN_YR
        puptakez <<- puptakez * NDAYS_IN_YR
        nmax <<- nmax * NDAYS_IN_YR
        pmax <<- pmax * NDAYS_IN_YR
        p_atm_deposition <<- p_atm_deposition * NDAYS_IN_YR
        p_rate_par_weather <<- p_rate_par_weather * NDAYS_IN_YR
        max_p_biochemical <<- max_p_biochemical * NDAYS_IN_YR
        rate_sorb_ssorb <<- rate_sorb_ssorb * NDAYS_IN_YR
        rate_ssorb_occ <<- rate_ssorb_occ * NDAYS_IN_YR
    } else {
        rateuptake <<- rateuptake / NDAYS_IN_YR
        prateuptake <<- prateuptake / NDAYS_IN_YR
        rateloss <<- rateloss / NDAYS_IN_YR
        prateloss <<- prateloss / NDAYS_IN_YR
        retransmob <<- retransmob / NDAYS_IN_YR
        pfdecay <<- pfdecay / NDAYS_IN_YR
        fdecaydry <<- fdecaydry / NDAYS_IN_YR
        crdecay <<- crdecay / NDAYS_IN_YR
        prdecay <<- prdecay / NDAYS_IN_YR
        rdecaydry <<- rdecaydry / NDAYS_IN_YR
        bdecay <<- bdecay / NDAYS_IN_YR
        wdecay <<- wdecay / NDAYS_IN_YR
        sapturnover <<- sapturnover / NDAYS_IN_YR
        kdec1 <<- kdec1 / NDAYS_IN_YR
        kdec2 <<- kdec2 / NDAYS_IN_YR
        kdec3 <<- kdec3 / NDAYS_IN_YR
        kdec4 <<- kdec4 / NDAYS_IN_YR
        kdec5 <<- kdec5 / NDAYS_IN_YR
        kdec6 <<- kdec6 / NDAYS_IN_YR
        kdec7 <<- kdec7 / NDAYS_IN_YR
        nuptakez <<- nuptakez / NDAYS_IN_YR
        puptakez <<- puptakez / NDAYS_IN_YR
        nmax <<- nmax / NDAYS_IN_YR
        pmax <<- pmax / NDAYS_IN_YR
        p_atm_deposition <<- p_atm_deposition / NDAYS_IN_YR
        p_rate_par_weather <<- p_rate_par_weather / NDAYS_IN_YR
        max_p_biochemical <<- max_p_biochemical / NDAYS_IN_YR
        rate_sorb_ssorb <<- rate_sorb_ssorb / NDAYS_IN_YR
        rate_ssorb_occ <<- rate_ssorb_occ / NDAYS_IN_YR
    }
}

day_end_calculations <- function(days_in_year, init) {
    # Calculate derived values from state variables.

    # Parameters:
    # -----------
    # day : integer
    #     day of simulation
    #
    # INIT : logical
    #     logical defining whether it is the first day of the simulation

    # update N:C and P:C of plant pool
    if (float_eq(shoot, 0.0)) {
        shootnc <<- 0.0
        shootpc <<- 0.0
    } else {
        shootnc <<- shootn / shoot
        shootpc <<- shootp / shoot
    }

    # Explicitly set the shoot N:C
    if (ncycle == FALSE)
        shootnc <<- prescribed_leaf_NC

    if (pcycle == FALSE)
        shootpc <<- prescribed_leaf_PC

    if (float_eq(root, 0.0)) {
        rootnc <<- 0.0
        rootpc <<- 0.0
    } else {
        rootnc <<- max(0.0, rootn / root)
        rootpc <<- max(0.0, rootp / root)
    }

    # total plant, soil & litter nitrogen
    soiln <<- inorgn + activesoiln + slowsoiln + passivesoiln
    litternag <<- structsurfn + metabsurfn
    litternbg <<- structsoiln + metabsoiln
    littern <<- litternag + litternbg
    plantn <<- shootn + rootn + crootn + branchn + stemn
    totaln <<- plantn + littern + soiln

    # total plant, soil & litter phosphorus
    inorgp <<- inorglabp + inorgsorbp + inorgssorbp + inorgoccp + inorgparp
    soilp <<- inorgp + activesoilp + slowsoilp + passivesoilp
    litterpag <<- structsurfp + metabsurfp
    litterpbg <<- structsoilp + metabsoilp
    litterp <<- litterpag + litterpbg
    plantp <<- shootp + rootp + crootp + branchp + stemp
    totalp <<- plantp + litterp + soilp # + inorgssorbp + inorgoccp + inorgparp

    # total plant, soil, litter and system carbon
    soilc <<- activesoil + slowsoil + passivesoil
    littercag <<- structsurf + metabsurf
    littercbg <<- structsoil + metabsoil
    litterc <<- littercag + littercbg
    plantc <<- root + croot + shoot + stem + branch
    totalc <<- soilc + litterc + plantc

    # optional constant passive pool
    if (passiveconst) {
        passivesoil <<- passivesoilz
        passivesoiln <<- passivesoilnz
        passivesoilp <<- passivesoilpz
    }

    if (init == FALSE)
        # Required so max leaf & root N:C can depend on Age
        age <<- age + 1.0 / days_in_year
}

get_soil_fracs <- function(soil_type) {
     # Based on Table 2 in Cosby et al 1984, page 2.
     # Fractions of silt, sand and clay (in that order)

    fsoil <- rep(0, 3)

    if (soil_type == "sand") {
        fsoil[1] <- 0.05
        fsoil[2] <- 0.92
        fsoil[3] <- 0.03
    } else if (soil_type == "loamy_sand") {
        fsoil[1] <- 0.12
        fsoil[2] <- 0.82
        fsoil[3] <- 0.06
    } else if (soil_type == "sandy_loam") {
        fsoil[1] <- 0.32
        fsoil[2] <- 0.58
        fsoil[3] <- 0.1
    } else if (soil_type == "loam") {
        fsoil[1] <- 0.39
        fsoil[2] <- 0.43
        fsoil[3] <- 0.18
    } else if (soil_type == "silty_loam") {
        fsoil[1] <- 0.7
        fsoil[2] <- 0.17
        fsoil[3] <- 0.13
    } else if (soil_type == "sandy_clay_loam") {
        fsoil[1] <- 0.15
        fsoil[2] <- 0.58
        fsoil[3] <- 0.27
    } else if (soil_type == "clay_loam") {
        fsoil[1] <- 0.34
        fsoil[2] <- 0.32
        fsoil[3] <- 0.34
    } else if (soil_type == "silty_clay_loam") {
        fsoil[1] <- 0.56
        fsoil[2] <- 0.1
        fsoil[3] <- 0.34
    } else if (soil_type == "sandy_clay") {
        fsoil[1] <- 0.06
        fsoil[2] <- 0.52
        fsoil[3] <- 0.42
    } else if (soil_type == "silty_clay") {
        fsoil[1] <- 0.47
        fsoil[2] <- 0.06
        fsoil[3] <- 0.47
    } else if (soil_type == "clay") {
        fsoil[1] <- 0.2
        fsoil[2] <- 0.22
        fsoil[3] <- 0.58
    } else {
        stop(paste("Could not understand soil type:" + soil_type))
    }
    return (fsoil)
}

get_soil_params <- function(soil_type) {
    # For a given soil type, get the parameters for the soil
    # moisture availability based on Landsberg and Waring, with updated
    # parameters from Landsberg and Sands (2011), pg 190, Table 7.1
    #
    # Table also has values from Saxton for soil texture, perhaps makes more
    # sense to use those than Cosby? Investigate?
    #
    # Reference
    # ---------
    # Landsberg and Sands (2011) Physiological ecology of forest production.
    # Landsberg and Waring (1997) Forest Ecology & Management, 95, 209-228.
    # Clapp & Hornberger (1978) Water Resources Research, 14, 601–604.

    if (soil_type == "clay") {
        c_theta <- 0.4
        n_theta <- 3.0
    } else if (soil_type == "clay_loam") {
        c_theta <- 0.5
        n_theta <- 5.0
    } else if (soil_type == "loam") {
        c_theta <- 0.55
        n_theta <- 6.0
    } else if (soil_type == "loamy_sand") {
        c_theta <- 0.65
        n_theta <- 8.0
    } else if (soil_type == "sand") {
        c_theta <- 0.7
        n_theta <- 9.0
    } else if (soil_type == "sandy_clay") {
        c_theta <- 0.45
        n_theta <- 4.0
    } else if (soil_type == "sandy_clay_loam") {
        c_theta <- 0.525
        n_theta <- 5.5
    } else if (soil_type == "sandy_loam") {
        c_theta <- 0.6
        n_theta <- 7.0
    } else if (soil_type == "silt") {
        c_theta <- 0.625
        n_theta <- 7.5
    } else if (soil_type == "silty_clay") {
        c_theta <- 0.425
        n_theta <- 3.5
    } else if (soil_type == "silty_clay_loam") {
        c_theta <- 0.475
        n_theta <- 4.5
    } else if (soil_type == "silty_loam") {
        c_theta <- 0.575
        n_theta <- 6.5
    } else {
        stop(paste("There are no parameters for your soil type:", soil_type))
    }

    return (c(c_theta, n_theta))
}

calc_soil_params <- function(fsoil) {
    # Cosby parameters for use within the Clapp Hornberger soil hydraulics
    # scheme are calculated based on the texture components of the soil.
    #
    # NB: Cosby et al were ambiguous in their paper as to what log base to
    # use.  The correct implementation is base 10, as below.
    #
    # Parameters:
    # ----------
    # fsoil : list
    #     fraction of silt, sand, and clay (in that order
    #
    # Returns:
    # --------
    # theta_fc : float
    #     volumetric soil water concentration at field capacity
    # theta_wp : float
    #     volumetric soil water concentration at the wilting point
    #
    # soil suction of 3.364m and 152.9m, or equivalent of -0.033 & -1.5 MPa
    pressure_head_wilt <- -152.9
    pressure_head_crit <- -3.364

    # *Note* subtle unit change to be consistent with fractions as opposed
    # to percentages of sand, silt, clay, e.g. I've changed the slope in
    # the "b" Clapp paramter from 0.157 to 15.7
    #
    # Also Cosby is unclear about which log base were used. 'Generally' now
    # assumed that logarithms to the base 10

    # Clapp Hornberger exponent [-]
    b <- 3.1 + 15.7 * fsoil[CLAY] - 0.3 * fsoil[SAND]

    # soil matric potential at saturation, taking inverse of log (base10)
    # units = m

    psi_sat <- CM_2_M * -(10.0 ^ (1.54 - 0.95 * fsoil[SAND] + 0.63 * fsoil[SILT]))
    psi_sat_mpa <- psi_sat * METER_OF_HEAD_TO_MPA

    # volumetric soil moisture concentrations at the saturation point
    theta_sp <- 0.505 - 0.037 * fsoil[CLAY] - 0.142 * fsoil[SAND]

    # volumetric soil moisture concentrations at the wilting point
    # assumed to equal suction of -1.5 MPa or a depth of water of 152.9 m

    theta_wp <- theta_sp * ((psi_sat / pressure_head_wilt) ^ (1.0 / b))

    # volumetric soil moisture concentrations at field capacity assumed to
    # equal a suction of -0.0033 MPa or a depth of water of 3.364 m

    theta_fc <- theta_sp * ((psi_sat / pressure_head_crit) ^ (1.0 / b))

    return (c(theta_fc, theta_wp, theta_sp, b, psi_sat_mpa))

}

zbrent <- function(func, x1, x2, tol, root_biomass,
                   surf_biomass, rooted_layers, top_lyr_thickness, root_reach) {

    #Using Brent’s method, find the root of a function func known to lie
    #between x1 and x2. The root, returned as zbrent, will be refined
    #until its accuracy is tol.

    #Numerical Recipies in C, chapter 9.3

    a <- x1
    b <- x2
    fa <- func(a, root_biomass, surf_biomass, rooted_layers, top_lyr_thickness, root_reach)
    fb <- func(b, root_biomass, surf_biomass, rooted_layers, top_lyr_thickness, root_reach)
    ITMAX <- 100   # Maximum allowed number of iterations.
    EPS <- 3.0e-8  # Machine floating-point precision.

    if (fb * fa > 0.0) {
        stop("ERROR: Root must be bracketed in ZBRENT")
	}
	fc <- fb
	for (iter in 1:ITMAX) {
        if (fb * fc > 0.0) {
            c <- a
            fc <- fa
            d <- b - a
            e <- b - a
        }
        if (abs(fc) < abs(fb)) {
            a <- b
            b <- c
            c <- a
            fa <- fb
            fb <- fc
            fc <- fa
        }
    	tol1 <- 2.0 * EPS * abs(b) + 0.5 * tol
    	xm <- 0.5 * (c-b)
    	if (abs(xm) <= tol1 || fb == 0.0)
            return (b)
    	if (abs(e) >= tol1 && abs(fa) > abs(fb)) {
            s <- fb / fa
            if (a == c) {
                p <- 2.0 * xm * s
                q <- 1.0 - s
            } else {
                q <- fa / fc
                r <- fb / fc
                p <- s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0))
                q <- (q - 1.0) * (r - 1.0) * (s - 1.0)
            }
            if (p > 0.0)  q <- -q
            p <- abs(p)
            min1 <- 3.0 * xm * q - abs(tol1 * q)
            min2 <- abs(e*q)
            if (2.0 * p < min(min1, min2)) {
                e <- d
                d <- p / q
    		} else {
                d <- xm
                e <- d
    		}
    	} else {
            d <- xm
            e <- d
    	}
    	a <- b
    	fa <- fb
    	if (abs(d) > tol1) {
            b <- b + d
        } else {
            b <- b + (if (xm > 0.0) abs(tol1) else -abs(tol1))
        }
    	fb <- func(b, root_biomass, surf_biomass, rooted_layers,
                   top_lyr_thickness, root_reach)
    }

    stop("Maximum number of iterations exceeded in ZBRENT")
}

saxton_field_capacity <- function(xval, potA, potB, dummy1, dummy2, dummy3) {
    # field capacity calculations for saxton eqns

    air_entry_swp <- 10.0 # (kPa)
    MPa_2_kPa <- 1000.0

    # calculate the soil water potential (MPa)
    swp <- -0.001 * potA * (xval ^ potB)
    return (MPa_2_kPa * swp + air_entry_swp)
}

calc_saxton_stuff <- function(fsoil) {
    # Calculate the key parameters of the Saxton equations:
    # cond1, 2, 3 and potA, B
    #
    # NB: Currently I'm assuming a single texture for the entire root zone,
    #     we could set it by layer obviously and we probably should...
    #
    # Reference:
    # ---------
    # * Saxton, K. E., Rawls, W. J., Romberger, J. S. & Papendick, R. I.
    #   (1986) Estimating generalized soil-water characteristics from texture.
    #   Soil Science Society of America Journal, 90, 1031-1036.
    mult1 <- 100.0
    mult2 <- 2.778E-6
    mult3 <- 1000.0
    A <- -4.396
    B <- -0.0715
    CC <- -4.880E-4
    D <- -4.285E-5
    E <- -3.140
    F <- -2.22E-3
    G <- -3.484E-5
    H <- 0.332
    J <- -7.251E-4
    K <- 0.1276
    P <- 12.012
    Q <- -7.551E-2
    R <- -3.895
    T <- 3.671e-2
    U <- -0.1103
    V <- 8.7546E-4
    x1 <- 0.1
    x2 <- 0.7
    tol <- 0.0001
    dummy <- 0.0
    sand <- fsoil[SAND] * 100.0
    clay <- fsoil[CLAY] * 100.0

    # As we aren't currently changing texture by layer, the loop is redundant,
    # but it is probably best to leave it under the assumption we change
    # this later

    for (i in 1:n_layers) {
        potA[i] <<- exp(A + B * clay + CC * sand * sand + D * sand * sand * clay) * 100.0
        potB[i] <<- E + F * clay * clay + G * sand * sand * clay
        cond1[i] <<- mult2
        cond2[i] <<- P + Q * sand
        cond3[i] <<- R + T * sand + U * clay + V * clay * clay

        porosity[i] <<- H + J * sand + K * log10(clay)

        # field capacity is water content at which SWP = -10 kPa
        field_capacity[i] <<- zbrent(saxton_field_capacity, x1, x2, tol,
                                      potA[i], potB[i],
                                      dummy, dummy, dummy)
    }
}

calc_soil_conductivity <- function(water_frac, cond1, cond2, cond3) {
    # Soil hydraulic conductivity (m s-1 ) per soil layer based on
    # Saxton et al. (1986) equations. Used in the soil drainage integrator

    # Avoid floating-underflow error
    if (water_frac < 0.05) {
        scond <- 1E-30
    } else {
        scond <- cond1 * exp(cond2 + cond3 / water_frac)
    }
    return (scond)
}

calc_soil_root_resistance <- function() {

    # head of pressure (MPa/m)
    head <- 0.009807
    root_xsec_area <- pi * root_radius * root_radius

    # Store each layers resistance, used in LWP calculatons
    rsum <- 0.0
    for (i in 1:rooted_layers) {

        # converts from ms-1 to m2 s-1 MPa-1
        Lsoil <- soil_conduct[i] / heads

        if (Lsoil < 1e-35) {
            # prevent floating point error
            soilR[i] <<- 1e35
        } else {
            rs <- sqrt(1.0 / (root_length[i] * M_PI))
            arg1 <- log(rs / root_radius)
            arg2 <- 2.0 * pi * root_length[i] * thickness[i] * Lsoil

            # soil resistance, convert from MPa s m2 m-3 to MPa s m2 mmol-1
            soilR1 <- (arg1 / arg2) * 1E-6 * 18. * 0.001

            # Need to combine resistances in parallel, but we only want the
            # soil term as the root component is part of the plant resistance
            rsum <- rsum + 1.0 / soilR1

            # second component of below ground resistance related to root
            # hydraulics
            soilR2 <- root_resist / (root_mass[i] * thickness[i])
            # soilR[i] <- soilR1 + soilR2 # MPa s m2 mmol-1
            soilR[i] <<- soilR1 + soilR2 # MPa s m2 mmol-1
        }
    }

    total_soil_resist <<- 1.0 / rsum
}

calc_soil_water_potential <- function() {
    # Calculate the SWP (MPa) in each soil layer based on algorithms from
    # Saxton et al. (1986). We are not updating the water fraction in each
    # layer

    for (i in 1:rooted_layers) {
        if (water_frac[i] > 0.0) {
            swp[i] <- -0.001 * potA[i] * (water_frac[i] ^ potB[i])
        } else {
            swp[i] <- -9999.0
        }
    }
}

calc_water_uptake_per_layer <- function() {
    # Determine which layer water is extracted from. This is achieved by
    # roughly estimating the maximum rate of water supply from each rooted
    # soil layer, using SWP and hydraulic resistance of each layer. Actual
    # water from each layer is determined using the estimated value as a
    # weighted factor.

    total_est_evap <- 0.0
    weighted_swp <<- 0.0

    # Estimate max transpiration from gradient-gravity / soil resistance
    for (i in 1:rooted_layers) {
        est_evap[i] <<- max(0.0, (swp[i] - min_lwp) / soilR[i])
        total_est_evap <- total_est_evap + est_evap[i]
    }

    if (total_est_evap > 0.0) {
        # fraction of water taken from layer
        for (i in 1:rooted_layers) {
            weighted_swp <<- weighted_swp + swp[i] * est_evap[i]
            fraction_uptake[i] <<- est_evap[i] / total_est_evap
        }
        weighted_swp <<- weighted_swp / total_est_evap
    } else {
        # No water was evaporated
        for (i in 1:rooted_layers) {
            fraction_uptake[i] <<- 1.0 / rooted_layers
        }
    }

    if (fraction_uptake[1] > 1 || fraction_uptake[1] < 0) {
        stop("Problem with the uptake fraction")
    }
}

initialise_soils_sub_daily <- function() {
    # Initialise soil water state & parameters
    #     - We have two options: (i) simple two-layer bucket, or
    #     - multi-layer soil model.
    #
    # Currently, I'm just using the same soil type for all layers in the
    # multi-layer model, once it all works we can add something to allow
    # texture to be changes by layer/depth.
    #
    # NB. get_soil_fracs, calc_soil_params, get_soil_params live in
    #     water_balance.c

    # site params not known, so derive them based on Cosby et al
    if (calc_sw_params) {
        fsoil_top <- get_soil_fracs(topsoil_type)
        fsoil_root <- get_soil_fracs(rootsoil_type)

        # top soil
        res__ <- calc_soil_params(fsoil_top)
        theta_fc_topsoil <<- res__[1]
        theta_wp_topsoil <<- res__[2]
        theta_sp_topsoil <<- res__[3]
        b_topsoil <<- res__[4]
        psi_sat_topsoil <<- res__[5]

        # Plant available water in top soil (mm)
        wcapac_topsoil <<- topsoil_depth * (theta_fc_topsoil - theta_wp_topsoil)
        # Root zone
        res__ <- calc_soil_params(fsoil_root)
        theta_fc_root <<- res__[1]
        theta_wp_root <<- res__[2]
        theta_sp_root <<- res__[3]
        b_root <<- res__[4]
        psi_sat_root <<- res__[5]

        # Plant available water in rooting zone (mm)
        wcapac_root <<- rooting_depth * (theta_fc_root - theta_wp_root)
    }


    # calculate Landsberg and Waring SW modifier parameters if not
    #    specified by the user based on a site calibration
    if (ctheta_topsoil < -900.0 && ntheta_topsoil < -900.0 && ctheta_root < -900.0 && ntheta_root < -900.0) {
        res__ <- get_soil_params(topsoil_type)
        ctheta_topsoil <<- res__[1]
        ntheta_topsoil <<- res__[2]
        res__ <- get_soil_params(rootsoil_type)
        ctheta_root <<- res__[1]
        ntheta_root <<- res__[2]
    }

    # Set up all the hydraulics stuff
    if (water_balance == HYDRAULICS) {
        calc_saxton_stuff(fsoil_root)

        for (i in 1:wetting) {
            wetting_bot[i] <<- 0.0
            wetting_top[i] <<- 0.0
        }

        # saturate the top layer
        wetting_bot[1] <<- thickness[1]

        # Initalise SW fraction - we should read this from param file
        initial_water <<- 0.0
        for (i in 1:n_layers) {
            water_frac[i] <<- 0.4
            initial_water <<- initial_water + 1E3 * (water_frac[i] * thickness[i])
        }

        # The loop needs to be outside the func as we need to be able to
        # calculate the soil conductance per layer and call this via
        # the integration func when we update the soil water balance

        for (i in 1:n_layers) {
            soil_conduct[i] <<- calc_soil_conductivity(water_frac[i], cond1[i], cond2[i], cond3[i])
        }

        calc_soil_root_resistance()
        calc_soil_water_potential()

        # Calculate the weighted soil-water-potential
        calc_water_uptake_per_layer()
    }
}

initialise_soils_day <- function() {
    # Initialise soil water state & parameters

    # site params not known, so derive them based on Cosby et al

    if (calc_sw_params) {
        fsoil_top <- get_soil_fracs(topsoil_type)
        fsoil_root <- get_soil_fracs(rootsoil_type)

        # top soil
        res__ <- calc_soil_params(fsoil_top)
        theta_fc_topsoil <<- res__[1]
        theta_wp_topsoil <<- res__[2]
        theta_sp_topsoil <<- res__[3]
        b_topsoil <<- res__[4]
        psi_sat_topsoil <<- res__[5]

        # Plant available water in top soil (mm)
        wcapac_topsoil <<- topsoil_depth * (theta_fc_topsoil - theta_wp_topsoil)

        # Root zone
        res__ <- calc_soil_params(fsoil_root)
        theta_fc_root <<- res__[1]
        theta_wp_root <<- res__[2]
        theta_sp_root <<- res__[3]
        b_root <<- res__[4]
        psi_sat_root <<- res__[5]

        # Plant available water in rooting zone (mm)
        wcapac_root <<- rooting_depth * (theta_fc_root - theta_wp_root)
    }


        # calculate Landsberg and Waring SW modifier parameters if not
        # specified by the user based on a site calibration
    if (ctheta_topsoil < -900.0 && ntheta_topsoil  < -900.0 && ctheta_root < -900.0 && ntheta_root < -900.0) {
        res__ <- get_soil_params(topsoil_type)
        ctheta_topsoil <<- res__[1]
        ntheta_topsoil <<- res__[2]
        res__ <- get_soil_params(rootsoil_type)
        ctheta_root <<- res__[1]
        ntheta_root <<- res__[2]
    }
}

figure_out_years_with_disturbances <- function() {
    yrs <- 0
    cnt <- 0
    if (burn_specific_yr < -900.0) {
        year_of_disturbance <- burn_specific_yr
        yrs[1] <- burn_specific_yr
    } else {
        yrs_till_event <- time_till_next_disturbance()
        year <- ma$year[day_idx]
        # year_of_disturbance = year + yrs_till_event
        year_of_disturbance <- 1996

        # figure out the years of the disturbance events
        prjday <- 0

        for (nyr in 1:(num_years - 1)) {
            year <- ma$year[day_idx]
            if (is_leap_year(year))
                prjday <- prjday + 366
            else
                prjday <- prjday + 365

            if (year == year_of_disturbance) {
                yrs_till_event <- time_till_next_disturbance()
                if (cnt == 0) {
                    yrs[1] <- year_of_disturbance
                    cnt <- cnt + 1
                } else {
                    cnt <- cnt + 1
                    yrs <- append(yrs, year_of_disturbance)
                }

                # See if there is another event?
                year_of_disturbance <- year + yrs_till_event
            }
        }
    }

    return (list(yrs, cnt))
}

calc_ini_grass_pheno_stuff <- function(project_day) {
    # Series of constraints based on temp && precip need to be
    # pre-calculated for grasses to determine leaf on/off

    # Save this as we need to loop over the data once to pre-calculate
    # everything

    project_day_save <- project_day
    tmax_ann <- 0.0
    Tmin_avg <- 0.0
    tmin_ann <- 70.0
    tavg_ann <- 0.0
    ppt_sum <- 0.0

    for (d in 1:num_days) {
        tair <- ma$tair[project_day]
        tmax <- ma$tmax[project_day]
        tmin <- ma$tmin[project_day]
        ppt_sum <- ppt_sum + ma$rain[project_day]
        Tmin_avg <- Tmin_avg + ma$tmin[project_day]

        if (tmax > tmax_ann)
           tmax_ann <- tmax

        if (tmin < tmin_ann)
           tmin_ann <- tmin

        tavg_ann <- tavg_ann + tair
        project_day <- project_day + 1
    }
    Tmin_avg <- Tmin_avg / num_days

    # reset date index
    project_day <- project_day_save

    Trange <- tmax_ann - tmin_ann
    tavg_ann <- tavg_ann / num_days

    # Cool or warm grassland Definitions are from Botta, Table 1, pg 712.
    # But grass temp thresholds are from Foley et al.

    if (Trange > 20.0 || tmin_ann < 5.0) {
        # cool
        grass_temp_threshold <- 0.0
    } else if (Trange <= 20.0 || tmin_ann >= 5.0) {
        # warm
        grass_temp_threshold <- 5.0
    } else {
        stop("Problem grass thresholds\n")
    }

    # 92% of tmax_ann is the threshold used in grass offset below
    # Note this has to be done below the range calcs as they use the tmax

    tmax_ann <- tmax_ann * 0.92

    # Based on White et al. 1997 this is a threshold for grasses so they
    # have enough accumulated rain. It is essentially a fudge for soil
    # moisture availability && the 15% is somewhat arbitary
    ppt_sum_crit <- ppt_sum * 0.15

    return(c(grass_temp_threshold, tmax_ann, Tmin_avg, ppt_sum_crit))
}

calculate_leafon_off <- function(
    grass_temp_threshold, tmax_ann,
    Tmin_avg, ppt_sum_crit,
    project_day, gdd_thresh) {

    # Tmax
    accumulated_ncd <- 0.0
    accum_gdd <- 0.0
    # Tmin_boxcar
    drop_leaves <- FALSE


    leaf_on_found <- FALSE
    leaf_off_found <- FALSE

    if (num_days == 366) {
        nov_doy <- 307
    } else {
        nov_doy <- 306
    }

    ppt_sum <- 0.0
    for (d in 1:num_days) {
        Tmean <- ma$tair[project_day]
        Tday <- ma$tday[project_day]
        Tsoil <- ma$tsoil[project_day]
        # Tmax = ma$tmax[project_day]
        ppt_sum <- ppt_sum + ma$rain[project_day]

        # Calculate ppt total from the next 7 days
        if (d <= 358) {
            st <- project_day + 1
            en <- project_day + 8
            ppt_sum_next <- 0.0
            for (dd in st:(en - 1)) {
                ppt_sum_next <- ppt_sum_next + ma$rain[dd]
            }
        } else {
            # i.e. end of year, didn't find this so have no effect
            ppt_sum_next <- 0.0
        }

        # Calculate ppt total from the previous 30 days
        if (project_day <= 30) {
            ppt_sum_prev <- 0.0
        } else {
            st <- project_day - 30
            en <- project_day
            ppt_sum_prev <- 0.0
            for (dd in st:(en - 1)) {
                ppt_sum_prev <- ppt_sum_prev + ma$rain[dd]
            }
        }

        if (d <= 362) {
            Tsoil_next_3days <- ((ma$tsoil[project_day] +
                                 ma$tsoil[project_day + 1] +
                                 ma$tsoil[project_day + 2]) / 3.0)

            #Tmin_boxcar = ((ma$tmin[project_day-1] +
                            # ma$tmin[project_day] +
                            # ma$tmin[project_day+1]) / 3.0)

        } else {
            # i.e. end of year, didn't find this so have no effect
            Tsoil_next_3days <- 999.9
        }

        # Sum the daily mean air temperature above 5degC starting on Jan 1
        accum_gdd <- accum_gdd + calc_gdd(Tmean)

        # Calculate leaf on
        if (alloc_model == GRASSES) {
            if (leaf_on_found == FALSE &&
                accum_gdd >= gdd_thresh &&
                ppt_sum >= ppt_sum_crit) {

                leaf_on <- d
                leaf_on_found <- TRUE
            }
        } else {
            if (leaf_on_found == FALSE && accum_gdd >= gdd_thresh) {
                  leaf_on <- d
                  leaf_on_found <- TRUE
            }
        }

        # Calculate leaf off
        if (alloc_model == GRASSES) {
            if (leaf_off_found == FALSE) {

                # Leaf drop constraint is based on Foley et al. 1996
                #
                # The 243 is just a safe guard to make sure we avoid
                # predicting offset in late spring (White et al. 1997).
                #
                # This Tmean is the mean daytime temp, but I wonder if it
                # should be the full 24 daytime temp mean?

                # if (d >= 243) {
                #     printf("%d %f %f\n", d, Tmean, grass_temp_threshold)
                # }
                if (d > 243 && Tday <= grass_temp_threshold) {
                    leaf_off_found <- TRUE
                    leaf_off <- d
                }


                # Leaf drop constraint is based on White et al. 1997
                #
                #  - test for hot && dry conditions.
                #
                # if (ppt_sum_prev < 11.4 &&
                #     ppt_sum_next < 9.7 &&
                #     Tmax > tmax_ann &&
                #     d > 243) {
                #
                #     *leaf_off_found = TRUE
                #     *leaf_off = d
                #
                #     - test for cold offset condition
                # } else if (d > 243 && Tmin_boxcar <= Tmin_avg) {
                #
                #     *leaf_off_found = TRUE
                #     *leaf_off = d
                #
                # }
            }
        } else {
            if (leaf_off_found == FALSE && accum_gdd >= gdd_thresh) {
                # I am prescribing that no leaves can fall off before doy=180
                # Had issue with KSCO simulations where the photoperiod was
                # less than the threshold very soon after leaf out.
                if (d > 183) {
                    drop_leaves <- leaf_drop(sday_length[d-1], Tsoil, Tsoil_next_3days)
                    if (drop_leaves) {
                        leaf_off_found <- TRUE
                        leaf_off <- d
                    }
                }
            }
        }
        # Calculated NCD from fixed date following Murray et al 1989.
        if (d + 1 >= nov_doy) {
            accumulated_ncd <- accumulated_ncd + calc_ncd(Tmean)
        }
        project_day <- project_day + 1
    }


     # updated stored param, note this will be written out if the user
     #   dumps the current state, which makes sense as we may want pass the
     #   stat between spinup and a simulation
    previous_ncd <<- accumulated_ncd

    return (c(leaf_on, leaf_off, leaf_on_found, leaf_off_found))
}

calculate_days_left_in_growing_season <- function(leaf_on, leaf_off, len_groloss) {
    # Calculate 2 arrays to signify the days left of growing period
    # an days left before all the leaves fall off. In both cases these will
    # be 2 lists, with 0.0 outside of the growing period and a series of
    # numbers e.g. day 46 to 0 for growing_days. 0.5 is subtracted from the
    # doy to get round the issue of approximating an integral with discrete
    # time steps trapezoidal type solution

    for (doy in 1:num_days) {

        if (doy > leaf_off - len_groloss && doy <= leaf_off) {
            remaining_days[doy] <<- (doy - 0.5) - leaf_off + len_groloss
        } else {
            remaining_days[doy] <<- 0.0
        }

        if (doy > leaf_on && doy <= len_groloss + leaf_on) {
            growing_days[doy] <<- len_groloss + leaf_on - (doy - 0.5)
        } else {
            growing_days[doy] <<- 0.0
        }

        if (doy > leaf_on && doy < leaf_off) {
            leaf_out_days[doy] <<- 1.0
        } else {
            leaf_out_days[doy] <<- 0.0
        }
    }
}

calculate_growing_season_fluxes <- function(len_groloss) {

    denominator <- len_groloss * len_groloss

    # C allocation rates across growing season
    lrate <<- 2.0 * c_to_alloc_shoot / denominator
    wrate <<- 2.0 * c_to_alloc_stem / denominator
    brate <<- 2.0 * c_to_alloc_branch / denominator
    crate <<- 2.0 * c_to_alloc_croot / denominator

    # N allocation rates across growing season
    lnrate <<- 2.0 * n_to_alloc_shoot / denominator
    bnrate <<- 2.0 * n_to_alloc_branch / denominator
    wnimrate <<- 2.0 * n_to_alloc_stemimm / denominator
    wnmobrate <<- 2.0 * n_to_alloc_stemmob / denominator
    cnrate <<- 2.0 * n_to_alloc_croot / denominator

    # P allocation rates across growing season
    lprate <<- 2.0 * p_to_alloc_shoot / denominator
    bprate <<- 2.0 * p_to_alloc_branch / denominator
    wpimrate <<- 2.0 * p_to_alloc_stemimm / denominator
    wpmobrate <<- 2.0 * p_to_alloc_stemmob / denominator
    cprate <<- 2.0 * p_to_alloc_croot / denominator
}

phenology <- function() {
    # There are two phenology schemes currently implemented, one which should
    # generally be applicable for deciduous broadleaf forests && one for
    # short-lived grasses.
    #
    # The tree scheme is based on Botta et al. who calibrated their model against
    # EO data. The grasses scheme is really a combination of approaches, I've
    # tried to keep it is as simple as possible...user beware :)
    #
    # There are some key issues:
    #  - leaf drop for trees won't work where the daylength isn't applicable, i.e
    #    the tropics would be my guess.
    #  - the grass phenology drop takes no account of available water. It of
    #    course should, but this would need a complete re-write of the logic.
    #    Perhaps someone else can do that.
    #
    # One potential thing to look at if you ever get bored is:
    #
    #  - Look at Caldararu et al 2014, Phenology as a strategy for carbon
    #    optimality: a global model. Seems interesting at a quick glance.
    #
    # The distribution of C&N is pre-calculated here using a ramping function
    # based on leaf out/off dates.
    #
    # Finally, no account has been taken for the southern hemisphere! This won't
    # work there.
    #
    # References:
    # -----------
    # * Botta, A. et al. (2000) GCB, 6, 709-725.
    # * Foley, J. et al. (1996) GBC, 10, 603-628.
    # * Krinner, G. et al. (2005) GBC, 19, GB1015
    # * White, M. A. et al. (1997) GBC, 11, 217-234.

    # (days) Leaf flush params following Botta.
    pa <- -68.0

    # (days) Leaf flush params following Botta.
    pb <- 638.0

    # (1/days) Leaf flush params following Botta.
    pc <- -0.01

    leaf_on <- 0
    leaf_off <- 0
    len_groloss <- 0.0
    project_day <- day_idx
    gdd_thresh <- -999.9


    # Krinner et al. 2005, page 26, alternatively Foley et al. 1996 suggests
    # the same value = 100 for both pathways

    if (alloc_model == GRASSES) {
        if (ps_pathway == C3)
            gdd_thresh <- 185.
        else if (ps_pathway == C4)
            gdd_thresh <- 400.
        res__ <- calc_ini_grass_pheno_stuff(project_day)
        grass_temp_threshold <- res__[1]
        tmax_ann <- res__[2]
        Tmin_avg <- res__[3]
        ppt_sum_crit <- res__[4]
    } else {
        gdd_thresh <- gdd_chill_thresh(pa, pb, pc, previous_ncd)
    }


    res__ <- calculate_leafon_off(grass_temp_threshold, tmax_ann, Tmin_avg, ppt_sum_crit, project_day, gdd_thresh)
    leaf_on <- res__[1]
    leaf_off <- res__[2]
    leaf_on_found <- res__[3]
    leaf_off_found <- res__[4]

    # No leaf drop found, try a warmer temperature i.e. 5 instead of 0,
    # if this doesn't work there really is an issue (or there is no leaf
    # drop and we have an evergreen grass...

    if (leaf_off_found == FALSE) {
        grass_temp_threshold <- 5.0
        calculate_leafon_off(c, ma, p, s, grass_temp_threshold, tmax_ann,
                             Tmin_avg, ppt_sum_crit, project_day,
                             leaf_on, leaf_off, leaf_on_found,
                             leaf_off_found, gdd_thresh)
        res__ <- calculate_leafon_off(grass_temp_threshold, tmax_ann, Tmin_avg, ppt_sum_crit, project_day, gdd_thresh)
        leaf_on <- res__[1]
        leaf_off <- res__[2]
        leaf_on_found <- res__[3]
        leaf_off_found <- res__[4]
    }

    # if widening the temperature threshold didn't produce a suitable leaf
    # drop date we will follow biome-bgc and assume the leaves fall on the
    # last day

    if (leaf_off_found == FALSE) {
        leaf_off <- 365
    }



    if (leaf_on_found == FALSE) {
        stop("Problem in phenology leaf *ON* not found")
    }


    # Length of time taken for new growth from storage to be allocated.
    # This is either some site-specific calibration or the midpoint of the
    # length of the growing season. The litterfall takes place over an
    # identical period. Dividing by a larger number would increase the
    # rate the C&N is allocated.

    growing_seas_len <<- leaf_off - leaf_on
    if (store_transfer_len < -900) {
        len_groloss <- floor(growing_seas_len / 2.0)
    } else {
        len_groloss <- store_transfer_len
    }

    calculate_days_left_in_growing_season(leaf_on, leaf_off, len_groloss)
    calculate_growing_season_fluxes(len_groloss)

    # printf("%d %d\n", leaf_on, leaf_off)
}

zero_stuff <- function() {
    shoot <<- 0.0
    shootn <<- 0.0
    shootp <<- 0.0
    shootnc <<- 0.0
    shootpc <<- 0.0
    lai <<- 0.0
    cstore <<- 0.0
    nstore <<- 0.0
    pstore <<- 0.0
    anpp <<- 0.0

    if (deciduous_model) {
        avg_alleaf <<- 0.0
        avg_alroot <<- 0.0
        avg_alcroot <<- 0.0
        avg_albranch  <<- 0.0
        avg_alstem <<- 0.0
    }
}

unpack_met_data <- function(hod, day_length) {
    # unpack met forcing
    if (sub_daily) {
        rain <<- ma$rain[hour_idx]
        wind <<- ma$wind[hour_idx]
        press <<- ma$pres[hour_idx] * KPA_2_PA
        vpd <<- ma$vpd[hour_idx] * KPA_2_PA
        tair <<- ma$tair[hour_idx]
        tsoil <<- ma$tsoil[hour_idx]
        par <<- ma$par[hour_idx]
        sw_rad <<- ma$par[hour_idx] * PAR_2_SW # W m-2
        Ca <<- ma$co2[hour_idx]

        # NDEP is per 30 min so need to sum 30 min data
        if (hod == 0) {
            ndep <<- ma$ndep[hour_idx]
            nfix <<- ma$nfix[hour_idx]
            pdep <<- ma$pdep[hour_idx]
        } else {
            ndep <<- ndep + ma$ndep[hour_idx]
            nfix <<- nfix + ma$nfix[hour_idx]
            pdep <<- pdep + ma$pdep[hour_idx]
        }
    } else {
        Ca <<- ma$co2[day_idx]
        tair <<- ma$tair[day_idx]
        tair_am <<- ma$tam[day_idx]
        tair_pm <<- ma$tpm[day_idx]
        par <<- ma$par_am[day_idx] + ma$par_pm[day_idx]

        # Conversion factor for PAR to SW rad
        c1 <- MJ_TO_J * J_2_UMOL / (day_length * 60.0 * 60.0) * PAR_2_SW
        c2 <- MJ_TO_J * J_2_UMOL / (day_length / 2.0 * 60.0 * 60.0) * PAR_2_SW
        sw_rad <<- par * c1
        sw_rad_am <<- ma$par_am[day_idx] * c2
        sw_rad_pm <<- ma$par_pm[day_idx] * c2
        rain <<- ma$rain[day_idx]
        vpd_am <<- ma$vpd_am[day_idx] * KPA_2_PA
        vpd_pm <<- ma$vpd_pm[day_idx] * KPA_2_PA
        wind_am <<- ma$wind_am[day_idx]
        wind_pm <<- ma$wind_pm[day_idx]
        press <<- ma$pres[day_idx] * KPA_2_PA
        ndep <<- ma$ndep[day_idx]
        nfix <<- ma$nfix[day_idx]
        pdep <<- ma$pdep[day_idx]
        tsoil <<- ma$tsoil[day_idx]
        Tk_am <<- ma$tam[day_idx] + DEG_TO_KELVIN
        Tk_pm <<- ma$tpm[day_idx] + DEG_TO_KELVIN

    }

    # N deposition + biological N fixation
    ninflow <<- ndep + nfix

    # P deposition to fluxess
    p_atm_dep <<- p_atm_deposition
}

calculate_litterfall <- function(doy) {

    # Leaf/root litter rates are higher during dry periods and therefore is
    # dependent on soil water content
    fdecay <- decay_in_dry_soils(pfdecay, fdecaydry)
    rdecay <- decay_in_dry_soils(prdecay, rdecaydry)

    # litter N:C ratios, roots and shoot
    ncflit <- shootnc * (1.0 - fretrans)
    ncrlit <- rootnc * (1.0 - rretrans)

    # litter P:C ratios, roots and shoot
    pcflit <- shootpc * (1.0 - fretransp)
    pcrlit <- rootpc * (1.0 - rretrans)

    # C litter production
    deadroots <<- rdecay * root
    deadcroots <<- crdecay * croot
    deadstems <<- wdecay * stem
    deadbranch <<- bdecay * branch
    deadsapwood <<- (wdecay + sapturnover) * sapwood


    if (deciduous_model) {
        deadleaves <<- lrate * remaining_days[doy]
    } else {
        deadleaves <<- fdecay * shoot
    }

    # N litter production
    deadleafn <<- deadleaves * ncflit

    # P litter production
    deadleafp <<- deadleaves * pcflit

    # Assuming fraction is retranslocated before senescence, i.e. a fracion
    #    of nutrients is stored within the plant
    deadrootn <<- deadroots * ncrlit
    deadcrootn <<- crdecay * crootn * (1.0 - cretrans)
    deadbranchn <<- bdecay * branchn * (1.0 - bretrans)

    deadrootp <<- deadroots * pcrlit
    deadcrootp <<- crdecay * crootp * (1.0 - cretrans)
    deadbranchp <<- bdecay * branchp * (1.0 - bretrans)

    # N in stemwood litter - only mobile n is retranslocated
    deadstemn <<- wdecay * (stemnimm + stemnmob *
                    (1.0 - wretrans))

    # P in stemwood litter - only mobile p is retranslocated
    deadstemp <<- wdecay * (stempimm + stempmob *
                   (1.0 - wretrans))

    # Animal grazing?

    # Daily...
    if (grazing == 1) {
        daily_grazing_calc(fdecay)

    # annually
    } else if (grazing == 2 && disturbance_doy == doy) {
        annual_grazing_calc()

    # no grazing
    } else {
        ceaten <<- 0.0
        neaten <<- 0.0
        peaten <<- 0.0
    }

    return (c(fdecay, rdecay))
}

check_for_fire <- function(year, distrubance_yrs, num_disturbance_yrs) {
    # Check if the current year has a fire, if so "burn" and then
    # return an indicator to tell the main code to reset the stress stream
    fire_found <- FALSE

    for (nyr in 1:num_disturbance_yrs) {
        if (year == distrubance_yrs[nyr]) {
            fire_found <- TRUE
        }
    }

    return (fire_found)
}

fire <- function() {
    # Fire...
    #
    # * 100 percent of aboveground biomass
    # * 100 percent of surface litter
    # * 50 percent of N volatilized to the atmosphere
    # * 50 percent of N returned to inorgn pool
    # * 100 percent of N returned to inorgmin pool  Wang 2015 Chemosphere DOI:10.1016/j.chemosphere.2014.01.084
    # * Coarse roots are not damaged by fire!
    #
    # vaguely following ...
    # http:#treephys.oxfordjournals.org/content/24/7/765.full.pdf

    totaln <- branchn + shootn + stemn + structsurfn
    inorgn <<- inorgn + totaln / 2.0

    totalp <- branchp + shootp + stemp + structsurfp
    inorglabp <<- inorglabp + totalp / 2.0

    # re-establish everything with C/N ~ 25 and C/P ~ 2500
    if (alloc_model == GRASSES) {
        branch <<- 0.0
        branchn <<- 0.0
        branchp <<- 0.0
        sapwood <<- 0.0
        stem <<- 0.0
        stemn <<- 0.0
        stemp <<- 0.0
        stemnimm <<- 0.0
        stemnmob <<- 0.0
        stempimm <<- 0.0
        stempmob <<- 0.0
    } else {
        branch <<- 0.001
        branchn <<- 0.00004
        branchp <<- 0.0000004
        sapwood <<- 0.001
        stem <<- 0.001
        stemn <<- 0.00004
        stemnimm <<- 0.00004
        stemnmob <<- 0.0
        stemp <<- 0.0000004
        stempimm <<- 0.0000004
        stempmob <<- 0.0
    }

    age <<- 0.0
    metabsurf <<- 0.0
    metabsurfn <<- 0.0
    metabsurfp <<- 0.0
    prev_sma <<- 1.0
    root <<- 0.001
    rootn <<- 0.00004
    rootp <<- 0.0000004
    shoot <<- 0.01

    lai <<- psla * M2_AS_HA / KG_AS_TONNES / cfracts * shoot
    shootn <<- 0.004
    shootp <<- 0.00004
    structsurf <<- 0.001
    structsurfn <<- 0.00004
    structsurfp <<- 0.0000004

    # reset litter flows
    deadroots <<- 0.0
    deadstems <<- 0.0
    deadbranch <<- 0.0
    deadsapwood <<- 0.0
    deadleafn <<- 0.0
    deadrootn <<- 0.0
    deadbranchn <<- 0.0
    deadstemn <<- 0.0
    deadleafp <<- 0.0
    deadrootp <<- 0.0
    deadbranchp <<- 0.0
    deadstemp <<- 0.0

    # update N:C, P:C of plant pools
    if (float_eq(shoot, 0.0)) {
        shootnc <<- 0.0
        shootpc <<- 0.0
    } else {
        shootnc <<- shootn / shoot
        shootpc <<- shootp / shoot
    }

    if (ncycle == FALSE)
        shootnc <<- prescribed_leaf_NC

    if (pcycle == FALSE)
        shootpc <<- prescribed_leaf_PC

    if (float_eq(root, 0.0)) {
        rootnc <<- 0.0
        rootpc <<- 0.0
    } else {
        rootnc <<- max(0.0, rootn / root)
        rootpc <<- max(0.0, rootp / root)
    }
}

hurricane_f <- function() {
    # Specifically for the florida simulations - reduce LAI by 40%

    # Reduce LAI by 40%
    lai <<- lai - lai * 0.4

    # adjust C in the foliage
    orig_shoot_c <- shoot
    shoot <<- lai / (psla * M2_AS_HA / KG_AS_TONNES / cfracts)
    lost_c <- orig_shoot_c - shoot
    lost_n <- shootnc * lost_c
    lost_p <- shootpc * lost_c
    shootn <<- shootn - lost_n
    shootp <<- shootp - lost_p

    # Drop straight to floor, no retranslocation

    # C structural
    if (float_eq(lost_c, 0.0)) {
        nc_leaf_litter <- 0.0
        pc_leaf_litter <- 0.0
    } else {
        nc_leaf_litter <- lost_n / lost_c
        pc_leaf_litter <- lost_p / lost_c
    }

    if (float_eq(nc_leaf_litter, 0.0)) {
        # catch divide by zero if we have no leaves
        lnleaf <- 0.0
    } else {
        lnleaf <- ligshoot / cfracts / nc_leaf_litter
    }

    fmleaf <- max(0.0, 0.85 - (0.018 * lnleaf))
    surf_struct_litter <<- surf_struct_litter + lost_c * (1.0 - fmleaf)

    # C metabolic
    surf_metab_litter <<- surf_metab_litter + lost_c * fmleaf

    # N structural
    if (float_eq(surf_struct_litter, 0.0)) {
        n_surf_struct_litter <<- n_surf_struct_litter + 0.0
    } else {
        n_surf_struct_litter <<- n_surf_struct_litter + (lost_n * surf_struct_litter *
                                    structrat / surf_struct_litter)
    }

    # P structural
    if (float_eq(surf_struct_litter, 0.0)) {
        p_surf_struct_litter <<- p_surf_struct_litter + 0.0
    } else {
        p_surf_struct_litter <<- p_surf_struct_litter + (lost_p * ssurf_struct_litter *
                                    structratp / surf_struct_litter)
    }

    # N metabolic pools
    n_surf_metab_litter <<- n_surf_metab_litter + lost_n - n_surf_struct_litter

    # P metabolic pools
    p_surf_metab_litter <<- p_surf_metab_litter + lost_p - p_surf_struct_litter

    # structsurf += lost_c
    # structsurfn += lost_n
}

calculate_soil_water_fac <- function() {
    # Estimate a relative water availability factor [0..1]
    #
    # A drying soil results in physiological stress that can induce stomatal
    # closure and reduce transpiration. Further, N mineralisation depends on
    # top soil moisture.
    #
    # qs = 0.2 in SDGVM
    #
    # References:
    # -----------
    # * Landsberg and Waring (1997) Forest Ecology and Management, 95, 209-228.
    #   See  Figure 2.
    # * Egea et al. (2011) Agricultural Forest Meteorology, 151, 1370-1384.
    #
    # But similarly see:
    # * van Genuchten (1981) Soil Sci. Soc. Am. J, 44, 892--898.
    # * Wang and Leuning (1998) Ag Forest Met, 91, 89-111.
    #
    # * Pepper et al. (2008) Functional Change Biology, 35, 493-508
    #
    # Returns:
    # --------
    # wtfac_topsoil : float
    #     water availability factor for the top soil [0,1]
    # wtfac_root : float
    #     water availability factor for the root zone [0,1]

    #psi_swp_topsoil

    #if (water_balance == HYDRAULICS) {
    #    continue
    #    # should put non-stomatal limitation stuff in here.

    if (water_balance == BUCKET && sw_stress_model == 0) {
        # JULES type model, see Egea et al. (2011)
        wtfac_topsoil <<- calc_beta(pawater_topsoil, topsoil_depth,
                                     theta_fc_topsoil, theta_wp_topsoil,
                                     qs)

        wtfac_root <<- calc_beta(pawater_root, rooting_depth,
                                     theta_fc_root, theta_wp_root,
                                     qs)

    } else if (water_balance == BUCKET && sw_stress_model == 1) {
        # Landsberg and Waring, (1997)
        moisture_ratio_topsoil <- pawater_topsoil / wcapac_topsoil
        moisture_ratio_root <- pawater_root / wcapac_root

        wtfac_topsoil <<- calc_sw_modifier(moisture_ratio_topsoil,
                                            ctheta_topsoil,
                                            ntheta_topsoil)
        wtfac_root <<- calc_sw_modifier(moisture_ratio_root, ctheta_root,
                                         ntheta_root)

    } else if (water_balance == BUCKET && sw_stress_model == 2) {
        # Zhou et al.(2013) Agricultural & Forest Met. 182–183, 204–214
        # Assuming that overnight 􏰀pre-dawn leaf water potential =
        # pre-dawn soil water potential.

        #fprintf(stderr, "Zhou model not implemented\n")
        #exit(EXIT_FAILURE)

        # Hardwiring this for testing. Values taken from Table, 1 in
        # De Kauwe et al. 2015, Biogeosciences
        b <- 0.82
        sf <- 1.9
        psi_f <- -1.85

        wtfac_topsoil <<- exp(b * predawn_swp)
        wtfac_root <<- exp(b * predawn_swp)

        #wtfac_topsoil_ns = (1.0 + exp(sf * psi_f)) / \
        #                      (1.0 + exp(sf * (psi_f - predawn_swp)))
        #wtfac_root_ns = (1.0 + exp(sf * psi_f)) / \
        #                      (1.0 + exp(sf * (psi_f - predawn_swp)))

        #
        # wtfac_topsoil = exp(g1_b * psi_s_topsoil)
        # wtfac_root = exp(g1_b * psi_s_root)
        #
        # ! SW modifier for Vcmax (non-stomatal limitation)
        # wtfac_topsoil_ns = (1.0 + exp(vcmax_sf * vcmax_psi_f)) / \
        #                       (1.0 + exp(vcmax_sf * \
        #                                 (vcmax_psi_f - psi_s_topsoil)))
        # wtfac_root_ns = (1.0 + exp(vcmax_sf * vcmax_psi_f)) / \
        #                     (1.0 + exp(vcmax_sf * \
        #                                 (vcmax_psi_f - psi_s_root)))
    }
}

calculate_jmax_and_vcmax_with_p <- function(Tk, N0, P0, mt) {
    # Calculate the maximum RuBP regeneration rate for light-saturated
    # leaves at the top of the canopy (Jmax) and the maximum rate of
    # rubisco-mediated carboxylation at the top of the canopy (Vcmax).
    #
    # Parameters:
    # ----------
    # Tk : float
    # air temperature (Kelvin)
    # N0 : float
    # leaf N   (g N m-2)
    # P0 : float
    # leaf P   (g P m-2)
    #
    # Returns:
    # --------
    # jmax : float (umol/m2/sec)
    # the maximum rate of electron transport at 25 degC
    # vcmax : float (umol/m2/sec)
    # the maximum rate of electron transport at 25 degC

    vcmax__ <- 0.0
    jmax__ <- 0.0

    if (modeljm == 0) {
        jmax__ <- jmax
        vcmax__ <- pvcmax
    } else if (modeljm == 1) {
        if (aci_relationship == WALKER) {
            # Walker et al. 2014 global synthesis relationship
            # the maximum rate of electron transport at 25 degC
            log_vcmax <- 3.946 + 0.921 * log(N0) + 0.121 * log(P0) + 0.282 * log(N0) * log(P0)
            vcmax25 <- exp(log_vcmax)

            # the maximum rate of electron transport at 25 degC
            log_jmax <- 1.246 + 0.886 * log_vcmax + 0.089 * log(P0)
            jmax25 <- exp(log_jmax)
        } else if (aci_relationship == ELLSWORTH) {
            # Ellsworth et al. 2015 PCE EucFACE relationship without TPU limitation
            # need to convert SLA from m2 kg-1 to m2 g-1
            jmax25n <- jmaxna * N0 + jmaxnb
            jmax25p <- jmaxpa * P0 + jmaxpb
            jmax25 <- min(jmax25n, jmax25p)

            # need to convert SLA from m2 kg-1 to m2 g-1
            vcmax25n <- vcmaxna * N0 + vcmaxnb
            vcmax25p <- vcmaxpa * P0 + vcmaxpb
            vcmax25 <- min(vcmax25n, vcmax25p)
        }
        # Temperature-dependent relationship ,
        # this response is well-behaved for TLEAF < 0.0
        jmax__ <- peaked_arrh(mt, jmax25, eaj, Tk, delsj, edj)

        vcmax__ <- arrh(mt, vcmax25, eav, Tk)

    } else if (modeljm == 2) {
        vcmax25 <- vcmaxna * N0 + vcmaxnb
        vcmax__ <- arrh(mt, vcmax25, eav, Tk)

        jmax25 <- jv_slope * vcmax25 - jv_intercept
        jmax__ <- peaked_arrh(mt, jmax25, eaj, Tk, delsj, edj)
    } else if (modeljm == 3) {
        # the maximum rate of electron transport at 25 degC
        jmax25 <- jmax

        # this response is well-behaved for TLEAF < 0.0
        jmax__ <- peaked_arrh(mt, jmax25, eaj, Tk,
                            delsj, edj)

        # the maximum rate of electron transport at 25 degC
        vcmax25 <- pvcmax
        vcmax__ <- arrh(mt, vcmax25, eav, Tk)
    }

    # reduce photosynthetic capacity with moisture stress
    jmax__ <- jmax__ * wtfac_root
    vcmax__ <- vcmax__ * wtfac_root
    #  Function allowing Jmax/Vcmax to be forced linearly to zero at low T
    jmax__ <- adj_for_low_temp(jmax__, Tk)
    vcmax__ <- adj_for_low_temp(vcmax__, Tk)

    return (c(jmax__, vcmax__))
}

calculate_jmax_and_vcmax <- function(Tk, N0, mt) {
    # Calculate the maximum RuBP regeneration rate for light-saturated
    # leaves at the top of the canopy (Jmax) and the maximum rate of
    # rubisco-mediated carboxylation at the top of the canopy (Vcmax).
    #
    # Parameters:
    # ----------
    # Tk : float
    #     air temperature (Kelvin)
    # N0 : float
    #     leaf N   (g N m-2)
    #
    #
    # Returns:
    # --------
    # jmax : float (umol/m2/sec)
    #     the maximum rate of electron transport at 25 degC
    # vcmax : float (umol/m2/sec)
    #     the maximum rate of electron transport at 25 degC

    vcmax__ <- 0.0
    jmax__ <- 0.0

    if (modeljm == 0) {
        jmax__ <- jmax
        vcmax__ <- pvcmax
    } else if (modeljm == 1) {

        # current unit for sla: m2 kg-1 convert into m2 g-1 for Walker relationship
        if (aci_relationship == WALKER) {
            log_vcmax <- 1.993 + 2.555 * log(N0) - 0.372 * log(psla / 1000.0) + 0.422 * log(N0) * log(psla / 1000.0)
            vcmax25 <- exp(log_vcmax)
            log_jmax <- 1.197 + 0.847 * log_vcmax
            jmax25 <- exp(log_jmax)
        } else if (aci_relationship == ELLSWORTH) {
            jmax25 <- jmaxna * N0 + jmaxnb
            vcmax25 <- vcmaxna * N0 + vcmaxnb
        } else if (aci_relationship == BASELINE) {
            jmax25 <- jmaxna * N0 + jmaxnb
            vcmax25 <- vcmaxna * N0 + vcmaxnb
        }


        vcmax__ <- arrh(mt, vcmax25, eav, Tk)

        jmax__ <- peaked_arrh(mt, jmax25, eaj, Tk, delsj, edj)

    } else if (modeljm == 2) {
        vcmax25 <- vcmaxna * N0 + vcmaxnb
        vcmax__ <- arrh(mt, vcmax25, eav, Tk)

        jmax25 <- jv_slope * vcmax25 - jv_intercept
        jmax__ <- peaked_arrh(mt, jmax25, eaj, Tk, delsj, edj)

    } else if (modeljm == 3) {
        # the maximum rate of electron transport at 25 degC
        jmax25 <- jmax

        # this response is well-behaved for TLEAF < 0.0
        jmax__ <- peaked_arrh(mt, jmax25, eaj, Tk, delsj, edj)

        # the maximum rate of electron transport at 25 degC
        vcmax25 <- vcmax
        vcmax__ <- arrh(mt, vcmax25, eav, Tk)

    }

    # reduce photosynthetic capacity with moisture stress
    jmax__ <- jmax__ * wtfac_root
    vcmax__ <- vcmax__ * wtfac_root
    #  Function allowing Jmax/Vcmax to be forced linearly to zero at low T
    jmax__ <- adj_for_low_temp(jmax__, Tk)
    vcmax__ <- adj_for_low_temp(vcmax__, Tk)
    return (c(jmax__, vcmax__))

}

mate_C3_photosynthesis <- function(daylen, ncontent, pcontent) {
    # MATE simulates big-leaf C3 photosynthesis (GPP) based on Sands (1995),
    # accounting for diurnal variations in irradiance and temp (am [sunrise-noon],
    # pm[noon to sunset]).
    #
    # MATE is connected to G'DAY via LAI and leaf N, P content.
    #
    # Plant autotrophic respiration is calculated via carbon-use efficiency
    # (CUE=NPP/GPP).
    #
    # References:
    # -----------
    # * Medlyn, B. E. et al (2011) Global Change Biology, 17, 2134-2144.
    # * McMurtrie, R. E. et al. (2008) Functional Change Biology, 35, 521-34.
    # * Sands, P. J. (1995) Australian Journal of Plant Physiology, 22, 601-14.
    #
    # Rubisco kinetic parameter values are from:
    # * Bernacchi et al. (2001) PCE, 24, 253-259.
    # * Medlyn et al. (2002) PCE, 25, 1167-1179, see pg. 1170.

    mt <- measurement_temp + DEG_TO_KELVIN

    # Calculate mate params & account for temperature dependencies
    N0 <- calculate_top_of_canopy_n(ncontent)   #Unit: g N m-2

    if (pcycle == TRUE) {
        P0 <- calculate_top_of_canopy_p(pcontent)   #Unit: g P m-2
    } else {
        P0 <- 0.0
    }

    gamma_star_am <- calculate_co2_compensation_point(Tk_am, mt)
    gamma_star_pm <- calculate_co2_compensation_point(Tk_pm, mt)

    Km_am <- calculate_michaelis_menten_parameter(Tk_am, mt)
    Km_pm <- calculate_michaelis_menten_parameter(Tk_pm, mt)

    if (pcycle == TRUE) {
        res__ <- calculate_jmax_and_vcmax_with_p(Tk_am, N0, P0, mt)
        jmax_am <- res__[1]
        vcmax_am <- res__[2]
        res__ <- calculate_jmax_and_vcmax_with_p(Tk_pm, N0, P0, mt)
        jmax_pm <- res__[1]
        vcmax_pm <- res__[2]
    } else {
        res__ <- calculate_jmax_and_vcmax(Tk_am, N0, mt)
        jmax_am <- res__[1]
        vcmax_am <- res__[2]
        res__ <- calculate_jmax_and_vcmax(Tk_pm, N0, mt)
        jmax_pm <- res__[1]
        vcmax_pm <- res__[2]
    }

    svcmax <<- (vcmax_am + vcmax_pm) / 2.0

    #fprintf(stderr, "jmax %f vcmax %f\n", jmax_am, vcmax_am)

    ci_am <- calculate_ci(vpd_am, Ca)
    ci_pm <- calculate_ci(vpd_pm, Ca)

    # quantum efficiency calculated for C3 plants
    alpha_am <- calculate_quantum_efficiency(ci_am, gamma_star_am)
    alpha_pm <- calculate_quantum_efficiency(ci_pm, gamma_star_pm)

    # Rubisco carboxylation limited rate of photosynthesis
    ac_am <- assim(ci_am, gamma_star_am, vcmax_am, Km_am)
    ac_pm <- assim(ci_pm, gamma_star_pm, vcmax_pm, Km_pm)

    # Light-limited rate of photosynthesis allowed by RuBP regeneration
    aj_am <- assim(ci_am, gamma_star_am, jmax_am / 4.0, 2.0 * gamma_star_am)
    aj_pm <- assim(ci_pm, gamma_star_pm, jmax_pm / 4.0, 2.0 * gamma_star_pm)

    if (pcycle == TRUE) {
        if (triose_p == TRUE) {
            # Triose-phosphates limited rate of photosynthesis
            ap_am <- assim_p(P0)
            ap_pm <- assim_p(P0)

            # light-saturated photosynthesis rate at the top of the canopy
            asat_am <- min(aj_am, ac_am, ap_am)
            asat_pm <- min(aj_pm, ac_pm, ap_pm)
        } else {
            asat_am <- min(aj_am, ac_am)
            asat_pm <- min(aj_pm, ac_pm)
        }
    } else {
        asat_am <- min(aj_am, ac_am)
        asat_pm <- min(aj_pm, ac_pm)
    }

    # Covert PAR units (umol PAR MJ-1)
    conv <- MJ_TO_J * J_2_UMOL
    par <<- par * conv

    # LUE (umol C umol-1 PAR)  note conversion in epsilon
    lue_am <- epsilon(asat_am, par, alpha_am, daylen)
    lue_pm <- epsilon(asat_pm, par, alpha_pm, daylen)

    # use average to simulate canopy photosynthesis
    lue_avg <- (lue_am + lue_pm) / 2.0

    # absorbed photosynthetically active radiation (umol m-2 s-1)
    if (float_eq(lai, 0.0))
        apar <<- 0.0
    else
        apar <<- par * fipar

    # convert umol m-2 d-1 gC m-2 d-1
    conv <- UMOL_TO_MOL * MOL_C_TO_GRAMS_C
    gpp_gCm2 <<- apar * lue_avg * conv
    gpp_am <<- (apar / 2.0) * lue_am * conv
    gpp_pm <<- (apar / 2.0) * lue_pm * conv
    npp_gCm2 <<- gpp_gCm2 * cue

    #fprintf(stderr, "gpp_gCm2 %f, apar %f, lue_avg %f, asat %f, aj %f, ac %f\n",
    #        gpp_gCm2, apar, lue_avg, asat_am, aj_am, ac_am)


    # g C m-2 to tonnes hectare-1 day-1
    gpp <<- gpp_gCm2 * G_AS_TONNES / M2_AS_HA
    npp <<- npp_gCm2 * G_AS_TONNES / M2_AS_HA

    # save apar in MJ m-2 d-1
    apar <<- apar * UMOL_2_JOL * J_TO_MJ
}

calculate_vcmax_parameter <- function(Tk, N0, mt) {
    # Calculate the maximum rate of rubisco-mediated carboxylation at the
    # top of the canopy
    #
    # References:
    # -----------
    # * Massad, R-S., Tuzet, A. and Bethenod, O. (2007) The effect of
    #   temperature on C4-type leaf photosynthesis parameters. Plant, Cell and
    #   Environment, 30, 1191-1204.
    #
    # http://www.cesm.ucar.edu/models/cesm1.0/clm/CLM4_Tech_Note.pdf
    # Table 8.2 has PFT values...
    #
    # Parameters:
    # ----------
    # Tk : float
    #     air temperature (kelvin)
    # N0 : float
    #     leaf N
    # P0 : float
    #  leaf P
    #
    # Returns:
    # -------
    # vcmax : float, list [am, pm]
    #     maximum rate of Rubisco activity
    # Massad et al. 2007
    # Ea = self.params.eav
    Hd <- self.params.edv
    delS <- self.params.delsv
    Ea <- 67294.0
    Hd <- 144568.0
    delS <- 472.0

    # the maximum rate of electron transport at 25 degC
    vcmax25 <- vcmaxna * N0 + vcmaxnb
    vcmax <- peaked_arrh(mt, vcmax25, Ea, Tk, delS, Hd)

    # reduce photosynthetic capacity with moisture stress
    vcmax <- vcmax * wtfac_root

    # Function allowing Jmax/Vcmax to be forced linearly to zero at low T
    vcmax <- adj_for_low_temp(vcmax, Tk)

    return (c(vcmax, vcmax25))
}

mate_C4_photosynthesis <- function(daylen, ncontent, pcontent) {
    # MATE simulates big-leaf C3 photosynthesis (GPP) based on Sands (1995),
    # accounting for diurnal variations in irradiance and temp (am [sunrise-noon],
    # pm[noon to sunset]).
    #
    # MATE is connected to G'DAY via LAI and leaf N content.
    #
    # Plant autotrophic respiration is calculated via carbon-use efficiency
    # (CUE=NPP/GPP).
    #
    # References:
    # -----------
    # * Collatz, G, J., Ribas-Carbo, M. and Berry, J. A. (1992) Coupled
    #   Photosynthesis-Stomatal Conductance Model for Leaves of C4 plants.
    #   Aust. J. Plant Physiol., 19, 519-38.
    # * von Caemmerer, S. (2000) Biochemical Models of Leaf Photosynthesis. Chp 4.
    #   Modelling C4 photosynthesis. CSIRO PUBLISHING, Australia. pg 91-122.
    # * Sands, P. J. (1995) Australian Journal of Plant Physiology, 22, 601-14.
    #
    # Temperature dependancies:
    # * Massad, R-S., Tuzet, A. and Bethenod, O. (2007) The effect of temperature
    #   on C4-type leaf photosynthesis parameters. Plant, Cell and Environment,
    #   30, 1191-1204.
    #
    # Intrinsic Quantum efficiency (mol mol-1), no Ci or temp dependancey
    # in c4 plants see:
    # * Ehleringer, J. R., 1978, Oecologia, 31, 255-267 or Collatz 1998.
    # * Value taken from Table 1, Collatz et al.1998 Oecologia, 114, 441-454.

    mt <- measurement_temp + DEG_TO_KELVIN

    # curvature parameter, transition between light-limited and
    # carboxylation limited flux. Collatz table 2
    beta1 <- 0.83

    # curvature parameter, co-limitaiton between flux determined by
    # Rubisco and light and CO2 limited flux. Collatz table 2
    beta2 <- 0.93

    # initial slope of photosynthetic CO2 response (mol m-2 s-1),
    # Collatz table 2
    kslope <- 0.7

    # Calculate mate params & account for temperature dependencies
    N0 <- calculate_top_of_canopy_n(ncontent)
    P0 <- calculate_top_of_canopy_p(pcontent)

    ci_am <- calculate_ci(vpd_am, Ca)
    ci_pm <- calculate_ci(vpd_pm, Ca)

    # Temp dependancies from Massad et al. 2007
    res__ <- calculate_vcmax_parameter(Tk_am, N0, mt)
    vcmax_am <- res__[1]
    vcmax25_am <- res__[2]
    res__ <- calculate_vcmax_parameter(Tk_pm, N0, mt)
    vcmax_pm <- res__[1]
    vcmax25_pm <- res__[2]

    # Covert solar irradiance to PAR (umol PAR MJ-1)
    conv <- MJ_TO_J * J_2_UMOL
    par <<- conv
    par_per_sec <- par / (60.0 * 60.0 * daylen)

    # Rubisco and light-limited capacity (Appendix, 2B)
    M_am <- quadratic(beta1, -(vcmax_am + alpha_c4 * par_per_sec),
                    (vcmax_am * alpha_c4 * par_per_sec))
    M_pm <- quadratic(beta1, -(vcmax_pm + alpha_c4 * par_per_sec),
                    (vcmax_pm * alpha_c4 * par_per_sec))

    # The limitation of the overall rate by M and CO2 limited flux
    A_am <- quadratic(beta2, -(M_am + kslope * ci_am),
                    (M_am * kslope * ci_am))
    A_pm <- quadratic(beta2, -(M_pm + kslope * ci_pm),
                    (M_pm * kslope * ci_pm))

    # These respiration terms are just for assimilation calculations,
    # autotrophic respiration is stil assumed to be half of GPP
    Rd_am <- calc_respiration(Tk_am, vcmax25_am)
    Rd_pm <- calc_respiration(Tk_pm, vcmax25_pm)

    # Net (saturated) photosynthetic rate, not sure if this makes sense.
    asat_am <- A_am - Rd_am
    asat_pm <- A_pm - Rd_pm

    # LUE (umol C umol-1 PAR)
    lue_am <- epsilon(asat_am, par, alpha_c4, daylen)
    lue_pm <- epsilon(asat_pm, par, alpha_c4, daylen)

    # use average to simulate canopy photosynthesis
    lue_avg <- (lue_am + lue_pm) / 2.0

    # absorbed photosynthetically active radiation (umol m-2 s-1)
    if (float_eq(lai, 0.0)) {
        apar <<- 0.0
    } else {
        apar <<- par * fipar
    }

    # convert umol m-2 d-1 gC m-2 d-1
    conv <- UMOL_TO_MOL * MOL_C_TO_GRAMS_C
    gpp_gCm2 <<- apar * lue_avg * conv
    gpp_am <<- (apar / 2.0) * lue_am * conv
    gpp_pm <<- (apar / 2.0) * lue_pm * conv
    npp_gCm2 <<- gpp_gCm2 * cue

    # g C m-2 to tonnes hectare-1 day-1
    gpp <<- gpp_gCm2 * G_AS_TONNES / M2_AS_HA
    npp <<- npp_gCm2 * G_AS_TONNES / M2_AS_HA

    # save apar in MJ m-2 d-1
    apar <<- apar * UMOL_2_JOL * J_TO_MJ
}

calc_autotrophic_respiration <- function() {
    # Autotrophic respiration is the sum of the growth component (Rg)
    # and the the temperature-dependent maintenance respiration (Rm) of
    # leaves (Rml), fine roots (Rmr) and wood (Rmw)

    k <- 0.0548

    # respiration rate (gC gN-1 d-1) on a 10degC base
    rk <- resp_coeff * k

    if (ncycle == FALSE) {
        shootn__ <- (shoot * 0.03) * TONNES_HA_2_G_M2
        rootn__ <- (root * 0.02) * TONNES_HA_2_G_M2
        stemn__ <- (stem * 0.003) * TONNES_HA_2_G_M2
    } else {
        shootn__ <- shootn * TONNES_HA_2_G_M2
        rootn__ <- rootn * TONNES_HA_2_G_M2
        stemn__ <- stemn * TONNES_HA_2_G_M2
    }

    if (pcycle == FALSE) {
        shootp__ <- (shoot * 0.003) * TONNES_HA_2_G_M2
        rootp__ <- (root * 0.002) * TONNES_HA_2_G_M2
        stemp__ <- (stem * 0.00003) * TONNES_HA_2_G_M2
    } else {
        shootp__ <- shootp * TONNES_HA_2_G_M2
        rootp__ <- rootp * TONNES_HA_2_G_M2
        stemp__ <- stemp * TONNES_HA_2_G_M2
    }

    # Maintenance respiration, the cost of metabolic processes in
    # living tissues, differs according to tissue N
    #Rml = rk * shootn * lloyd_and_taylor(tair)

    # leaf dark respiration ~ leaf P, N, vcmax, and TWQ
    # where TWQ is the mean temperature of the warmest quarter
    Rml <- 1.2636 + (0.0728 * shootn__) + (0.015 * shootp__) + (0.0095 * svcmax) - (0.0358 * twq)

    Rmw <- 0.2 * rk * stemn__ * lloyd_and_taylor(tair) # should really be sapwoo
    Rmr <- rk * rootn__ * lloyd_and_taylor(tsoil)
    Rm <- (Rml + Rmw + Rmr) * GRAM_C_2_TONNES_HA

    # After maintenance respiration is subtracted from GPP, 25% of the
    # remainder is taken as growth respiration, the cost of producing new
    # tissues
    Rg <- max(0.0, (gpp - Rm) * 0.25)
    auto_resp <<- Rm + Rg

    # Should be revisited, but it occurs to me that during spinup
    # the tissue initialisation could be greater than the incoming gpp and so
    # we might never grow if we respire all out C. Clearly were we using a
    # storage pool this wouldn't be such a drama. For now bound it...
    #
    # De Lucia et al. Global Change Biology (2007) 13, 1157–1167:
    # CUE varied from 0.23 to 0.83 from a literature survey
    cue <- (gpp - auto_resp) / gpp
    if (cue < 0.2) {
        auto_resp <<- gpp * 0.8
    } else if (cue > 0.8) {
        auto_resp <<- gpp * 0.2
    }
}

carbon_daily_production <- function(daylen) {
    # Calculate GPP, NPP and plant respiration at the daily timestep
    #
    # Parameters:
    # -----------
    # daylen : float
    #     daytime length (hrs)
    #
    # References:
    # -----------
    # * Jackson, J. E. and Palmer, J. W. (1981) Annals of Botany, 47, 561-565.

    if (lai > 0.0) {
        # average leaf nitrogen content (g N m-2 leaf)
        leafn <- (shootnc * cfracts / psla * KG_AS_G)
        # leafn <- (0.035 * cfracts / sla * KG_AS_G)
        # average leaf phosphorus content (g P m-2 leaf)
        leafp <- (shootpc * cfracts / psla * KG_AS_G)

        # total nitrogen content of the canopy
        ncontent <- leafn * lai
        # total phosphorus content of the canopy
        pcontent <- leafp * lai

    } else {
        ncontent <- 0.0
        pcontent <- 0.0
    }

    # When canopy is not closed, canopy light interception is reduced
    #     - calculate the fractional ground cover
    if (lai < lai_closed) {
        # discontinuous canopies
        fc <- lai / lai_closed
    } else {
        fc <- 1.0
    }

    # fIPAR - the fraction of intercepted PAR = IPAR/PAR incident at the
    #    top of the canopy, accounting for partial closure based on Jackson
    #    and Palmer (1979).
    if (lai > 0.0) {
        fipar <<- ((1.0 - exp(-kext * lai / fc)) * fc)
    } else {
        fipar <<- 0.0
    }

    if (water_stress) {
        # Calculate the soil moisture availability factors [0,1] in the
        #    topsoil and the entire root zone
        calculate_soil_water_fac()
    } else {
        # really this should only be a debugging option!
        wtfac_topsoil <<- 1.0
        wtfac_root <<- 1.0
    }
    # Estimate photosynthesis
    if (assim_model == BEWDY){
        stop()
    } else if (assim_model == MATE) {
        if (ps_pathway == C3) {
            mate_C3_photosynthesis(daylen, ncontent, pcontent)
        } else {
            mate_C4_photosynthesis(daylen, ncontent, pcontent)
        }
    } else {
        stop("Unknown photosynthesis model'")
    }


    # Calculate plant respiration
    if (respiration_model == FIXED) {
        # Plant respiration assuming carbon-use efficiency.
        auto_resp <<- gpp * (1.0 - cue)
    } else if (respiration_model == VARY) {
        calc_autotrophic_respiration()
    }

    npp <<- max(0.0, gpp - auto_resp)
    npp_gCm2 <<- npp * TONNES_HA_2_G_M2
}

calc_interception <- function() {
    # Estimate canopy interception.
    #
    # 1. At the day scale using a simple model from Landsberg
    # 2. At the sub-daily time scale using the logic from CABLE, but ignoring
    #    canopy drainage e.g. a Rutter type model for the moment.
    #
    # Parameters:
    # -------
    # rain : float
    #     rainfall [mm d-1]
    #
    # References:
    # ----------
    # * Wang (2011)
    # * Landsberg and Sands

    if (sub_daily) {

        # Max canopy intercept (mm): BATS-type canopy saturation
        # proportional to LAI
        canopy_capacity <- 0.1 * lai

        # Calculate canopy intercepted rainfall
        if (rain > 0.0) {
            max_interception <- min(rain, canopy_capacity - canopy_store)
            interception <- min(0.0, max_interception)
            if (tair < 0.0) {
                interception <- 0.0
            }
        } else {
            interception <- 0.0
        }

        # Define canopy throughfall
        throughfall <- rain - interception

        # Add canopy interception to canopy storage term
        canopy_store <<- canopy_store + interception

        # Calculate canopy water storage excess
        if (canopy_store > canopy_capacity) {
            canopy_spill <- canopy_store - canopy_capacity
        } else {
            canopy_spill <- 0.0
        }

        # Move excess canopy water to throughfall
        throughfall <- throughfall + canopy_spill

        # Update canopy storage term
        canopy_store <<- canopy_store - canopy_spill

        # remove canopy evap flux
        if (canopy_store > canopy_evap) {
            canopy_store <<- canopy_store - canopy_evap
        } else {
            # reduce evaporation to water available
            canopy_evap <- canopy_store
            canopy_store <<- 0.0
        }

    } else {

        if (rain > 0.0) {
            # throughfall  = MAX(0.0, rain * rfmult - lai * wetloss)
            # canopy_evap = rain - *throughfall
            # interception = 0.0


            canopy_evap <- (rain * intercep_frac *
                            min(1.0, lai / max_intercep_lai))
            throughfall <- rain - canopy_evap
            interception <- 0.0


        } else {
            canopy_evap <- 0.0
            throughfall <- 0.0
            interception <- 0.0
        }
        canopy_store <<- 0.0
    }

    return (c(throughfall, interception, canopy_evap))
}

update_water_storage <- function(throughfall, interception, canopy_evap, transpiration, soil_evap) {
    # Calculate top soil, root zone plant available water & runoff.
    # NB. et, transpiration & soil evap may all be adjusted in
    # if we don't have sufficient water

    # This is used to account for transpiration losses from the top layer.
    transpiration_topsoil <- wtfac_topsoil * fractup_soil * transpiration

    # Top soil layer
    previous <- pawater_topsoil
    delta_topsoil <- 0.0
    topsoil_loss <- transpiration_topsoil + soil_evap
    pawater_topsoil <<- pawater_topsoil + throughfall - topsoil_loss

    # We have attempted to evap more water than we have
    if (pawater_topsoil < 0.0) {
        # make the layer completely dry
        pawater_topsoil <<- 0.0

        # ** if there was any water in the layer before we over-evaporated
        # ** then use this to do some of the evaporation required

    if (previous > 0.0) {
        soil_evap <- previous
        transpiration_topsoil <- previous - soil_evap
        topsoil_loss <- transpiration_topsoil + soil_evap
    } else {
        soil_evap <- 0.0
        transpiration_topsoil <- 0.0
        topsoil_loss <- 0.0
    }

    delta_topsoil <- previous - pawater_topsoil

    # ** We have more water than the layer can hold, so set the layer to the
    # ** maximum

    } else if (pawater_topsoil > wcapac_topsoil) {
        pawater_topsoil <<- wcapac_topsoil
        delta_topsoil <- previous - pawater_topsoil

        # ** We have enough water to meet demands

    } else {
        delta_topsoil <- previous - pawater_topsoil
    }

    drainage <- throughfall - topsoil_loss + delta_topsoil

    # ** Root zone
    # ** - this is the layer we are actually taking all the water out of.
    # **   it really encompasses the topsoil so as well, so we need to have
    # **   the soil evpaoration here as well, although we aren't adjusting
    # **   that if water isn't available as we've already calculated that
    # **   above based on the top soil layer. Ditto the transpiration taken
    # **   from the top soil layer.


    previous <- pawater_root
    transpiration_root <- transpiration - transpiration_topsoil
    pawater_root <<- pawater_root + drainage - transpiration_root

    # Default is we have no runoff
    runoff <- 0.0

    # We attempted to extract more water than the rootzone holds
    if (pawater_root < 0.0) {

        # make the layer completely dry
        pawater_root <<- 0.0

        # ** if there was any water in the layer before we over-evaporated
        # ** then use this to do some of the evaporation required

    if (previous > 0.0) {
      transpiration_root <- previous
    } else {
      transpiration_root <- 0.0
    }

    # We have more water than the rootzone can hold runoff
    } else if (pawater_root > wcapac_root) {
        runoff <- pawater_root - wcapac_root
        pawater_root <<- wcapac_root
    }

    # Update transpiration & et accounting for the actual available water
    transpiration <- transpiration_topsoil + transpiration_root
    et <- transpiration + soil_evap + canopy_evap

    delta_sw_store <<- pawater_root - previous

    # calculated at the end of the day for sub_daily
    if (!sub_daily) {
        if (water_stress) {
          # Calculate the soil moisture availability factors [0,1] in the
          # topsoil and the entire root zone
          calculate_soil_water_fac()
        } else {
          # really this should only be a debugging option!
          wtfac_topsoil <<- 1.0
          wtfac_root <<- 1.0
        }
    }

    return (c(transpiration, soil_evap, et, runoff))
}

update_daily_water_struct <- function(day_soil_evap, day_transpiration, day_et, day_interception, day_thoughfall, day_canopy_evap, day_runoff) {

    # add half hour fluxes to day total store
    soil_evap <<- day_soil_evap
    transpiration <<- day_transpiration
    et <<- day_et
    interception <<- day_interception
    throughfall <<- day_thoughfall
    canopy_evap <<- day_canopy_evap
    runoff <<- day_runoff
}

calculate_water_balance <- function(daylen, trans_leaf, omega_leaf, rnet_leaf) {
    # Calculate the water balance (including all water fluxes).
    #
    # Parameters:
    # ----------
    # control : structure
    # control structure
    # fluxes : structure
    # fluxes structure
    # met : structure
    # meteorological drivers structure
    # params : structure
    # parameters structure
    # day : int
    # project day. (Dummy argument, only passed for daily model)
    # daylen : double
    # length of day in hours. (Dummy argument, only passed for daily model)
    # trans_leaf : double
    # leaf transpiration (Dummy argument, only passed for sub-daily model)
    # omega_leaf : double
    # decoupling coefficient (Dummy argument, only passed for sub-daily model)
    # rnet_leaf : double
    # total canopy rnet (Dummy argument, only passed for sub-daily model)

    SEC_2_DAY <- 60.0 * 60.0 * daylen
    DAY_2_SEC <- 1.0 / SEC_2_DAY

    # don't need to work out the canopy evap
    res__ <- calc_interception()
    throughfall <- res__[1]
    interception <- res__[2]
    canopy_evap <- res__[3]

    net_rad_am <- calc_net_radiation(sw_rad_am, tair_am)
    net_rad_pm <- calc_net_radiation(sw_rad_pm, tair_pm)
    soil_evap__ <- soil_evap * MOLE_WATER_2_G_WATER * G_TO_KG * (60.0 * 60.0 * daylen)

    # gC m-2 day-1 umol m-2 s-1
    conv <- GRAMS_C_TO_MOL_C * MOL_TO_UMOL * DAY_2_SEC
    gpp_am__ <- gpp_am * conv
    gpp_pm__ <- gpp_pm * conv

    res__ <- penman_canopy_wrapper(press, vpd_am, tair_am, wind_am, net_rad_am, Ca, gpp_am__)
    ga_am <- res__[1]
    gs_am <- res__[2]
    transpiration_am <- res__[3]
    LE_am <- res__[4]
    omega_am <- res__[5]

    res__ <- penman_canopy_wrapper(press, vpd_pm, tair_pm, wind_pm, net_rad_pm, Ca, gpp_pm__)
    ga_pm <- res__[1]
    gs_pm <- res__[2]
    transpiration_pm <- res__[3]
    LE_pm <- res__[4]
    omega_pm <- res__[5]

    # mol m-2 s-1 to mm/day
    conv <- MOLE_WATER_2_G_WATER * G_TO_KG * SEC_2_DAY
    transpiration <- (transpiration_am + transpiration_pm) * conv

    omega <<- (omega_am + omega_pm) / 2.0

    # output in mol H20 m-2 s-1
    gs_mol_m2_sec <<- gs_am + gs_pm
    ga_mol_m2_sec <<- ga_am + ga_pm


    # ** NB. et, transpiration & soil evap may all be adjusted in
    # ** update_water_storage if we don't have sufficient water

    et <- transpiration + soil_evap__ + canopy_evap

    res__ <- update_water_storage(throughfall, interception, canopy_evap, transpiration, soil_evap__)
    transpiration <- res__[1]
    soil_evap__ <- res__[2]
    et <- res__[3]
    runoff <- res__[4]

    update_daily_water_struct(soil_evap__, transpiration, et, interception, throughfall, canopy_evap, runoff)
}

carbon_allocation <- function(npitfac, doy) {
    # C distribution - allocate available C through system
    #
    # Parameters:
    # -----------
    # npitfac : float
    #     leaf N:C as a fraction of 'Ncmaxfyoung' (max 1.0)

    if (deciduous_model) {
        days_left <- growing_days[doy]
        cpleaf <<- lrate * days_left
        cpbranch <<- brate * days_left
        cpstem <<- wrate * days_left
        cproot <<- c_to_alloc_root * 1.0 / num_days
        cpcroot <<- crate * days_left
    } else {
        cpleaf <<- npp * alleaf
        cproot <<- npp * alroot
        cpcroot <<- npp * alcroot
        cpbranch <<- npp * albranch
        cpstem <<- npp * alstem
    }

    # evaluate SLA of new foliage accounting for variation in SLA
    # with tree and leaf age (Sands and Landsberg, 2002). Assume
    # SLA of new foliage is linearly related to leaf N:C ratio
    # via nitfac. Based on date from two E.globulus stands in SW Aus, see
    # Corbeels et al (2005) Ecological Modelling, 187, 449-474.
    # (m2 onesided/kg DW)
    # This needs to be updated to consider P effect

    psla <<- slazero + npitfac * (slamax - slazero)

    if (deciduous_model) {
        if (shoot == 0.0) {
            lai <<- 0.0
        } else if (leaf_out_days[doy] > 0.0) {
            lai <<- lai + (cpleaf *
                      (psla * M2_AS_HA / (KG_AS_TONNES * cfracts)) -
                      (deadleaves + ceaten) * lai / shoot)
        } else {
            lai <<- 0.0
        }
    } else {
        # update leaf area [m2 m-2]
        if (shoot ==0.0) {
            lai <<- 0.0
        } else {
            lai <<- lai + (cpleaf *
                      (psla * M2_AS_HA / (KG_AS_TONNES * cfracts)) -
                      (deadleaves + ceaten) * lai / shoot)
        }
    }

    if (fixed_lai) {
        lai <<- fix_lai
    }

    #fprintf(stderr, "lai %f\n", lai)
}

calculate_cnp_wood_ratios <- function(npitfac, nitfac, pitfac) {
    # Estimate the N:C and P:C ratio in the branch and stem. Option to vary
    # the N:C and P:C ratio of the stem following Jeffreys (1999) or keep it a fixed
    # fraction
    #
    # Parameters:
    # -----------
    # npitfac: float
    #    min of nitfac and pitfac
    # nitfac: float
    #    leaf N:C as a fraction of Ncmaxyoung
    # pitfac : float
    #    leaf P:C as a fraction of Pcmaxyoung
    #
    # Returns:
    # --------
    # ncbnew : float
    #     N:C ratio of branch
    # nccnew : double
    #     N:C ratio of coarse root
    # ncwimm : float
    #     N:C ratio of immobile stem
    # ncwnew : float
    #     N:C ratio of mobile stem
    # pcbnew : float
    #     P:C ratio of branch
    # pccnew : double
    #     P:C ratio of coarse root
    # pcwimm : float
    #     P:C ratio of immobile stem
    # pcwnew : float
    #     P:C ratio of mobile stem
    #
    # References:
    # ----------
    # * Jeffreys, M. P. (1999) Dynamics of stemwood nitrogen in Pinus radiata
    #   with modelled implications for forest productivity under elevated
    #   atmospheric carbon dioxide. PhD.

    # calculate N:C ratios
    if (nitfac < npitfac) {
        # n:c ratio of new branch wood
        ncbnew <- ncbnew + nitfac * (ncbnew - ncbnewz)

        # n:c ratio of coarse root
        nccnew <- nccnew + nitfac * (nccnew - nccnewz)

        # fixed N:C in the stemwood
        if (fixed_stem_nc) {
            # n:c ratio of stemwood - immobile pool and new ring
            ncwimm <- ncwimm + nitfac * (ncwimm - ncwimmz)

            # New stem ring N:C at critical leaf N:C (mobile)
            ncwnew <- ncwnew + nitfac * (ncwnew - ncwnewz)

        # vary stem N:C based on reln with foliage, see Jeffreys PhD thesis.
        # Jeffreys 1999 showed that N:C ratio of new wood increases with
        # foliar N:C ratio,modelled here based on evidence as a linear
        # function.
        } else {
            ncwimm <- max(0.0, (0.0282 * shootnc + 0.000234) * fhw)

            # New stem ring N:C at critical leaf N:C (mobile)
            ncwnew <- max(0.0, 0.162 * shootnc - 0.00143)
        }
    } else {
        # n:c ratio of new branch wood
        ncbnew <- ncbnew + npitfac * (ncbnew - ncbnewz)

        # n:c ratio of coarse root
        nccnew <- nccnew + npitfac * (nccnew - nccnewz)

        # fixed N:C in the stemwood
        if (fixed_stem_nc) {
            # n:c ratio of stemwood - immobile pool and new ring
            ncwimm <- ncwimm + npitfac * (ncwimm - ncwimmz)

            # New stem ring N:C at critical leaf N:C (mobile)
            ncwnew <- ncwnew + npitfac * (ncwnew - ncwnewz)

        # vary stem N:C based on reln with foliage, see Jeffreys PhD thesis.
        # Jeffreys 1999 showed that N:C ratio of new wood increases with
        # foliar N:C ratio,modelled here based on evidence as a linear
        # function.
        } else {
            ncwimm <- max(0.0, (0.0282 * shootnc + 0.000234) * fhw)

            # New stem ring N:C at critical leaf N:C (mobile)
            ncwnew <- max(0.0, 0.162 * shootnc - 0.00143)
        }
    }

    # calculate P:C ratios
    if (pitfac < npitfac) {
        # p:c ratio of new branch wood
        pcbnew <- pcbnew + pitfac * (pcbnew - pcbnewz)

        # p:c ratio of coarse root
        pccnew <- pccnew + pitfac * (pccnew - pccnewz)

        # fixed P:C in the stemwood
        if (fixed_stem_pc) {
            # p:c ratio of stemwood - immobile pool and new ring
            pcwimm <- pcwimm + pitfac * (pcwimm - pcwimmz)

            # New stem ring P:C at critical leaf P:C (mobile)
            pcwnew <- pcwnew + pitfac * (pcwnew - pcwnewz)

            # vary stem P:C based on reln with foliage,
            # equation based on data from Attiwill 1978 - 1980 paper series.
        } else {
            pcwimm <- max(0.0, -0.0016 * shootpc + 0.000003)

            # New stem ring P:C at critical leaf P:C (mobile),
            # equation based on data from Attiwill 1978 - 1980 paper series
            pcwnew <- max(0.0, -0.0022 * shootpc + 0.000009)
        }
    } else {
        # p:c ratio of new branch wood
        pcbnew <- pcbnew + npitfac * (pcbnew - pcbnewz)

        # p:c ratio of coarse root
        pccnew <- pccnew + npitfac * (pccnew - pccnewz)

        # fixed P:C in the stemwood
        if (fixed_stem_pc) {
            # p:c ratio of stemwood - immobile pool and new ring
            pcwimm <- pcwimm + npitfac * (pcwimm - pcwimmz)

            # New stem ring P:C at critical leaf P:C (mobile)
            pcwnew <- pcwnew + npitfac * (pcwnew - pcwnewz)

        # vary stem P:C based on reln with foliage,
        # equation based on data from Attiwill 1978 - 1980 paper series.
        } else {
            pcwimm <- max(0.0, -0.0016 * shootpc + 0.000003)

            # New stem ring P:C at critical leaf P:C (mobile),
            # equation based on data from Attiwill 1978 - 1980 paper series
            pcwnew <- max(0.0, -0.0022 * shootpc + 0.000009)
        }
    }

    return (c(ncbnew, nccnew, ncwimm, ncwnew, pcbnew, pccnew, pcwimm, pcwnew))
}

cut_back_production <- function(tot, xcbnew, xccnew, xcwimm, xcwnew, doy) {

    # default is we don't need to recalculate the water balance,
    # however if we cut back on NPP due to available N and P below then we do
    # need to do this
    recalc_wb <- FALSE

    # Need to readjust the LAI for the reduced growth as this will
    # have already been increased. First we need to figure out how
    # much we have increased LAI by, important it is done here
    # before cpleaf is reduced!
    if (float_eq(shoot, 0.0)) {
        lai_inc <- 0.0
    } else {
        lai_inc <- (cpleaf *
                   (psla * M2_AS_HA / (KG_AS_TONNES * cfracts)) -
                   (deadleaves + ceaten) * lai / shoot)
    }

    npp <- npp * tot / (npstemimm + npstemmob + npbranch + npcroot)

    # need to adjust growth values accordingly as well
    cpleaf <<- npp * alleaf
    cproot <<- npp * alroot
    cpcroot <<- npp * alcroot
    cpbranch <<- npp * albranch
    cpstem <<- npp * alstem



    if (pcycle) {
        ppbranch <<- npp * albranch * xcbnew
        ppstemimm <<- npp * alstem * xcwimm
        ppstemmob <<- npp * alstem * (xcwnew - xcwimm)
        ppcroot <<- npp * alcroot * xccnew
    } else {
        npbranch <<- npp * albranch * xcbnew
        npstemimm <<- npp * alstem * xcwimm
        npstemmob <<- npp * alstem * (xcwnew - xcwimm)
        npcroot <<- npp * alcroot * xccnew
    }

    # Save WUE before cut back
    wue <<- gpp_gCm2 / transpiration

    # Also need to recalculate GPP and thus Ra and return a flag
    # so that we know to recalculate the water balance.
    gpp <<- npp / cue
    conv <- G_AS_TONNES / M2_AS_HA
    gpp_gCm2 <<- gpp / conv
    gpp_am <<- gpp_gCm2 / 2.0
    gpp_pm <<- gpp_gCm2 / 2.0

    # New respiration flux
    auto_resp <<- gpp - npp
    recalc_wb <- TRUE

    # Now reduce LAI for down-regulated growth.
    if (deciduous_model) {
        if (float_eq(shoot, 0.0)) {
            lai <<- 0.0
        } else if (leaf_out_days[doy] > 0.0) {
            lai <<- lai - lai_inc
            lai <<- lai + (cpleaf *
                       (psla * M2_AS_HA /
                       (KG_AS_TONNES * cfracts)) -
                       (deadleaves + ceaten) * lai / shoot)
        } else {
            lai <<- 0.0
        }
    } else {
        # update leaf area [m2 m-2]
        if (float_eq(shoot, 0.0)) {
            lai <<- 0.0
        } else {
            lai <<- lai - lai_inc
            lai <<- lai + (cpleaf *
                       (psla * M2_AS_HA /
                       (KG_AS_TONNES * cfracts)) -
                       (deadleaves + ceaten) * lai / shoot)
        }
    }

    return (recalc_wb)
}

np_allocation <- function(ncbnew, nccnew, ncwimm, ncwnew, pcbnew,
                  pccnew, pcwimm, pcwnew, fdecay, rdecay, doy) {
    # Nitrogen and phosphorus distribution - allocate available N and
    # P (mineral) through system. N and P is first allocated to the woody
    # component, surplus N and P is then allocated to the shoot and roots
    # with flexible ratios.
    #
    # References:
    # -----------
    # McMurtrie, R. E. et al (2000) Plant and Soil, 224, 135-152.
    #
    # Parameters:
    # -----------
    # ncbnew : float
    #     N:C ratio of branch
    # ncwimm : float
    #     N:C ratio of immobile stem
    # ncwnew : float
    #     N:C ratio of mobile stem
    # pcbnew : float
    #     P:C ratio of branch
    # pcwimm : float
    #     P:C ratio of immobile stem
    # pcwnew : float
    #     P:C ratio of mobile stem
    # fdecay : float
    #     foliage decay rate
    # rdecay : float
    #     fine root decay rate

    depth_guess <- 1.0

    # default is we don't need to recalculate the water balance,
    # however if we cut back on NPP due to available N and P below then we do
    # need to do this
    recalc_wb <- FALSE

    # N and P retranslocated proportion from dying plant tissue and stored within
    # the plant
    retrans <<- nitrogen_retrans(fdecay, rdecay, doy)
    retransp <<- phosphorus_retrans(fdecay, rdecay, doy)
    nuptake <<- calculate_nuptake()
    puptake <<- calculate_puptake()

    #  Ross's Root Model.
    if (model_optroot) {

        # convert t ha-1 day-1 to gN m-2 year-1
        nsupply <- (calculate_nuptake() * TONNES_HA_2_G_M2 * DAYS_IN_YRS)

        # covnert t ha-1 to kg DM m-2
        rtot <- root * TONNES_HA_2_KG_M2 / cfracts
        # nuptake_old = nuptake

        res__ <- calc_opt_root_depth(d0x, r0, topsoil_depth * MM_TO_M, rtot, nsupply, depth_guess)
        root_depth <<- res__[1]
        nuptake <<- res__[2]
        rabove <<- res__[3]

        # covert nuptake from gN m-2 year-1  to t ha-1 day-1
        nuptake <<- nuptake * G_M2_2_TONNES_HA * YRS_IN_DAYS

        # covert from kg DM N m-2 to t ha-1
        deadroots <<- prdecay * rabove * cfracts * KG_M2_2_TONNES_HA
        deadrootn <<- rootnc * (1.0 - rretrans) * deadroots

    }

    # Mineralised nitrogen lost from the system by volatilisation/leaching
    nloss <<- rateloss * inorgn

    # Mineralised P lost from the system by leaching
    if (inorgsorbp > 0.0) {
        ploss <<- prateloss * inorglabp
    } else {
        ploss <<- 0.0
    }

    # total nitrogen/phosphorus to allocate
    ntot <- max(0.0, nuptake + retrans)
    ptot <- max(0.0, puptake + retransp)

    if (deciduous_model) {
        # allocate N to pools with fixed N:C ratios

        # N flux into new ring (immobile component structrual components)
        npstemimm <<- wnimrate * growing_days[doy]

        # N flux into new ring (mobile component can be retrans for new
        # woody tissue)
        npstemmob <<- wnmobrate * growing_days[doy]
        nproot <<- n_to_alloc_root / num_days
        npcroot <<- cnrate * growing_days[doy]
        npleaf <<- lnrate * growing_days[doy]
        npbranch <<- bnrate * growing_days[doy]

        # allocate P to pools with fixed P:C ratios
        ppstemimm <<- wpimrate * growing_days[doy]
        ppstemmob <<- wpmobrate * growing_days[doy]
        pproot <<- p_to_alloc_root / num_days
        ppcroot <<- cprate * growing_days[doy]
        ppleaf <<- lprate * growing_days[doy]
        ppbranch <<- bprate * growing_days[doy]

    } else {
        # allocate N to pools with fixed N:C ratios

        # N flux into new ring (immobile component  structural components)
        npstemimm <<- npp * alstem * ncwimm

        # N flux into new ring (mobile component  can be retrans for new
        # woody tissue)
        npstemmob <<- npp * alstem * (ncwnew - ncwimm)
        npbranch <<- npp * albranch * ncbnew
        npcroot <<- npp * alcroot * nccnew

        # allocate P to pools with fixed P:C ratios
        ppstemimm <<- npp * alstem * pcwimm
        ppstemmob <<- npp * alstem * (pcwnew - pcwimm)
        ppbranch <<- npp * albranch * pcbnew
        ppcroot <<- npp * alcroot * pccnew

        # If we have allocated more N than we have avail, cut back C prodn
        arg <- npstemimm + npstemmob + npbranch + npcroot
        if (arg > ntot && fixleafnc == FALSE && fixed_lai && ncycle) {
            recalc_wb <- cut_back_production(ntot, ncbnew, nccnew, ncwimm, ncwnew, doy)
            print("in cut back N \n")
        }

        # If we have allocated more P than we have avail, cut back C prodn
        arg <- ppstemimm + ppstemmob + ppbranch + ppcroot
        if (arg > ptot && fixleafpc == FALSE && fixed_lai && pcycle) {
            recalc_wb <- cut_back_production(ptot, pcbnew, pccnew, pcwimm, pcwnew, doy)
        }

        # Nitrogen reallocation to flexible-ratio pools
        ntot <- ntot - npbranch + npstemimm + npstemmob + npcroot
        ntot <- max(0.0, ntot)

        # allocate remaining N to flexible-ratio pools
        npleaf <<- ntot * alleaf / (alleaf + alroot * ncrfac)
        nproot <<- ntot - npleaf

        # Phosphorus reallocation to flexible-ratio pools
        ptot <- ptot - ppbranch + ppstemimm + ppstemmob + ppcroot
        ptot <- max(0.0, ptot)

        # allocate remaining P to flexible-ratio pools
        ppleaf <<- ptot * alleaf / (alleaf + alroot * pcrfac)
        pproot <<- ptot - ppleaf
    }

    return (recalc_wb)
}

calc_root_exudation <- function() {
    # Rhizodeposition (root_exc) is assumed to be a fraction of the
    # current root growth rate (cproot), which increases with increasing
    # N stress of the plant.

    if (float_eq(shoot, 0.0) || float_eq(shootn, 0.0)) {
        # nothing happens during leaf off period
        CN_leaf <- 0.0
        frac_to_rexc <- 0.0
    } else {

        if (deciduous_model) {
            # broadleaf
            CN_ref <- 25.0
        } else {
            # conifer
            CN_ref <- 42.0
        }

        # The fraction of growth allocated to rhizodeposition, constrained
        # to solutions lower than 0.5

        CN_leaf <- 1.0 / shootnc
        arg <- max(0.0, (CN_leaf - CN_ref) / CN_ref)
        frac_to_rexc <- min(0.5, a0rhizo + a1rhizo * arg)
    }

    # Rhizodeposition
    root_exc <<- frac_to_rexc * cproot
    if (float_eq(cproot, 0.0)) {
        root_exn <<- 0.0
    } else {
        # N flux associated with rhizodeposition is based on the assumption
        # that the CN ratio of rhizodeposition is equal to that of fine root growth

        root_exn <<- root_exc * (nproot / cproot)
    }

    # Need to remove exudation C & N fluxes from fine root growth fluxes so
    # that things balance.

    cproot <<- cproot - root_exc
    nproot <<- nproot - root_exn
}

calculate_cnp_store <- function() {
    # Calculate labile C, N and P stores from which growth is allocated in the
    # following year.

    cstore <<- cstore + npp
    nstore <<- nstore + nuptake + retrans
    pstore <<- pstore + puptake + retransp
    anpp <<- anpp + npp
}

update_plant_state <- function(fdecay, rdecay, doy) {
    # Daily change in C content
    #
    # Parameters:
    # -----------
    # fdecay : float
    #     foliage decay rate
    # rdecay : float
    #     fine root decay rate

    # Carbon pools

    shoot <<- shoot + cpleaf - deadleaves - ceaten
    root <<- root + cproot - deadroots
    croot <<- croot + cpcroot - deadcroots
    branch <<- branch + cpbranch - deadbranch
    stem <<- stem + cpstem - deadstems

    # fprintf(stderr, "cproot %f\n", cproot)
    # fprintf(stderr, "deadroots %f\n", deadroots)

    # annoying but can't see an easier way with the code as it is.
    # If we are modelling grases, i.e. no stem them without this
    # the sapwood will end up being reduced to a silly number as
    # deadsapwood will keep being removed from the pool, even though there
    # is no wood.
    if (float_eq(stem, 0.01)) {
        sapwood <<- 0.01
    } else if (stem < 0.01) {
        sapwood <<- 0.01
    } else {
        sapwood <<- sapwood + cpstem - deadsapwood
    }

    # Nitrogen and Phosphorus pools

    if (deciduous_model) {
        shootn <<- shootn + (npleaf - (lnrate * remaining_days[doy]) - neaten)
        shootp <<- shootp + (ppleaf - (lprate * remaining_days[doy]) - peaten)
    } else {
        shootn <<- shootn + npleaf - fdecay * shootn - neaten
        shootp <<- shootp + ppleaf - fdecay * shootp - peaten
    }

    branchn <<- branchn + npbranch - bdecay * branchn
    rootn <<- rootn + nproot - rdecay * rootn
    crootn <<- crootn + npcroot - crdecay * crootn
    stemnimm <<- stemnimm + npstemimm - wdecay * stemnimm
    stemnmob <<- stemnmob + (npstemmob - wdecay * stemnmob - retransmob * stemnmob)
    stemn <<- stemnimm + stemnmob

    branchp <<- branchp + ppbranch - bdecay * branchp

    rootp <<- rootp + pproot - rdecay * rootp

    # fprintf(stderr, "nuptake %f\n", nuptake*100000)
    # fprintf(stderr, "puptake %f\n", puptake*100000)
    # fprintf(stderr, "nproot %f\n", nproot)
    # fprintf(stderr, "pproot %f\n", pproot)
    # fprintf(stderr, "rootc %f\n", root)
    # fprintf(stderr, "rootn %f\n", rootn)
    # fprintf(stderr, "rootp %f\n", rootp)
    # fprintf(stderr, "ncrfac calc %f\n", (rootn/root)/(shootn/shoot))
    # fprintf(stderr, "pcrfac calc %f\n", (rootp/root)/(shootp/shoot))

    crootp <<- crootp + ppcroot - crdecay * crootp
    stempimm <<- stempimm + ppstemimm - wdecay * stempimm

    stempmob <<- stempmob + (ppstemmob - wdecay * stempmob - retransmob * stempmob)

    stemp <<- stempimm + stempmob


    if (deciduous_model == FALSE) {
        # =============================
        #  Enforce maximum N:C and P:C ratios.
        # =============================

        # If foliage or root N/C exceeds its max, then N uptake is cut back
        # Similarly, of foliage or root P/C exceeds max, then P uptake is cut back

        # maximum leaf n:c and p:c ratios is function of stand age
        #  - switch off age effect by setting ncmaxfyoung = ncmaxfold
        #  - switch off age effect by setting pcmaxfyoung = pcmaxfold
        age_effect <- (age - ageyoung) / (ageold - ageyoung)
        ncmaxf <- ncmaxfyoung - (ncmaxfyoung - ncmaxfold) * age_effect
        pcmaxf <- pcmaxfyoung - (pcmaxfyoung - pcmaxfold) * age_effect

        if (ncmaxf < ncmaxfold)
            ncmaxf <- ncmaxfold

        if (ncmaxf > ncmaxfyoung)
            ncmaxf <- ncmaxfyoung

        if (pcmaxf < pcmaxfold)
            pcmaxf <- pcmaxfold

        if (pcmaxf > pcmaxfyoung)
            pcmaxf <- pcmaxfyoung

        extrasn <- 0.0
        if (lai > 0.0) {

            if (shootn > (shoot * ncmaxf)) {
                extrasn <- shootn - shoot * ncmaxf

                # Ensure N uptake cannot be reduced below zero.
                if (extrasn >  nuptake)
                    extrasn <- nuptake

                shootn <<- shootn - extrasn
                # nuptake -= extrasn
            }
        }

        extrasp <- 0.0
        if (lai > 0.0) {

            if (shootp > (shoot * pcmaxf)) {
                extrasp <- shootp - shoot * pcmaxf

                # Ensure P uptake cannot be reduced below zero.
                if (extrasp >  puptake)
                  extrasp <- puptake

                shootp <<- shootp - extrasp
                # puptake -= extrasp
            }
        }

        # if root N:C ratio exceeds its max, then nitrogen uptake is cut
        # back. n.b. new ring n/c max is already set because it is related
        # to leaf n:c

        # max root n:c
        ncmaxr <- ncmaxf * ncrfac
        extrarn <- 0.0
        if (rootn > (root * ncmaxr)) {
            extrarn <- rootn - root * ncmaxr

            # Ensure N uptake cannot be reduced below zero.
            if ((extrasn + extrarn) > nuptake)
                extrarn <- nuptake - extrasn

            rootn <<- rootn - extrarn
            nuptake <<- nuptake - (extrarn + extrasn)
        }

        # max root p:c
        pcmaxr <- pcmaxf * pcrfac
        extrarp <- 0.0
        if (rootp > (root * pcmaxr)) {
            extrarp <- rootp - root * pcmaxr

            # Ensure P uptake cannot be reduced below zero.
            if ((extrasp + extrarp) > puptake)
                extrarp <- puptake - extrasp

            rootp <<- rootp - extrarp
            puptake <<- puptake - (extrarp + extrasp)
        }
    }

    # Update deciduous storage pools
    if (deciduous_model)
        calculate_cnp_store()
}

precision_control <- function() {
    # Detect very low values in state variables and force to zero to
    # avoid rounding and overflow errors

    tolerance <- 1E-10

    # C, N & P state variables
    if (shoot < tolerance) {
        deadleaves <<- deadleaves + shoot
        deadleafn <<- deadleafn + shootn
        deadleafp <<- deadleafp + shootp
        shoot <<- 0.0
        shootn <<- 0.0
        shootp <<- 0.0
    }

    if (branch < tolerance) {
        deadbranch <<- deadbranch + branch
        deadbranchn <<- deadbranchn + branchn
        deadbranchp <<- deadbranchp + branchp
        branch <<- 0.0
        branchn <<- 0.0
        branchp <<- 0.0
    }

    if (root < tolerance) {
        deadrootn <<- deadrootn + rootn
        deadrootp <<- deadrootp + rootp
        deadroots <<- deadroots + root
        root <<- 0.0
        rootn <<- 0.0
        rootp <<- 0.0
    }

    if (croot < tolerance) {
        deadcrootn <<- deadcrootn + crootn
        deadcrootp <<- deadcrootp + crootp
        deadcroots <<- deadcroots + croot
        croot <<- 0.0
        crootn <<- 0.0
        crootp <<- 0.0
    }

    # Not setting these to zero as this just leads to errors with desert
    # regrowth...instead seeding them to a small value with a CN~25 and CP~300.

    if (stem < tolerance) {
        deadstems <<- deadstems + stem
        deadstemn <<- deadstemn + stemn
        deadstemp <<- deadstemp + stemp
        stem <<- 0.001
        stemn <<- 0.00004
        stemp <<- 0.000003
        stemnimm <<- 0.00004
        stemnmob <<- 0.0
        stempimm <<- 0.000003
        stempmob <<- 0.0
    }

    # need separate one as this will become very small if there is no
    # mobile stem N/P
    if (stemnmob < tolerance) {
        deadstemn <<- deadstemn + stemnmob
        stemnmob <<- 0.0
    }

    if (stemnimm < tolerance) {
        deadstemn <<- deadstemn + stemnimm
        stemnimm <<- 0.00004
    }

    if (stempmob < tolerance) {
        deadstemp <<- deadstemp + stempmob
        stempmob <<- 0.0
    }

    if (stempimm < tolerance) {
        deadstemp <<- deadstemp + stempimm
        stempimm <<- 0.000003
    }
}

update_water_storage_recalwb <- function() {
    # Calculate root and top soil plant available water and runoff.
    # Soil drainage is estimated using a "leaky-bucket" approach with two
    # soil layers. In reality this is a combined drainage and runoff
    # calculation, i.e. "outflow". There is no drainage out of the "bucket"
    # soil.
    # Returns:
    # --------
    # outflow : float
    # outflow [mm d-1]

    # This is used to account for transpiration losses from the top layer.
    transpiration_topsoil <- (wtfac_topsoil * fractup_soil * transpiration)

    # Top soil layer
    previous <- pawater_topsoil
    delta_topsoil <- 0.0
    topsoil_loss <- transpiration_topsoil + soil_evap
    pawater_topsoil <- pawater_topsoil + throughfall - topsoil_loss

    # We have attempted to evap more water than we have
    if (pawater_topsoil < 0.0) {

        # make the layer completely dry
        pawater_topsoil <<- 0.0

        # if there was any water in the layer before we over-evaporated
        # then use this to do some of the evaporation required

        if (previous > 0.0) {
            # soil_evap = previous / 2.0
            # transpiration_topsoil = previous / 2.0
            soil_evap <<- previous
            transpiration_topsoil <- previous - soil_evap
            topsoil_loss <- transpiration_topsoil + soil_evap
        } else {
            soil_evap <<- 0.0
            transpiration_topsoil <- 0.0
            topsoil_loss <- 0.0
        }
        delta_topsoil <- previous - pawater_topsoil

        # We have more water than the layer can hold, so set the layer to the
        # maximum
    } else if (pawater_topsoil > wcapac_topsoil) {
        pawater_topsoil <<- wcapac_topsoil
        delta_topsoil <- previous - pawater_topsoil

        # We have enough water to meet demands
    } else {
        delta_topsoil <- previous - pawater_topsoil
    }

    drainage <- throughfall - topsoil_loss + delta_topsoil

    # Root zone
    # - this is the layer we are actually taking all the water out of.
    #   it really encompasses the topsoil so as well, so we need to have
    #   the soil evpaoration here as well, although we aren't adjusting
    #   that if water isn't available as we've already calculated that
    #   above based on the top soil layer. Ditto the transpiration taken
    #   from the top soil layer.

    previous <- pawater_root
    transpiration_root <- transpiration - transpiration_topsoil
    pawater_root <- pawater_root + drainage - transpiration_root

    # Default is we have no runoff
    runoff <<- 0.0

    # We attempted to extract more water than the rootzone holds
    if (pawater_root < 0.0) {

        # make the layer completely dry
        pawater_root <<- 0.0

        # if there was any water in the layer before we over-evaporated
        # then use this to do some of the evaporation required
        if (previous > 0.0) {
            transpiration_root <- previous
        } else {
            transpiration_root <- 0.0
        }

    # We have more water than the rootzone can hold runoff
    } else if (pawater_root > wcapac_root) {
        runoff <<- pawater_root - wcapac_root
        pawater_root <<- wcapac_root
    }

    # Update transpiration & et accounting for the actual available water
    transpiration <<- transpiration_topsoil + transpiration_root
    et <<- transpiration + soil_evap + canopy_evap

    delta_sw_store <<- pawater_root - previous

    # calculated at the end of the day for sub_daily
    if (water_stress) {
        # Calculate the soil moisture availability factors [0,1] in the
        # topsoil and the entire root zone
        calculate_soil_water_fac()
    } else {
        # really this should only be a debugging option!
        wtfac_topsoil <<- 1.0
        wtfac_root <<- 1.0
    }
}

zero_carbon_day_fluxes <- function() {
    gpp_gCm2 <<- 0.0
    npp_gCm2 <<- 0.0
    gpp <<- 0.0
    npp <<- 0.0
    auto_resp <<- 0.0
    apar <<- 0.0
}

zero_water_day_fluxes <- function() {
    et <<- 0.0
    soil_evap <<- 0.0
    transpiration <<- 0.0
    interception <<- 0.0
    canopy_evap <<- 0.0
    throughfall <<- 0.0
    runoff <<- 0.0
    gs_mol_m2_sec <<- 0.0
    day_ppt <<- 0.0
}

unpack_solar_geometry <- function() {

    # This geometry calculations are suprisingly intensive which is a waste
    # during spinup, so we are now doing this once and then we are just
    # accessing the 30-min value from the array position

    # calculate_solar_geometry(cw, p, doy, hod)
    # get_diffuse_frac(cw, doy, sw_rad)
    cos_zenith <<- cz_store[hour_idx]
    elevation <<- ele_store[hour_idx]
    diffuse_frac <<- df_store[hour_idx]
}

calculate_absorbed_radiation <- function(par) {
    # Calculate absorded irradiance of sunlit and shaded fractions of
    # the canopy. The total irradiance absorbed by the canopy and the
    # sunlit/shaded components are all expressed on a ground-area basis!
    #
    # NB:  sin_beta == cos_zenith
    #
    # References:
    # -----------
    # * De Pury & Farquhar (1997) PCE, 20, 537-557.
    #
    # but see also:
    # * Wang and Leuning (1998) AFm, 91, 89-111.
    # * Dai et al. (2004) Journal of Climate, 17, 2281-2299.

    rho_cd <- 0.036    # canopy reflection coeffcient for diffuse PAR
    rho_cb <- 0.029     # canopy reflection coeffcient for direct PAR
    omega <- 0.15                # leaf scattering coefficient of PAR
    k_dash_b <- 0.46 / cos_zenith      # beam & scat PAR ext coef
    k_dash_d <- 0.719     # diffuse & scattered PAR extinction coeff

    # Ross-Goudriaan function is the ratio of the projected area of leaves
    # in the direction perpendicular to the direction of incident solar
    # radiation and the actual leaf area. See Sellers (1985), eqn 13/
    # note this is taken from CABLE code (Kowalczyk '06, eqn 28/29)

    psi1 <- 0.5 - 0.633 * lad
    psi2 <- 0.877 * (1.0 - 2.0 * psi1)
    Gross <- psi1 + psi2 * cos_zenith

    # beam extinction coefficient for black leaves
    kb <<- Gross / cos_zenith

    # Direct-beam irradiance absorbed by sunlit leaves - de P & F, eqn 20b
    Ib <- par * cwdirect_frac
    beam <- Ib * (1.0 - omega) * (1.0 - exp(kb * lai))

    # Diffuse irradiance absorbed by sunlit leaves - de P & F, eqn 20c
    Id <- par * diffuse_frac
    arg1 <- Id * (1.0 - rho_cd)
    arg2 <- 1.0 - exp(-(k_dash_d + kb) * lai)
    arg3 <- k_dash_d / (k_dash_d + kb)
    shaded <- arg1 * arg2 * arg3

    # Scattered-beam irradiance abs. by sunlit leaves - de P & F, eqn 20d
    arg1 <- (1.0 - rho_cb) * (1.0 - exp(-(k_dash_b + kb) * lai))
    arg2 <- k_dash_b / (k_dash_b + kb)
    arg3 <- (1.0 - omega) * (1.0 - exp(-2.0 * kb * lai)) / 2.0
    scattered <- Ib * (arg1 * arg2 - arg3)

    # Total irradiance absorbed by the canopy (Ic) - de P & F, eqn 13
    arg1 <- (1.0 - rho_cb) * Ib * (1.0 - exp(-k_dash_b * lai))
    arg2 <- (1.0 - rho_cd) * Id * (1.0 - exp(-k_dash_d * lai))
    total_canopy_irradiance <- arg1 + arg2

    # Irradiance absorbed by the sunlit fraction of the canopy
    apar_leaf[SUNLIT] <<- beam + scattered + shaded

    # Irradiance absorbed by the shaded fraction of the canopy
    apar_leaf[SHADED] <<- total_canopy_irradiance - apar_leaf[SUNLIT]

    # Calculate sunlit &shdaded LAI of the canopy - de P * F eqn 18
    lai_leaf[SUNLIT] <<- (1.0 - exp(-kb * lai)) / kb
    lai_leaf[SHADED] <<- lai - lai_leaf[SUNLIT]
}

calculate_top_of_canopy_leafn <- function() {
    # Calculate the N at the top of the canopy (g N m-2), N0.
    #
    # References:
    # -----------
    # * Chen et al 93, Oecologia, 93,63-69.

    # leaf mass per area (g C m-2 leaf)
    LMA <- 1.0 / psla * cfracts * KG_AS_G

    if (lai > 0.0) {
        # the total amount of nitrogen in the canopy
        Ntot <- shootnc * LMA * lai

        # top of canopy leaf N (gN m-2)
        N0 <<- Ntot * kn / (1.0 - exp(-kn * lai))
    } else {
        N0 <<- 0.0
    }
}

calculate_top_of_canopy_leafp <- function() {
    # Calculate the P at the top of the canopy (g P m-2), P0, based on
    # calculate_top_of_canopy_leafp relationship

    # leaf mass per area (g C m-2 leaf)
    LMA <- 1.0 / psla * cfracts * KG_AS_G

    if (lai > 0.0) {
        # the total amount of phosphorus in the canopy
        Ptot <- shootpc * LMA * lai

        # top of canopy leaf N (gN m-2)
        P0 <<- Ptot * kp / (1.0 - exp(-kp * lai))
    } else {
        P0 <<- 0.0
    }
}

calc_leaf_to_canopy_scalar <- function() {
    # Calculate scalar to transform leaf Vcmax and Jmax values to big leaf
    # values. Following Wang & Leuning, as long as sunlit and shaded
    # leaves are treated seperately, values of parameters in the coupled
    # model for the two big leaves can be closely approximated by
    # integrating values for individual leaves.
    #
    # - Inserting eqn C6 & C7 into B5
    #
    # per unit ground area
    #
    # Parameters:
    # ----------
    # canopy_wk : structure
    #     various canopy values: in this case the sunlit or shaded LAI &
    #     cos_zenith angle.
    # scalar_sun : float
    #     scalar for sunlit leaves, values returned in unit ground area
    #     (returned)
    # scalar_sha : float
    #     scalar for shaded leaves, values returned in unit ground area
    #     (returned)
    #
    # References:
    # ----------
    # * Wang and Leuning (1998) AFm, 91, 89-111  particularly the Appendix.

    lai_sun <- lai_leaf[SUNLIT]
    lai_sha <- lai_leaf[SHADED]

    cscalar[SUNLIT] <<- (1.0 - exp(-(kb + kn) * lai_sun)) / (kb + kn)
    cscalar[SHADED] <<- (1.0 - exp(-kn * lai_sha)) / kn - cscalar[SUNLIT]
}

initialise_leaf_surface <- function() {
    # initialise values of Tleaf, Cs, dleaf at the leaf surface
    tleaf[ileaf] <<- tair
    dleaf <<- vpd
    Cs <<- Ca
}

calculate_jmaxt_vcmaxt <- function(tleaf) {
    # Calculate the potential electron transport rate (Jmax) and the
    # maximum Rubisco activity (Vcmax) at the leaf temperature.
    #
    # For Jmax  peaked arrhenius is well behaved for tleaf < 0.0
    #
    # Parameters:
    # ----------
    # tleaf : float
    #     air temperature (deg C)
    # jmax : float
    #     the potential electron transport rate at the leaf temperature
    #     (umol m-2 s-1)
    # vcmax : float
    #     the maximum Rubisco activity at the leaf temperature (umol m-2 s-1)

    lower_bound <- 0.0
    upper_bound <- 10.0
    tref <- measurement_temp
    cscalar <- cscalar[ileaf]
    conv <- MMOL_2_MOL * 31.0

    if (modeljm == 0) {
        if (ileaf == SUNLIT) {
            jmax <- jmax * cscalar
            vcmax <- pvcmax * cscalar
        } else {
            jmax <- jmax * cscalar
            vcmax <- psvcmax * cscalar
        }
    } else if (modeljm == 1) {
        vcmax25 <- (vcmaxna * N0 + vcmaxnb) * cscalar

        if (pcycle == TRUE) {
            jmax25n <- (jmaxna * N0 + jmaxnb) * cscalar
            # Ellsworth et al. 2015 Plant, Cell and Environment
            jmax25p <- (400.99 * (P0 * conv) + 88.56) * cscalar
            jmax25 <- min(jmax25n, jmax25p)
        } else {
            jmax25 <- (jmaxna * N0 + jmaxnb) * cscalar
        }
        vcmax <- arrhenius(vcmax25, eav, tleaf, tref)
        jmax <- peaked_arrhenius(jmax25, eaj, tleaf, tref, delsj, edj)
    } else if (modeljm == 2) {
        # NB when using the fixed JV reln, we only apply scalar to Vcmax
        vcmax25 <- (vcmaxna * N0 + vcmaxnb) * cscalar
        jmax25 <- (jv_slope * vcmax25 - jv_intercept)

        vcmax <- arrhenius(vcmax25, eav, tleaf, tref)
        jmax <- peaked_arrhenius(jmax25, eaj, tleaf, tref, delsj, edj)
    } else if (modeljm == 3) {
        jmax25 <- jmax * cscalar
        vcmax25 <- pvcmax * cscalar
        vcmax <- arrhenius(vcmax25, eav, tleaf, tref)
        jmax <- peaked_arrhenius(jmax25, eaj, tleaf, tref, delsj, edj)
    } else {
        stop("You haven't set Jmax/Vcmax model: modeljm \n")
    }

    # reduce photosynthetic capacity with moisture stress
    if (water_balance == BUCKET) {
        jmax <- jmax * wtfac_root
        vcmax <- vcmax * wtfac_root
    } # Should add non-stomal limitation here

    # Jmax/Vcmax forced linearly to zero at low T
    if (tleaf < lower_bound) {
        jmax <- 0.0
        vcmax <- 0.0
    } else if (tleaf < upper_bound) {
        jmax <- jmax * (tleaf - lower_bound) / (upper_bound - lower_bound)
        vcmax <- vcmax * (tleaf - lower_bound) / (upper_bound - lower_bound)
    }

    return (c(jmax, vcmax))
}

photosynthesis_C3 <- function() {
    # Calculate photosynthesis following Farquhar & von Caemmerer, this is an
    # implementation of the routinue in MAESTRA
    #
    # References:
    # -----------
    # * GD Farquhar, S Von Caemmerer (1982) Modelling of photosynthetic
    #   response to environmental conditions. Encyclopedia of plant
    #   physiology 12, 549-587.
    # * Medlyn, B. E. et al (2011) Global Change Biology, 17, 2134-2144.
    # * Medlyn et al. (2002) PCE, 25, 1167-1179, see pg. 1170.

    # Rd0 = 0.92  Dark respiration rate make a paramater!
    qudratic_error <- FALSE
    g0_zero <- 1E-09 # numerical issues, don't use zero

    # unpack some stuff
    idx <- ileaf
    par <- apar_leaf[idx]
    tleaf <- tleaf[idx]

    # Calculate photosynthetic parameters from leaf temperature.
    gamma_star <- calc_co2_compensation_point(tleaf)
    km <- calculate_michaelis_menten(tleaf)
    res__ <- calculate_jmaxt_vcmaxt(tleaf)
    jmax <- res__[1]
    vcmax <- res__[2]

    # leaf respiration in the light, Collatz et al. 1991
    rd <- 0.015 * vcmax
    #rd <- calc_leaf_day_respiration(tleaf, Rd0)

    # actual electron transport rate
    qudratic_error <- FALSE
    large_root <- FALSE
    res__ <- quad(theta, -(alpha_j * par + jmax),
             alpha_j * par * jmax, large_root, qudratic_error)
    J <- res__[1]

    # RuBP regeneration rate
    Vj <- J / 4.0

    # Deal with extreme cases
    if (jmax <= 0.0 || vcmax <= 0.0 || is.nan(J)) {
        an_leaf[idx] <<- -rd
        gsc_leaf[idx] <<- g0_zero
    } else {
        # Hardwiring this for Medlyn gs model for the moment, till I figure
        # out the best structure

        # For the medlyn model this is already in conductance to CO2, so the
        # 1.6 from the corrigendum to Medlyn et al 2011 is missing here
        dleaf_kpa <- dleaf * PA_2_KPA
        if (dleaf_kpa < 0.05) {
            dleaf_kpa <- 0.05
        }

        # This is calculated by SPA hydraulics so we don't need to account for
        # water stress on gs.
        if (water_balance == HYDRAULICS) {
            gs_over_a <- (1.0 + g1 / sqrt(dleaf_kpa)) / Cs
        } else {
            gs_over_a <- (1.0 + (g1 * wtfac_root) / sqrt(dleaf_kpa)) / Cs
        }
        g0 <- g0_zero

        # Solution when Rubisco activity is limiting
        A <- g0 + gs_over_a * (vcmax - rd)
        B <- ( (1.0 - Cs * gs_over_a) * (vcmax - rd) + g0 * (km - Cs) -
               gs_over_a * (vcmax * gamma_star + km * rd) )
        C <- ( -(1.0 - Cs * gs_over_a) * (vcmax * gamma_star + km * rd) -
               g0 * km * Cs )

        # intercellular CO2 concentration
        qudratic_error <- FALSE
        large_root <- TRUE
        res__ <- quad(A, B, C, large_root, qudratic_error)
        Ci <- res__[1]
        qudratic_error <- res__[2]

        if (qudratic_error || Ci <= 0.0 || Ci > Cs) {
            Ac <- 0.0
        } else {
            Ac <- vcmax * (Ci - gamma_star) / (Ci + km)
        }

        # Solution when electron transport rate is limiting
        A <- g0 + gs_over_a * (Vj - rd)
        B <- ( (1. - Cs * gs_over_a) * (Vj - rd) + g0 *
              (2. * gamma_star - Cs) - gs_over_a *
              (Vj * gamma_star + 2.0 * gamma_star * rd) )
        C <- ( -(1.0 - Cs * gs_over_a) * gamma_star * (Vj + 2.0 * rd) -
               g0 * 2.0 * gamma_star * Cs )

        # Intercellular CO2 concentration
        qudratic_error <- FALSE
        large_root <- TRUE
        res__ <- quad(A, B, C, large_root, qudratic_error)
        Ci <- res__[1]

        Aj <- Vj * (Ci - gamma_star) / (Ci + 2.0 * gamma_star)

        # Below light compensation point?
        if (Aj - rd < 1E-6) {
            Ci <- Cs
            Aj <- Vj * (Ci - gamma_star) / (Ci + 2.0 * gamma_star)
        }
        an_leaf[idx] <<- min(Ac, Aj) - rd
        rd_leaf[idx] <<- rd
        gsc_leaf[idx] <<- max(g0, g0 + gs_over_a * an_leaf[idx])
    }

    # Pack calculated values into a temporary array as we may need to
    # recalculate A if water is limiting, i.e. the Emax case below
    if (water_balance == HYDRAULICS) {
        ts_Cs <<- Cs
        ts_vcmax <<- vcmax
        ts_km <<- km
        ts_gamma_star <<- gamma_star
        ts_rd <<- rd
        ts_Vj <<- Vj
    }
}

photosynthesis_C3_emax <- function() {
    # Calculate photosynthesis as above but for here we are resolving Ci and
    # A for a given gs (Jarvis style) to get the Emax solution.
    #
    # This follows MAESPA code.

    g0_zero <- 1E-09 # numerical issues, don't use zero
    qudratic_error <- FALSE

    # Unpack calculated properties from first photosynthesis solution
    idx <- ileaf
    Cs <- ts_Cs
    vcmax <- ts_vcmax
    km <- ts_km
    gamma_star <- ts_gamma_star
    rd <- ts_rd
    Vj <- ts_Vj

    # A very low minimum  for numerical stability.
    if (gsc_leaf[idx] < g0_zero) {
        gsc_leaf[idx] <<- g0_zero
    }
    gs <- gsc_leaf[idx]

    # Solution when Rubisco rate is limiting
    # A = 1.0 / gs
    # B = (0.0 - vcmax) / gs - Cs - km
    # C = vcmax * (Cs - gamma_star)

    # From MAESTRA, not sure of the reason for the difference.
    A <- 1.0 / gs
    B <- (rd - vcmax) / gs - Cs - km
    C <- vcmax * (Cs - gamma_star) - rd * (Cs + km)

    qudratic_error <- FALSE
    large_root <- FALSE
    res__ <- quad(A, B, C, large_root, qudratic_error)
    Ac <- res__[1]
    qudratic_error <- res__[2]
    if (qudratic_error) {
        Ac <- 0.0
    }

    # Solution when electron transport rate is limiting
    A <- 1.0 / gs
    B <- (rd - Vj) / gs - Cs - 2.0 * gamma_star
    C <- Vj * (Cs - gamma_star) - rd * (Cs + 2.0 * gamma_star)

    qudratic_error <- FALSE
    large_root <- FALSE
    res__ <- quad(A, B, C, large_root, qudratic_error)
    Aj <- res__[1]
    qudratic_error <- res__[2]
    if (qudratic_error) {
        Aj <- 0.0
    }

    an_leaf[idx] <<- min(Ac, Aj)
}

calculate_emax <- function() {
    # Assumption that during the day transpiration cannot exceed a maximum
    # value, Emax. At this point we've reached a leaf water potential minimum.
    # Once this point is reached transpiration, gs and A are reclulated
    #
    # Reference:
    # * Duursma et al. 2008, Tree Physiology 28, 265–276
    #

    # plant component of the leaf-specific hydraulic conductance
    # (mmol m–2 s–1 MPa–1)
    kp <- 2.0
    idx <- ileaf
    stressed <- FALSE

    # Need to work out the direct/diffuse weighting to adjust the
    # stressed gs/transpiration calculation.
    #if (ileaf == 0) {
    #    frac = 1.0 - diffuse_frac
    #} else {
    #    frac = diffuse_frac
    #}

    # Hydraulic conductance of the entire soil-to-leaf pathway
    ktot <- 1.0 / (total_soil_resist + 1.0 / kp)

    # Maximum transpiration rate (mmol m-2 s-1)
    emax_leaf <- ktot * (weighted_swp - min_lwp)

    # Leaf transpiration (mmol m-2 s-1), i.e. ignoring boundary layer effects!
    etest <- MOL_2_MMOL * (vpd / press) * gsc_leaf[idx] * GSVGSC

    # leaf water potential (MPa)
    lwp_leaf[idx] <<- calc_lwp(ktot, etest)

    if (etest > emax_leaf) {
        stressed <- TRUE

        # Calculate gs (mol m-2 s-1) given emax_leaf
        gsv <- MMOL_2_MOL * emax_leaf / (vpd / press)
        gsc_leaf[idx] <<- gsv / GSVGSC


        # Need to make sure transpiration solution is consistent, force
        # Tleaf to Tair as we aren't solving this
        trans_leaf[idx] <<- emax_leaf * MMOL_2_MOL

        tleaf[idx] <<- tair
        tleaf_new <<- tair
        rnet_leaf[idx] <<- calc_leaf_net_rad(vpd, apar_leaf[idx] * PAR_2_SW)

        # Minimum leaf water potential reached so recalculate LWP (MPa)
        lwp_leaf[idx] <<- calc_lwp(ktot, emax_leaf)

        # Re-solve An for the new gs
        photosynthesis_C3_emax()

        # Need to calculate an effective beta to use in soil decomposition
        #fwsoil_leaf[idx] = emax_leaf / etest
        fwsoil_leaf[idx] <<- exp(g1 * predawn_swp)
    } else {
        fwsoil_leaf[idx] <<- 1.0

    }
    return (stressed)
}

penman_leaf_wrapper <- function(tleaf, rnet, gsc) {
    # Calculates transpiration by leaves using the Penman-Monteith
    #
    # Parameters:
    # ----------
    # press : float
    #     atmospheric pressure (Pa)
    # rnet : float
    #     net radiation (J m-2 s-1)
    # vpd : float
    #     vapour pressure deficit of air (Pa)
    # tair : float
    #     air temperature (deg C)
    # transpiration : float
    #     transpiration (mol H2O m-2 s-1) (returned)
    # LE : float
    #     latent heat flux, W m-2 (returned)

    # Radiation conductance (mol m-2 s-1)
    gradn <- calc_radiation_conductance(tair)

    # Boundary layer conductance for heat - single sided, forced
    # convection (mol m-2 s-1)
    gbhu <- calc_bdn_layer_forced_conduct(tair, press, wind, leaf_width)

    # Boundary layer conductance for heat - single sided, free convection
    gbhf <- calc_bdn_layer_free_conduct(tair, tleaf, press, leaf_width)

    # Total boundary layer conductance for heat
    gbh <- gbhu + gbhf

    # Total conductance for heat - two-sided
    gh <- 2.0 * (gbh + gradn)

    gbv <- GBVGBH * gbh
    gsv <- GSVGSC * gsc

    # Total leaf conductance to water vapour
    gv <- (gbv * gsv) / (gbv + gsv)
    gbc <- gbh / GBHGBC

    lambda <- calc_latent_heat_of_vapourisation(tair)
    gamma <- calc_pyschrometric_constant(press, lambda)
    slope <- calc_slope_of_sat_vapour_pressure_curve(tair)

    res__ <- penman_monteith(press, vpd, rnet, slope, lambda, gamma, gh, gv)
    gh <- res__[1]
    gv <- res__[2]
    transpiration <- res__[3]
    LE <- res__[4]

    # Calculate decoupling coefficient (McNaughton and Jarvis 1986)
    epsilon <- slope / gamma
    omega <- (1.0 + epsilon) / (1.0 + epsilon + gbv / gsv)
    return (c(transpiration, LE, gbc, gh, gv, omega))
}

solve_leaf_energy_balance <- function() {
    # Wrapper to solve conductances, transpiration and calculate a new
    # leaf temperautre, vpd and Cs at the leaf surface.
    #
    # - The logic broadly follows MAESTRA code, with some restructuring.
    #
    # References
    # ----------
    # * Wang & Leuning (1998) Agricultural & Forest Meterorology, 91, 89-111.

    idx <- ileaf
    sw_rad <- apar_leaf[idx] * PAR_2_SW # W m-2
    rnet_leaf[idx] <<- calc_leaf_net_rad(tair, vpd, sw_rad)
    res__ <- penman_leaf_wrapper(tleaf[idx], rnet_leaf[idx], gsc_leaf[idx])
    transpiration <- res__[1]
    LE <- res__[2]
    gbc <- res__[3]
    gh <- res__[4]
    gv <- res__[5]
    omega <- res__[6]

    # store in structure
    trans_leaf[idx] <<- transpiration
    omega_leaf[idx] <<- omega

    # calculate new Cs, dleaf & tleaf
    Tdiff <- (rnet_leaf[idx] - LE) / (CP * MASS_AIR * gh)
    tleaf_new <<- tair + Tdiff / 4.0
    Cs <<- Ca - an_leaf[idx] / gbc
    dleaf <<- trans_leaf[idx] * press / gv
}

zero_hourly_fluxes <- function() {
    # sunlit / shaded loop
    for (i in 1:NUM_LEAVES) {
        an_leaf[i] <<- 0.0
        rd_leaf[i] <<- 0.0
        gsc_leaf[i] <<- 0.0
        trans_leaf[i] <<- 0.0
        rnet_leaf[i] <<- 0.0
        apar_leaf[i] <<- 0.0
        omega_leaf[i] <<- 0.0
    }
}

scale_leaf_to_canopy <- function() {
    an_canopy <<- an_leaf[SUNLIT] + an_leaf[SHADED]
    rd_canopy <<- rd_leaf[SUNLIT] + rd_leaf[SHADED]
    gsc_canopy <<- gsc_leaf[SUNLIT] + gsc_leaf[SHADED]
    apar_canopy <<- apar_leaf[SUNLIT] + apar_leaf[SHADED]
    trans_canopy <<- trans_leaf[SUNLIT] + trans_leaf[SHADED]
    omega_canopy <<- (omega_leaf[SUNLIT] + omega_leaf[SHADED]) / 2.0
    rnet_canopy <<- rnet_leaf[SUNLIT] + rnet_leaf[SHADED]

    if (water_balance == HYDRAULICS) {
        lwp_canopy <<- (lwp_leaf[SUNLIT] + lwp_leaf[SHADED]) / 2.0

        beta <- (fwsoil_leaf[SUNLIT] + fwsoil_leaf[SHADED]) / 2.0
        wtfac_topsoil <<- beta
        wtfac_root <<- beta
    }
}

sum_hourly_carbon_fluxes <- function() {

    # umol m-2 s-1 gC m-2 30 min-1
    gpp_gCm2 <<- gpp_gCm2 + an_canopy * UMOL_TO_MOL * MOL_C_TO_GRAMS_C * SEC_2_HLFHR
    npp_gCm2 <<- gpp_gCm2 * cue
    gpp <<- gpp_gCm2 * GRAM_C_2_TONNES_HA
    npp <<- npp_gCm2 * GRAM_C_2_TONNES_HA
    auto_resp <<- gpp - npp

    # umol m-2 s-1  J m-2 s-1  MJ m-2 30 min-1
    apar <<- apar + apar_canopy * UMOL_2_JOL * J_TO_MJ * SEC_2_HLFHR
    gs_mol_m2_sec <<- gs_mol_m2_sec + gsc_canopy
}

zero_water_movement <- function() {
    # Losses and gains of water both from PPT and between layers need to be
    # zero'd at the start of each timestep
    for (i in range(1:n_layers)) {
        water_loss[i] <<- 0.0
        water_gain[i] <<- 0.0
        ppt_gain[i] <<- 0.0
        fraction_uptake[i] <<- 0.0
        est_evap[i] <<- 0.0
    }
}

calc_wetting_layers <- function(soil_evap, surface_water) {
    # Tracks surface wetting and drying in the top soil layer and so the
    # thickness of the uppermost dry layer and thus soil evaporation

    seconds_per_step <- 1800.0
    dmin <- 0.001
    airspace <- porosity[0]

    # soil LE should be withdrawn from the wetting layer with the
    # smallest depth..
    ar1 <- 0
    min_val <- 9999.9

    for (i in 1:wetting) {
        if (wetting_bot[i] > 0.0 && wetting_bot[i] < min_val) {
            ar1 <- i
            min_val <- wetting_bot[i]
        }
    }

    # Need to make this into a negative flux of energy from the soil's
    # perspective to match remaining logic here, i.e. negative leaving the
    # surface
    if (soil_evap > 0.0)
        soil_evap <- soil_evap * -1.0

    # Calulate the net change in wetting in the top zone
    netc <- (soil_evap * MM_TO_M) / airspace + (surface_water * MM_TO_M) / airspace

    # wetting
    if (netc > 0.0) {

        # resaturate the layer if top is dry and recharge is greater
        # than dry_thick


        if ((netc > wetting_top[ar1]) && (wetting_top[ar1] > 0.0)) {

            # extra water to deepen wetting layer
            diff <- netc - wetting_top[ar1]
            wetting_top[ar1] <<- 0.0
            if (ar1 > 0) {
                # Not in primary layer (primary layer can't extend deeper)
                wetting_bot[ar1] <<- wetting_bot[ar1] + diff
            }
            dry_thick <<- dmin
        } else {

            if (wetting_top[ar1] == 0.0) {

                # surface is already wet, so extend depth of this wet zone
                if (ar1 > 0) {
                    # not in primary lay (primary layer can't extend deeper)
                    wetting_bot[ar1] <<- wetting_bot[ar1] + netc
                    if (wetting_bot[ar1] >= wetting_top[ar1-1]) {
                        # Layers are conterminous..
                        wetting_top[ar1-1] <<- wetting_top[ar1]
                        wetting_top[ar1] <<- 0.     # remove layer
                        wetting_bot[ar1] <<- 0.     # remove layer
                    }
                }
            } else {

                # Create a new wetting zone
                wetting_top[ar1+1] <<- 0.0
                wetting_bot[ar1+1] <<- netc
            }
            dry_thick <<- dmin
        }

    # Drying
    } else {
        # Drying increases the depth to top of wet soil layers
        wetting_top[ar1] <<- wetting_top[ar1] - netc

        # Wetting layer is dried out.
        if (wetting_top[ar1] > wetting_bot[ar1]) {
            # How much more drying is there?
            diff <- wetting_top[ar1] - wetting_bot[ar1]
            wetting_top[ar1] <<- 0.0
            wetting_bot[ar1] <<- 0.0
            ar2 <- ar1 - 1

            # Move to deeper wetting layer
            if (ar2 > 0) {
                # dry out deeper layer
                wetting_top[ar2] <<- wetting_top[ar2] + diff
                dry_thick <<- max(dmin, wetting_top[ar2])
            # no deeper layer
            } else {
                # layer 1 is dry
                dry_thick <<- thickness[0]
            }
        } else {
            dry_thick <<- max(dmin, wetting_top[ar1])
        }
    }

    if (dry_thick == 0.0) {
        stop("Problem in dry_thick\n")
    }
}

extract_water_from_layers <- function(soil_evap, transpiration) {

    # Extract soil evaporation and transpiration from the soil profile

    # Is soil evap taken from first or second layer?
    if (dry_thick < thickness[0]) {
        # The dry zone does not extend beneath the top layer
        rr <- 0
    } else {
        # The dry zone does extend beneath the top layer
        rr <- 1
    }

    if (soil_evap > 0.0) {
        water_loss[rr] <<- water_loss[rr] + soil_evap * MM_TO_M
    } # ignoring water gain due to due formation...


    # Determing water loss from each layer due to transpiration
    for (i in 1:rooted_layers) {
        water_loss[i] <<- water_loss[i] + (transpiration * MM_TO_M) * fraction_uptake[i]
    }
}

soil_water_store <- function(time_dummy, y, dydt,
                      unsat, drain_layer, cond1,
                      cond2, cond3) {

    # numerical lib vectors are index 1..n, so we need to index the return
    # from the odeint func with 1, not 0
    index <- 1

    drainage <- calc_soil_conductivity(y[index], cond1, cond2, cond3)

    # Convert units, soil conductivity is in m s-1
    drainage <- drainage * SEC_2_HLFHR

    # gravitational drainage above field_capacity
    if (y[index] <= drain_layer) {
        drainage <- 0.0
    }

    # layer below cannot accept more water than unsat
    if (drainage > unsat) {
        drainage <- unsat
    }

    # waterloss from this layer
    dydt[index] <- -drainage

    return (dydt)
}

rkck <- function(y, dydx, n, x, h, yout,
	      yerr, aa, bb, cc, dd, ee,
	      derivs) {
	a2 <- 0.2
    a3 <- 0.3
    a4 <- 0.6
    a5 <- 1.0
    a6 <- 0.875
    b21 <- 0.2
    b31 <- 3.0/40.0
    b32 <- 9.0/40.0
    b41 <- 0.3
    b42 <- -0.9
    b43 <- 1.2
    b51 <- -11.0/54.0
    b52 <- 2.5
    b53 <- -70.0/27.0
    b54 <- 35.0/27.0
    b61 <- 1631.0/55296.0
    b62 <- 175.0/512.0
    b63 <- 575.0/13824.0
    b64 <- 44275.0/110592.0
    b65 <- 253.0/4096.0
    c1 <- 37.0/378.0
    c3 <- 250.0/621.0
    c4 <- 125.0/594.0
    c6 <- 512.0/1771.0
    dc5 <- -277.0/14336.0
	dc1 <- c1 - 2825.0/27648.0
    dc3 <- c3 - 18575.0/48384.0
    dc4 <- c4 - 13525.0/55296.0
    dc6 <- c6 - 0.25

	#       *ak2, *ak3, *ak4, *ak5, *ak6, *ytemp
	# ak2=dvector(1,n)
	# ak3=dvector(1,n)
	# ak4=dvector(1,n)
	# ak5=dvector(1,n)
	# ak6=dvector(1,n)
	# ytemp=dvector(1,n)

	for (i in 1:n)
		ytemp[i] <<- y[i]+b21*h*dydx[i]
	ak2 <<- derivs(x+a2*h,ytemp,ak2, aa, bb, cc, dd, ee)
	for (i in 1:n)
		ytemp[i] <<- y[i]+h*(b31*dydx[i]+b32*ak2[i])
	ak3 <<- derivs(x+a3*h,ytemp,ak3, aa, bb, cc, dd, ee)
	for (i in 1:n)
		ytemp[i] <<- y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i])
	ak4 <<- derivs(x+a4*h,ytemp,ak4, aa, bb, cc, dd, ee)
	for (i in 1:n)
		ytemp[i] <<- y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i])
	ak5 <<- derivs(x+a5*h,ytemp,ak5, aa, bb, cc, dd, ee)
	for (i in 1:n)
		ytemp[i] <<- y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i])
	ak6 <<- derivs(x+a6*h,ytemp,ak6, aa, bb, cc, dd, ee)
	for (i in 1:n)
		yout[i] <- y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i])
	for (i in 1:n)
		yerr[i] <- h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i])

    return (c(yout, yerr))
	# free_dvector(ytemp, 1, n)
	# free_dvector(ak6, 1, n)
	# free_dvector(ak5, 1, n)
	# free_dvector(ak4, 1, n)
	# free_dvector(ak3, 1, n)
	# free_dvector(ak2, 1, n)
}

rkqs <- function(y, dydx, n, x, htry, eps,
	      yscal, hdid, hnext,
		  aa, bb, cc, dd, ee,
	      derivs) {
	# *yerr,*ytemp

	# yerr=dvector(1,n)
	# ytemp=dvector(1,n)

	h <- htry
	while (TRUE) {
		res__ <- rkck(y, dydx, n, x, h, ytemp, yerr, aa, bb, cc, dd, ee, derivs)
        ytemp <<- res__[1]
        yerr <<- res__[2]
		errmax <- 0.0
		for (i in 1:n) {
            errmax <- max(errmax, abs(yerr[i] / yscal[i]))
        }
		errmax <- errmax / eps
		if (errmax > 1.0) {
			h <- SAFETY * h * (errmax ^ PSHRNK)
			if (h < 0.1 * h) h <- h * 0.1
			xnew <- x + h
			if (xnew == x) {
                nrerror("stepsize underflow in rkqs")
            }
			next
		} else {
			if (errmax > ERRCON) {
                hnext <- SAFETY * h * (errmax ^ PGROW)
            } else {
                hnext <- 5.0 * h
            }
            hdid <- h
			x <- x + hdid
			for (i in 1:n) {
                y[i] <- ytemp[i]
            }
			break
		}
	}
    return (c(y, x, hdid, hnext))
	# free_dvector(ytemp,1,n)
	# free_dvector(yerr,1,n)
}

odeint <- function(ystart, nvar, x1, x2, eps,
            h1, hmin, aa, bb, cc, dd, ee,
            derivs, rkqs) {
	#  *yscal=NULL, *y=NULL, *dydx=NULL, *xp=NULL
	# **yp=NULL

	# initialising this within the func, which means this isn't generic
	# kmax = 100

	# xp = dvector(1, kmax)
	# yp = dmatrix(1,nvar,1,kmax)
	# yscal=dvector(1,nvar)
	# y=dvector(1,nvar)
	# dydx=dvector(1,nvar)


    dxsav <- (x2 - x1) / 20.0
    x <- x1

	h <- SIGN(h1, x2 - x1)
	nok <- 0
    nbad <- 0
    kount <- 0
    hdid <- 0
    hnext <- 0

    # printf("* %lf\n", ystart[1])
	for (i in 1:nvar) {
        y[i] <<- ystart[i]
    }
	if (kmax > 0) {
        xsav <- x - dxsav * 2.0
    }
	for (nstp in 1:MAXSTP) {
		dydx <- derivs(x, y, dydx, aa, bb, cc, dd, ee)
		for (i in 1:nvar) {
			yscal[i] <<- abs(y[i]) + abs(dydx[i] * h) + TINY
        }
		if (kmax > 0 && kount < kmax - 1 && abs(x - xsav) > abs(dxsav)) {
            kount <- kount + 1
			xp[kount] <<- x
			for (i in 1:nvar) {
                yp[i][kount] <<- y[i]
            }
			xsav <- x
		}

		if ((x + h - x2) * (x + h - x1) > 0.0) {
            h <- x2 - x
        }
		res__ <- rkqs(y, dydx, nvar, x, h, eps, yscal, hdid, hnext,
                aa, bb, cc, dd, ee, derivs)
        y <- res__[1]
        x <- res__[2]
        hdid <- res__[3]
        hnext <- res__[4]
		if (hdid == h) {
            nok <- nok + 1
        } else {
            nbad <- nbad + 1
        }

		if ((x - x2) * (x2 - x1) >= 0.0) {
            for (i in 1:nvar) {
                ystart[i] <- y[i]
            }

			if (kmax) {
                kount <- kount + 1
				xp[kount] <<- x
                for (i in 1:nvar) {
                    yp[i][kount] <<- y[i]
                }
			}

			# free_dvector(dydx,1,nvar)
			# free_dvector(y,1,nvar)
			# free_dvector(yscal,1,nvar)
            # free_dvector(xp,1,kmax)
            # free_dmatrix(yp,1,nvar,1,kmax)

            # printf("** %lf\n", ystart[1])
			return (c(ystart, nok, nbad))
		}
		if (abs(hnext) <= hmin) {
            nrerror("Step size too small in odeint")
        }
		h <- hnext
	}
	nrerror("Too many steps in routine odeint")
}

calc_soil_balance <- function(soil_layer, water_lost) {
    # Integrator for soil gravitational drainage

    N <- 1
    eps <- 1.0e-4        # precision
    h1 <- .001           # first guess at integrator size
    hmin <- 0.0          # minimum value of the integrator step
    x1 <- 1.0             # initial time
    x2 <- 2.0             # final time

    # value affecting the max time interval at which variables should b calc
    soilpor <- porosity[soil_layer]
    # *ystart <- NULL
    # ystart <- dvector(1,N)

    for (i in 1:N) {
        ystart[i] <<- 0.0
    }

    # unsaturated volume of layer below (m3 m-2)
    unsat <- max(0.0, (porosity[soil_layer + 1] -
                      water_frac[soil_layer + 1]) *
                      thickness[soil_layer + 1] / thickness[soil_layer])

    # soil water capacity of the current layer
    drain_layer <- field_capacity[soil_layer]
    liquid <- water_frac[soil_layer]


    # initial conditions i.e. is there liquid water and more
    # water than layer can hold

    if (liquid > 0.0 && liquid > drain_layer) {

        # ystart is a vector 1..N, so need to index from 1 not 0
        ystart[1] <<- water_frac[soil_layer]

        # Runge-Kunte ODE integrator used to estimate soil gravitational
        # drainage during each time-step
        res__ <- odeint(ystart, N, x1, x2, eps, h1, hmin, unsat,
               drain_layer, cond1[soil_layer], cond2[soil_layer],
               cond3[soil_layer], soil_water_store, rkqs)
        ystart <<- res__[1]
        nok <- res__[2]
        nbad <- res__[3]

        # ystart is a vector 1..N, so need to index from 1
        new_water_frac <- ystart[1]

        # convert from water fraction to absolute amount (m)
        change <- (water_frac[soil_layer] - new_water_frac) *
                   thickness[soil_layer]

        # update soil layer below with drained liquid
        if (soil_layer + 1 < n_layers) {
            water_gain[soil_layer + 1] <<- water_gain[soil_layer + 1] + change
        } else {
            # We are draining through the bottom soil layer, add to runoff
            water_lost <- water_lost + change
        }
        water_loss[soil_layer] <<- water_loss[soil_layer] + change

    }

    if (water_loss[soil_layer] < 0.0) {
        print(sprintf("waterloss probem in soil_balance: %d %f\n",
                soil_layer, water_loss[soil_layer]))
    }

    return (water_lost)
}

update_soil_water_storage <- function(soil_evap, transpiration) {
    # Update the soil water storage at the end of the timestep

    root_zone_total <- 0.0
    for (i in 1:n_layers) {

        # water content of soil layer (m)
        water_content <- water_frac[i] * thickness[i]

        needed <- water_content + water_gain[i] +
                  ppt_gain[i] - water_loss[i]

        # Is soil evap taken from first or second layer?
        if (dry_thick < thickness[0]) {
            # The dry zone does not extend beneath the top layer
            rr <- 0
        } else {
            # The dry zone does extend beneath the top layer
            rr <- 1
        }

        # Correction for potential to over-evaporate if using Emax drought
        # stress correction. This stops that happening.
        if (i == rr) {
            if (needed < 0.0) {
                prev_soil_evap <- soil_evap
                soil_evap <- max(0.0, soil_evap + (needed * M_TO_MM))
                if (soil_evap > 0.0) {
                    taken <- (prev_soil_evap - soil_evap) * MM_TO_M
                    needed <- needed - taken
                    water_loss[i] <<- water_loss[i] + taken
                }
            }
        } else {
            if (needed < 0.0) {
                prev_trans <- transpiration
                transpiration <- max(0.0, transpiration + (needed * M_TO_MM))
                taken <- (prev_trans - transpiration) * MM_TO_M
                water_loss[i] <<- water_loss[i] + taken
            }
        }

        # NB water gain here is drainage from the layer above
        water_content <- max(0.0, water_content +
                                  water_gain[i] +
                                  ppt_gain[i] -
                                  water_loss[i])

        # Determine volumetric water content water content of layer (m3 m-3)
        water_frac[i] <<- water_content / thickness[i]

        # update old GDAY effective two-layer buckets
        # - this is just for outputting, these aren't used.
        if (i == 0) {
            pawater_topsoil <<- water_content * M_TO_MM
        } else {
            root_zone_total <- root_zone_total + water_content * M_TO_MM
        }
    }
    pawater_root <<- root_zone_total

    return (c(soil_evap, transpiration))
}

sum_hourly_water_fluxes <- function(soil_evap_hlf_hr,
                             transpiration_hlf_hr, et_hlf_hr,
                             interception_hlf_hr,
                             thoughfall_hlf_hr,
                             canopy_evap_hlf_hr,
                             runoff_hlf_hr,
                             omega_hlf_hr,
                             rain_hlf_hr) {

    # add half hour fluxes to day total store
    soil_evap <<- soil_evap + soil_evap_hlf_hr
    transpiration <<- transpiration + transpiration_hlf_hr
    et <<- et + et_hlf_hr
    interception <<- interception + interception_hlf_hr
    throughfall <<- throughfall + thoughfall_hlf_hr
    canopy_evap <<- canopy_evap + canopy_evap_hlf_hr
    runoff <<- runoff + runoff_hlf_hr
    omega <<- omega + omega_hlf_hr # average at the end of hour loop
    day_ppt <<- day_ppt + rain_hlf_hr
}

calculate_water_balance_sub_daily <- function(daylen, trans, omega_leaf, rnet_leaf) {
    # Calculate the water balance (including all water fluxes).
    # - we are using all the hydraulics instead
    #
    # Parameters:
    # ----------
    # control : structure
    #     control structure
    # fluxes : structure
    #     fluxes structure
    # met : structure
    #     meteorological drivers structure
    # params : structure
    #     parameters structure
    # day : int
    #     project day. (Dummy argument, only passed for daily model)
    # daylen : double
    #     length of day in hours. (Dummy argument, only passed for daily
    #     model)
    # trans : double
    #     transpiration (Dummy argument, only passed for sub-daily
    #     model)
    # omega_leaf : double
    #     decoupling coefficient (Dummy argument, only passed for sub-daily
    #     model)
    # rnet_leaf : double
    #     total canopy rnet (Dummy argument, only passed for sub-daily model)

    # Water drained through the bottom soil layer
    water_lost <- 0.0

    if (water_balance == HYDRAULICS) {

        zero_water_movement()

        # calculate potential canopy evap rate, this may be reduced later
        # depending on canopy water storage
        canopy_evap <- calc_canopy_evaporation(rnet_leaf)

        # mol m-2 s-1 to mm d-1
        conv <- MOLE_WATER_2_G_WATER * G_TO_KG * SEC_2_HLFHR
        canopy_evap <- canopy_evap * conv

        # We could now replace this interception bit with the Rutter scheme?
        res__ <- calc_interception()
        surface_water <- res__[1]
        interception <- res__[2]
        canopy_evap <- res__[3]

        net_rad <- calc_net_radiation(sw_rad, tair)
        soil_evap <- calc_soil_evaporation(net_rad)
        soil_evap <- soil_evap * MOLE_WATER_2_G_WATER * G_TO_KG * SEC_2_HLFHR


        # mol m-2 s-1 to mm/30 min
        transpiration <- trans * MOLE_WATER_2_G_WATER * G_TO_KG * SEC_2_HLFHR

        et <- transpiration + soil_evap + canopy_evap

        # The loop needs to be outside the func as we need to be able to
        # calculate the soil conductance per layer and call this via
        # the integration func when we update the soil water balance

        for (i in 1:n_layers) {
            soil_conduct[i] <<- calc_soil_conductivity(water_frac[i],
                                                        cond1[i],
                                                        cond2[i],
                                                        cond3[i])
        }

        calc_soil_water_potential()
        calc_soil_root_resistance()

        # If we have leaves we are transpiring
        if (lai > 0.0) {
            calc_water_uptake_per_layer()
        }

        # Calculates the thickness of the top dry layer and determines water
        # lost in upper layers due to evaporation
        calc_wetting_layers(soil_evap, surface_water)
        extract_water_from_layers(soil_evap, transpiration)

        # determines water movement between soil layers due drainage
        # down the profile

        for (i in 1:n_layers) {
            water_lost <- calc_soil_balance(i, water_lost)
        }

        # how much surface water infiltrantes the first soil layer in the
        # current time step? Water which does not infiltrate in a single step
        # is considered runoff

        runoff <- calc_infiltration(surface_water)
        runoff <- runoff + water_lost * M_TO_MM

        # Don't see point of calculating these again
        # Find SWP & soil resistance without updating waterfrac yet
        calc_soil_water_potential()
        calc_soil_root_resistance()

        res__ <- update_soil_water_storage(soil_evap, transpiration)
        soil_evap <- res__[1]
        transpiration <- res__[2]
        et <- transpiration + soil_evap + canopy_evap
    } else {

        # Simple soil water bucket appoximation

        # calculate potential canopy evap rate, this may be reduced later
        # depending on canopy water storage
        canopy_evap <- calc_canopy_evaporation(rnet_leaf)

        # mol m-2 s-1 to mm/day
        conv <- MOLE_WATER_2_G_WATER * G_TO_KG * SEC_2_HLFHR
        canopy_evap <- canopy_evap * conv
        res__ <- calc_interception()
        surface_water <- res__[1]
        interception <- res__[2]
        canopy_evap <- res__[3]

        net_rad <- calc_net_radiation(sw_rad, tair)
        soil_evap <- calc_soil_evaporation(net_rad)
        soil_evap <- soil_evap * MOLE_WATER_2_G_WATER * G_TO_KG * SEC_2_HLFHR

        # mol m-2 s-1 to mm/30 min
        transpiration <- trans * MOLE_WATER_2_G_WATER * G_TO_KG * SEC_2_HLFHR

        # NB. et, transpiration & soil evap may all be adjusted in
        # update_water_storage if we don't have sufficient water

        et <- transpiration + soil_evap + canopy_evap


        res__ <- update_water_storage(surface_water, interception, canopy_evap, transpiration, soil_evap)
        transpiration <- res__[1]
        soil_evap <- res__[2]
        et <- res__[3]
        runoff <- res__[4]

    }

    sum_hourly_water_fluxes(soil_evap, transpiration, et, interception,
                            surface_water, canopy_evap, runoff, omega_leaf, rain)
}

write_subdaily_outputs_ascii <- function(year, doy, hod) {
    # Write sub-daily canopy fluxes - very basic for now

    # time stuff
    ofp_sd[nrow(ofp_sd) + 1, ] <<- c(year, doy, hod, an_canopy, rd_canopy, gsc_canopy, apar_canopy, trans_canopy, tleaf_new)
}

check_water_balance <- function (previous_sw, current_sw, previous_cs,
                                 current_cs, year, doy) {

    delta_sw <- current_sw - previous_sw
    delta_cs <- current_cs - previous_cs

    day_wbal <<- day_ppt - (runoff + et + delta_cs + delta_sw)

    print(sprintf(
        "%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
        year, doy, day_ppt, runoff,
        et, delta_sw, delta_cs, previous_sw, current_sw,
        predawn_swp, midday_lwp, day_wbal
    ))
}

canopy <- function() {
    # Canopy module consists of two parts:
    # (1) a radiation sub-model to calculate apar of sunlit/shaded leaves
    #     - this is all handled in radiation.c
    # (2) a coupled model of stomatal conductance, photosynthesis and
    #     the leaf energy balance to solve the leaf temperature and partition
    #     absorbed net radiation between sensible and latent heat.
    # - The canopy is represented by a single layer with two big leaves
    #   (sunlit & shaded).
    #
    # - The logic broadly follows MAESTRA code, with some restructuring.
    #
    # References
    # ----------
    # * Wang & Leuning (1998) Agricultural & Forest Meterorology, 91, 89-111.
    # * Dai et al. (2004) Journal of Climate, 17, 2281-2299.
    # * De Pury & Farquhar (1997) PCE, 20, 537-557.
    iter <- 0
    sitermax <- 100
    dummy <- 0
    debug <- TRUE
    stressed <- FALSE
    dummy2 <- 0.0

    # loop through the day
    zero_carbon_day_fluxes()
    zero_water_day_fluxes()
    previous_sw <- pawater_topsoil + pawater_root
    previous_cs <- canopy_store
    sunlight_hrs <- 0
    doy <- ma$doy[hour_idx]
    year <- ma$year[hour_idx]

    for (hod in 1:num_hlf_hrs) {
        unpack_met_data(hod, dummy2)

        # calculates diffuse frac from half-hourly incident radiation
        unpack_solar_geometry()

        # Is the sun up?
        if (elevation > 0.0 && par > 20.0) {
            calculate_absorbed_radiation(par)
            calculate_top_of_canopy_leafn()
            calculate_top_of_canopy_leafp()
            calc_leaf_to_canopy_scalar()

            # sunlit / shaded loop

            for (i in 1:NUM_LEAVES) {
                ileaf <<- i
                # initialise values of Tleaf, Cs, dleaf at the leaf surface
                initialise_leaf_surface()

                # Leaf temperature loop
                while (TRUE) {

                    if (ps_pathway == C3) {
                        photosynthesis_C3()
                    } else {
                        # Nothing implemented
                        stop("C4 photosynthesis not implemented\n")
                    }

                    if (an_leaf[ileaf] > 1E-04) {

                        if (ater_balance == HYDRAULICS) {
                            # Ensure transpiration does not exceed Emax, if it
                            # does we recalculate gs and An
                            stressed <- calculate_emax()
                        }

                        if (stressed == FALSE) {
                            # Calculate new Cs, dleaf, Tleaf
                            solve_leaf_energy_balance()
                        }
                    } else {
                        break
                    }

                    if (iter >= itermax) {
                        stop("No convergence in canopy loop:\n")
                    } else if (abs(tleaf[ileaf] - tleaf_new) < 0.02) {
                        break
                    }

                    # Update temperature & do another iteration
                    tleaf[ileaf] <<- tleaf_new
                    iter <- iter + 1
                } # end of leaf temperature loop
            } # end of sunlit/shaded leaf loop
        } else {
            zero_hourly_fluxes()

            # set tleaf to tair during the night
            tleaf[SUNLIT] <<- tair
            tleaf[SHADED] <<- tair


            # pre-dawn soil water potential (MPa), clearly one should link this
            # the actual sun-rise :). Here 10 = 5 am, 10 is num_half_hr

            if (water_balance == HYDRAULICS && hod == 10) {
                predawn_swp <<- weighted_swp
                # _calc_soil_water_potential(c, p, s)=
            }

        }
        scale_leaf_to_canopy()
        if (water_balance == HYDRAULICS && hod == 24) {
            midday_lwp <<- lwp_canopy
        }
        sum_hourly_carbon_fluxes()

        calculate_water_balance_sub_daily(dummy, trans_canopy, omega_canopy, rnet_canopy)

        if (print_options == SUBDAILY && spin_up == FALSE) {
            write_subdaily_outputs_ascii(year, doy, hod)
        }
        hour_idx <<- hour_idx + 1
        sunlight_hrs <- sunlight_hrs + 1
    } # end of hour loop

    # work out average omega for the day over sunlight hours
    omega <<- omega / sunlight_hrs

    if (water_stress) {
        # Calculate the soil moisture availability factors [0,1] in the
        # topsoil and the entire root zone
        calculate_soil_water_fac()

        #printf("%lf %.10lf\n", wtfac_root, saved_swp)
    } else {
        wtfac_topsoil <<- 1.0
        wtfac_root <<- 1.0
    }

    current_sw <- pawater_topsoil + pawater_root
    current_cs <- canopy_store
    check_water_balance(previous_sw, current_sw, previous_cs, current_cs, year, doy)
}

calc_day_growth <- function(day_length, doy, fdecay, rdecay) {
    dummy <- 0.0

    # Store the previous days soil water store
    previous_topsoil_store <- pawater_topsoil
    previous_rootzone_store <- pawater_root

    previous_sw <- pawater_topsoil + pawater_root
    previous_cs <- canopy_store
    year <- ma$year[day_idx]

    if (sub_daily) {
        # calculate 30 min two-leaf GPP/NPP, respiration and water fluxes
        canopy()
    } else {
        # calculate daily GPP/NPP, respiration and update water balance
        carbon_daily_production(day_length)
        calculate_water_balance(day_length, dummy, dummy, dummy)


        current_sw <- pawater_topsoil + pawater_root
        current_cs <- canopy_store
        day_ppt <<- rain

        # check_water_balance(c, f, s, previous_sw, current_sw, previous_cs,
        #                     current_cs, year, doy)
    }

    #  leaf N:C as a fraction of Ncmaxyoung, i.e. the max N:C ratio of
    # foliage in young stand, and leaf P:C as a fraction of Pcmaxyoung
    nitfac <- min(1.0, shootnc / ncmaxfyoung)
    pitfac <- min(1.0, shootpc / pcmaxfyoung)

    # checking for pcycle control parameter
    if (pcycle == TRUE) {
        npitfac <- min(nitfac, pitfac)
    } else {
        npitfac <- nitfac
    }

    # figure out the C allocation fractions
    if (deciduous_model){
        # Allocation is annually for deciduous "tree" model, but we need to
        # keep a check on stresses during the growing season and the LAI
        # figure out limitations during leaf growth period. This also
        # applies for deciduous grasses, need to do the growth stress
        # calc for grasses here too.
        if (leaf_out_days[doy] > 0.0) {

            calc_carbon_allocation_fracs(npitfac)

            # store the days allocation fraction, we average these at the
            # end of the year (for the growing season)
            avg_alleaf <<- avg_alleaf + alleaf
            avg_albranch <<- avg_albranch + albranch
            avg_alstem <<- avg_alstem + alstem
            avg_alroot <<- avg_alroot + alroot
            avg_alcroot <<- avg_alcroot + alcroot

        }
    } else {
        # daily allocation...
        calc_carbon_allocation_fracs(npitfac)
    }

    # Distribute new C, N and P through the system
    carbon_allocation(npitfac, doy)

    res__ <- calculate_cnp_wood_ratios(npitfac, nitfac, pitfac)
    ncbnew <- res__[1]
    nccnew <- res__[2]
    ncwimm <- res__[3]
    ncwnew <- res__[4]
    pcbnew <- res__[5]
    pccnew <- res__[6]
    pcwimm <- res__[7]
    pcwnew <- res__[8]


    recalc_wb <- np_allocation(ncbnew, nccnew, ncwimm, ncwnew,
                              pcbnew, pccnew, pcwimm, pcwnew,
                              fdecay, rdecay, doy)

    if (exudation && alloc_model != GRASSES) {
        calc_root_exudation()
    }

    # If we didn't have enough N available to satisfy wood demand, NPP
    # is down-regulated and thus so is GPP. We also need to recalculate the
    # water balance given the lower GPP.
    if (recalc_wb) {
        pawater_topsoil <<- previous_topsoil_store
        pawater_root <<- previous_rootzone_store

        if (sub_daily) {
            # reduce transpiration to match cut back GPP
            # -there isn't an obvious way to make this work at the 30 min
            #  timestep, so invert T from WUE assumption and use that
            #  to recalculate the end day water balance

            transpiration <<- gpp_gCm2 / wue
            update_water_storage_recalwb()

        } else {
            calculate_water_balance(day_length, dummy, dummy, dummy)
        }
    }

    update_plant_state(fdecay, rdecay, doy)

    precision_control()
}

calculate_decay_rates <- function() {
    # Model decay rates - decomposition rates have a strong temperature
    # and moisture dependency. Note same temperature is assumed for all 3
    # SOM pools, found by Knorr et al (2005) to be untrue. N mineralisation
    # depends on top soil moisture (most variable) (Connell et al. 1995)
    #
    # References:
    # -----------
    # Knorr et al. (2005) Nature, 433, 298-301.
    # Connell et al. (1995) Biol. Fert. Soils, 20, 213-220.

    # abiotic decomposition factor - impact of soil moisture
    # and soil temperature on microbial activity
    adfac <- wtfac_root * tfac_soil_decomp

    # Effect of soil texture (silt + clay content) on active SOM turnover
    # higher turnover for sandy soils
    soil_text <- 1.0 - (0.75 * finesoil)

    # Impact of lignin content
    lignin_cont_leaf <- exp(-3.0 * ligshoot)
    lignin_cont_root <- exp(-3.0 * ligroot)

    # decay rate of surface structural pool
    decayrate[1] <<- kdec1 * lignin_cont_leaf * adfac

    # decay rate of surface metabolic pool
    decayrate[2] <<- kdec2 * adfac

    # decay rate of soil structural pool
    decayrate[3] <<- kdec3 * lignin_cont_root * adfac

    # decay rate of soil metabolic pool
    decayrate[4] <<- kdec4 * adfac

    # decay rate of active pool
    decayrate[5] <<- kdec5 * soil_text * adfac

    # decay rate of slow pool
    decayrate[6] <<- kdec6 * adfac

    # decay rate of passive pool
    decayrate[7] <<- kdec7 * adfac
}

ratio_of_litternc_to_live_leafnc <- function() {
    # ratio of litter N:C to live leaf N:C
    #
    # Returns:
    # --------
    # nc_leaf_litter : float
    #     N:C ratio of litter to foliage

    if (use_eff_nc){
        nc_leaf_litter <- liteffnc * (1.0 - fretrans)
    } else {
        if (float_eq(deadleaves, 0.0)){
            nc_leaf_litter <- 0.0
        } else {
            nc_leaf_litter <- deadleafn / deadleaves
        }
    }

    return (nc_leaf_litter)
}

ratio_of_litternc_to_live_rootnc <- function() {
    # ratio of litter N:C to live root N:C
    #
    # Returns:
    # --------
    # nc_root_litter : float
    #     N:C ratio of litter to live root

    if (use_eff_nc){
        nc_root_litter <- liteffnc * ncrfac *  (1.0 - rretrans)
    } else {
        if (float_eq(deadroots, 0.0)){
            nc_root_litter <- 0.0
        } else{
            nc_root_litter <- deadrootn / deadroots
        }
    }

    return (nc_root_litter)
}

calc_ligin_nratio_leaves <- function() {
    # Estimate Lignin/N ratio, as this dictates the how plant litter is
    # seperated between metabolic and structural pools.
    #
    # Returns:
    # --------
    # lnleaf : float
    #     lignin:N ratio of leaf


    nc_leaf_litter <- ratio_of_litternc_to_live_leafnc()

    if (float_eq(nc_leaf_litter, 0.0)) {
        # catch divide by zero if we have no leaves
        lnleaf <- 0.0
    } else {
        lnleaf <- ligshoot / cfracts / nc_leaf_litter
    }

    return (lnleaf)

}

calc_ligin_nratio_fine_roots <- function() {
    # Estimate Lignin/N ratio, as this dictates the how plant litter is
    # seperated between metabolic and structural pools.
    #
    # Returns:
    # --------
    # lnroot : float
    #     lignin:N ratio of fine root

    nc_root_litter <- ratio_of_litternc_to_live_rootnc()


    if (float_eq(nc_root_litter, 0.0)) {
        # catch divide by zero if we have no roots
        lnroot <- 0.0
    } else {
        lnroot <- ligroot / cfracts / nc_root_litter
    }

    return (lnroot)
}

flux_from_grazers <- function() {
    #  Input from faeces
    if (grazing) {
        fmfaeces <<- metafract(ligfaeces * faecescn / cfracts)
        faecesc <<- ceaten * fracfaeces
    } else {
        fmfaeces <<- 0.0
        faecesc <<- 0.0
    }
}

partition_plant_litter <- function() {
    # Partition litter from the plant (surface) and roots into metabolic
    # and structural pools

    # Surface (leaves, branches, stem) Litter

    # ...to the structural pool
    leaf_material <- deadleaves * (1.0 - fmleaf)
    wood_material <- deadbranch + deadstems
    faeces_material <- faecesc * (1.0 - fmfaeces)
    surf_struct_litter <<- leaf_material + wood_material + faeces_material

    # ...to the metabolic pool
    surf_metab_litter <<- deadleaves * fmleaf + faecesc * fmfaeces

    # Root Litter

    # ...to the structural pool
    soil_struct_litter <<- deadroots * (1.0 - fmroot) + deadcroots

    # ...to the metabolic pool
    soil_metab_litter <<- deadroots * fmroot
}

cfluxes_from_structural_pool <- function() {

    # Send structural c fluxes to other SOM pools

    structout_surf <- structsurf * decayrate[1]
    structout_soil <- structsoil * decayrate[3]

    # C flux surface structural pool  slow pool
    surf_struct_to_slow <<- structout_surf * ligshoot * 0.7

    # C flux surface structural pool active pool
    surf_struct_to_active <<- structout_surf * (1.0 - ligshoot) * 0.55

    # C flux soil structural pool slow pool
    soil_struct_to_slow <<- structout_soil * ligroot * 0.7

    # soil structural pool active pool
    soil_struct_to_active <<- structout_soil * (1.0 - ligroot) * 0.45

    # Respiration fluxes

    # CO2 lost during transfer of structural C to the slow pool
    co2_to_air[1] <<- (structout_surf *
                        (ligshoot * 0.3 + (1.0 - ligshoot) * 0.45))

    # CO2 lost during transfer structural C  to the active pool
    co2_to_air[2] <<- (structout_soil *
                        (ligroot * 0.3 + (1.0 - ligroot) * 0.55))

}

cfluxes_from_metabolic_pool <- function() {
    # Send C from metabolic pools to other SOM pools

    # C flux surface metabolic pool  active pool
    surf_metab_to_active <<- metabsurf * decayrate[2] * 0.45

    # C flux soil metabolic pool  active pool
    soil_metab_to_active <<- metabsoil * decayrate[4] * 0.45

    # Respiration fluxes
    co2_to_air[3] <<- metabsurf * decayrate[2] * 0.55
    co2_to_air[4] <<- metabsoil * decayrate[4] * 0.55
}

cfluxes_from_active_pool <- function(frac_microb_resp) {
    # Send C fluxes from active pool to other SOM pools

    activeout <- activesoil * decayrate[5]

    # C flux active pool slow pool
    active_to_slow <<- activeout * (1.0 - frac_microb_resp - 0.004)

    # (Parton 1993)
    # active_to_slow <<- (activeout * (1.0 - frac_microb_resp - 0.003 -
    #                      0.032 * Claysoil))


    # C flux active pool passive pool
    active_to_passive <<- activeout * 0.004

    # Respiration fluxes
    co2_to_air[5] <<- activeout * frac_microb_resp
}

cfluxes_from_slow_pool <- function() {
    # Send C fluxes from slow pool to other SOM pools

    slowout <- slowsoil * decayrate[6]

    # C flux slow pool active pool
    slow_to_active <<- slowout * 0.42

    # slow pool passive pool
    slow_to_passive <<- slowout * 0.03

    # Respiration fluxes
    co2_to_air[6] <<- slowout * 0.55
}

cfluxes_from_passive_pool <- function() {

    # C flux passive pool  active pool
    passive_to_active <<- passivesoil * decayrate[7] * 0.45

    # Respiration fluxes
    co2_to_air[7] <<- passivesoil * decayrate[7] * 0.55
}

calculate_soil_respiration <- function() {
    # calculate the total soil respiration (heterotrophic) flux, i.e.
    # the amount of CO2 released back to the atmosphere

    # total CO2 production
    hetero_resp <<- (co2_to_air[1] + co2_to_air[2] +
                      co2_to_air[3] + co2_to_air[4] + co2_to_air[5] +
                      co2_to_air[6] + co2_to_air[7])

    # insert following line so value of respiration obeys c conservation if
    # assuming a fixed passive pool
    if (passiveconst == TRUE) {
        hetero_resp <<- (hetero_resp + active_to_passive +
                          slow_to_passive - passivesoil *
                          decayrate[7])
    }
}

precision_control_soil_c <- function() {
    # Detect very low values in state variables and force to zero to
    # avoid rounding and overflow errors

    tolerance <- 1E-08

    # C & N state variables
    if (metabsurf < tolerance) {
        excess <- metabsurf
        surf_metab_to_active <<- excess * 0.45
        co2_to_air[3] <<- excess * 0.55
        metabsurf <<- 0.0
    }

    if (metabsoil < tolerance) {
        excess <- metabsoil
        soil_metab_to_active <<- excess * 0.45
        co2_to_air[4] <<- excess * 0.55
        metabsoil <<- 0.0
    }
}

calculate_cpools <- function() {
    # Calculate new soil carbon pools.

    # Update pools
    structsurf <<- structsurf + (surf_struct_litter -
                     (surf_struct_to_slow + surf_struct_to_active +
                      co2_to_air[1]))

    structsoil <<- structsoil + (soil_struct_litter -
                     (soil_struct_to_slow + soil_struct_to_active +
                      co2_to_air[2]))

    metabsurf <<- metabsurf + (surf_metab_litter -
                     (surf_metab_to_active + co2_to_air[3]))

    metabsoil <<- metabsoil + (soil_metab_litter -
                     (soil_metab_to_active + co2_to_air[4]))

    # store the C SOM fluxes for Nitrogen/Phosphorus calculations
    c_into_active <<- (surf_struct_to_active + soil_struct_to_active +
                        surf_metab_to_active + soil_metab_to_active +
                        slow_to_active + passive_to_active)

    c_into_slow <<- (surf_struct_to_slow + soil_struct_to_slow +
                      active_to_slow)

    c_into_passive <<- active_to_passive + slow_to_passive

    activesoil <<- activesoil + (c_into_active -
                      (active_to_slow + active_to_passive +
                       co2_to_air[5]))

    slowsoil <<- slowsoil + (c_into_slow -
                    (slow_to_active + slow_to_passive +
                     co2_to_air[6]))

    passivesoil <<- passivesoil + (c_into_passive -
                        (passive_to_active + co2_to_air[7]))

    # When nothing is being added to the metabolic pools, there is the
    # potential scenario with the way the model works for tiny bits to be
    # removed with each timestep. Effectively with time this value which is
    # zero can end up becoming zero but to a silly decimal place

    precision_control_soil_c()
}

calc_root_exudation_uptake_of_C <- function() {
    # The amount of C which enters the active pool varies according to the
    # CUE of SOM in response to root exudation (REXCUE). REXCUE determines
    # the fraction of REXC that enters the active pool as C. The remaining
    # flux is respired.


    # REXCUE determines which fraction of REXC enters the active pool as C
    # (delta_Cact). The remaining fraction of REXC is respired as CO2.

    active_CN <- activesoil / activesoiln

    if (root_exu_CUE < -0.5) {
        # flexible cue
        #  - The constraint of 0.3<=REXCUE<=0.6 is based on observations of
        #    the physical limits of microbes

        if (float_eq(root_exc, 0.0)) {
            rex_NC <- 0.0
        } else {
            rex_NC <- root_exn / root_exc
        }
        rexc_cue <<- max(0.3, min(0.6, rex_NC * active_CN))
    } else {
        rexc_cue <<- root_exu_CUE
    }

    C_to_active_pool <- root_exc * rexc_cue
    activesoil <<- activesoil + C_to_active_pool

    # Update respiration fluxes.


    # CUE of microbial rhizodeposition uptake is constant, so the fraction
    # of the rhizodeposition will be used for immediate respiration

    co2_released_exud <<- (1.0 - rexc_cue) * root_exc
    hetero_resp <<- hetero_resp + co2_released_exud
}

calculate_csoil_flows <- function(tsoil, doy) {
    # Fraction of C lost due to microbial respiration
    frac_microb_resp <- 0.85 - (0.68 * finesoil)

    # need to store grazing flag. Allows us to switch on the annual
    # grazing event, but turn it off for every other day of the year.
    cntrl_grazing <- grazing
    if (grazing == 2 && disturbance_doy == doy+1) {
        grazing <<- TRUE
    }

    tfac_soil_decomp <<- calc_soil_temp_factor(tsoil)

    # calculate model decay rates
    calculate_decay_rates()

    # plant litter inputs to the metabolic and structural pools determined
    # by ratio of lignin/N ratio

    lnleaf <- calc_ligin_nratio_leaves()
    lnroot <- calc_ligin_nratio_fine_roots()
    fmleaf <- metafract(lnleaf)
    fmroot <- metafract(lnroot)

    # input from faeces
    flux_from_grazers()
    partition_plant_litter()
    cfluxes_from_structural_pool()
    cfluxes_from_metabolic_pool()
    cfluxes_from_active_pool(frac_microb_resp)
    cfluxes_from_slow_pool()
    cfluxes_from_passive_pool()
    calculate_soil_respiration()

    # update the C pools
    calculate_cpools()

    # calculate NEP
    nep <<- npp - hetero_resp - ceaten * (1.0 - fracfaeces)

    # save fluxes for NCEAS output
    co2_rel_from_surf_struct_litter <<- co2_to_air[1]
    co2_rel_from_soil_struct_litter <<- co2_to_air[2]
    co2_rel_from_surf_metab_litter <<- co2_to_air[3]
    co2_rel_from_soil_metab_litter <<- co2_to_air[4]
    co2_rel_from_active_pool <<- co2_to_air[5]
    co2_rel_from_slow_pool <<- co2_to_air[6]
    co2_rel_from_passive_pool <<- co2_to_air[7]

    # switch off grazing if this was just activated as an annual event
    grazing <<- cntrl_grazing

    if (exudation && alloc_model != GRASSES) {
        calc_root_exudation_uptake_of_C()
    }
}

grazer_inputs <- function() {
    # Grazer inputs from faeces and urine, flux detd by faeces c:n
    if (grazing) {
        faecesn <<- faecesc / faecescn
    } else {
        faecesn <<- 0.0
    }

    # make sure faecesn <= total n input to soil from grazing
    arg <- neaten * fractosoil
    if (faecesn > arg)
        faecesn <<- neaten * fractosoil

    # urine=total-faeces
    if (grazing) {
        nurine <<- neaten * fractosoil - faecesn
    } else {
        nurine <<- 0.0
    }

    if (nurine < 0.0)
        nurine <<- 0.0
}

n_inputs_from_plant_litter <- function() {
    # inputs from plant litter.
    #
    # surface and soil pools are independent. Structural input flux n:c can
    # be either constant or a fixed fraction of metabolic input flux.
    #
    # Returns:
    # --------
    # nsurf : float
    #     N input from surface pool
    # nsoil : float
    #     N input from soil pool

    # surface and soil inputs (faeces n goes to abovgrd litter pools)
    nsurf <- deadleafn + deadbranchn + deadstemn + faecesn
    nsoil <- deadrootn + deadcrootn

    return (c(nsurf, nsoil))
}

partition_plant_litter_n <- function(nsurf, nsoil) {
    # Partition litter N from the plant (surface) and roots into metabolic
    # and structural pools
    #
    # Parameters:
    # -----------
    # nsurf : float
    #     N input from surface pool
    # nsoil : float
    #     N input from soil pool

    # constant structural input n:c as per century
    if (strfloat == 1) {

        # structural input n:c is a fraction of metabolic
        c_surf_struct_litter <- (surf_struct_litter * structrat +
                                surf_metab_litter)

        if (float_eq(c_surf_struct_litter, 0.0)) {
             n_surf_struct_litter <<- 0.0
        } else {
             n_surf_struct_litter <<- (nsurf * surf_struct_litter *
                                        structrat / c_surf_struct_litter)
        }


        c_soil_struct_litter <- (soil_struct_litter * structrat +
                                soil_metab_litter)

        if (float_eq(c_soil_struct_litter, 0.0)) {
            n_soil_struct_litter <<- 0.0
        } else {
            n_soil_struct_litter <<- (nsurf * soil_struct_litter *
                                       structrat / c_soil_struct_litter)
        }

    } else {

        # n flux surface structural pool
        n_surf_struct_litter <<- surf_struct_litter / structcn

        # n flux soil structural pool
        n_soil_struct_litter <<- soil_struct_litter / structcn

        # if not enough N for structural, all available N goes to structural
        if (n_surf_struct_litter > nsurf)
             n_surf_struct_litter <<- nsurf
        if (n_soil_struct_litter > nsoil)
            n_soil_struct_litter <<- nsoil
    }


    # remaining N goes to metabolic pools
    n_surf_metab_litter <<- nsurf - n_surf_struct_litter
    n_soil_metab_litter <<- nsoil - n_soil_struct_litter
}

nfluxes_from_structural_pools <- function() {
    # from structural pool
    structout_surf <- structsurfn * decayrate[1]
    structout_soil <- structsoiln * decayrate[3]

    sigwt <- structout_surf / (ligshoot * 0.7 + (1.0 - ligshoot) * 0.55)

    # N flux from surface structural pool slow pool
    n_surf_struct_to_slow <<- sigwt * ligshoot * 0.7

    # N flux surface structural pool active pool
    n_surf_struct_to_active <<- sigwt * (1.0 - ligshoot) * 0.55

    sigwt <- structout_soil / (ligroot * 0.7 + (1. - ligroot) * 0.45)


    # N flux from soil structural pool slow pool
    n_soil_struct_to_slow <<- sigwt * ligroot * 0.7

    # N flux from soil structural pool active pool
    n_soil_struct_to_active <<- sigwt * (1.0 - ligroot) * 0.45
}

nfluxes_from_metabolic_pool <- function() {

    # N flux surface metabolic pool active pool
    n_surf_metab_to_active <<- metabsurfn * decayrate[2]

    # N flux soil metabolic pool  active pool
    n_soil_metab_to_active <<- metabsoiln * decayrate[4]
}

nfluxes_from_active_pool <- function(frac_microb_resp) {
    # N fluxes from active pool
    activeout <- activesoiln * decayrate[5]
    sigwt <- activeout / (1.0 - frac_microb_resp)

    # N flux active pool slow pool
    n_active_to_slow <<- sigwt * (1.0 - frac_microb_resp - 0.004)

    # N flux active pool passive pool
    n_active_to_passive <<- sigwt * 0.004
}

nfluxes_from_slow_pool <- function() {
    # N fluxes from slow pools

    slowout <- slowsoiln * decayrate[6]
    sigwt <- slowout / 0.45

    # C flux slow pool active pool
    n_slow_to_active <<- sigwt * 0.42

    # slow pool passive pool
    n_slow_to_passive <<- sigwt * 0.03
}

nfluxes_from_passive_pool <- function() {
    # N fluxes from passive pool

    # C flux passive pool active pool
    n_passive_to_active <<- passivesoiln * decayrate[7]
}

calculate_n_mineralisation <- function() {
    # N gross mineralisation rate is given by the excess of N outflows
    # over inflows. Nitrogen mineralisation is the process by which organic
    # N is converted to plant available inorganic N, i.e. microbes decompose
    # organic N from organic matter to ammonia (NH3) and ammonium (NH4),
    # called ammonification.
    #
    # Returns:
    # --------
    # value : float
    #     Gross N mineralisation

    ngross <<- (n_surf_struct_to_slow + n_surf_struct_to_active +
                  n_soil_struct_to_slow + n_soil_struct_to_active +
                  n_surf_metab_to_active + n_soil_metab_to_active +
                  n_active_to_slow + n_active_to_passive +
                  n_slow_to_active + n_slow_to_passive +
                  n_passive_to_active)
}

calculate_n_immobilisation <- function() {
    # N immobilised in new soil organic matter, the reverse of
    # mineralisation. Micro-organisms in the soil compete with plants for N.
    # Immobilisation is the process by which nitrate and ammonium are taken up
    # by the soil organisms and thus become unavailable to the plant
    # organic N).
    #
    # When C:N ratio is high the microorganisms need more nitrogen from
    # the soil to decompose the carbon in organic materials. This N will be
    # immobilised until these microorganisms die and the nitrogen is
    # released.
    #
    # General equation for new soil N:C ratio vs Nmin, expressed as linear
    # equation passing through point Nmin0, actncmin (etc). Values can be
    # Nmin0=0, Actnc0=Actncmin
    #
    # if Nmin < Nmincrit:
    #     New soil N:C = soil N:C (when Nmin=0) + slope * Nmin
    #
    # if Nmin > Nmincrit
    #     New soil N:C = max soil N:C
    #
    # NB N:C ratio of new passive SOM can change even if assume Passiveconst
    #
    # Returns:
    # --------
    # nimob : float
    #     N immobilsed

    # N:C new SOM - active, slow and passive
    active_nc_slope <- calculate_nc_slope(actncmax, actncmin)
    slow_nc_slope <- calculate_nc_slope(slowncmax, slowncmin)
    passive_nc_slope <- calculate_nc_slope(passncmax, passncmin)

    # convert units
    nmin <- nmin0 / M2_AS_HA * G_AS_TONNES

    arg1 <- (passncmin - passive_nc_slope * nmin) * c_into_passive
    arg2 <- (slowncmin - slow_nc_slope * nmin) * c_into_slow
    arg3 <- c_into_active * (actncmin - active_nc_slope * nmin)
    numer1 <- arg1 + arg2 + arg3

    arg1 <- c_into_passive * passncmax
    arg2 <- c_into_slow * slowncmax
    arg3 <- c_into_active * actncmax
    numer2 <- arg1 + arg2 + arg3

    arg1 <- c_into_passive * passive_nc_slope
    arg2 <- c_into_slow * slow_nc_slope
    arg3 <- c_into_active * active_nc_slope
    denom <- arg1 + arg2 + arg3

    # evaluate N immobilisation in new SOM
    nimmob <- numer1 + denom * inorgn
    if (nimmob > numer2)
        nimmob <- numer2

    return (c(nimmob, active_nc_slope, slow_nc_slope, passive_nc_slope))
}

calc_n_net_mineralisation <- function() {
    # N Net mineralisation from microbial activity
    nmineralisation <<- ngross - nimmob + nlittrelease
}

calc_root_exudation_uptake_of_N <- function() {
    # When N mineralisation is large enough to allow a small amount of N
    # immobilisation, the amount of N which enters the active pool is
    # calculated according to REXC divided by the CN of the active pool. When
    # exudation enters the active pool, the CN ratio of the exudates drops
    # from REXC/REXN to the CN of the active pool. Which is consistent with
    # the CENTURY framework, where C flows between pools lead to either
    # mineralisation (N gain) or immobilisation (N loss) due to differences
    # in the CN ratio of the outgoing and incoming pools.
    #
    # The amount of N added to the active pool is independent of the CUE of
    # the microbial pool in response to root exudation (REXCUE).

    N_available <- inorgn + (ninflow + nurine + nmineralisation -
                               nloss - nuptake)

    active_NC <- activesoiln / activesoil
    delta_Nact <- root_exc * rexc_cue * active_NC

    # Demand for N from exudation to meet the C:N ratio of the active pool,
    # given the amount of N you add.

    N_miss <- delta_Nact - root_exn

    if (N_miss <= 0.0) {

        # Root exudation includes more N than is needed by the microbes, the
        # excess is mineralised

        nmineralisation <<- nmineralisation - N_miss
        N_to_active_pool <- root_exn + N_miss
    } else {

        # Not enough N in the soil to meet demand, so we are providing all
        # the N we have, which means that the C:N ratio of the active pool
        # changes.

        if (N_miss > N_available) {
            N_to_active_pool <- root_exn + N_available
            nmineralisation <<- nmineralisation - N_available
        } else {

            # Enough N to meet demand, so takes N from the mineralisation
            # and the active pool maintains the same C:N ratio.

            N_to_active_pool <- root_exn + N_miss
            nmineralisation <<- nmineralisation - N_miss
        }
    }

    # update active pool
    activesoiln <<- activesoiln + N_to_active_pool
}

adjust_residence_time_of_slow_pool <- function() {
    # Priming simulations the residence time of the slow pool is flexible,
    # as the flux out of the active pool (factive) increases the residence
    # time of the slow pool decreases.

    # total flux out of the factive pool
    factive <<- (active_to_slow + active_to_passive +
                  co2_to_air[5] + co2_released_exud)

    if (float_eq(factive, 0.0)) {
        # Need to correct units of rate constant
        rt_slow_pool <- 1.0 / (kdec6 * NDAYS_IN_YR)
    } else {
        rt_slow_pool <- (1.0 / prime_y) /
                        max(0.3, (factive / (factive + prime_z)))

        # GDAY uses decay rates rather than residence times...
        kdec6 <<- 1.0 / rt_slow_pool

        # rate constant needs to be per day inside GDAY
        kdec6 <<- kdec6 / NDAYS_IN_YR

    }

    # Save for outputting purposes only
    rtslow <<- rt_slow_pool
}

precision_control_soil_n <- function() {
    # Detect very low values in state variables and force to zero to
    # avoid rounding and overflow errors

    tolerance <- 1E-08

    if (metabsurfn < tolerance) {
        excess <- metabsurfn
        n_surf_metab_to_active <<- excess
        metabsurfn <<- 0.0
    }

    if (metabsoiln < tolerance) {
        excess <- metabsoiln
        n_soil_metab_to_active <<- excess
        metabsoiln <<- 0.0
    }
}

calculate_npools <- function(active_nc_slope, slow_nc_slope, passive_nc_slope) {
    # Update N pools in the soil
    #
    # Parameters
    # ----------
    # active_nc_slope : float
    #     active NC slope
    # slow_nc_slope: float
    #     slow NC slope
    # passive_nc_slope : float
    #     passive NC slope

    # net N release implied by separation of litter into structural
    # & metabolic. The following pools only fix or release N at their
    # limiting n:c values.

    # N released or fixed from the N inorganic pool is incremented with
    # each call to nc_limit and stored in nlittrelease
    nlittrelease <<- 0.0

    structsurfn <<- structsurfn + (n_surf_struct_litter -
                        (n_surf_struct_to_slow +
                         n_surf_struct_to_active))

    structsoiln <<- structsoiln + (n_soil_struct_litter -
                       (n_soil_struct_to_slow + n_soil_struct_to_active))

    if (strfloat == 0) {
        structsurfn <<- structsurfn + nc_limit(structsurf, structsurfn,
                                   1.0/structcn, 1.0/structcn)
        structsoiln <<- structsoiln + nc_limit(structsoil, structsoiln,
                                   1.0/structcn, 1.0/structcn)
    }

    metabsurfn <<- metabsurfn + n_surf_metab_litter - n_surf_metab_to_active
    metabsurfn <<- metabsurfn + nc_limit(metabsurf, metabsurfn,1.0/25.0, 1.0/10.0)
    metabsoiln <<- metabsoiln + (n_soil_metab_litter - n_soil_metab_to_active)
    metabsoiln <<- metabsoiln + nc_limit(metabsoil, metabsoiln, 1.0/25.0,
                              1.0/10.0)

    # When nothing is being added to the metabolic pools, there is the
    # potential scenario with the way the model works for tiny bits to be
    # removed with each timestep. Effectively with time this value which is
    # zero can end up becoming zero but to a silly decimal place
    precision_control_soil_n()

    # Update SOM pools
    n_into_active <- (n_surf_struct_to_active + n_soil_struct_to_active +
                     n_surf_metab_to_active + n_soil_metab_to_active +
                     n_slow_to_active + n_passive_to_active)

    n_out_of_active <- n_active_to_slow + n_active_to_passive

    n_into_slow <- (n_surf_struct_to_slow + n_soil_struct_to_slow +
                   n_active_to_slow)

    n_out_of_slow <- n_slow_to_active + n_slow_to_passive
    n_into_passive <- n_active_to_passive + n_slow_to_passive
    n_out_of_passive <- n_passive_to_active

    # N:C of the SOM pools increases linearly btw prescribed min and max
    # values as the Nconc of the soil increases.
    arg <- inorgn - nmin0 / M2_AS_HA * G_AS_TONNES

    # active
    active_nc <- actncmin + active_nc_slope * arg
    if (active_nc > actncmax)
        active_nc <- actncmax

    # release N to Inorganic pool or fix N from the Inorganic pool in order
    # to normalise the N:C ratio of a net flux
    fixn <- nc_flux(c_into_active, n_into_active, active_nc)
    activesoiln <<- activesoiln + n_into_active + fixn - n_out_of_active

    # slow
    slow_nc <- slowncmin + slow_nc_slope * arg
    if (slow_nc > slowncmax)
        slow_nc <- slowncmax

    # release N to Inorganic pool or fix N from the Inorganic pool in order
    # to normalise the N:C ratio of a net flux
    fixn <- nc_flux(c_into_slow, n_into_slow, slow_nc)
    slowsoiln <<- slowsoiln + n_into_slow + fixn - n_out_of_slow

    # passive, update passive pool only if passiveconst=0
    pass_nc <- passncmin + passive_nc_slope * arg
    if (pass_nc > passncmax)
        pass_nc <- passncmax

    # release N to Inorganic pool or fix N from the Inorganic pool in order
    # to normalise the N:C ratio of a net flux
    fixn <- nc_flux(c_into_passive, n_into_passive, pass_nc)
    passivesoiln <<- passivesoiln + n_into_passive + fixn - n_out_of_passive

    # Daily increment of soil inorganic N pool, diff btw in and effluxes
    # (grazer urine n goes directly into inorganic pool) nb inorgn may be
    # unstable if rateuptake is large
    inorgn <<- inorgn + (ninflow + nurine + nmineralisation -
                  nloss - nuptake)

    # fprintf(stderr, "inorgn=%f\n", inorgn)

    # nmineralisation = ngross - nimmob + nlittrelease
}

calculate_nsoil_flows <- function(doy) {

    # need to store grazing flag. Allows us to switch on the annual
    # grazing event, but turn it off for every other day of the year.
    cntrl_grazing <- grazing
    if (grazing == 2 && disturbance_doy == doy + 1) {
        grazing <<- TRUE
    }

    # Fraction of C lost due to microbial respiration
    frac_microb_resp <- 0.85 - (0.68 * finesoil)

    grazer_inputs()
    res__ <- n_inputs_from_plant_litter()
    nsurf <- res__[1]
    nsoil <- res__[2]
    partition_plant_litter_n(nsurf, nsoil)

    # SOM nitrogen effluxes.  These are assumed to have the source n:c
    # ratio prior to the increase of N:C due to co2 evolution.
    nfluxes_from_structural_pools()
    nfluxes_from_metabolic_pool()
    nfluxes_from_active_pool(frac_microb_resp)
    nfluxes_from_slow_pool()
    nfluxes_from_passive_pool()

    # gross N mineralisation
    calculate_n_mineralisation()

    # calculate N immobilisation
    res__ <- calculate_n_immobilisation()
    nimmob <<- res__[1]
    active_nc_slope <- res__[2]
    slow_nc_slope <- res__[3]
    passive_nc_slope <- res__[4]

    # calculate N net mineralisation
    calc_n_net_mineralisation()

    if (exudation && alloc_model != GRASSES) {
        calc_root_exudation_uptake_of_N()
    }

    if (adjust_rtslow) {
        adjust_residence_time_of_slow_pool()
    } else {
        # Need to correct units of rate constant
        rtslow <<- 1.0 / (kdec6 * NDAYS_IN_YR)
    }

    # Update model soil N pools
    calculate_npools(active_nc_slope, slow_nc_slope, passive_nc_slope)

    # switch off grazing if this was just activated as an annual event
    grazing <<- cntrl_grazing
}

grazer_inputs_p <- function() {
    # Grazer inputs from faeces and urine, flux detd by faeces c:p

    if (grazing) {
        faecesp <<- faecesc / faecescp
    } else {
        faecesp <<- 0.0
    }

    # make sure faecesp <= total p input to soil from grazing
    arg <- peaten * fractosoilp
    if (faecesp > arg)
        faecesp <<- peaten * fractosoilp

    # urine=total-faeces
    if (grazing)
        purine <<- peaten * fractosoilp - faecesp
    else
        purine <<- 0.0

    if (purine < 0.0)
        purine <<- 0.0
}

p_inputs_from_plant_litter <- function() {
    # inputs from plant litter.
    #
    # surface and soil pools are independent. Structural input flux p:c can
    # be either constant or a fixed fraction of metabolic input flux.
    #
    # Returns:
    # --------
    # psurf : float
    # P input from surface pool
    # psoil : float
    # P input from soil pool

    # surface and soil inputs (faeces p goes to abovgrd litter pools)
    psurf <- deadleafp + deadbranchp + deadstemp + faecesp
    psoil <- deadrootp + deadcrootp

    return (c(psurf, psoil))
}

partition_plant_litter_p <- function(psurf, psoil) {
    # Partition litter P from the plant (surface) and roots into metabolic
    # and structural pools
    #
    # Parameters:
    # -----------
    # psurf : float
    # P input from surface pool
    # psoil : float
    # P input from soil pool

    # constant structural input p:c as per century
    if (strpfloat == 1) {
        # structural input p:c is a fraction of metabolic
        c_surf_struct_litter <- (surf_struct_litter * structratp +
                                surf_metab_litter)

        if (float_eq(c_surf_struct_litter, 0.0)) {
            p_surf_struct_litter <<- 0.0
        } else {
            p_surf_struct_litter <<- (psurf * surf_struct_litter *
                                    structratp / c_surf_struct_litter)
        }


        c_soil_struct_litter <- (soil_struct_litter * structratp +
                                soil_metab_litter)

        if (float_eq(c_soil_struct_litter, 0.0)) {
            p_soil_struct_litter <<- 0.0
        } else {
            p_soil_struct_litter <<- (psurf * soil_struct_litter *
                                    structratp / c_soil_struct_litter)
        }

    } else {

        # p flux surface structural pool
        p_surf_struct_litter <<- surf_struct_litter / structcp

        # p flux soil structural pool
        p_soil_struct_litter <<- soil_struct_litter / structcp

        # if not enough P for structural, all available P goes to structural
        if (p_surf_struct_litter > psurf)
            p_surf_struct_litter <<- psurf
        if (p_soil_struct_litter > psoil)
            p_soil_struct_litter <<- psoil
    }

    # remaining P goes to metabolic pools
    p_surf_metab_litter <<- psurf - p_surf_struct_litter
    p_soil_metab_litter <<- psoil - p_soil_struct_litter
}

calculate_p_immobilisation <- function() {
    # P immobilised in new soil organic matter, the reverse of
    # mineralisation. Micro-organisms in the soil compete with plants for P.
    # Immobilisation occurs when plant available P forms are consumed by microbes, turning
    # the P into organic P forms that are not available to plants.
    #
    # When C:P ratio is high the microorganisms need more P from
    # the soil to decompose the carbon in organic materials. This P will be
    # immobilised until these microorganisms die and the P is
    # released.
    #
    # General equation for new soil P:C ratio vs Pmin, expressed as linear
    # equation passing through point Pmin0, actpcmin (etc). Values can be
    # Pmin0=0, Actpc0=Actpcmin
    #
    # if Pmin < Pmincrit:
    # New soil P:C = soil P:C (when Pmin=0) + slope * Pmin
    #
    # if Pmin > Pmincrit
    # New soil P:C = max soil P:C
    #
    # NB P:C ratio of new passive SOM can change even if assume Passiveconst
    #
    # Returns:
    # --------
    # pimob : float
    # P immobilsed

    # P:C new SOM - active, slow and passive
    active_pc_slope <- calculate_pc_slope(actpcmax, actpcmin)
    slow_pc_slope <- calculate_pc_slope(slowpcmax, slowpcmin)
    passive_pc_slope <- calculate_pc_slope(passpcmax, passpcmin)

    # convert units
    pmin <- pmin0 / M2_AS_HA * G_AS_TONNES

    arg1 <- (passpcmin - passive_pc_slope * pmin) * c_into_passive
    arg2 <- (slowpcmin - slow_pc_slope * pmin) * c_into_slow
    arg3 <- c_into_active * (actpcmin - active_pc_slope * pmin)
    numer1 <- arg1 + arg2 + arg3

    arg1 <- c_into_passive * passpcmax
    arg2 <- c_into_slow * slowpcmax
    arg3 <- c_into_active * actpcmax
    numer2 <- arg1 + arg2 + arg3

    arg1 <- c_into_passive * passive_pc_slope
    arg2 <- c_into_slow * slow_pc_slope
    arg3 <- c_into_active * active_pc_slope
    denom <- arg1 + arg2 + arg3

    # evaluate P immobilisation in new SOM
    pimmob <- numer1 + denom * inorglabp
    if (pimmob > numer2)
        pimmob <- numer2

    return (c(pimmob, active_pc_slope, slow_pc_slope, passive_pc_slope))
}

pfluxes_from_structural_pools <- function() {
    # from structural pool
    structout_surf <- structsurfp * decayrate[1]
    structout_soil <- structsoilp * decayrate[3]

    sigwt <- structout_surf / (ligshoot * 0.7 + (1.0 - ligshoot) * 0.55)

    # P flux from surface structural pool slow pool
    p_surf_struct_to_slow <<- sigwt * ligshoot * 0.7

    # P flux surface structural pool active pool
    p_surf_struct_to_active <<- sigwt * (1.0 - ligshoot) * 0.55

    sigwt <- structout_soil / (ligroot * 0.7 + (1. - ligroot) * 0.45)


    # P flux from soil structural pool slow pool
    p_soil_struct_to_slow <<- sigwt * ligroot * 0.7

    # N flux from soil structural pool active pool
    p_soil_struct_to_active <<- sigwt * (1.0 - ligroot) * 0.45
}

pfluxes_from_metabolic_pool <- function() {

    # P flux surface metabolic pool active pool
    p_surf_metab_to_active <<- metabsurfp * decayrate[2]

    # P flux soil metabolic pool  active pool
    p_soil_metab_to_active <<- metabsoilp * decayrate[4]
}


pfluxes_from_active_pool <- function(frac_microb_resp) {
    # P fluxes from active pool
    activeout <- activesoilp * decayrate[5]
    sigwt <- activeout / (1.0 - frac_microb_resp)

    # P flux active pool slow pool
    p_active_to_slow <<- sigwt * (1.0 - frac_microb_resp - 0.004)

    # P flux active pool passive pool
    p_active_to_passive <<- sigwt * 0.004
}

pfluxes_from_slow_pool <- function() {
    # P fluxes from slow pools

    slowout <- slowsoilp * decayrate[6]
    sigwt <- slowout / 0.45

    # fprintf(stderr, "slowsoilp %f\n", slowsoilp)

    # P flux slow pool active pool
    p_slow_to_active <<- sigwt * 0.42

    # slow pool passive pool
    p_slow_to_passive <<- sigwt * 0.03
}

pfluxes_from_passive_pool <- function() {
    # P fluxes from passive pool

    # P flux passive pool active pool
    p_passive_to_active <<- passivesoilp * decayrate[7]
}

calculate_p_parent_fluxes <- function() {
    # Calculate weathering of parent P materials, i.e.
    # the fluxes enterring into mineral P pool
    # Fluxes in = out so that parent P pool is a constant pool


    # parent material weathering
    p_par_to_min <<- p_rate_par_weather * inorgparp
}

calculate_p_mineralisation <- function() {
    # P gross mineralisation rate is given by the excess of P outflows
    # over inflows. P mineralisation is the process by which organic P is
    # converted to plant available inorganic P, i.e. the microbial conversion
    # of organic P to H2PO4- or HPO42- forms of plant available P known as
    # orthophosphates.
    #
    # Returns:
    # --------
    # value : float
    # Gross P mineralisation
    # Unit: t/ha/d

    pgross <<- (p_surf_struct_to_slow + p_surf_struct_to_active +
                  p_soil_struct_to_slow + p_soil_struct_to_active +
                  p_surf_metab_to_active + p_soil_metab_to_active +
                  p_active_to_slow + p_active_to_passive +
                  p_slow_to_active + p_slow_to_passive +
                  p_passive_to_active)
}

calc_p_net_mineralisation <- function() {
    # P Net mineralisation from microbial activity,
    # excluding the (- p_sorb_to_ssorb + p_ssorb_to_sorb activity)
    pmineralisation <<- pgross - pimmob + plittrelease
}

calculate_p_biochemical_mineralisation <- function() {
    # Calculate biochemical P mineralisation rate for slow SOM pool only
    #
    # Ref: Wang et al., 2007, GB1018.
    #
    # Parameters
    # ----------
    # n_cost_of_p: float
    # N cost of plant root P uptake [g N (g P)-1]
    #
    # crit_n_cost_of_p: float
    # The critical value of N cost of root P uptake above which
    # phosphatase production starts [g N (g P)-1]
    #
    # max_p_biochemical: float
    # Maximum rate of biochemical P mineralisation [g P m-2 y-1]
    #
    # biochemical_p_constant: float
    # Michaelis-Menten constant for biochemical P mineralisation [g N (g P)-1]
    #
    # c_gain_of_p: float
    # Carbon per unit of P uptake (NPP/Puptake) [g C (g P)-1]
    #
    # c_gain_of_n: float
    # Carbon per unit of N uptake (NPP/Nuptake) [g C (g N)-1]
    #
    # Return
    # ----------
    # p_slow_biochemical: float
    # biochemical P mineralisation rate in slow SOM pool
    # Unit [t/ha]


    # Calculate C gain per unit P
    if (puptake > 0.0) {
        c_gain_of_p <- npp / puptake
    } else {
        c_gain_of_p <- 0.0
    }

    # Calculate C gain per unit N
    if (nuptake > 0.0) {
        c_gain_of_n <- npp / nuptake
    } else {
        c_gain_of_n <- 0.0
    }

    # Calculate n cost of p
    if (c_gain_of_n > 0.0) {
        n_cost_of_p <- c_gain_of_p / c_gain_of_n
    } else {
        n_cost_of_p <- 0.0
    }

    # Phosphatase production start when N cost of P is greater than critical
    # N cost of P value, and when carbon uptake is more strongly limited by
    # the uptake of P than by the uptake of N
    if (c_gain_of_p > c_gain_of_n && n_cost_of_p > crit_n_cost_of_p) {
        p_slow_biochemical <<- max_p_biochemical *
                                ((n_cost_of_p - crit_n_cost_of_p) /
                                (n_cost_of_p - crit_n_cost_of_p +
                                biochemical_p_constant))
    } else {
        p_slow_biochemical <<- 0.0
    }
}

calculate_p_ssorb_to_sorb <- function() {
    # calculate P transfer from strongly sorbed P pool to
    # sorbed P pool
    #
    # Parameters
    # ----------
    # phtextint: float
    #     intercept for the texture equation of strongly sorbed P depends upon
    #     pH input
    #
    # phtextslope: float
    #     slope value used in determining effect of sand content on ssorb P
    #     flow to mineral P
    #
    # psecmn: float
    #     controls the flow from secondary to mineral P
    #
    # Returns:
    # ----------
    # p_ssorb_to_min: float
    #     flux rate of p strongly sorbed pool to p sorbed pool

    cntrl_text_p <- text_effect_p

    if (cntrl_text_p == 1) {

        dely <- phtextmax - phtextmin
        delx <- phmax - phmin

        xslope <- dely/delx
        yint <- phtextmin - xslope * phmin

        if (soilph < phmin) {
            phtextint <- phtextmin
        } else if (soilph > phmax) {
            phtextint <- phtextmax
        } else {
            phtextint <- xslope * soilph + yint
        }

        p_ssorb_to_min <<- max(0.0, (phtextint + phtextslope *
                                     (1.0 - finesoil)) * inorgssorbp)

    } else {
        if (inorgssorbp > 0.0) {
            p_ssorb_to_min <<- psecmnp * inorgssorbp
        } else {
            p_ssorb_to_min <<- 0.0
        }

    }
}

calculate_p_sorb_to_ssorb <- function() {

    # P flux from sorbed pool to strongly sorbed P pool
    if (inorgsorbp > 0.0) {
        p_min_to_ssorb <<- rate_sorb_ssorb * inorgsorbp
    } else {
        p_min_to_ssorb <<- 0.0
    }
}

calculate_p_ssorb_to_occ <- function() {
    # P flux from strongly sorbed pool to occluded P pool
    if (inorgssorbp > 0.0) {
        p_ssorb_to_occ <<- rate_ssorb_occ * inorgssorbp
    } else {
        p_ssorb_to_occ <<- 0.0
    }
}

soil_sorption_parameters <- function(soil_order) {
    # Parameterize Smax and Ks parameters based on soil order
    # Ref. Yang et al., 2016, GRL, Table S2

    if (soil_order == "andisol") {
        smax <<- 10
    } else if (soil_order == "gelisol") {
        smax <<- 5
    } else if (soil_order == "histosol") {
        smax <<- 5
    } else if (soil_order == "entisol") {
        smax <<- 5
    } else if (soil_order == "inceptisol") {
        smax <<- 5
    } else if (soil_order == "aridsol") {
        smax <<- 7
    } else if (soil_order == "vertisol") {
        smax <<- 7
    } else if (soil_order == "mollisol") {
        smax <<- 7
    } else if (soil_order == "alfisol") {
        smax <<- 7
    } else if (soil_order == "spodosol") {
        smax <<- 9
    } else if (soil_order == "ultisol") {
        smax <<- 9
    } else if (soil_order == "oxisol") {
        smax <<- 9
    } else {
        stop("Could not understand soil order")
    }
}

calculate_p_min_fluxes <- function() {
    # Calculate the mineral P fluxes (in and out)
    min_frac_p_available_to_plant <- 0.4
    max_frac_p_available_to_plant <- 0.8
    mineral_n_with_max_p <- 0.02              # Unit [t N ha-1]

    # Note: pmineralisation can be negative, and therefore, during spinup,
    # when sorbp stock is low, the current approach resulted in negative sorbp
    # in some cases, but it does not affect the final equilibrated state, so
    # leave as is until better method is found
    tot_in <- p_par_to_min + pmineralisation +
             purine + p_slow_biochemical +
             p_ssorb_to_min

    if (inorglabp > 0) {
        tot_out <- puptake + ploss + p_min_to_ssorb
    } else {
        puptake <<- 0.0
        ploss <<- 0.0
        p_min_to_ssorb <<- 0.0
        tot_out <- puptake + ploss + p_min_to_ssorb
    }

    # Use soil order to obtain smax and ks values
    soil_sorption_parameters(soil_order)

    # Calculate lab P dynamics
    numer <- smax * ks
    denom1 <- (inorglabp + ks) * (inorglabp + ks)
    p_lab_in <<- tot_in / (1.0 + numer / denom1)
    p_lab_out <<- tot_out / (1.0 + numer / denom1)

    # calculate sorb P dynamics
    denom2 <- (inorglabp + ks) * (inorglabp + ks) + numer
    p_sorb_in <<- tot_in * (numer / denom2)
    p_sorb_out <<- tot_out * (numer / denom2)

    # calculating fraction of labile P available for plant uptake
    p_lab_avail <<- max(min_frac_p_available_to_plant,
                     min(min_frac_p_available_to_plant + inorgn *
                    (max_frac_p_available_to_plant - min_frac_p_available_to_plant) /
                     mineral_n_with_max_p, max_frac_p_available_to_plant))
}

precision_control_soil_p <- function() {
    # Detect very low values in state variables and force to zero to
    # avoid rounding and overflow errors

    tolerance <- 1E-08

    if (metabsurfp < tolerance) {
        excess <- metabsurfp
        p_surf_metab_to_active <<- excess
        metabsurfp <<- 0.0
    }

    if (metabsoilp < tolerance) {
        excess <- metabsoilp
        p_soil_metab_to_active <<- excess
        metabsoilp <<- 0.0
    }
}

calculate_ppools <- function(active_pc_slope, slow_pc_slope, passive_pc_slope) {
    # Update P pools in the soil
    #
    # Parameters
    # ----------
    # active_pc_slope : float
    #     active PC slope
    # slow_pc_slope: float
    #     slow PC slope
    # passive_pc_slope : float
    #     passive PC slope

    # net P release implied by separation of litter into structural
    # & metabolic. The following pools only fix or release P at their
    # limiting p:c values.

    # P released or fixed from the P inorganic labile pool is incremented with
    # each call to pc_limit and stored in plittrelease
    plittrelease <<- 0.0

    structsurfp <<- structsurfp + (p_surf_struct_litter -
                      (p_surf_struct_to_slow +
                       p_surf_struct_to_active))

    structsoilp <<- structsoilp + (p_soil_struct_litter -
                      (p_soil_struct_to_slow +
                       p_soil_struct_to_active))

    if (strpfloat == 0) {
        structsurfp <<- structsurfp + pc_limit(structsurf, structsurfp,
                               1.0/structcp, 1.0/structcp)
        structsoilp <<- structsoilp + pc_limit(structsoil, structsoilp,
                               1.0/structcp, 1.0/structcp)
    }

    # pcmin & pcmax from Parton 1989 fig 2
    metabsurfp <<- metabsurfp + p_surf_metab_litter - p_surf_metab_to_active
    metabsurfp <<- metabsurfp + pc_limit(metabsurf, metabsurfp,
                              1.0/150.0, 1.0/80.0)

    # pcmin & pcmax from Parton 1989 fig 2
    metabsoilp <<- metabsoilp + (p_soil_metab_litter - p_soil_metab_to_active)
    metabsoilp <<- metabsoilp + pc_limit(metabsoil, metabsoilp,
                              1.0/150.0, 1.0/80.0)


    # When nothing is being added to the metabolic pools, there is the
    # potential scenario with the way the model works for tiny bits to be
    # removed with each timestep. Effectively with time this value which is
    # zero can end up becoming zero but to a silly decimal place
    precision_control_soil_p()

    # Update SOM pools
    p_into_active <- (p_surf_struct_to_active + p_soil_struct_to_active +
    p_surf_metab_to_active + p_soil_metab_to_active +
    p_slow_to_active + p_passive_to_active)

    p_out_of_active <- p_active_to_slow + p_active_to_passive

    p_into_slow <- (p_surf_struct_to_slow + p_soil_struct_to_slow +
                   p_active_to_slow)

    p_out_of_slow <- p_slow_to_active + p_slow_to_passive + p_slow_biochemical
    p_into_passive <- p_active_to_passive + p_slow_to_passive
    p_out_of_passive <- p_passive_to_active

    # P:C of the SOM pools increases linearly btw prescribed min and max
    # values as the Pconc of the soil increases.
    arg <- inorglabp - pmin0 / M2_AS_HA * G_AS_TONNES

    # active
    active_pc <- actpcmin + active_pc_slope * arg
    if (active_pc > actpcmax)
        active_pc <- actpcmax

    # release P to Inorganic labile pool or fix P from the Inorganic pool in order
    # to normalise the P:C ratio of a net flux
    fixp <- pc_flux(c_into_active, p_into_active, active_pc)
    activesoilp <<- activesoilp + p_into_active + fixp - p_out_of_active

    # slow
    slow_pc <- slowpcmin + slow_pc_slope * arg
    if (slow_pc > slowpcmax)
        slow_pc <- slowpcmax

    # release P to Inorganic pool or fix P from the Inorganic pool in order
    # to normalise the P:C ratio of a net flux
    fixp <- pc_flux(c_into_slow, p_into_slow, slow_pc)
    slowsoilp <<- slowsoilp + p_into_slow + fixp - p_out_of_slow

    # passive, update passive pool only if passiveconst=0
    pass_pc <- passpcmin + passive_pc_slope * arg
    if (pass_pc > passpcmax)
        pass_pc <- passpcmax

    # release P to Inorganic pool or fix P from the Inorganic pool in order
    # to normalise the P:C ratio of a net flux
    # fixp <- pc_flux(c_into_passive, p_into_passive, pass_pc)
    passivesoilp <<- passivesoilp + p_into_passive + fixp - p_out_of_passive

    #fprintf(stderr, "inorglabp 1 %f\n", inorglabp)

    # Daily increment of soil inorganic labile and sorbed P pool
    inorglabp <<- inorglabp + p_lab_in - p_lab_out
    inorgsorbp <<- inorgsorbp + p_sorb_in - p_sorb_out

    #fprintf(stderr, "psorb calc %f\n", (9 * inorglabp)/(0.0012+inorglabp))

    # Daily increment of soil inorganic available P pool (lab + sorb)
    inorgavlp <<- inorglabp + inorgsorbp

    # Daily increment of soil inorganic secondary P pool (strongly sorbed)
    inorgssorbp <<- inorgssorbp + p_min_to_ssorb - p_ssorb_to_occ - p_ssorb_to_min

    # Daily increment of soil inorganic occluded P pool
    inorgoccp <<- inorgoccp + p_ssorb_to_occ

    # Daily increment of soil inorganic parent P pool
    inorgparp <<- inorgparp + p_atm_dep - p_par_to_min

    #fprintf(stderr, "inorglabp %f\n", inorglabp)
    #fprintf(stderr, "inorgsorbp %f\n", inorgsorbp)
    #fprintf(stderr, "inorgssorbp %f\n", inorgssorbp)
    #fprintf(stderr, "inorgoccp %f\n", inorgoccp)
    #fprintf(stderr, "inorgparp %f\n", inorgparp)
}

calculate_psoil_flows <- function(doy) {
    # need to store grazing flag. Allows us to switch on the annual
    # grazing event, but turn it off for every other day of the year.
    cntrl_grazing <- grazing
    if (grazing == 2 && disturbance_doy == doy+1) {
        grazing <<- TRUE
    }

    # Fraction of C lost due to microbial respiration
    frac_microb_resp <- 0.85 - (0.68 * finesoil)

    grazer_inputs_p()
    res__ <- p_inputs_from_plant_litter()
    psurf <- res__[1]
    psoil <- res__[2]
    partition_plant_litter_p(psurf, psoil)

    # SOM phosphorus effluxes.  These are assumed to have the source p:c
    # ratio prior to the increase of P:C due to co2 evolution.
    pfluxes_from_structural_pools()
    pfluxes_from_metabolic_pool()
    pfluxes_from_active_pool(frac_microb_resp)
    pfluxes_from_slow_pool()
    pfluxes_from_passive_pool()

    # calculate P parent influxe to mineral P
    calculate_p_parent_fluxes()

    # gross P mineralisation
    calculate_p_mineralisation()

    # calculate P immobilisation
    res__ <- calculate_p_immobilisation()
    pimmob <<- res__[1]
    active_pc_slope <- res__[2]
    slow_pc_slope <- res__[3]
    passive_pc_slope <- res__[4]

    # calculate P net mineralisation, excluding biochemical mineralisation
    calc_p_net_mineralisation()

    # calculate P biochemical mineralisation
    calculate_p_biochemical_mineralisation()

    # SIM phosphorus dynamics
    calculate_p_ssorb_to_sorb()
    calculate_p_ssorb_to_occ()
    calculate_p_sorb_to_ssorb()

    # calculate P lab and sorb fluxes from gross P flux
    calculate_p_min_fluxes()


    # Update model soil P pools
    calculate_ppools(active_pc_slope, slow_pc_slope, passive_pc_slope)

    # switch off grazing if this was just activated as an annual event
    grazing <<- cntrl_grazing
}

reset_all_n_pools_and_fluxes <- function() {
    # If the N-Cycle is turned off the way I am implementing this is to
    # do all the calculations and then reset everything at the end. This is
    # a waste of resources but saves on multiple IF statements.

    # State

    shootn <<- 0.0
    rootn <<- 0.0
    crootn <<- 0.0
    branchn <<- 0.0
    stemnimm <<- 0.0
    stemnmob <<- 0.0
    structsurfn <<- 0.0
    metabsurfn <<- 0.0
    structsoiln <<- 0.0
    metabsoiln <<- 0.0
    activesoiln <<- 0.0
    slowsoiln <<- 0.0
    passivesoiln <<- 0.0
    inorgn <<- 0.0
    stemn <<- 0.0
    stemnimm <<- 0.0
    stemnmob <<- 0.0
    nstore <<- 0.0

    # Fluxes

    nuptake <<- 0.0
    nloss <<- 0.0
    npassive <<- 0.0
    ngross <<- 0.0
    nimmob <<- 0.0
    nlittrelease <<- 0.0
    nmineralisation <<- 0.0
    npleaf <<- 0.0
    nproot <<- 0.0
    npcroot <<- 0.0
    npbranch <<- 0.0
    npstemimm <<- 0.0
    npstemmob <<- 0.0
    deadleafn <<- 0.0
    deadrootn <<- 0.0
    deadcrootn <<- 0.0
    deadbranchn <<- 0.0
    deadstemn <<- 0.0
    neaten <<- 0.0
    nurine <<- 0.0
    leafretransn <<- 0.0
    n_surf_struct_litter <<- 0.0
    n_surf_metab_litter <<- 0.0
    n_soil_struct_litter <<- 0.0
    n_soil_metab_litter <<- 0.0
    n_surf_struct_to_slow <<- 0.0
    n_soil_struct_to_slow <<- 0.0
    n_surf_struct_to_active <<- 0.0
    n_soil_struct_to_active <<- 0.0
    n_surf_metab_to_active <<- 0.0
    n_surf_metab_to_active <<- 0.0
    n_active_to_slow <<- 0.0
    n_active_to_passive <<- 0.0
    n_slow_to_active <<- 0.0
    n_slow_to_passive <<- 0.0
    n_passive_to_active <<- 0.0
}

reset_all_p_pools_and_fluxes <- function() {
    # If the P-Cycle is turned off the way I am implementing this is to
    # do all the calculations and then reset everything at the end. This is
    # a waste of resources but saves on multiple IF statements.

    # State

    shootp <<- 0.0
    rootp <<- 0.0
    crootp <<- 0.0
    branchp <<- 0.0
    stempimm <<- 0.0
    stempmob <<- 0.0
    structsurfp <<- 0.0
    metabsurfp <<- 0.0
    structsoilp <<- 0.0
    metabsoilp <<- 0.0
    activesoilp <<- 0.0
    slowsoilp <<- 0.0
    passivesoilp <<- 0.0
    inorgp <<- 0.0
    inorgavlp <<- 0.0
    inorglabp <<- 0.0
    inorgsorbp <<- 0.0
    inorgssorbp <<- 0.0
    inorgoccp <<- 0.0
    inorgparp <<- 0.0
    stemp <<- 0.0
    stempimm <<- 0.0
    stempmob <<- 0.0
    pstore <<- 0.0

    # Fluxes

    puptake <<- 0.0
    ploss <<- 0.0
    ppassive <<- 0.0
    pgross <<- 0.0
    pimmob <<- 0.0
    plittrelease <<- 0.0
    pmineralisation <<- 0.0
    ppleaf <<- 0.0
    pproot <<- 0.0
    ppcroot <<- 0.0
    ppbranch <<- 0.0
    ppstemimm <<- 0.0
    ppstemmob <<- 0.0
    deadleafp <<- 0.0
    deadrootp <<- 0.0
    deadcrootp <<- 0.0
    deadbranchp <<- 0.0
    deadstemp <<- 0.0
    peaten <<- 0.0
    purine <<- 0.0
    leafretransp <<- 0.0
    p_surf_struct_litter <<- 0.0
    p_surf_metab_litter <<- 0.0
    p_soil_struct_litter <<- 0.0
    p_soil_metab_litter <<- 0.0
    p_surf_struct_to_slow <<- 0.0
    p_soil_struct_to_slow <<- 0.0
    p_surf_struct_to_active <<- 0.0
    p_soil_struct_to_active <<- 0.0
    p_surf_metab_to_active <<- 0.0
    p_surf_metab_to_active <<- 0.0
    p_active_to_slow <<- 0.0
    p_active_to_passive <<- 0.0
    p_slow_to_active <<- 0.0
    p_slow_to_passive <<- 0.0
    p_slow_biochemical <<- 0.0
    p_passive_to_active <<- 0.0
    p_lab_in <<- 0.0
    p_lab_out <<- 0.0
    p_sorb_in <<- 0.0
    p_sorb_out <<- 0.0
    p_min_to_ssorb <<- 0.0
    p_ssorb_to_min <<- 0.0
    p_ssorb_to_occ <<- 0.0
    p_par_to_min <<- 0.0
    p_atm_dep <<- 0.0
}

write_daily_outputs_ascii <- function(year, doy) {
    # Write daily state and fluxes headers to an output CSV file. Note we
    # are not writing anything useful like units as there is a wrapper
    # script to translate the outputs to a nice CSV file with input met
    # data, units and nice header information.
    ofp[nrow(ofp) + 1,] <<- c(
        year, doy,
        wtfac_root, wtfac_topsoil, pawater_root,  pawater_topsoil,
        shoot, lai, branch, stem, root, croot,
        shootn, branchn, stemn, rootn, crootn,
        shootp, branchp, stemp, rootp, crootp,
        cstore, nstore, pstore,
        soilc, soiln, soilp, inorgn,
        inorgp, inorgavlp, inorglabp, inorgsorbp, inorgssorbp, inorgoccp, inorgparp,
        litterc, littercag, littercbg, litternag, litternbg,
        litterpag, litterpbg,
        activesoil, slowsoil, passivesoil,
        activesoiln, slowsoiln, passivesoiln, activesoilp, slowsoilp, passivesoilp,
        et, transpiration, soil_evap, canopy_evap,
        runoff, gs_mol_m2_sec, ga_mol_m2_sec,
        deadleaves, deadbranch, deadstems, deadroots, deadcroots,
        deadleafn, deadbranchn, deadstemn, deadrootn, deadcrootn,
        deadleafp, deadbranchp, deadstemp, deadrootp, deadcrootp,
        nep, gpp, npp, hetero_resp, auto_resp, apar,
        cpleaf, cpbranch, cpstem, cproot, cpcroot,
        npleaf, npbranch, npstemimm, npstemmob, nproot, npcroot,
        ppleaf, ppbranch, ppstemimm, ppstemmob, pproot, ppcroot,
        nuptake, ngross, nmineralisation, nloss,
        puptake, pgross, pmineralisation, ploss,
        leafretransn,
        leafretransp
    )
}

calculate_average_alloc_fractions <- function(days) {
    avg_alleaf <<- avg_alleaf / days
    avg_alroot <<- avg_alroot / days
    avg_alcroot <<- avg_alcroot / days
    avg_albranch <<- avg_albranch / days
    avg_alstem <<- avg_alstem / days

    alleaf <<- avg_alleaf
    alroot <<- avg_alroot
    alcroot <<- avg_alcroot
    albranch <<- avg_albranch
    alstem <<- avg_alstem

    # Because we are taking the average growing season fracs the total may
    # end up being just under 1, due to rounding. So put the missing
    # bit into the leaves - arbitary decision there

    excess <- 1.0 - alleaf - alroot - alcroot - albranch - alstem
    alleaf <<- alleaf + excess
}

update_roots <- function() {
    # Given the amount of roots grown by GDAY predict the assoicated rooting
    # distribution accross soil layers
    # - These assumptions come from Mat's SPA model.
    #
    soil_layers <- vector("numeric", n_layers)
    C_2_BIOMASS <- 2.0
    min_biomass <- 10.0  # To exend at least a layer g C m-2
    x1 <- 0.1            # lower bound for brent search
    x2 <- 10.0           # upper bound for brent search
    tol <- 0.0001        # tolerance for brent search

    # Enforcing a minimum fine root mass, otherwise during spinup this can go
    # wrong.
    fine_root_min <- 50.0
    if (root < fine_root_min) {
        fine_root <- fine_root_min
    } else {
        fine_root <- root
    }

    root_biomass <- max(min_biomass, fine_root * TONNES_HA_2_G_M2 * C_2_BIOMASS)
    # root_biomass <- MAX(min_biomass,  305.0 * C_2_BIOMASS)

    root_cross_sec_area <- pi * root_radius * root_radius   # (m2)
    root_depth <- max_depth * root_biomass / (root_k + root_biomass)

    rooted_layers <<- 0

    for (i in 1:n_layers) {
        if (layer_depth[i] > root_depth) {
            rooted_layers <<- i
            break
        }
    }

    # how for into the soil do the reach extend?
    root_reach <- layer_depth[rooted_layers]

    # Enforce 50 % of root mass in the top 1/4 of the rooted layers.
    mult <- min(1.0 / thickness[0],
               max(2.0, 11.0 * exp(-0.006 * root_biomass)))

    # assume surface root density is related to total root biomass by a
    # scalar dependent on root biomass

    surf_biomass <- root_biomass * mult

    if (rooted_layers > 1) {
        # determine slope of root distribution given rooting depth
        # and ratio of root mass to surface root density

        slope <- zbrent(calc_root_dist, x1, x2, tol, root_biomass,
                       surf_biomass, rooted_layers, thickness[0],
                       root_reach)

        prev <- 1.0 / slope
        cumulative_depth <- 0.0

        for (i in 1:rooted_layers) {
            cumulative_depth <- cumulative_depth + thickness[i]
            curr <- 1.0 / slope * exp(-slope * cumulative_depth)
            root_mass[i] <<- (prev - curr) * surf_biomass

            # (m m-3 soil)
            root_length[i] <<- root_mass[i] / (root_density *
                                                   root_cross_sec_area)
            prev <- curr
        }
    } else {
        root_mass[0] <<- root_biomass
    }
}
