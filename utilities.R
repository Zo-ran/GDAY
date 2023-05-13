is_leap_year <- function (yr) {
    return (yr %% 400 == 0 || (yr %% 100 != 0 && yr %% 4 == 0))
}

round_to_value <- function(number, roundto) {
    return (round(number / roundto) * roundto)
}

DEG2RAD <- function(DEG) {
    return (DEG * pi / 180.0)
}

RAD2DEG <- function(RAD) {
    return (180.0 * RAD / pi)
}

float_eq <- function(a, b) {
    # Are two floats approximately equal...?
    #
    # Reference:
    # ----------
    # D. E. Knuth. The Art of Computer Programming. Sec. 4.2.2 pp. 217-8.

    return (abs(a - b) <= EPSILON * abs(a))
}

alloc_goal_seek <- function (simulated, target, alloc_max, sensitivity) {
     # Sensitivity parameter characterises how allocation fraction respond
     #   when the leaf:sapwood area ratio departs from the target value
     #   If sensitivity close to 0 then the simulated leaf:sapwood area ratio
     #   will closely track the target value
    frac <- 0.5 + 0.5 * (1.0 - simulated / target) / sensitivity

    return (max(0.0, alloc_max * min(1.0, frac)))
}

time_till_next_disturbance <- function() {
    # calculate the number of years until a disturbance event occurs
    # assuming a return interval of X years
    #
    # - section 3.4.1 D. Knuth, The Art of Computer Programming.
    #
    # Parameters
    # ----------
    # return_interval : int/float
    #     interval disturbance return at in years
    # rate = 1.0 / return_interval

    #return int(-log(1.0 - random.random()) / rate)

    return (11)
}

calc_warmest_quarter_temp <- function(yr) {
    # calculate mean temperature of the warmest quarter
    # Ref: Atkin et al. (2015) New Phytologist
    # Table 6 best model for area-based broadleaved tree leaf respiration
    #
    # Parameters:
    # ----------
    # TWQ: float
    #      mean temperature of the warmest quarter
    # q1 - q4: float
    #      quarterly-based mean air temperature
    t1 <- 0
    t2 <- 0
    t3 <- 0
    t4 <- 0

    if (is_leap_year(yr)) {
        for (doy in 1:366) {
            if (doy >= 1 && doy <= 91) {
                t1 <- t1 + ma$tair[doy]
            } else if (doy >= 92 && doy <= 182) {
                t2 <- t2 + ma$tair[doy]
            } else if (doy >= 183 && doy <= 274) {
                t3 <- t3 + ma$tair[doy]
            } else if (doy >= 275 && doy <= 366) {
                t4 <- t4 + ma$tair[doy]
            }
        }
    } else {
        for (doy in 1:365) {
            if (doy >= 1 && doy <= 90) {
                t1 <- t1 + ma$tair[doy]
            } else if (doy >= 91 && doy <= 181) {
                t2 <- t2 + ma$tair[doy]
            } else if (doy >= 182 && doy <= 273) {
                t3 <- t3 + ma$tair[doy]
            } else if (doy >= 274 && doy <= 365) {
                t4 <- t4 + ma$tair[doy]
            }
        }
    }

    # compute quarterly mean temperature
    if (is_leap_year(yr)) {
       q1 <- t1 / 91.0
    } else {
       q1 <- t1 / 90.0
    }

    q2 <- t2 / 91.0
    q3 <- t3 / 92.0
    q4 <- t4 / 92.0

    TWQ <- max(q1, q2, q3, q4)

    twq <<- TWQ
}

calc_day_length <- function(doy, num_days, latitude) {

    # Daylength in hours
    #
    # Eqns come from Leuning A4, A5 and A6, pg. 1196
    #
    # Reference:
    # ----------
    # Leuning et al (1995) Plant, Cell and Environment, 18, 1183-1200.
    #
    # Parameters:
    # -----------
    # doy : int
    #     day of year, 1=jan 1
    # yr_days : int
    #     number of days in a year, 365 or 366
    # latitude : float
    #     latitude [degrees]
    #
    # Returns:
    # --------
    # dayl : float
    #     daylength [hrs]

    deg2rad <- pi / 180.0
    latr <- latitude * deg2rad
    sindec <- -sin(23.5 * deg2rad) * cos(2.0 * pi * (doy + 10.0) / num_days)
    a <- sin(latr) * sindec
    b <- cos(latr) * cos(asin(sindec))

    return (12.0 * (1.0 + (2.0 / pi) * asin(a / b)))
}

calculate_daylength <- function(num_days, latitude) {
    # wrapper to put the day length into an array
    for (i in 1:num_days) {
        day_length[i] <<- calc_day_length(i, num_days, latitude)
    }
}

gdd_chill_thresh <- function(pa, pb, pc, ncd) {
    # Leaf out has a chilling requirement, num chill days reduces the GDD
    # requirement
    return (pa + pb * exp(pc * ncd))
}

calc_gdd <- function(Tavg) {
    # calculate the number of growing degree days, hypothesis is that
    # leaves appear after a threshold has been passed.
    Tbase <- 5.0 # degC

    return (max(0.0, Tavg - Tbase))
}

leaf_drop <- function(daylen, Tsoil, Tsoil_next_3days) {
    # Thresholds to drop leaves come from White et al.
    # Note 655 minutes = 10.916 hrs.
    #
    # - Dependance on daylength means that this is only valid outside of the
    #   tropics.
    #
    # References:
    # -----------
    # White, M. A. et al. (1997) GBC, 11, 217-234.

    if ((daylen <= 10.9166667 && Tsoil <= 11.15) || Tsoil_next_3days < 2.0) {
        return (TRUE)
    } else {
        return (FALSE)
    }
}

calc_ncd <- function(Tmean) {
        # Calculate the number of chilling days from fixed dates (1 Nov),
        # following Murray et al. 1989, same as Botta does. """

    if (Tmean < 5.0) {
        return (1.0)
    } else {
        return (0.0)
    }
}

decay_in_dry_soils <- function(decay_rate, decay_rate_dry) {
    # Decay rates (e.g. leaf litterfall) can increase in dry soil, adjust
    # decay param. This is based on field measurements by F. J. Hingston
    # (unpublished) cited in Corbeels.
    #
    # Parameters:
    # -----------
    # decay_rate : float
    #     default model parameter decay rate [tonnes C/ha/day]
    # decay_rate_dry : float
    #     default model parameter dry deacy rate [tonnes C/ha/day]
    #
    # Returns:
    # --------
    # decay_rate : float
    #     adjusted deacy rate if the soil is dry [tonnes C/ha/day]
    #
    # Reference:
    # ----------
    # Corbeels et al. (2005) Ecological Modelling, 187, 449-474.

    # turn into fraction...
    smc_root <- pawater_root / wcapac_root

    new_decay_rate <- (decay_rate_dry - (decay_rate_dry - decay_rate) *
                     (smc_root - watdecaydry) /
                     (watdecaywet - watdecaydry))

    if (new_decay_rate < decay_rate)
        new_decay_rate <- decay_rate

    if (new_decay_rate > decay_rate_dry)
        new_decay_rate <- decay_rate_dry

    return (new_decay_rate)
}

daily_grazing_calc <- function(fdecay) {
    # daily grass grazing...
    #
    # Parameters:
    # -----------
    # fdecay : float
    #     foliage decay rate
    #
    # Returns:
    # --------
    # ceaten : float
    #     C consumed by grazers [tonnes C/ha/day]
    # neaten : float
    #     N consumed by grazers [tonnes N/ha/day]
    # peaten : float
    #     P consumed by grazers [tonnes P/ha/day]

    ceaten <<- fdecay * fracteaten / (1.0 - fracteaten) * shoot
    neaten <<- fdecay * fracteaten / (1.0 - fracteaten) * shootn
    peaten <<- fdecay * fracteaten / (1.0 - fracteaten) * shootp

}

annual_grazing_calc <- function() {
    # Annual grass grazing...single one off event
    #
    #
    # Returns:
    # --------
    # ceaten : float
    #     C consumed by grazers [tonnes C/ha/day]
    # neaten : float
    #     N consumed by grazers [tonnes N/ha/day]
    # peaten : float
    #     P consumed by grazers [tonnes P/ha/day]

    ceaten <<- shoot * fracteaten
    neaten <<- shootn * fracteaten
    peaten <<- shootp * fracteaten
}

calc_beta <- function(paw, depth, fc, wp, exponent) {
    # Soil water modifier, standard JULES/CABLE type approach
    #
    # equation 16 in Egea
    #
    # Note: we don't need to subtract the wp in the denominator here
    #       because our plant available water (paw) isn't bounded by
    #       the wilting point, it reaches zero
    #
    # Reference:
    # ----------
    # * Egea et al. (2011) Agricultural and Forest Meteorology.

    theta <- paw / depth
    beta <- (theta / (fc - wp)) ^ exponent
    if (beta > fc) {
        beta <- 1.0
    } else if (beta <= wp) {
        beta <- 0.0
    }

    return (beta)
}

calc_sw_modifier <- function(theta, c_theta, n_theta) {
    # Soil water modifier, equation 2 in Landsberg and Waring.
    # Note: "The values of c_theta and n_theta are, nevertheless, chosen
    #       without specific empirical justification" :)
    #
    # Reference:
    # ----------
    # * Landsberg and Waring (1997) Forest Ecology and Management 95, 209-228.
    return (1.0  / (1.0 + ((1.0 - theta) / c_theta) ^ n_theta))
}

calculate_top_of_canopy_n <- function(ncontent)  {
    # Calculate the canopy N at the top of the canopy (g N m-2), N0.
    # Assuming an exponentially decreasing N distribution within the canopy:
    #
    # Note: swapped kext with kn
    #
    # Returns:
    # -------
    # N0 : float (g N m-2)
    #     Top of the canopy N
    #
    # References:
    # -----------
    # * Chen et al 93, Oecologia, 93,63-69.

    if (lai > 0.0) {
        # calculation for canopy N content at the top of the canopy
        N0 <- ncontent * kn / (1.0 - exp(-kn * lai))
    } else {
        N0 <- 0.0
    }

    return (N0)
}

calculate_top_of_canopy_p <- function(pcontent)  {
    # Calculate the canopy P at the top of the canopy (g P m-2), P0.
    # Assuming an exponentially decreasing P distribution within the canopy:
    #
    # Note: swapped kext with kp
    #
    # Returns:
    # -------
    # P0 : float (g P m-2)
    # Top of the canopy P


    if (lai > 0.0) {
        # calculation for canopy P content at the top of the canopy
        P0 <- pcontent * kp / (1.0 - exp(-kp * lai))
    } else {
        P0 <- 0.0
    }

    return (P0)
}

arrh <- function(mt, k25, Ea, Tk) {
    # Temperature dependence of kinetic parameters is described by an
    # Arrhenius function
    #
    # Parameters:
    # ----------
    # k25 : float
    #     rate parameter value at 25 degC
    # Ea : float
    #     activation energy for the parameter [J mol-1]
    # Tk : float
    #     leaf temperature [deg K]
    #
    # Returns:
    # -------
    # kt : float
    #     temperature dependence on parameter
    #
    # References:
    # -----------
    # * Medlyn et al. 2002, PCE, 25, 1167-1179.
    return (k25 * exp((Ea * (Tk - mt)) / (mt * RGAS * Tk)))
}

calculate_co2_compensation_point <- function(Tk, mt) {
    # CO2 compensation point in the absence of mitochondrial respiration
    # Rate of photosynthesis matches the rate of respiration and the net CO2
    # assimilation is zero.
    #
    # Parameters:
    # ----------
    # Tk : float
    #     air temperature (Kelvin)
    #
    # Returns:
    # -------
    # gamma_star : float
    #     CO2 compensation point in the abscence of mitochondrial respiration
    return (arrh(mt, gamstar25, eag, Tk))
}

calculate_michaelis_menten_parameter <- function(Tk, mt) {
    # Effective Michaelis-Menten coefficent of Rubisco activity
    #
    # Parameters:
    # ----------
    # Tk : float
    #     air temperature (Kelvin)
    #
    # Returns:
    # -------
    # Km : float
    #     Effective Michaelis-Menten constant for Rubisco catalytic activity
    #
    # References:
    # -----------
    # Rubisco kinetic parameter values are from:
    # * Bernacchi et al. (2001) PCE, 24, 253-259.
    # * Medlyn et al. (2002) PCE, 25, 1167-1179, see pg. 1170.

    # Michaelis-Menten coefficents for carboxylation by Rubisco
    Kc <- arrh(mt, kc25, eac, Tk)

    # Michaelis-Menten coefficents for oxygenation by Rubisco
    Ko <- arrh(mt, ko25, eao, Tk)

    # return effective Michaelis-Menten coefficient for CO2
    return ( Kc * (1.0 + oi / Ko) )

}

peaked_arrh <- function(mt, k25, Ea, Tk, deltaS, Hd) {
    # Temperature dependancy approximated by peaked Arrhenius eqn,
    # accounting for the rate of inhibition at higher temperatures.
    #
    # Parameters:
    # ----------
    # k25 : float
    #     rate parameter value at 25 degC
    # Ea : float
    #     activation energy for the parameter [J mol-1]
    # Tk : float
    #     leaf temperature [deg K]
    # deltaS : float
    #     entropy factor [J mol-1 K-1)
    # Hd : float
    #     describes rate of decrease about the optimum temp [J mol-1]
    #
    # Returns:
    # -------
    # kt : float
    #     temperature dependence on parameter
    #
    # References:
    # -----------
    # * Medlyn et al. 2002, PCE, 25, 1167-1179.

    arg1 <- arrh(mt, k25, Ea, Tk)
    arg2 <- 1.0 + exp((mt * deltaS - Hd) / (mt * RGAS))
    arg3 <- 1.0 + exp((Tk * deltaS - Hd) / (Tk * RGAS))


    return (arg1 * arg2 / arg3)
}

adj_for_low_temp <- function(param, Tk) {
    # Function allowing Jmax/Vcmax to be forced linearly to zero at low T
    #
    # Parameters:
    # ----------
    # Tk : float
    #     air temperature (Kelvin)
    lower_bound <- 0.0
    upper_bound <- 10.0

    Tc <- Tk - DEG_TO_KELVIN

    if (Tc < lower_bound)
        param <- 0.0
    else if (Tc < upper_bound)
        param <- param * (Tc - lower_bound) / (upper_bound - lower_bound)

    return (param)
}

calculate_ci <- function(vpd, Ca) {
    # Calculate the intercellular (Ci) concentration
    #
    # Formed by substituting gs = g0 + 1.6 * (1 + (g1/sqrt(D))) * A/Ca into
    # A = gs / 1.6 * (Ca - Ci) and assuming intercept (g0) = 0.
    #
    # Parameters:
    # ----------
    # vpd : float
    #     vapour pressure deficit [Pa]
    # Ca : float
    #     ambient co2 concentration
    #
    # Returns:
    # -------
    # ci:ca : float
    #     ratio of intercellular to atmospheric CO2 concentration
    #
    # References:
    # -----------
    # * Medlyn, B. E. et al (2011) Global Change Biology, 17, 2134-2144.

    ci <- 0.0

    if (gs_model == MEDLYN) {
        g1w <- g1 * wtfac_root
        cica <- g1w / (g1w + sqrt(vpd * PA_2_KPA))
        ci <- cica * Ca
    } else {
        stop("Only Belindas gs model is implemented")
    }

    return (ci)
}

assim <- function(ci, gamma_star, a1, a2) {
    # Morning and afternoon calcultion of photosynthesis with the
    # limitation defined by the variables passed as a1 and a2, i.e. if we
    # are calculating vcmax or jmax limited.
    #
    # Parameters:
    # ----------
    # ci : float
    #     intercellular CO2 concentration.
    # gamma_star : float
    #     CO2 compensation point in the abscence of mitochondrial respiration
    # a1 : float
    #     variable depends on whether the calculation is light or rubisco
    #     limited.
    # a2 : float
    #     variable depends on whether the calculation is light or rubisco
    #     limited.
    #
    # Returns:
    # -------
    # assimilation_rate : float
    #     assimilation rate assuming either light or rubisco limitation.
    if (ci < gamma_star) {
        return (0.0)
    } else {
        return (a1 * (ci - gamma_star) / (a2 + ci))
    }
}

calculate_quantum_efficiency <- function(ci, gamma_star) {
    # Quantum efficiency for AM/PM periods replacing Sands 1996
    # temperature dependancy function with eqn. from Medlyn, 2000 which is
    # based on McMurtrie and Wang 1993.
    #
    # Parameters:
    # ----------
    # ci : float
    #     intercellular CO2 concentration.
    # gamma_star : float [am/pm]
    #     CO2 compensation point in the abscence of mitochondrial respiration
    #
    # Returns:
    # -------
    # alpha : float
    #     Quantum efficiency
    #
    # References:
    # -----------
    # * Medlyn et al. (2000) Can. J. For. Res, 30, 873-888
    # * McMurtrie and Wang (1993) PCE, 16, 1-13.

    return (assim(ci, gamma_star, alpha_j / 4.0, 2.0 * gamma_star))
}

assim_p <- function(P0) {
    # // Calculate photosynthesis assimilation rate based on P limitation
    # // Ref: Ellsworth et al., 2015, Plant, Cell and Environment, able 2
    # //
    # // Returns:
    # // -------
    # // ap: float
    # // assimilation rate assuming P limitation

    conv <- MMOL_2_MOL * 31.0
    tp <- 6.51 + 14.64 * (P0 * conv)
    ap <- 3.0 * tp

    return(ap)
}

epsilon <- function(asat, par, alpha, daylen) {
    # Canopy scale LUE using method from Sands 1995, 1996.
    #
    # Sands derived daily canopy LUE from Asat by modelling the light response
    # of photosysnthesis as a non-rectangular hyperbola with a curvature
    # (theta) and a quantum efficiency (alpha).
    #
    # Assumptions of the approach are:
    #  - horizontally uniform canopy
    #  - PAR varies sinusoidally during daylight hours
    #  - extinction coefficient is constant all day
    #  - Asat and incident radiation decline through the canopy following
    #    Beer's Law.
    #  - leaf transmission is assumed to be zero.
    #
    # * Numerical integration of "g" is simplified to 6 intervals.
    #
    # Parameters:
    # ----------
    # asat : float
    #     Light-saturated photosynthetic rate at the top of the canopy
    # par : float
    #     photosyntetically active radiation (umol m-2 d-1)
    # theta : float
    #     curvature of photosynthetic light response curve
    # alpha : float
    #     quantum yield of photosynthesis (mol mol-1)
    #
    # Returns:
    # -------
    # lue : float
    #     integrated light use efficiency over the canopy (umol C umol-1 PAR)
    #
    # Notes:
    # ------
    # NB. I've removed solar irradiance to PAR conversion. Sands had
    # gamma = 2000000 to convert from SW radiation in MJ m-2 day-1 to
    # umol PAR on the basis that 1 MJ m-2 = 2.08 mol m-2 & mol to umol = 1E6.
    # We are passing PAR in umol m-2 d-1, thus avoiding the above.
    #
    # References:
    # -----------
    # See assumptions above...
    # * Sands, P. J. (1995) Australian Journal of Plant Physiology,
    #   22, 601-14.

    # subintervals scalar, i.e. 6 intervals
    delta <- 0.16666666667

    # number of seconds of daylight
    h <- daylen * SECS_IN_HOUR

    if (asat > 0.0) {
        # normalised daily irradiance
        q <- pi * kext * alpha * par / (2.0 * h * asat)
        integral_g <- 0.0
        for (i in seq(1, 12, by=2)) {
            sinx <- sin(pi * i / 24)
            arg1 <- sinx
            arg2 <- 1.0 + q * sinx
            arg3 <- (sqrt((1.0 + q * sinx) ^ 2 - 4.0 * theta * q * sinx))
            integral_g <- integral_g + arg1 / (arg2 + arg3)
        }
        integral_g <- integral_g * delta
        lue <- alpha * integral_g * pi
    } else {
        lue <- 0.0
    }

    return (lue)
}

quadratic <- function(a, b, c) {
    # minimilist quadratic solution
    #
    # Parameters:
    # ----------
    # a : float
    #     co-efficient
    # b : float
    #     co-efficient
    # c : float
    #     co-efficient
    #
    # Returns:
    # -------
    # root : float

    # discriminant
    d <- b ^ 2.0 - 4.0 * a * c

    # Negative quadratic equation
    root <- (-b - sqrt(d)) / (2.0 * a)

    return (root)
}

calc_respiration <- function(Tk, vcmax25) {
    # Mitochondrial respiration may occur in the mesophyll as well as in the
    # bundle sheath. As rubisco may more readily refix CO2 released in the
    # bundle sheath, Rd is described by its mesophyll and bundle-sheath
    # components: Rd = Rm + Rs
    #
    # Parameters:
    # ----------
    # Tk : float
    #     air temperature (kelvin)
    # vcmax25 : float, list
    #
    # Returns:
    # -------
    # Rd : float, list [am, pm]
    #     (respiration in the light) 'day' respiration (umol m-2 s-1)
    #
    # References:
    # -----------
    # Tjoelker et al (2001) GCB, 7, 223-230.
    Tref <- 25.0

    # ratio between respiration rate at one temperature and the respiration
    # rate at a temperature 10 deg C lower
    Q10 <- 2.0

    # scaling constant to Vcmax25, value = 0.015 after Collatz et al 1991.
    # Agricultural and Forest Meteorology, 54, 107-136. But this if for C3
    # Value used in JULES for C4 is 0.025, using that one, see Clark et al.
    # 2011, Geosci. Model Dev, 4, 701-722.
    fdr <- 0.025

    # specific respiration at a reference temperature (25 deg C)
    Rd25 <- fdr * vcmax25

    return (Rd25 * (Q10 ^ (((Tk - DEG_TO_KELVIN) - Tref) / 10.0)))
}

lloyd_and_taylor <- function(temp) {
    # Modified Arrhenius equation (Lloyd & Taylor, 1994)
    # The modification introduced by Lloyd & Taylor (1994) represents a
    # decline in the parameter for activation energy with temperature.
    #
    # Parameters:
    # -----------
    # temp : float
    #      temp deg C

    return (exp(308.56 * ((1.0 / 56.02) - (1.0 / (temp + 46.02)))))
}

calc_net_radiation <- function(sw_rad, tair) {

    # Net loss of long-wave radn, Monteith & Unsworth '90, pg 52, eqn 4.17
    net_lw <- 107.0 - 0.3 * tair            # W m-2

    # Net radiation recieved by a surf, Monteith & Unsw '90, pg 54 eqn 4.21
    #     - note the minus net_lw is correct as eqn 4.17 is reversed in
    #       eqn 4.21, i.e Lu-Ld vs. Ld-Lu
    #     - NB: this formula only really holds for cloudless skies!
    #     - Bounding to zero, as we can't have negative soil evaporation, but you
    #       can have negative net radiation.
    #     - units: W m-2

    net_rad <- max(0.0, (1.0 - albedo) * sw_rad - net_lw)

    return (net_rad)
}

calc_stomatal_conductance <- function(vpd, Ca, gpp) {
    #     Calculate stomatal conductance using Belinda's model. For the medlyn
    #     model this is already in conductance to CO2, so the 1.6 from the
    #     corrigendum to Medlyn et al 2011 is missing here
    #
    # References:
    # -----------
    # For conversion factor for conductance see...
    # * Jones (1992) Plants and microclimate, pg 56 + Appendix 3
    # * Diaz et al (2007) Forest Ecology and Management, 244, 32-40.
    #
    # Stomatal Model:
    # * Medlyn et al. (2011) Global Change Biology, 17, 2134-2144.
    # **Note** Corrigendum   Global Change Biology, 18, 3476.
    #
    # Parameters:
    # -----------
    # g1 : float
    #     slope
    # wtfac : float
    #     water availability scaler [0,1]
    # vpd : float
    #     vapour pressure deficit (Pa)
    # Ca : float
    #     atmospheric co2 [umol mol-1]
    # gpp : float
    #     photosynthesis at the canopy scale (umol m-2 s-1)
    #
    # Returns:
    # --------
    # gs : float
    #     stomatal conductance (mol CO2 m-2 s-1)

    g1__ <- g1 * wtfac_root
    g0__ <- 0.0  # g0
    gs_over_a <- (1.0 + g1__ / sqrt(vpd * PA_2_KPA)) / Ca
    gsc <- max(g0__, g0__ + gs_over_a * gpp)

    # mol m-2 s-1
    return (gsc)

}

canopy_boundary_layer_conduct <- function(canht, wind, press, tair) {
    # Canopy boundary layer conductance, ga (from Jones 1992 p 68)
    #
    #
    # Notes:
    # ------
    # 'Estimates of ga for pine canopies from LAI of 3 to 6 vary from
    # 3.5 to 1.1 mol m-2 s-1  (Kelliher et al., 1993 Juang et al., 2007).'
    # Drake et al, 2010, 17, pg. 1526.
    #
    # References:
    # ------------
    # * Jones 1992, pg. 67-8.
    # * Monteith and Unsworth (1990), pg. 248. Note this in the inverted form
    #   of what is in Monteith (ga = 1 / ra)
    # * Allen et al. (1989) pg. 651.
    # * Gash et al. (1999) Ag forest met, 94, 149-158.
    #
    # Parameters:
    # -----------
    # params : p
    #     parameters structure
    # canht : float
    #     canopy height (m)
    # wind : float
    #     wind speed (m s-1)
    # press : float
    #     atmospheric pressure (Pa)
    # tair : float
    #     air temperature (deg C)
    #
    # Returns:
    # --------
    # ga : float
    #     canopy boundary layer conductance (mol m-2 s-1)

    # z0m roughness length governing momentum transfer [m]
    vk <- 0.41

    # Convert from mm s-1 to mol m-2 s-1
    cmolar <- press / (RGAS * (tair + DEG_TO_KELVIN))

    # roughness length for momentum
    z0m <- dz0v_dh * canht

    #  z0h roughness length governing transfer of heat and vapour [m]
    # *Heat tranfer typically less efficent than momentum transfer. There is
    #  a lot of variability in values quoted for the ratio of these two...
    #  JULES uses 0.1, Campbell and Norman '98 say z0h = z0m / 5. Garratt
    #  and Hicks, 1973/ Stewart et al '94 say z0h = z0m / 7. Therefore for
    #  the default I am following Monteith and Unsworth, by setting the
    #  ratio to be 1, the code below is identical to that on page 249,
    #  eqn 15.7

    z0h <- z0h_z0m * z0m

    #s zero plan displacement height [m]
    d <- displace_ratio * canht

    arg1 <- (vk * vk) * wind
    arg2 <- log((canht - d) / z0m)
    arg3 <- log((canht - d) / z0h)

    ga <- (arg1 / (arg2 * arg3)) * cmolar

    return (ga)
}

calc_sat_water_vapour_press <- function(tac) {
    # Calculate saturated water vapour pressure (Pa) at
    # temperature TAC (Celsius). From Jones 1992 p 110 (note error in
    # a - wrong units)
    return (613.75 * exp(17.502 * tac / (240.97 + tac)))
}

calc_slope_of_sat_vapour_pressure_curve <- function(tair) {
    # Constant slope in Penman-Monteith equation
    #
    # Parameters:
    # -----------
    # tavg : float
    #     average daytime temperature
    #
    # Returns:
    # --------
    # slope : float
    #     slope of saturation vapour pressure curve [Pa K-1]

    # Const slope in Penman-Monteith equation  (Pa K-1)
    arg1 <- calc_sat_water_vapour_press(tair + 0.1)
    arg2 <- calc_sat_water_vapour_press(tair)
    slope <- (arg1 - arg2) / 0.1

    return (slope)
}


calc_pyschrometric_constant <- function(press, lambda) {
    # Psychrometric constant ratio of specific heat of moist air at
    # a constant pressure to latent heat of vaporisation.
    #
    # Parameters:
    # -----------
    # press : float
    #     air pressure (Pa)
    # lambda : float
    #      latent heat of water vaporization (J mol-1)
    #
    #
    # Returns:
    # --------
    # gamma : float
    #     pyschrometric constant [Pa K-1]

    return ( CP * MASS_AIR * press / lambda )

}

calc_latent_heat_of_vapourisation <- function(tair) {
    # Latent heat of water vapour at air temperature
    #
    # Returns:
    # -----------
    # lambda : float
    #     latent heat of water vaporization [J mol-1]
    return ( (H2OLV0 - 2.365E3 * tair) * H2OMW )

}

penman_monteith <- function(press, vpd, rnet, slope, lambda, gamma, gh, gv) {
    # Calculates transpiration using the Penman-Monteith
    #
    # Parameters:
    # ----------
    # press : float
    #     atmospheric pressure (Pa)
    # vpd : float
    #     vapour pressure deficit of air (Pa)
    # rnet : float
    #     net radiation (J m-2 s-1)
    # slope : float
    #     slope of VPD/T curve, Pa K-1
    # lambda : flot
    #     latent heat of water at air T, J mol-1
    # gamma : float
    #     psychrometric constant, J mol-1
    # gh : float
    #     boundary layer conductance to heat (free & forced & radiative
    #     components), mol m-2 s-1
    # gv : float
    #     conductance to water vapour (stomatal & bdry layer components),
    #     mol m-2 s-1
    # transpiration : float
    #     transpiration (mol H2O m-2 s-1 returned)
    # LE : float
    #     latent heat flux (W m-2 returned)
    # omega : float
    #     decoupling coefficient (unitless returned)
    #
    # References:
    # ------------
    # * Medlyn et al. (2007), Tree Physiology, 27, 1687-1699.

    if (gv > 0.0) {
        arg1 <- slope * rnet + vpd * gh * CP * MASS_AIR
        arg2 <- slope + gamma * gh / gv
        LE <- arg1 / arg2 # W m-2
        transpiration <- LE / lambda # mol H20 m-2 s-1
    } else {
        transpiration <- 0.0
        LE <- 0
    }

    # Should not be negative - not sure gv>0.0 catches it as g0 = 1E-09?
    transpiration <- max(0.0, transpiration)

    return (c(gh, gv, transpiration, LE))
}

penman_canopy_wrapper <- function(press, vpd, tair, wind, rnet, ca, gpp) {
    # Calculates transpiration at the canopy scale (or big leaf) using the
    # Penman-Monteith
    #
    # Parameters:
    # ----------
    # parms : structure
    #     parameters
    # state : structure
    #     state variables
    # press : float
    #     atmospheric pressure (Pa)
    # vpd : float
    #     vapour pressure deficit of air (Pa)
    # tair : float
    #     air temperature (deg C)
    # wind : float
    #     wind speed (m s-1)
    # rnet : float
    #     net radiation (J m-2 s-1)
    # ca : float
    #     ambient CO2 concentration (umol mol-1)
    # gpp : float
    #     gross primary productivity (umol m-2 s-1)
    # ga : float
    #     canopy scale boundary layer conductance (mol m-2 s-1 returned)
    # gsv : float
    #     stomatal conductance to H20 (mol m-2 s-1 returned)
    # transpiration : float
    #     transpiration (mol H2O m-2 s-1 returned)
    # LE : float
    #     latent heat flux (W m-2 returned)

    # stomtal conductance to CO2
    gsc <- calc_stomatal_conductance(vpd, ca, gpp)

    # stomtal conductance to H2O
    gsv <- GSVGSC * gsc

    ga <- canopy_boundary_layer_conduct(canht, wind, press, tair)

    # Total leaf conductance to water vapour
    gv <- 1.0 / (1.0 / gsv + 1.0 / ga)

    lambda <- calc_latent_heat_of_vapourisation(tair)
    gamma <- calc_pyschrometric_constant(press, lambda)
    slope <- calc_slope_of_sat_vapour_pressure_curve(tair)

    res__ <- penman_monteith(press, vpd, rnet, slope, lambda, gamma, ga, gv)
    ga <- res__[1]
    gv <- res__[2]
    transpiration <- res__[3]
    LE <- res__[4]


    # Calculate decoupling coefficient (McNaughton and Jarvis 1986)
    epsilon <- slope / gamma
    omega <- (1.0 + epsilon) / (1.0 + epsilon + ga / gsv)
    return (c(ga, gsv, transpiration, LE, omega))
}

nitrogen_retrans <- function(fdecay, rdecay, doy) {
    # Nitrogen retranslocated from senesced plant matter.
    # Constant rate of n translocated from mobile pool
    #
    # Parameters:
    # -----------
    # fdecay : float
    #     foliage decay rate
    # rdecay : float
    #     fine root decay rate
    #
    # Returns:
    # --------
    # N retrans : float
    #     N retranslocated plant matter

    if (deciduous_model) {
        leafretransn <- fretrans * lnrate * remaining_days[doy]
    } else {
        leafretransn <- fretrans * fdecay * shootn
    }

    rootretransn <- rretrans * rdecay * rootn
    crootretransn <- cretrans * crdecay * crootn
    branchretransn <- bretrans * bdecay * branchn
    stemretransn <- (wretrans * wdecay * stemnmob + retransmob * stemnmob)

    # store for NCEAS output
    leafretransn <<- leafretransn

    return (leafretransn + rootretransn + crootretransn + branchretransn + stemretransn)
}

phosphorus_retrans <- function(fdecay, rdecay, doy) {
    # Phosphorus retranslocated from senesced plant matter.
    # Constant rate of p translocated from mobile pool
    #
    # Parameters:
    # -----------
    # fdecay : float
    # foliage decay rate
    # rdecay : float
    # fine root decay rate
    #
    # Returns:
    # --------
    # P retrans : float
    # P retranslocated plant matter

    if (deciduous_model) {
        leafretransp <- fretransp * lprate * remaining_days[doy]
    } else {
        leafretransp <- fretransp * fdecay * shootp
    }

    rootretransp <- rretrans * rdecay * rootp
    crootretransp <- cretrans * crdecay * crootp
    branchretransp <- bretrans * bdecay * branchp
    stemretransp <- (wretrans * wdecay * stempmob + retransmob * stempmob)

    # store for NCEAS output
    leafretransp <<- leafretransp

    return (leafretransp + rootretransp + crootretransp + branchretransp + stemretransp)
}

calculate_nuptake <- function() {
    # N uptake depends on the rate at which soil mineral N is made
    # available to the plants.
    #
    # Returns:
    # --------
    # nuptake : float
    #     N uptake
    #
    # References:
    # -----------
    # * Dewar and McMurtrie, 1996, Tree Physiology, 16, 161-171.
    # * Raich et al. 1991, Ecological Applications, 1, 399-429.

    if (nuptake_model == 0) {
        # Constant N uptake
        nuptake <- nuptakez
    } else if (nuptake_model == 1) {
        # evaluate nuptake : proportional to dynamic inorganic N pool
        nuptake <- rateuptake * inorgn
    } else if (nuptake_model == 2) {
        # N uptake is a saturating function on root biomass following
        # Dewar and McMurtrie, 1996.

        # supply rate of available mineral N
        U0 <- rateuptake * inorgn
        Kr <- kr
        nuptake <- max(U0 * root / (root + Kr), 0.0)

        # fprintf(stderr, "inorgn %f\n", inorgn)
        # fprintf(stderr, "nuptake %f\n", nuptake)

        # Make minimum uptake rate supply rate for deciduous_model cases
        # otherwise it is possible when growing from scratch we don't have
        # enough root mass to obtain N at the annual time step
        # I don't see an obvious better solution?
    } else {
        stop("Unknown N uptake option\n")
    }

    return (nuptake)
}


calculate_puptake <- function() {
    # P uptake depends on the rate at which soil mineral P is made
    # available to the plants.
    #
    # Returns:
    # --------
    # puptake : float
    # P uptake

    if (puptake_model == 0) {
        # Constant P uptake
        puptake <- puptakez
    } else if (puptake_model == 1) {
        # evaluate puptake : proportional to lab P pool that is
        # available to plant uptake
        puptake <- prateuptake * inorglabp * p_lab_avail
    } else if (puptake_model == 2) {
        # P uptake is a saturating function on root biomass, as N

        # supply rate of available mineral P
        if (inorgsorbp > 0.0) {
            U0 <- prateuptake * inorglabp * p_lab_avail
        } else {
            U0 <- min((p_par_to_min + pmineralisation +
                     purine + p_slow_biochemical),
                     (prateuptake * inorglabp * p_lab_avail))
        }

        Kr <- krp
        puptake <- max(U0 * root / (root + Kr), 0.0)
    } else {
        stop("Unknown P uptake option\n")
    }

    return (puptake)
}

rtot <- function(dmax, r0, d0) {
    # Estimate the total root biomass per unit ground area, i.e. the
    # integral of root mass per unit soil volume, R(z), over depth z from the
    # soil surface to the maximim rooting depth, dmax. (Eqn 8, in McM 2012)
    #
    # Parameters:
    # -----------
    # dmax : float
    #     Rooting depth [m]
    # r0 : float
    #     Root C at half-max N uptake.
    # d0 : float
    #     Length scale for exponential decline of Umax(z)
    #
    # Returns:
    # --------
    # rtot : float
    #     Total root C mass given a rooting depth

    return (r0 * (2.0 * d0 * (exp(dmax / (2.0 * d0)) - 1.0) - dmax))
}

rtot_wrapper <- function(dmax, rtoti, r0, d0) {
    # Wrapper method that calls rtot. Need to subtract rtoti because we
    # are optimising the depth (Dmax) but the rtot estimate has to match the
    # rtot from GDAY, i.e. diff = 0.0
    #
    # Parameters:
    # -----------
    # rtoti : float
    #     Initial fine root root C mass [from G'DAY]
    # rtot : function
    #     call rtot function to estimate a value of rtot given a rooting
    #     depth iteration
    #
    # Returns:
    # --------
    # val  : float
    #     A optimised rooting depth iteration

    return (rtot(dmax, r0, d0) - rtoti)
}

rtot_derivative <- function(dmax, rtoti, r0, d0) {
    # Derivative of maximum root depth equation, rtot
    #
    # Parameters:
    # -----------
    # dmax : float
    #     Rooting depth [m]
    # rtoti : float
    #     Initial fine root root C mass [from G'DAY]
    # r0 : float
    #     Root C at half-max N uptake.
    # d0 : float
    #     Length scale for exponential decline of Umax(z)
    #
    # Returns:
    # --------
    # val : float
    #     derivative of rtot

    return (r0 * (exp(0.5 * dmax / d0) - 1.0) )
}

newton <- function(func, fprime, x0, arg1, arg2, arg3) {
    # Newton-Raphson: finds a zero of the func, given an inital guess
    #
    # Parameters
    # ----------
    # f : function
    #     The function whose zero is wanted.
    # x0 : float
    #     An initial guess.
    # fprime : function
    #     The derivative of the function
    # args : tuple, optional
    #     Extra arguments to be used in the function call.
    # tol : float, optional
    #     The allowable error of the zero value.
    # maxiter : int, optional
    #     Maximum number of iterations.
    #
    # Returns
    # -------
    # val : float
    #     Estimated location where function is zero
    maxiter <- 250
    tol <- 1E-6
    for (iter in 1:maxiter) {

        x <- (x0 - (func(x0, arg1, arg2, arg3) / fprime(x0, arg1, arg2, arg3)))
        if (abs(x - x0) < tol) {
            return (x)
        }
        x0 <- x
    }
    stop(stderr, "Minimum not found!!\n")
}

estimate_max_root_depth <- function(rtoti, depth_guess, r0, d0) {
    # Determing the maximum rooting depth through solving Eqn. B6. for
    # rooting depth
    #
    # Parameters:
    # -----------
    # rtoti : float
    #     Initial fine root root C mass [from G'DAY]
    # depth_guess : float
    #     initial starting guess at the root depth [m]
    #
    # Returns:
    # --------
    # rooting_depth : float
    #     optimised rooting depth [m]

    fPtr <- rtot_wrapper
    fprimePtr <- rtot_derivative

    root_depth <- newton(fPtr, fprimePtr, depth_guess, rtoti, r0, d0)

    return (root_depth)
}

calc_umax <- function(nsupply, top_soil_depth, d0) {
    # Calculate potential N uptake integrated over all soil depths
    #
    # Parameters
    # ----------
    # nsupply : float
    #     N supply rate to a specified soil depth (probably 30 cm)
    #
    # Returns
    # -------
    # Umax : float
    #     potential N uptake integrated over all soil depths
    return (nsupply / (1.0 - exp(-top_soil_depth / d0)))
}

calc_plant_nuptake <- function(root_depth, nsupply, d0, top_soil_depth) {
    # Plant N uptake (Utot) as a func of maximum rooting depth
    #
    # This is the alternative eqn from McM word document
    #
    # Parameters
    # ----------
    # root_depth : float
    #     max rooting depth [m]
    # z : float
    #     incremental depth provided by integration func
    # nsupply : float
    #     soil N supply rate to plants per day [N/m2]
    # top_soil_depth : float
    #     Depth of soil assumed by G'DAY model [m]
    #
    # Returns
    # -------
    # nuptake : float
    #     plant N uptake

    Umax <- calc_umax(nsupply, top_soil_depth, d0)
    arg <- 1.0 - exp(-root_depth / (2.0 * d0))

    return (Umax * (arg * arg))
}

calculate_root_mass_above_depth <- function(rtoti, root_depth,
                                       r0, d0, top_soil_depth) {
    # Estimate cumulative root mass above depth, 30 cm for the G'DAY model
    #
    # Parameters
    # ----------
    # rtoti : float
    #     Initial fine root root C mass [from G'DAY]
    # root_depth : float
    #     model rooting depth (m)
    #
    # Returns
    # -------
    # val : float
    #     cumulative root C mass above soil depth assumed by G'DAY model, 30cm

    arg1 <-- rtoti + 2.0 * r0 * d0 + root_depth * r0
    arg2 <-- 1.0 - exp(-top_soil_depth / (2.0 * d0))
    return (arg1 * arg2 - r0 * top_soil_depth)
}

calc_opt_root_depth <- function(d0, r0, top_soil_depth,
                         rtoti, nsupply, depth_guess) {
    # Parameters:
    # -----------
    # d0 : float
    #     Length scale for exponential decline of Umax(z)
    # r0 : float
    #     root C at half-maximum N uptake (kg C/m3)
    # top_soil_depth : float
    #     depth of soil assumed by G'DAY, note Ross comment about 30 cm (units=m)
    #     [email]
    # rtoti : float
    #     Initial fine root C mass   from G'DAY
    # nsupply : float
    #     daily net N mineralisation in top soil layer from G'DAY
    # depth_guess : float
    #     Initial guess at the rooting depth, used as the first point in the
    #     root depth optimisation scheme [m].
    #
    # Returns:
    # --------
    # root_depth : float
    #     rooting depth [m]
    # nuptake : float
    #     N uptake from roots [gN m-2 yr-1]
    # rabove : float

    # Determine maximum rooting depth for model for a value root C
    depth <- estimate_max_root_depth(rtoti, depth_guess, r0, d0)
    root_depth <- depth
    # Optimised plant N uptake
    nuptake <- calc_plant_nuptake(depth, nsupply, d0, top_soil_depth)

    # G'DAY requires root litter input to the top 30 cm of soil, so
    # return the roots above this depth
    print(depth)

    rabove <- calculate_root_mass_above_depth(rtoti, depth, r0, d0,
                                              top_soil_depth)

    return (c(root_depth, nuptake, rabove))
}

arrhenius <- function(k25, Ea, T, Tref) {
    # Temperature dependence of kinetic parameters is described by an
    # Arrhenius function
    #
    # Parameters:
    # ----------
    # k25 : float
    #     rate parameter value at 25 degC
    # Ea : float
    #     activation energy for the parameter [J mol-1]
    # T : float
    #     temperature [deg C]
    # Tref : float
    #     measurement temperature [deg C]
    #
    # Returns:
    # -------
    # kt : float
    #     temperature dependence on parameter
    #
    # References:
    # -----------
    # * Medlyn et al. 2002, PCE, 25, 1167-1179.
    Tk <- T + DEG_TO_KELVIN
    TrefK <- Tref + DEG_TO_KELVIN

    return (k25 * exp(Ea * (T - Tref) / (RGAS * Tk * TrefK)))
}

calc_co2_compensation_point <- function(tleaf) {
    # CO2 compensation point in the absence of non-photorespiratory
    # respiration.
    #
    # Parameters:
    # ----------
    # tleaf : float
    #     leaf temperature (deg C)
    #
    # Returns:
    # -------
    # gamma_star : float
    #     CO2 compensation point in the abscence of mitochondrial respiration
    #     (umol mol-1)
    #
    # References:
    # -----------
    # * Bernacchi et al 2001 PCE 24: 253-260

    return (arrhenius(gamstar25, eag, tleaf, measurement_temp))
}

calculate_michaelis_menten <- function (tleaf) {
    # Effective Michaelis-Menten coefficent of Rubisco activity
    #
    # Parameters:
    # ----------
    # tleaf : float
    #     leaf temperature (deg C)
    #
    # Returns:
    # -------
    # Km : float
    #     Effective Michaelis-Menten constant for Rubisco catalytic activity
    #     (umol mol-1)
    #
    # References:
    # -----------
    # Rubisco kinetic parameter values are from:
    # * Bernacchi et al. (2001) PCE, 24, 253-259.
    # * Medlyn et al. (2002) PCE, 25, 1167-1179, see pg. 1170.

    # Michaelis-Menten coefficents for carboxylation by Rubisco
    Kc <- arrhenius(kc25, eac, tleaf, measurement_temp)

    # Michaelis-Menten coefficents for oxygenation by Rubisco
    Ko <- arrhenius(ko25, eao, tleaf, measurement_temp)

    # return effective Michaelis-Menten coefficient for CO2
    Km <- Kc * (1.0 + oi / Ko)

    return (Km)
}

peaked_arrhenius <- function(k25, Ea, T, Tref, deltaS, Hd) {
    # Temperature dependancy approximated by peaked Arrhenius eqn,
    # accounting for the rate of inhibition at higher temperatures.
    #
    # Parameters:
    # ----------
    # k25 : float
    #     rate parameter value at 25 degC
    # Ea : float
    #     activation energy for the parameter [J mol-1]
    # T : float
    #     temperature [deg C]
    # Tref : float
    #     measurement temperature [deg C]
    # deltaS : float
    #     entropy factor [J mol-1 K-1)
    # Hd : float
    #     describes rate of decrease about the optimum temp [J mol-1]
    #
    # Returns:
    # -------
    # kt : float
    #     temperature dependence on parameter
    #
    # References:
    # -----------
    # * Medlyn et al. 2002, PCE, 25, 1167-1179.

    Tk <- T + DEG_TO_KELVIN
    TrefK <- Tref + DEG_TO_KELVIN

    arg1 <- arrhenius(k25, Ea, T, Tref)
    arg2 <- 1.0 + exp((deltaS * TrefK - Hd) / (RGAS * TrefK))
    arg3 <- 1.0 + exp((deltaS * Tk - Hd) / (RGAS * Tk))

    return (arg1 * arg2 / arg3)
}

quad <- function(a, b, c, large, error) {
    # quadratic solution
    #
    # Parameters:
    # ----------
    # a : float
    #     co-efficient
    # b : float
    #     co-efficient
    # c : float
    #     co-efficient
    #
    # Returns:
    # -------
    # root : float

    # discriminant */
    d <- (b * b) - 4.0 * a * c
    if (d < 0.0) {
        # fprintf(stderr, "imaginary root found\n")
        # exit(EXIT_FAILURE)
        error <- TRUE
    }

    if (large) {
        if (a == 0.0 && b > 0.0) {
            root <- -c / b
        } else if (a == 0.0 && b == 0.0) {
            root <- 0.0
            if (c != 0.0) {
                # fprintf(stderr, "Can't solve quadratic\n")
                # exit(EXIT_FAILURE)
                error <- TRUE
            }
        } else {
            root <- (-b + sqrt(d)) / (2.0 * a)
        }

    } else {
        if (a == 0.0 && b > 0.0) {
            root <- -c / b
        } else if (a == 0.0 && b == 0.0) {
            root <- 0.0
            if (c != 0.0) {
                # fprintf(stderr, "Can't solve quadratic\n")
                # exit(EXIT_FAILURE)
                error <- TRUE
            }
        } else {
            root <- (-b - sqrt(d)) / (2.0 * a)
        }
    }
    return (c(root, error))
}

calc_lwp <- function(kl, transpiration) {
    lwp <- weighted_swp - (transpiration / kl)
    if (lwp < -20.0) {
        lwp <- -20.0
    }

    return (lwp)
}

calc_leaf_net_rad <- function(tair, vpd, sw_rad) {
    # extinction coefficient for diffuse radiation and black leaves
    # (m2 ground m2 leaf)
    kd <- 0.8

    # isothermal net LW radiaiton at top of canopy, assuming emissivity of
    # the canopy is 1
    Tk <- tair + DEG_TO_KELVIN

    # Isothermal net radiation (Leuning et al. 1995, Appendix)
    ea <- calc_sat_water_vapour_press(tair) - vpd

    # catch for AWAP diurnal stuff until I better connect VPD and Tair
    if (ea < 0.0) {
        ea <- 0.0
    }

    # apparent emissivity for a hemisphere radiating at air temp eqn D4
    emissivity_atm <- 0.642 * ((ea / Tk) ^ (1.0 / 7.0))

    net_lw_rad <- (1.0 - emissivity_atm) * SIGMA * (Tk ^ 4.0)
    rnet <- leaf_abs * sw_rad - net_lw_rad * kd * exp(-kd * lai)

    return (rnet)
}

calc_radiation_conductance <- function(tair) {
    # Returns the 'radiation conductance' at given temperature.
    # Units: mol m-2 s-1
    #
    # References:
    # -----------
    # * Formula from Ying-Ping's version of Maestro, cf. Wang and Leuning
    #   1998, Table 1,
    # * See also Jones (1992) p. 108.
    # * And documented in Medlyn 2007, equation A3, although I think there
    #   is a mistake. It should be Tk**3 not Tk**4, see W & L.

    Tk <- tair + DEG_TO_KELVIN
    grad <- 4.0 * SIGMA * (Tk * Tk * Tk) * LEAF_EMISSIVITY / (CP * MASS_AIR)

    return (grad)
}

calc_bdn_layer_forced_conduct <- function(tair, press, wind, leaf_width) {
    # Boundary layer conductance for heat - single sided, forced convection
    # (mol m-2 s-1)
    # See Leuning et al (1995) PC&E 18:1183-1200 Eqn E1

    Tk <- tair + DEG_TO_KELVIN
    cmolar <- press / (RGAS * Tk)
    gbh <- 0.003 * sqrt(wind / leaf_width) * cmolar

    return (gbh)
}

calc_bdn_layer_free_conduct <- function(tair, tleaf, press, leaf_width) {
    # Boundary layer conductance for heat - single sided, free convection
    # (mol m-2 s-1)
    # See Leuning et al (1995) PC&E 18:1183-1200 Eqns E3 & E4

    Tk <- tair + DEG_TO_KELVIN
    cmolar <- press / (RGAS * Tk)
    leaf_width_cubed <- leaf_width * leaf_width * leaf_width

    if (float_eq((tleaf - tair), 0.0)) {
        gbh <- 0.0
    } else {
        grashof <- 1.6E8 * abs(tleaf - tair) * leaf_width_cubed
        gbh <- 0.5 * DHEAT * (grashof ^ 0.25) / leaf_width * cmolar
    }

    return (gbh)
}

calc_canopy_evaporation <- function(rnet) {
    # Use Penman eqn to calculate evaporation flux at the potential rate for
    # canopy evaporation
    #
    # units = (mm/day)
    #
    # Parameters:
    # -----------
    # tair : float
    #     temperature [degC]
    # net_rad : float
    #     net radiation [W m-2]
    # press : float
    #     air pressure [kPa]
    #
    # Returns:
    # --------
    # pot_evap : float
    #     evaporation [mm d-1]

    ga <- canopy_boundary_layer_conduct(canht, wind, press, tair)
    lambda <- calc_latent_heat_of_vapourisation(tair)
    gamma <- calc_pyschrometric_constant(press, lambda)
    slope <- calc_slope_of_sat_vapour_pressure_curve(tair)

    arg1 <- slope * rnet + vpd * ga * CP * MASS_AIR
    arg2 <- slope + gamma
    LE <- arg1 / arg2 # W m-2
    pot_evap <- LE / lambda # mol H20 m-2 s-1

    return (pot_evap)
}

calc_soil_evaporation <- function(net_rad) {
    # Use Penman eqn to calculate top soil evaporation flux at the
    # potential rate.
    #
    # Soil evaporation is dependent upon soil wetness and plant cover. The net
    # radiation term is scaled for the canopy cover passed to this func and
    # the impact of soil wetness is accounted for in the wtfac term. As the
    # soil dries the evaporation component reduces significantly.
    #
    # Key assumptions from Ritchie...
    #
    # * When plant provides shade for the soil surface, evaporation will not
    # be the same as bare soil evaporation. Wind speed, net radiation and VPD
    # will all belowered in proportion to the canopy density. Following
    # Ritchie role ofwind, VPD are assumed to be negligible and are therefore
    # ignored.
    #
    # These assumptions are based on work with crops and whether this holds
    # for tree shading where the height from the soil to the base of the
    # crown is larger is questionable.
    #
    # units = (mm/day)
    #
    # References:
    # -----------
    # * Ritchie, 1972, Water Resources Research, 8, 1204-1213.
    #
    # Parameters:
    # -----------
    # tair : float
    #     temperature [degC]
    # net_rad : float
    #     net radiation [W m-2]
    # press : float
    #     air pressure [kPa]
    #
    # Returns:
    # --------
    # soil_evap : float
    #     soil evaporation [mm d-1]


    lambda <- calc_latent_heat_of_vapourisation(tair)
    gamma <- calc_pyschrometric_constant(press, lambda)
    slope <- calc_slope_of_sat_vapour_pressure_curve(tair)

    # mol H20 m-2 s-1
    soil_evap <- ((slope / (slope + gamma)) * net_rad) / lambda

    # Surface radiation is reduced by overstory LAI cover. This empirical
    # fit comes from Ritchie (1972) and is formed by a fit between the LAI
    # of 5 crops types and the fraction of observed net radiation at the
    # surface. Whilst the LAI does cover a large range, nominal 0–6, there
    # are only 12 measurements and only three from LAI > 3. So this might
    # not hold as well for a forest canopy?
    # Ritchie 1972, Water Resources Research, 8, 1204-1213.

    if (lai > 0.0)
        soil_evap <- soil_evap * exp(-0.398 * lai)

    # reduce soil evaporation if top soil is dry
    soil_evap <- soil_evap * wtfac_topsoil

    return (soil_evap)
}

SIGN <- function(a, b) {
    if (b >= 0) {
        return (abs(a))
    } else {
        return (-abs(a))
    }
}

nrerror <- function(error_text) {
    # Numerical Recipes standard error handler
	print("Numerical Recipes run-time error...\n")
	print(error_text)
	print("...now exiting to system...\n")
	exit(1)
}

calc_infiltration <- function(surface_water) {
    # Takes surface_water and distrubutes it among top layers. Assumes
    # total infilatration in timestep.

    add <- surface_water * MM_TO_M

    for (i in 1:n_layers) {
        ppt_gain[i] <<- 0.0
    }

    runoff <- 0.0
    for (i in 1:n_layers) {

        # determine the available pore space in current soil layer
        wdiff <- max(0.0, (porosity[i] - water_frac[i]) *
                          thickness[i] - water_gain[i] +
                          water_loss[i])

        # is the input of water greater than available space?
        if (add > wdiff) {
            # if so fill and subtract from input and move on to the next layer
            ppt_gain[i] <<- wdiff
            add <- add - wdiff
        } else {
            # otherwise infiltate all in the current layer
            ppt_gain[i] <<- add
            add <- 0.0
        }

        # if we have added all available water we are done
        if (add <= 0.0) {
            break
        }
    }

    # if after all of this we have some water left assume it is runoff
    if (add >  0.0) {
       runoff <- add * M_TO_MM
    } else {
       runoff <- 0.0
    }

    return (runoff)
}

calc_soil_temp_factor <- function(tsoil) {
    # Soil-temperature activity factor (A9). Fit to Parton's fig 2a
    #
    # Parameters:
    # -----------
    # tsoil : double
    #     soil temperature (deg C)
    #
    # Returns:
    # --------
    # tfac : double
    #     soil temperature factor [degC]

    if (tsoil > 0.0) {
        tfac <- max(0.0, 0.0326 + 0.00351 * (tsoil ^ 1.652) -
                        ((tsoil / 41.748) ^ 7.19))
    } else {
        # negative number cannot be raised to a fractional power
        # number would need to be complex
        tfac <- 0.0
    }

    return (tfac)
}

metafract <- function(lig2n) {
    # Calculate what fraction of the litter will be partitioned to the
    # metabolic pool which is given by the lignin:N ratio.
    #
    # Parameters:
    # -----------
    # lig2n : float
    #     lignin to N ratio
    #
    # Returns:
    # --------
    # metabolic fraction : float
    #     partitioned fraction to metabolic pool [must be positive]

    # Original implementation based on Parton et al.
    return (max(0.0, 0.85 - (0.018 * lig2n)))
}

calculate_nc_slope <- function(ncmax, ncmin) {
    # Returns N:C ratio of the mineral pool slope
    #
    # based on fig 4 of Parton et al 1993. Standard slow pool C:N is different
    # to the values in Parton. Bill quotes 12-20, whereas McMurtrie et al '01
    # use 10-40.
    #
    # Parameters
    # ----------
    # ncmax : float
    #     SOM pools maximum N:C
    # ncmin: float
    #     SOM pools minimum N:C
    #
    # Returns:
    # --------
    # value : float
    #     SOM pool N:C ratio

    arg1 <- ncmax - ncmin
    arg2 <- nmincrit - nmin0
    conv <- M2_AS_HA / G_AS_TONNES

    return (arg1 / arg2 * conv)
}

nc_limit <- function(cpool, npool, ncmin, ncmax) {
    # Release N to 'Inorgn' pool or fix N from 'Inorgn', in order to keep
    # the  N:C ratio of a litter pool within the range 'ncmin' to 'ncmax'.
    #
    # Parameters:
    # -----------
    # cpool : float
    #     various C pool (state)
    # npool : float
    #     various N pool (state)
    # ncmin : float
    #     maximum N:C ratio
    # ncmax : float
    #     minimum N:C ratio
    #
    # Returns:
    # --------
    # fix/rel : float
    #     amount of N to be added/released from the inorganic pool
    nmax <- cpool * ncmax
    nmin <- cpool * ncmin

    if (npool > nmax) {
        # release
        rel <- npool - nmax
        nlittrelease <<- nlittrelease + rel
        return (-rel)
    } else if (npool < nmin) {
        # fix
        fix <- nmin - npool
        nlittrelease <<- nlittrelease - fix
        return (fix)
    } else {
        return (0.0)
    }
}

nc_flux <- function(cflux, nflux, nc_ratio) {
    # Release N to Inorganic pool or fix N from the Inorganic pool in order
    # to normalise the N:C ratio of a net flux
    #
    # Parameters:
    # -----------
    # cflux : float
    #     C flux into SOM pool
    # nflux : float
    #     N flux into SOM pool
    # nc_ratio : float
    #     preferred N:C ratio
    #
    # Returns:
    #     fix : float
    #     Returns the amount of N required to be fixed

    return ((cflux * nc_ratio) - nflux)
}

calculate_pc_slope <- function () {
    # Returns P:C ratio of the mineral pool slope
    # Need to check back for good relationships Currently using olde NC relationship
    # based on Parton et al., 1993
    #
    # Parameters
    # ----------
    # pcmax : float
    # SOM pools maximum P:C
    # pcmin: float
    # SOM pools minimum P:C
    #
    # Returns:
    # --------
    # value : float
    # SOM pool P:C ratio

    arg1 <- pcmax - pcmin
    arg2 <- pmincrit - pmin0
    conv <- M2_AS_HA / G_AS_TONNES

    return (arg1 / arg2 * conv)
}

pc_limit <- function(cpool, ppool, pcmin, pcmax) {
    # Release P to 'Inorglabp' pool or fix P from 'Inorglabp', in order to keep
    # the  P:C ratio of a litter pool within the range 'pcmin' to 'pcmax'.
    #
    # Parameters:
    # -----------
    # cpool : float
    #     various C pool (state)
    # ppool : float
    #     various P pool (state)
    # pcmin : float
    #     min P:C ratio
    # pcmax : float
    #     max P:C ratio
    #
    # Returns:
    # --------
    # fix/rel : float
    #     amount of P to be added/released from the inorganic pool

    pmax <- cpool * pcmax
    pmin <- cpool * pcmin

    # fprintf(stderr, "ppool %f\n", ppool*100000)
    # fprintf(stderr, "pmax %f\n", pmax*100000)
    # fprintf(stderr, "pmin %f\n", pmin*100000)

    if (ppool > pmax) {
        # release
        rel <- ppool - pmax
        plittrelease <<- plittrelease + rel
        return (-rel)
    } else if (ppool < pmin) {
        # fix
        fix <- pmin - ppool
        plittrelease <<- plittrelease - fix
        return (fix)
    } else {
        return (0.0)
    }
}

pc_flux <- function(cflux, pflux, pc_ratio) {
    # Release P to Inorganic pool or fix P from the Inorganic pool in order
    # to normalise the P:C ratio of a net flux
    #
    # Parameters:
    # -----------
    # cflux : float
    #     C flux into SOM pool
    # pflux : float
    #     P flux into SOM pool
    # pc_ratio : float
    #     preferred P:C ratio
    #
    # Returns:
    #     fix : float
    #     Returns the amount of P required to be fixed

    return ((cflux * pc_ratio) - pflux)
}

calculate_growth_stress_limitation <- function () {
    # Calculate level of stress due to nitrogen, phosphorus or water
    # availability

    nc_opt <- 0.04
    pc_opt <- 0.004

    # N limitation based on leaf NC ratio
    if (shootnc < nf_min) {
        nlim <- 0.0
    } else if (shootnc < nc_opt && shootnc > nf_min) {
        nlim <- 1.0 - ((nc_opt - shootnc) / (nc_opt - nf_min))
    } else {
        nlim <- 1.0
    }

    # Limitation by nutrients or water. Water constraint is
    # implicit, in that, water stress results in an increase of root mass,
    # which are assumed to spread horizontally within the rooting zone.
    # So in effect, building additional root mass doesnt alleviate the
    # water limitation within the model. However, it does more
    # accurately reflect an increase in root C production at a water
    # limited site. This implementation is also consistent with other
    # approaches, e.g. LPJ. In fact I dont see much evidence for models
    # that have a flexible bucket depth. Minimum constraint is limited to
    # 0.1, following Zaehle et al. 2010 (supp), eqn 18.

    current_limitation <- max(0.1, min(nlim,wtfac_root))

    if(pcycle == TRUE) {
        # P limitation based on leaf PC ratio
        if (shootpc < pf_min) {
            plim <- 0.0
        } else if (shootpc < pc_opt && shootpc > pf_min) {
            plim <- 1.0 - ((pc_opt - shootpc) / (pc_opt - pf_min))
        } else {
            plim <- 1.0
        }
        nutrient_lim <- min(nlim, plim)
        current_limitation <- max(0.1, min(nutrient_lim, wtfac_root))
    }

    return (current_limitation)
}

calc_root_dist <- function(slope, root_biomass, surf_biomass,
                      rooted_layers, top_lyr_thickness,
                      root_reach) {
    # This function is used in the in the zbrent numerical algorithm to
    # figure out the slope of the rooting distribution for a given depth

    one <- (1.0 - exp(-slope * rooted_layers * top_lyr_thickness)) / slope
    two <- root_biomass / surf_biomass
    arg1 <- (1.0 - exp(-slope * root_reach)) / slope
    arg2 <- root_biomass / surf_biomass

    return (arg1 - arg2)
}