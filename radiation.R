calculate_solar_noon <- function(et, longitude) {
    # Calculation solar noon - De Pury & Farquhar, '97: eqn A16
    #
    # Reference:
    # ----------
    # * De Pury & Farquhar (1997) PCE, 20, 537-557.
    #
    # Returns:
    # ---------
    # t0 - solar noon (hours).
    # */
    # double t0, Ls;
    #
    # all international standard meridians are multiples of 15deg east/west of greenwich
    Ls <- round_to_value(longitude, 15)
    t0 <- 12.0 + (4.0 * (Ls - longitude) - et) / 60.0

    return (t0)
}

calculate_hour_angle <- function(t, t0) {
    # Calculation solar noon - De Pury & Farquhar, '97: eqn A15

    # Reference:
    # ----------
    # * De Pury & Farquhar (1997) PCE, 20, 537-557.
    #
    # Returns:
    # ---------
    # h - hour angle (radians).

    return (pi * (t - t0) / 12.0)

}

day_angle <- function(doy) {
    # Calculation of day angle - De Pury & Farquhar, '97: eqn A18

    # Reference:
    # ----------
    # * De Pury & Farquhar (1997) PCE, 20, 537-557.
    # * J. W. Spencer (1971). Fourier series representation of the position of
    #   the sun.
    #
    # Returns:
    # ---------
    # gamma - day angle in radians.
    #
    return (2.0 * pi * (doy - 1.0) / 365.0)
}

calculate_solar_declination <- function(doy, gamma) {
    # Solar Declination Angle is a function of day of year and is indepenent
    # of location, varying between 23deg45' to -23deg45'
    #
    # Arguments:
    # ----------
    # doy : int
    #     day of year, 1=jan 1
    # gamma : double
    #     fractional year (radians)
    #
    # Returns:
    # --------
    # dec: float
    #     Solar Declination Angle [radians]
    #
    # Reference:
    # ----------
    # * De Pury & Farquhar (1997) PCE, 20, 537-557.
    # * Leuning et al (1995) Plant, Cell and Environment, 18, 1183-1200.
    # * J. W. Spencer (1971). Fourier series representation of the position of
    #   the sun.
    #  Solar Declination Angle (radians) A14 - De Pury & Farquhar
    decl <- -23.4 * (pi / 180) * cos(2.0 * pi * (doy + 10) / 365)

    return (decl)
}

calculate_eqn_of_time <- function(gamma) {
    # Equation of time - correction for the difference btw solar time
    # and the clock time.
    #
    # Arguments:
    # ----------
    # doy : int
    #     day of year
    # gamma : double
    #     fractional year (radians)
    #
    # References:
    # -----------
    # * De Pury & Farquhar (1997) PCE, 20, 537-557.
    # * Campbell, G. S. and Norman, J. M. (1998) Introduction to environmental
    #   biophysics. Pg 169.
    # * J. W. Spencer (1971). Fourier series representation of the position of
    #   the sun.
    # * Hughes, David W.; Yallop, B. D.; Hohenkerk, C. Y. (1989),
    #   "The Equation of Time", Monthly Notices of the Royal Astronomical
    #   Society 238: 1529–1535
    #

    # ** from Spencer '71. This better matches the de Pury worked example (pg 554)
    # ** The de Pury version is this essentially with the 229.18 already applied
    # ** It probably doesn't matter which is used, but there is some rounding
    # ** error below (radians)

    et <- 0.000075 + 0.001868 * cos(gamma) - 0.032077 * sin(gamma) - 0.014615 * cos(2.0 * gamma) - 0.04089 * sin(2.0 * gamma)

    # radians to minutes
    et <-  et * 229.18

    # radians to hours
    #et *= 24.0 / (2.0 * pi);

    # minutes - de Pury and Farquhar, 1997 - A17
    #et = (0.017 + 0.4281 * cos(gamma) - 7.351 * sin(gamma) - 3.349 * cos(2.0 * gamma) - 9.731 * sin(gamma));

    return (et)
}



calc_extra_terrestrial_rad <- function(doy, cos_zenith) {
    # Solar radiation incident outside the earth's atmosphere, e.g.
    # extra-terrestrial radiation. The value varies a little with the earths
    # orbit.
    #
    # Using formula from Spitters not Leuning!
    #
    # Arguments:
    # ----------
    # doy : double
    #     day of year
    # cos_zenith : double
    #     cosine of zenith angle (radians)
    #
    # Returns:
    # --------
    # So : float
    #     solar radiation normal to the sun's bean outside the Earth's atmosphere
    #     (J m-2 s-1)
    #
    # Reference:
    # ----------
    # * Spitters et al. (1986) AFM, 38, 217-229, equation 1.

    # Solar constant (J m-2 s-1)
    Sc <- 1370.0

    if (cos_zenith > 0.0) {
        # remember sin_beta = cos_zenith; trig funcs are cofuncs of each other
        # sin(x) = cos(90-x) and cos(x) = sin(90-x).
        So <- Sc * (1.0 + 0.033 * cos(doy / 365.0 * 2.0 * pi)) * cos_zenith
    } else {
        So <- 0.0
    }

    return (So)

}


estimate_clearness <- function(sw_rad, So) {
    # estimate atmospheric transmisivity - the amount of diffuse radiation
    # is a function of the amount of haze and/or clouds in the sky. Estimate
    # a proxy for this, i.e. the ratio between global solar radiation on a
    # horizontal surface at the ground and the extraterrestrial solar
    # radiation

    # catch possible divide by zero when zenith = 90.
    if (So <= 0.0) {
        tau <- 0.0
    } else {
        tau <- sw_rad / So
    }

    if (tau > 1.0) {
        tau <- 1.0
    } else if (tau < 0.0) {
        tau <- 0.0
    }

    return (tau)
}