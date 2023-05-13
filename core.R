run_sim <- function () {
    if (deciduous_model) {
        # Are we reading in last years average growing season?
        if (float_eq(avg_alleaf, 0.0) &&
            float_eq(avg_alstem, 0.0) &&
            float_eq(avg_albranch, 0.0) &&
            float_eq(avg_alleaf, 0.0) &&
            float_eq(avg_alroot, 0.0) &&
            float_eq(avg_alcroot, 0.0)) {
            npitfac <- 0.0
            calc_carbon_allocation_fracs(npitfac)
        } else {
            alleaf <<- avg_alleaf
            alstem <<- avg_alstem
            albranch <<- avg_albranch
            alroot <<- avg_alroot
            alcroot <<- avg_alcroot
        }
        allocate_stored_cnp()
    }
    if (print_options == SUBDAILY && spin_up == FALSE) {
        # open the 30 min outputs file and the daily output files
        if (output_ascii) {
            ofp_sd <<- write_output_subdaily_header()
            ofp <<- write_output_header()
        } else {
            stop("Nothing implemented for sub-daily binary\n")
        }
    } else if (print_options == DAILY && spin_up == FALSE) {
        # Daily outputs
        if (output_ascii) {
            ofp <<- write_output_header()
        } else {
            stop("Nothing implemented for daily binary\n")
            ofp_hdr <<- write_output_header()
        }
    } else if (print_options == END && spin_up == FALSE) {
        # Final state + param file
        # open_output_file(c, out_param_fname, &(ofp))
    }

    # Window size = root lifespan in days...
    # For deciduous species window size is set as the length of the
    # growing season in the main part of the code

    window_size <- as.integer(1.0 / prdecay * NDAYS_IN_YR)
    sma_obj <<- list(0, 0, window_size, rep(0, window_size), 0)
    names(sma_obj) <<- c("sma", "sum", "period", "values", "lv")
    if (prev_sma > -900) {
        for (i in 1:window_size) {
            SMA_ADD(prev_sma)
        }
    }
    # Set up SMA
    # - If we don't have any information about the N & water limitation, i.e.
    #   as would be the case with spin-up, assume that there is no limitation
    #   to begin with.
    if (prev_sma < -900)
        prev_sma <<- 1.0

    # Params are defined in per year, needs to be per day. Important this is
    # done here as rate constants elsewhere in the code are assumed to be in
    # units of days not years

    correct_rate_constants(FALSE)
    day_end_calculations(-99, TRUE)

    if (sub_daily) {
        initialise_soils_sub_daily()
    } else {
        initialise_soils_day()
    }

    if (water_balance == HYDRAULICS) {
        # Update the soil water storage
        root_zone_total <- 0.0
        for (i in 1:n_layers) {

            # water content of soil layer (m)
            water_content <- water_frac[i] * thickness[i]

            # update old GDAY effective two-layer buckets
            # - this is just for outputting, these aren't used.
            if (i == 1) {
                pawater_topsoil <<- water_content * M_TO_MM
            } else {
                root_zone_total <- root_zone_total + water_content * M_TO_MM
            }
        }
        pawater_root <<- root_zone_total
    } else {
        pawater_root <<- wcapac_root
        pawater_topsoil <<- wcapac_topsoil
    }

    if (fixed_lai) {
        lai <<- fix_lai
    } else {
        lai <<- max(0.01, (psla * M2_AS_HA / KG_AS_TONNES / cfracts * shoot))
    }

    if (disturbance) {
        res__ <- figure_out_years_with_disturbances()
        disturbance_yrs <- res__[[1]]
        num_disturbance_yrs <- res__[[2]]
    }


    # ======================
    #   Y E A R    L O O P
    # ======================
    day_idx <<- 1
    hour_idx <<- 1



    for (nyr in 1:num_years) {

        if (sub_daily) {
            year <- ma$year[hour_idx]
        } else {
            year <- ma$year[day_idx]
        }
        if (is_leap_year(year)) {
            num_days <<- 366
        } else {
            num_days <<- 365
        }

        calc_warmest_quarter_temp(year)

        calculate_daylength(num_days, latitude)

        if (deciduous_model) {
            phenology()

            # Change window size to length of growing season
            sma_obj <<- list(0, 0, growing_seas_len, rep(0, growing_seas_len), 0)
            names(sma_obj) <<- c("sma", "sum", "period", "values", "lv")
            if (prev_sma > -900) {
                for (i in 1:growing_seas_len) {
                    SMA_ADD(prev_sma)
                }
            }

            zero_stuff()
        }
        # ===================
        #   D A Y   L O O P
        # ===================
        for (doy in 1:num_days) {

            #if (year == 2001 && doy+1 == 230) {
            #   pdebug = TRUE;
            #}


            if (!sub_daily) {
                unpack_met_data(dummy, day_length[doy])
            }
            res__ <- calculate_litterfall(doy)
            fdecay <- res__[1]
            rdecay <- res__[2]

            if (disturbance && disturbance_doy == doy + 1) {
                # Fire Disturbance?
                fire_found <- FALSE
                fire_found <- check_for_fire(year, disturbance_yrs, num_disturbance_yrs)

                if (fire_found) {
                    fire()

                    # This will only work for evergreen, but that is fine
                    # this should be removed after KSCO is done

                    sma_obj <<- list(0, 0, window_size, rep(0, window_size), 0)
                    names(sma_obj) <<- c("sma", "sum", "period", "values", "lv")
                    if (prev_sma > -900) {
                        for (i in 1:window_size) {
                            SMA_ADD(prev_sma)
                        }
                    }
                }
            } else if (hurricane &&
                hurricane_yr == year &&
                hurricane_doy == doy) {

                # Hurricane?
                hurricane_f()
            }


            calc_day_growth(day_length[doy], doy, fdecay, rdecay)

            #printf("%d %f %f\n", doy, gpp*100, lai);
            calculate_csoil_flows(tsoil, doy)
            calculate_nsoil_flows(doy)

            if (pcycle == TRUE) {
                calculate_psoil_flows(doy)
            }

            # update stress SMA
            if (deciduous_model && leaf_out_days[doy] > 0.0) {
                # Allocation is annually for deciduous "tree" model, but we
                # need to keep a check on stresses during the growing season
                # and the LAI figure out limitations during leaf growth period.
                # This also applies for deciduous grasses, need to do the
                # growth stress calc for grasses here too.

                current_limitation <- calculate_growth_stress_limitation()
                SMA_ADD(current_limitation)
                prev_sma <<- sma_obj$sma
            } else if (deciduous_model == FALSE) {
                current_limitation <- calculate_growth_stress_limitation()
                SMA_ADD(current_limitation)
                prev_sma <<- sma_obj$sma
            }

            # if grazing took place need to reset "stress" running mean
            # calculation for grasses

            if (grazing == 2 && disturbance_doy == doy + 1) {
                sma_obj <<- list(0, 0, growing_seas_len, rep(0, growing_seas_len), 0)
                names(sma_obj) <<- c("sma", "sum", "period", "values", "lv")
            }

            # Turn off all N calculations
            if (ncycle == FALSE)
                reset_all_n_pools_and_fluxes()

            # Turn off all P calculations
            if (pcycle == FALSE)
                reset_all_p_pools_and_fluxes()

            # calculate C:N ratios and increment annual flux sum
            day_end_calculations(num_days, FALSE)

            if (print_options == SUBDAILY && spin_up == FALSE) {
                write_daily_outputs_ascii(year, doy)
            } else if (print_options == DAILY && spin_up == FALSE) {
                if(output_ascii) {
                    write_daily_outputs_ascii(year, doy)
                } else {
                    stop("Nothing implemented for daily binary\n")
                    # write_daily_outputs_binary(year, doy + 1)
                }
            }
            day_idx <<- day_idx + 1


            #printf("%d %d %f", (int)year, doy, water_frac[0] * thickness[0] * M_TO_MM);
            #printf("%d %d %f", (int)year, doy, water_frac[0]);
            #for (i = 1; i < n_layers; i++) {
            #
            #    #printf(" %f", water_frac[i] * thickness[i] * M_TO_MM);
            #    printf(" %f", water_frac[i]);
            #
            #}
            #printf("\n");
            #printf("%d %d %lf %lf %lf\n", (int)year, doy, saved_swp, wtfac_root, gpp*100);

            #printf("%d %d %lf %lf %lf %lf\n", (int)year, doy, gpp*100, transpiration, wtfac_root, saved_swp);
            #printf("%d %d %lf %lf %lf\n", (int)year, doy, gpp*100, transpiration, wtfac_root);


            # =======================
            #   E N D   O F   D A Y
            # =======================
        }


        # Allocate stored C,N and P for the following year
        if (deciduous_model) {
            calculate_average_alloc_fractions(growing_seas_len)
            allocate_stored_cnp()
        }

        # Adjust rooting distribution at the end of the year to account for
        # growth of new roots. It is debatable when this should be done. I've
        # picked the year end for computation reasons and probably because
        # plants wouldn't do this as dynamcially as on a daily basis. Probably
        if (water_balance == HYDRAULICS) {
            update_roots()
        }
    }
    # =========================
    #   E N D   O F   Y E A R
    # =========================
    correct_rate_constants(TRUE)

    if (print_options == DAILY && spin_up == FALSE) {
        write.csv(ofp, file = out_fname, row.names = FALSE)
    }

    if (print_options == END && spin_up == FALSE) {
        # write_final_state()
    }
}

spin_up_pools <- function () {

    #Spin up model plant & soil pools to equilibrium.

    # - Examine sequences of 50 years and check if C pools are changing
    #   by more than 0.005 units per 1000 yrs.
    #
    # References:
    # ----------
    # Adapted from...
    # * Murty, D and McMurtrie, R. E. (2000) Ecological Modelling, 134,
    #   185-205, specifically page 196.

    tol_c <- 5E-03
    tol_n <- 5E-03
    tol_p <- 5E-03
    prev_plantc <- 99999.9
    prev_soilc <- 99999.9
    prev_plantn <- 99999.9
    prev_soiln <- 99999.9
    prev_plantp <- 99999.9
    prev_soilp <- 99999.9

    # check for convergences in units of kg/m2
    conv <- TONNES_HA_2_KG_M2

    # If we are prescribing disturbance, first allow the forest to establish
    if (disturbance) {
        cntrl_flag <- disturbance
        disturbance <<- FALSE
        #  200 years (50 yrs x 4 cycles)
        for (i in 1:4) {
            run_sim() # run GDAY
        }
        disturbance <<- cntrl_flag
    }

    repeat {
        if (abs((prev_plantc) - (plantc)) < tol_c &&
            abs((prev_soilc) - (soilc)) < tol_c &&
            abs((prev_plantn) - (plantn)) < tol_n &&
            abs((prev_soiln) - (soiln)) < tol_n &&
            abs((prev_plantp) - (plantp)) < tol_p &&
            abs((prev_soilp) - (soilp)) < tol_p) {
            break
        } else {
            prev_plantc <- plantc
            prev_soilc <- soilc
            prev_plantn <- plantn
            prev_soiln <- soiln
            prev_plantp <- plantp
            prev_soilp <- soilp

            for (i in 1:20) {
                run_sim()
                # print(i)
            }
            if (pcycle) {
                # Have we reached a steady state?
                print(sprintf("Spinup: Leaf C - %f, Leaf CN - %f, Leaf CP - %f, Wood C - %f, Leaf N - %f, Leaf P - %f, Soil P - %f, LAI - %f\n",
                        shoot, shoot / shootn, shoot / shootp, stem, shootn, shootp, soilp, lai))
            } else {
                # Have we reached a steady state?
                print(sprintf("Spinup: Leaf C - %f, Leaf NC - %f,  Wood C - %f, LAI - %f\n",
                      shoot, 1.0 / shootnc, stem,  lai))
            }
        }
    }
}