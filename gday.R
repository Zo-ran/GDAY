gday <- function() {
    # potentially allocating 1 extra spot, but will be fine as we always
    # index by num_days
    day_length <<- vector("numeric", 366)

    # House keeping!
    if (water_balance == HYDRAULICS && sub_daily == FALSE) {
        stop("You can't run the hydraulics model with daily flag")
    }
    if (water_balance == HYDRAULICS) {
        allocate_numerical_libs_stuff()
        initialise_roots()
        setup_hydraulics_arrays()
    }

    ma <<- read.csv(met_fname, header=TRUE)
    num_years <<- max(ma$year) - min(ma$year) + 1
    total_num_days <<- dim(ma)[1]

    if (sub_daily) {
        fill_up_solar_arrays()
    }

    if (spin_up) {
        spin_up_pools()
    } else {
        run_sim()
    }
}
