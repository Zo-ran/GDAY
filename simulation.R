rm(list = ls())
source("constants.R")
source("control.R")
source("params.R")
source("fluxes.R")
source("state.R")
source("nrutil.R")
source("utilities.R")
source("radiation.R")
source("functions.R")
source("gday.R")
source("core.R")
# potentially allocating 1 extra spot, but will be fine as we always
# index by num_days

main <- function(SPIN_UP=TRUE, POST_INDUST=TRUE) {
    met_dir <- 'met_data/'
    out_dir <- 'outputs/'
    if (SPIN_UP) {
        met_fname <<- paste0(met_dir, 'EUC_met_data_amb_var_co2.csv')
        out_fname <<- paste0(out_dir, 'FACE_EUC_model_spinup_equilib.out')
        print_options <<- END
        spin_up <<- TRUE
        gday()
    }
    if (POST_INDUST) {
        met_fname <<- paste0(met_dir, 'EUC_met_data_amb_var_co2.csv')
        out_fname <<- paste0(out_dir, 'EUC_amb_equilib.csv')
        print_options <<- DAILY
        spin_up <<- FALSE
        gday()
    }
}


main(SPIN_UP=TRUE, POST_INDUST=TRUE)