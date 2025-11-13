#' ####################################################################### #
#' PROJECT: [Himalayan Summit Success - BRMS Models]
#'  DEPENDENCIES:
#'  - ModelData.csv in Exports directory
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PREAMBLE ================================================================
rm(list = ls()) # some may not like it, but it helps my workflow

## Packages ---------------------------------------------------------------
packages <- list(
    core = c("tidyverse"),
    viz = c("ggplot2"),
    stats = c("brms", "tidybayes"),
    utils = c("pbapply", "lubridate")
)

install_if_missing <- function(pkg_list) {
    new_pkgs <- unlist(pkg_list)[!unlist(pkg_list) %in% installed.packages()[, "Package"]]
    if (length(new_pkgs)) install.packages(new_pkgs, repos = "http://cran.us.r-project.org")
    invisible(lapply(unlist(pkg_list), library, character.only = TRUE))
}

# Install and load all packages
invisible(install_if_missing(packages))

## Directories ------------------------------------------------------------
### Define directories in relation to project directory
Dir.Base <- getwd() # identifying the current directory
Dir.Data <- file.path(Dir.Base, "Data") # folder path for data
Dir.Covariates <- file.path(Dir.Base, "Covariates") # folder path for covariates
Dir.Exports <- file.path(Dir.Base, "Exports") # folder path for exports
### create directories, if they don't exist yet
Dirs <- sapply(
    c(Dir.Data, Dir.Covariates, Dir.Exports),
    function(x) if (!dir.exists(x)) dir.create(x)
)
rm(Dirs) # removing temporary variable
ModelData_df <- read.csv(file.path(Dir.Exports, "ModelData.csv"))


# Optional: remove peaks with no variation
ModelData_df_nonconst <- ModelData_df |>
    dplyr::group_by(PKNAME) |>
    dplyr::filter(dplyr::n_distinct(Mortality_all) > 1) |>
    dplyr::ungroup()

# Fit the model
mortality_model <- brm(
    bf(
        Mortality_all ~ YEAR + (YEAR | PKNAME), # random intercepts & slopes per peak
        phi ~ 1, # constant precision
        zoi ~ YEAR + (YEAR | PKNAME) # model zero/one inflation similarly
    ),
    data = ModelData_df_nonconst,
    family = zero_one_inflated_beta(),
    chains = 4,
    cores = 4,
    iter = 1e4,
    warmup = 2e3,
    inits = "0",
    seed = 123,
    control = list(adapt_delta = 0.95, max_treedepth = 15)
)

library(tidyverse)

# Build prediction data frame
# newdata <- data.frame(
#     PKNAME = rep(ModelData_df_nonconst$PKNAME, each = length(min(ModelData_df_nonconst$YEAR):max(ModelData_df_nonconst$YEAR))),
#     YEAR = rep(min(ModelData_df_nonconst$YEAR):max(ModelData_df_nonconst$YEAR), length(ModelData_df_nonconst$PKNAME))
# )


newdata <- ModelData_df_nonconst %>%
    group_by(PKNAME) %>%
    summarise(YEAR = min(YEAR):max(YEAR))

# Conditional mortality (given >0)
mortality_draws <- mortality_model %>%
    add_epred_draws(newdata = newdata, re_formula = NULL, ndraws = 1e2) %>%
    rename(mort_rate = .epred)

# Probability of any mortality (zero/one hurdle)
hurdle_draws <- mortality_model %>%
    add_epred_draws(newdata = newdata, dpar = "zoi", re_formula = NULL, ndraws = 1e2) %>%
    rename(prob_mort = zoi)

combined_draws <- left_join(
    mortality_draws %>% select(PKNAME, YEAR, .draw, mort_rate),
    hurdle_draws %>% select(PKNAME, YEAR, .draw, prob_mort),
    by = c("PKNAME", "YEAR", ".draw")
)

# Find max values per peak for scaling
peak_max <- combined_draws %>%
    group_by(PKNAME) %>%
    summarise(max_mort_rate = max(mort_rate), max_prob = max(prob_mort))

# Join to combined_draws
combined_draws <- combined_draws %>%
    left_join(peak_max, by = "PKNAME") %>%
    mutate(mort_rate_scaled = mort_rate / max_mort_rate * max_prob) # scale to 0–max prob


ggplot(combined_draws, aes(x = YEAR)) +
    # Primary y-axis: probability of any mortality
    stat_lineribbon(
        aes(y = prob_mort),
        .width = c(.95, .80, .50),
        alpha = 0.4,
        fill = "steelblue"
    ) +
    # Secondary y-axis: conditional mortality (scaled)
    stat_lineribbon(
        aes(y = mort_rate_scaled),
        .width = c(.95, .80, .50),
        alpha = 0.3,
        fill = "darkred"
    ) +
    # Secondary axis transformation
    scale_y_continuous(
        name = "Probability of Any Mortality",
        sec.axis = sec_axis(
            ~ . * max(combined_draws$mort_rate) / max(combined_draws$prob_mort),
            name = "Conditional Mortality (given >0)"
        )
    ) +
    facet_wrap(~PKNAME, scales = "free_y") +
    labs(
        x = "Year",
        title = "Predicted Mortality per Peak Over Time",
        subtitle = "Blue = Probability of any mortality | Red = Conditional mortality given >0"
    ) +
    theme_minimal()


ggplot(combined_draws, aes(x = YEAR)) +
    # Probability of any mortality (primary, blue)
    stat_lineribbon(
        aes(y = prob_mort),
        .width = c(.95, .80, .50),
        alpha = 0.4,
        fill = "steelblue"
    ) +
    # Conditional mortality given >0 (red)
    stat_lineribbon(
        aes(y = mort_rate),
        .width = c(.95, .80, .50),
        alpha = 0.3,
        fill = "darkred"
    ) +
    facet_wrap(~PKNAME, scales = "free_y") +
    scale_y_continuous(
        name = "Rate / Probability",
        limits = c(0, 1)
    ) +
    labs(
        x = "Year",
        title = "Predicted Mortality per Peak Over Time",
        subtitle = "Blue = Probability of any mortality | Red = Conditional mortality given >0"
    ) +
    theme_minimal()
