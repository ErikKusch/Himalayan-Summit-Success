#' ####################################################################### #
#' PROJECT: [Himalayan Summit Success - Climate Data]
#' CONTENTS:
#'  A. DATA DOWNLOAD: Download and process climate data
#'  B. DATA EXTRACTION: Extract and prepare data for analysis
#'  C. ANALYSIS & VISUALIZATIONS: Statistical analysis and visualization
#'  DEPENDENCIES:
#'  - None
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PREAMBLE ================================================================
rm(list = ls()) # some may not like it, but it helps my workflow

## Packages ---------------------------------------------------------------
packages <- list(
    core = c("readr", "dplyr", "tidyr", "terra", "sf", "sp"),
    viz = c("ggplot2", "viridis", "cowplot", "mapview", "ggrepel", "tidyterra", "ggpubr"),
    spatial = c("rnaturalearth", "rnaturalearthdata"),
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

## KrigR Setup --------------------------------------------------------
if ("KrigR" %in% rownames(installed.packages()) == FALSE) {
    Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
    devtools::install_github("ErikKusch/KrigR")
}
library(KrigR)

## API Credentials --------------------------------------------------------
try(source("X - PersonalSettings.R"))
if (!exists("API_Key") | !exists("API_User")) { # CS API check: if CDS API credentials have not been specified elsewhere
    API_User <- readline(prompt = "Please enter your Climate Data Store API user number and hit ENTER.")
    API_Key <- readline(prompt = "Please enter your Climate Data Store API key number and hit ENTER.")
} # end of CDS API check
# NUMBER OF CORES
if (!exists("numberOfCores")) { # Core check: if number of cores for parallel processing has not been set yet
    numberOfCores <- as.numeric(readline(prompt = paste("How many cores do you want to allocate to these processes? Your machine has", parallel::detectCores())))
} # end of Core check

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

## Project / Analysis Settings --------------------------------------------
### Climate Variables Configuration -------
climate_vars <- data.frame(
    name = c("2m_temperature", "snow_cover", "snow_depth", "snowfall", "windspeed"),
    display = c("Air Temperature [K]", "Snow Cover [%]", "Snow Depth [m]", "Snowfall [m]", "Windspeed [m/s]"),
    subset = TRUE,
    stringsAsFactors = FALSE
)

# For backward compatibility
VNames_vec <- climate_vars$display
SubsetVariables_vec <- climate_vars$name[climate_vars$subset]

# A. DATA LOADING =========================================================
message("#### Loading Data from Disk ##################################")

## Expedition DATA -------------------------------------------------
Expeditions_df <- read.csv(file.path(Dir.Data, "CleanedExpeditions.csv"))

## Summits as Spatial Objects --------------------------------------
summits_df <- read_csv(file.path(Dir.Data, "selected_peaks_coordinates_counts.csv")) # load positions and names of summits
summits_df <- summits_df[!duplicated(summits_df$ID), ] # check for duplicates and eventually delete them
## making sp object of summits
summits_sp <- summits_df
coordinates(summits_sp) <- ~ LON + LAT
proj4string(summits_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
## making into sf object
summits_sf <- st_as_sf(summits_sp)
eightks_sf <- summits_sf[summits_sf$HEIGHTM >= 8000, ] # used for kriging demonstration
sevenks_sf <- summits_sf[summits_sf$HEIGHTM >= 7000, ] # used for time series extraction / analyses

## ERA5-Land Data Download -----------------------------------------
message("#### Download Data from CDS ##################################")
Variables_vec <- c(
    "2m_temperature",
    "skin_temperature",
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
    "snow_cover",
    "snow_density",
    "snow_depth",
    "snow_depth_water_equivalent",
    "snow_evaporation",
    "snowfall",
    "snowmelt",
    "temperature_of_snow_layer"
)
Years_vec <- 1951:2021

CDSData_ls <- lapply(Variables_vec, FUN = function(Var_Iter) {
    message(Var_Iter)
    if (file.exists(file.path(Dir.Data, paste0(Var_Iter, ".nc")))) {
        print("Already prepared")
        return(rast(file.path(Dir.Data, paste0(Var_Iter, ".nc"))))
    }

    Raw_rast <- CDownloadS(
        Variable = Var_Iter,
        CumulVar = ifelse(Var_Iter %in% "snowfall", TRUE, FALSE),
        DataSet = "reanalysis-era5-land",
        Type = NA,
        DateStart = paste0(Years_vec[1], "-01-01 00:00"),
        DateStop = paste0(tail(Years_vec, 1), "-12-31 23:00"),
        TResolution = "hour",
        TStep = 1,
        Extent = ext(c(78.83, 89.33, 25.143, 31.543)),
        Dir = Dir.Data,
        FileName = paste0(Var_Iter, "_Raw"),
        API_User = API_User,
        API_Key = API_Key,
        Cores = numberOfCores
    )

    Var_data <- KrigR:::Temporal.Aggr(Raw_rast,
        BaseResolution = "hour",
        BaseStep = 1,
        TResolution = "day",
        TStep = 1,
        FUN = mean,
        Cores = numberOfCores,
        TZone = "UTC"
    )
    terra::writeCDF(Var_data, file = file.path(Dir.Data, paste0(Var_Iter, ".nc")))
    Var_data
})
names(CDSData_ls) <- Variables_vec

## Windspeed Calculation -------------------------------------------
### Hourly -------
# if (!file.exists(file.path(Dir.Data, "windspeed_RAW.nc"))) {
#     RAW_ls <- lapply(c("10m_u_component_of_wind", "10m_v_component_of_wind"),
#         FUN = function(i) rast(file.path(Dir.Data, paste0(i, "_RAW.nc")))
#     )

#     Indices <- ceiling((1:terra::nlyr(RAW_ls[[1]])) / 2e4)
#     ru_ls <- terra::split(x = RAW_ls[[1]], f = Indices)
#     rv_ls <- terra::split(x = RAW_ls[[2]], f = Indices)

#     ret_ls <- pblapply(1:length(ru_ls), FUN = function(BASE_iter) {
#         sqrt(abs(ru_ls[[BASE_iter]])^2 + abs(rv_ls[[BASE_iter]])^2)
#     })
#     Windspeed <- do.call(c, ret_ls)
#     terraOptions(memmax = 10)
#     writeCDF(Windspeed[[1:nlyr(Windspeed) / 2]],
#         filename = file.path(Dir.Data, "windspeed_RAW.nc")
#     )
#     terraOptions(memfrac = .9)
# }

### Daily -------
if (!file.exists(file.path(Dir.Data, "windspeed.nc"))) {
    CDSData_ls$windspeed <- sqrt(abs(CDSData_ls$`10m_u_component_of_wind`)^2 +
        abs(CDSData_ls$`10m_v_component_of_wind`)^2)
    writeCDF(CDSData_ls$windspeed,
        filename = file.path(Dir.Data, "windspeed.nc")
    )
} else {
    CDSData_ls$windspeed <- rast(file.path(Dir.Data, "windspeed.nc"))
}
varnames(CDSData_ls$windspeed) <- "windspeed"

## Data Subsetting -------------------------------------------------
CDSData_ls <- CDSData_ls[SubsetVariables_vec]

## Season Definition -----------------------------------------------
message("#### Splitting Data into Seasons #############################")
CDS_PreMonsoon_ls <- lapply(CDSData_ls, FUN = function(x) {
    x[[format(terra::time(CDSData_ls[[1]]), "%m") %in% c("03", "04", "05")]] # March - May
})

CDS_PostMonsoon_ls <- lapply(CDSData_ls, FUN = function(x) {
    x[[format(terra::time(CDSData_ls[[1]]), "%m") %in% c("09", "10", "11")]] # September - November
})

# B. DATA EXTRACTION ======================================================
message("#### Extracting Data for Analyses ############################")

## Time-Series Extraction for Peaks --------------------------------
print("Creating time series for peaks...")
if (file.exists(file.path(Dir.Exports, "peaks_time_series.csv"))) {
    peaks_ts_df <- read.csv(file.path(Dir.Exports, "peaks_time_series.csv"))
} else {
    ### Extract daily time series -------
    print("                               ... extracting data")
    peak_ts_ls <- pblapply(names(CDSData_ls), FUN = function(VarName) {
        VarIter <- CDSData_ls[[VarName]]
        values_mat <- t(terra::extract(
            VarIter,
            sevenks_sf,
            method = "bilinear",
            fun = "mean"
        ))
        values_df <- as.data.frame(values_mat[-1, ])
        colnames(values_df) <- sevenks_sf$PKNAME
        values_df$Date <- time(VarIter)
        values_df$Variable <- VarName
        values_df$Month <- as.numeric(format(values_df$Date, "%m"))
        values_df$Season <- "Out Of Season"
        values_df$Season[values_df$Month %in% 3:5] <- "Pre"
        values_df$Season[values_df$Month %in% 9:11] <- "Post"
        values_df
    })
    names(peak_ts_ls) <- names(CDSData_ls)

    ### Calculate seasonal bounds -------
    print("                               ... calculating seasonal bounds")
    seasonal_bounds <- do.call(rbind, pblapply(names(CDSData_ls), FUN = function(VarName) {
        var_df <- peak_ts_ls[[VarName]]
        do.call(rbind, lapply(sevenks_sf$PKNAME, FUN = function(peak) {
            data.frame(
                Variable = VarName,
                Peak = peak,
                Season = rep(c("Pre", "Post"), each = 2),
                Bound = rep(c("Lower", "Upper"), 2),
                Value = c(
                    quantile(var_df[var_df$Season == "Pre", peak], probs = c(0.05, 0.95), na.rm = TRUE),
                    quantile(var_df[var_df$Season == "Post", peak], probs = c(0.05, 0.95), na.rm = TRUE)
                )
            )
        }))
    }))

    ### Create final time series dataframe -------
    print("                              ... making final data frame")
    peaks_ts_df <- do.call(rbind, pblapply(names(CDSData_ls), FUN = function(VarName) {
        var_df <- peak_ts_ls[[VarName]]
        var_long <- tidyr::pivot_longer(
            var_df,
            cols = sevenks_sf$PKNAME,
            names_to = "PeakID",
            values_to = "Value"
        )
        var_long <- var_long[var_long$Season != "Out Of Season", ]

        # Add extreme indicators
        var_long$Extreme <- "Normal"
        var_long$ExtremeRatioHigh <- var_long$ExtremeRatioLow <- NA
        for (peak in unique(var_long$PeakID)) {
            for (season in c("Pre", "Post")) {
                bounds <- seasonal_bounds[
                    seasonal_bounds$Variable == VarName &
                        seasonal_bounds$Peak == peak &
                        seasonal_bounds$Season == season,
                ]
                mask <- var_long$PeakID == peak & var_long$Season == season
                var_long$Extreme[mask & var_long$Value < bounds$Value[bounds$Bound == "Lower"]] <- "LOW"
                var_long$Extreme[mask & var_long$Value > bounds$Value[bounds$Bound == "Upper"]] <- "HIGH"
                var_long$ExtremeRatioHigh <- var_long$Value / bounds$Value[bounds$Bound == "Upper"]
                var_long$ExtremeRatioLow <- var_long$Value / bounds$Value[bounds$Bound == "Lower"]
            }
        }

        var_long
    }))

    ### Save time series data ---------------------------------------
    write.csv(
        peaks_ts_df,
        file = file.path(Dir.Exports, "peaks_time_series.csv"),
        row.names = FALSE
    )
}

## Predictability Metrics ------------------------------------------
message("Calculating predictability metrics...")

ar <- function(x, lag = 1) {
    # Remove NA values
    ts <- x[!is.na(x)]

    # Check if there is sufficient data for the specified lag
    if (length(ts) <= lag) {
        return(NA) # Not enough data for the specified lag
    }

    # Create lagged time series
    ts_lagged <- ts[1:(length(ts) - lag)]
    ts_original <- ts[(lag + 1):length(ts)]

    # Compute the correlation
    cor(ts_lagged, ts_original, use = "complete.obs")
}

if (file.exists(file.path(Dir.Exports, "peaks_predictability.csv"))) {
    pred_df <- read.csv(file.path(Dir.Exports, "peaks_predictability.csv"))
} else {
    pred_df <- do.call(rbind, pblapply(unique(peaks_ts_df$PeakID), FUN = function(peak) {
        peak_data <- peaks_ts_df[peaks_ts_df$PeakID == peak, ]

        do.call(rbind, lapply(unique(peak_data$Variable), FUN = function(var) {
            var_data <- peak_data[peak_data$Variable == var, ]
            # Add year column
            var_data$Year <- as.numeric(format(as.Date(var_data$Date), "%Y"))

            # Calculate metrics for each year
            do.call(rbind, lapply(unique(var_data$Year), FUN = function(yr) {
                yr_data <- var_data[var_data$Year == yr, ]

                # Calculate autocorrelations for different lags
                ar_values <- sapply(c(1, 2, 3, 5, 10), function(lag) {
                    c(
                        ar(yr_data$Value[yr_data$Season == "Pre"], lag = lag),
                        ar(yr_data$Value[yr_data$Season == "Post"], lag = lag)
                    )
                })

                # # make Extreme column character
                # yr_data$Extreme <- as.character(yr_data$Extreme)

                # Create data frame with all metrics
                data.frame(
                    PeakID = peak,
                    Variable = var,
                    Year = yr,
                    Season = c("Pre", "Post"),
                    CV = c(
                        sd(yr_data$Value[yr_data$Season == "Pre"], na.rm = TRUE) /
                            mean(yr_data$Value[yr_data$Season == "Pre"], na.rm = TRUE),
                        sd(yr_data$Value[yr_data$Season == "Post"], na.rm = TRUE) /
                            mean(yr_data$Value[yr_data$Season == "Post"], na.rm = TRUE)
                    ),
                    ExtremeFreq = c(
                        sum(yr_data$Extreme[yr_data$Season == "Pre"] != "Normal") /
                            sum(yr_data$Season == "Pre"),
                        sum(yr_data$Extreme[yr_data$Season == "Post"] != "Normal") /
                            sum(yr_data$Season == "Post")
                    ),
                    AR1 = ar_values[, 1],
                    AR2 = ar_values[, 2],
                    AR3 = ar_values[, 3],
                    AR5 = ar_values[, 4],
                    AR10 = ar_values[, 5]
                )
            }))
        }))
    }))

    write.csv(
        pred_df,
        file = file.path(Dir.Exports, "peaks_predictability.csv"),
        row.names = FALSE
    )
}

## Fusing with Expedition Data -------------------------------------
message("Fusing Weather and Climate Data with Expedition Data...")
if (file.exists(file.path(Dir.Exports, "ModelData.csv"))) {
    ModelData_df <- read.csv(file.path(Dir.Exports, "ModelData.csv"))
} else {
    Expeditions_df <- Expeditions_df[Expeditions_df$PKNAME %in% sevenks_sf$PKNAME, ] # reduce to only those summits for which we have coordinates, losing 31 rows of data in Expeditions out of a total of 1917
    Expeditions_df <- Expeditions_df[as.numeric(substr(Expeditions_df$TERMDATE, 1, 4)) > 1950, ] # climate data only available for 1951 onwards, losing a further 5 expeditions
    Expeditions_df <- Expeditions_df[as.numeric(substr(Expeditions_df$BCDATE, 1, 4)) < 2022, ] # climate data only available until end of 2021, losing a further 26 expeditions
    Expeditions_df <- Expeditions_df[format(as.Date(Expeditions_df$TERMDATE), "%m") %in% c("03", "04", "05", "09", "10", "11") & format(as.Date(Expeditions_df$BCDATE), "%m") %in% c("03", "04", "05", "09", "10", "11"), ] # make sure data falls into seasons

    ModelData_ls <- pblapply(1:nrow(Expeditions_df), FUN = function(expedition) {
        # expedition = 72
        # print(expedition)
        Iter_df <- Expeditions_df[expedition, ]
        weather_df <- peaks_ts_df[
            peaks_ts_df$PeakID == Iter_df$PKNAME &
                peaks_ts_df$Date >= as.Date(Iter_df$BCDATE) &
                peaks_ts_df$Date <= as.Date(Iter_df$TERMDATE),
        ]
        unique(weather_df$Season)

        bind_df <- do.call(cbind, lapply(unique(weather_df$Variable), FUN = function(var) {
            # var = unique(weather_df$Variable)[1]
            var_df <- weather_df[weather_df$Variable == var, ]
            ret_df <- data.frame(
                mean = mean(var_df$Value, na.rm = TRUE),
                sd = sd(var_df$Value, na.rm = TRUE),
                LOWn = sum(var_df$Extreme == "LOW", na.rm = TRUE),
                HIGHn = sum(var_df$Extreme == "HIGH", na.rm = TRUE),
                LOWFreq = sum(var_df$Extreme == "LOW", na.rm = TRUE) / nrow(var_df),
                HIGHFreq = sum(var_df$Extreme == "HIGH", na.rm = TRUE) / nrow(var_df),
                AR1 = ar(var_df$Value, lag = 1),
                AR2 = ar(var_df$Value, lag = 2),
                AR3 = ar(var_df$Value, lag = 3),
                AR5 = ar(var_df$Value, lag = 5),
                AR10 = ar(var_df$Value, lag = 10)
            )
            colnames(ret_df) <- paste(var, colnames(ret_df), sep = "_")
            ret_df
        }))

        cbind(Iter_df, bind_df)
    })

    ModelData_df <- do.call(rbind, ModelData_ls)
    # head(ModelData_df)
    write.csv(ModelData_df, file = file.path(Dir.Exports, "ModelData.csv"))
}

# C. ANALYSES & VISUALISATION =============================================
message("#### Analyses & Visualizations ###############################")
MainVars <- data.frame(
    VarName = c("2m_temperature", "windspeed", "snow_cover"),
    ClearName = c("Air Temperature [K]", "Wind Speed [m/s]", "Snow Cover [%]")
)
PeakGroups <- list(
    MainText = eightks_sf$PKNAME,
    Supplement = sevenks_sf$PKNAME[!(sevenks_sf$PKNAME %in% eightks_sf$PKNAME)]
)

## Figure 1 - Summit Bid Windows & Mortality ------------------------------
df_data <- ModelData_df
### Further data prepping -------
column_to_check <- "Mortality_all" # Here: Replace 'column_to_check' with the name of the column you want to check for NAs
df_data_combined_complete <- df_data[complete.cases(df_data[, column_to_check]), ] # Remove rows where the specified column contains NAs
df_data_combined_complete$year_0 <- df_data_combined_complete$YEAR - min(df_data_combined_complete$YEAR) # add year 0 as starting year
df_data_combined_complete$DIFF_BCtoSMT <- ifelse(df_data_combined_complete$DIFF_BCtoSMT < 0, 0, df_data_combined_complete$DIFF_BCtoSMT) # All values <0 are replaced by 0
column_to_check <- "DIFF_BCtoSMT"
df_data_combined_complete <- df_data[complete.cases(df_data[, column_to_check]), ] # Remove rows where the specified column contains NAs
df_filtered <- df_data_combined_complete %>% # remove BCDate < 0
    filter(DIFF_BCtoSMT >= 0)
df_filtered$year_0 <- df_filtered$YEAR - min(df_filtered$YEAR) # adding respective time stamp
df_filtered <- df_filtered %>%
    filter(HEIGHTM >= 7000)
df_filtered$year_0 <- df_filtered$YEAR - min(df_filtered$YEAR) # add year 0 as starting year
breaks <- seq(1051, 2021, by = 10)
df_filtered$YearBin <- cut(
    df_filtered$YEAR,
    breaks = breaks,
    include.lowest = TRUE,
    right = FALSE,
    labels = paste0(breaks[-length(breaks)], ":", breaks[-1])
)
df_filtered$FacetLabel <- paste(df_filtered$SEASON, df_filtered$PEAKID, sep = " – ")
densdf_SB <- df_filtered

### Summit bid window model -------
if (file.exists(file.path(Dir.Exports, "model_SB.RData"))) {
    load(file.path(Dir.Exports, "model_SB.RData"))
} else {
    model_SB <- brms::brm(
        DIFF_BCtoSMT ~ year_0 +
            +year_0:PEAKID +
            (1 | PEAKID),
        data = df_filtered,
        # family = "gaussian",
        family = poisson(),
        chains = 4,
        cores = 4,
        # cores = parallel::detectCores(),
        iter = 10000,
        warmup = 5000
    )
    save(model_SB, file = file.path(Dir.Exports, "model_SB.RData"))
}
summary(model_SB)
# plot(model_SB)
pp_check(model_SB, ndraws = 500)
# conditional_effects(model_SB, points = TRUE)
conditions <- data.frame(PEAKID = sort(unique(df_filtered$PEAKID)))
rownames(conditions) <- sort(unique(df_filtered$PEAKID))
conde_SB <- conditional_effects(model_SB,
    conditions = conditions,
    # method = "posterior_predict", # posterior predictive
    re_formula = NULL, # Include random effects
    effects = "year_0"
)
df_SB <- conde_SB[1]$year_0

### Mortality model -------
ids_to_remove <- c("SAIP", "PUTH", "GYAJ", "HIME", "LANG", "GANG", "CHAM") # Peaks with mortality = 0, i.e. no reported deaths

# Fiter dataset, i.e. deleting peaks with no mortality
df_filtered <- df_data_combined_complete %>%
    filter(!(PEAKID %in% ids_to_remove))
df_filtered$year_0 <- df_filtered$YEAR - min(df_filtered$YEAR) # adding respective time stamp
breaks <- seq(1051, 2021, by = 10)
df_filtered$YearBin <- cut(
    df_filtered$YEAR,
    breaks = breaks,
    include.lowest = TRUE,
    right = FALSE,
    labels = paste0(breaks[-length(breaks)], ":", breaks[-1])
)
df_filtered$FacetLabel <- paste(df_filtered$SEASON, df_filtered$PEAKID, sep = " – ")
densdf_MT <- df_filtered

if (file.exists(file.path(Dir.Exports, "model_Mt.RData"))) {
    load(file.path(Dir.Exports, "model_Mt.RData"))
} else {
    model_Mt <- brms::brm(
        Mortality_all ~ year_0 +
            year_0:PEAKID +
            (1 | PEAKID),
        data = df_filtered,
        # data = df_data_combined_complete,
        family = zero_one_inflated_beta(),
        chains = 4,
        cores = parallel::detectCores(),
        iter = 10000,
        warmup = 5000,
        init_r = 0.1
    )
    save(model_Mt, file = file.path(Dir.Exports, "model_Mt.RData"))
}
summary(model_Mt)
# plot(model_Mt)
pp_check(model_Mt, ndraws = 500)
# conditional_effects(model_Mt, points = TRUE)
# round(brms::ranef(model_Mt)[[1]], 2)
conditions <- data.frame(PEAKID = sort(unique(df_filtered$PEAKID)))
rownames(conditions) <- sort(unique(df_filtered$PEAKID))
conde_Mt <- conditional_effects(model_Mt,
    conditions = conditions,
    # method = "posterior_predict", # posterior predictive
    re_formula = NULL, # Include random effects
    effects = "year_0"
)
df_MT <- conde_Mt[1]$year_0

### Plotting -------
lapply(1:length(PeakGroups), FUN = function(PKNames) {
    # PKNames <- 1
    FName <- file.path(Dir.Exports, paste0("Figure1_", names(PeakGroups)[PKNames], ".png"))
    if (file.exists(file.path(Dir.Exports, paste0("Figure1_", names(PeakGroups)[PKNames], ".png")))) {
        return(FName)
    }

    IterIDs <- summits_df$ID[match(PeakGroups[[PKNames]], summits_df$PKNAME)]
    dfSB_iter <- df_SB[df_SB$PEAKID %in% IterIDs, ]
    dfMT_iter <- df_MT[df_MT$PEAKID %in% IterIDs, ]
    dfMT_iter$Source <- "MT"
    dfSB_iter$Source <- "SB"

    colnames(dfMT_iter)[2] <- "DataValue"
    colnames(dfSB_iter)[2] <- "DataValue"
    df_combined <- rbind(dfMT_iter, dfSB_iter)

    # SB_density <- ggplot(
    #     densdf_SB[densdf_SB$PEAKID %in% IterIDs, ],
    #     aes(
    #         x = DIFF_BCtoSMT,
    #         color = YearBin,
    #         fill = YearBin
    #     )
    # ) +
    #     geom_density(alpha = 0.1, adjust = 1) +
    #     # facet_wrap(~SEASON, scales = "free", ncol = 1) + # <-- ONE LINE OF FACETING
    #     scale_color_viridis_d(option = "C", direction = -1) +
    #     scale_fill_viridis_d(option = "C", direction = -1) +
    #     labs(
    #         x = "Summit Bid Window Estimates",
    #         y = "Density",
    #         color = "10-year bins",
    #         fill = "10-year bins",
    #     ) +
    #     theme_bw() +
    #     theme(
    #         legend.position = "bottom",
    #         legend.box = "horizontal",
    #         legend.direction = "horizontal",
    #         legend.justification = "center",
    #         legend.box.spacing = unit(0, "pt")
    #     ) +
    #     guides(
    #         color = guide_legend(nrow = 1, byrow = TRUE),
    #         fill  = guide_legend(nrow = 1, byrow = TRUE)
    #     )

    # Compute scaling factor
    scale_factor <- max(dfMT_iter$estimate__, na.rm = TRUE) / max(dfSB_iter$estimate__, na.rm = TRUE)

    BidMortality_gg <- ggplot() +
        # MT line + ribbon
        geom_ribbon(
            data = dfMT_iter,
            aes(x = year_0 + 1951, ymin = lower__, ymax = upper__, fill = "Mortality Estimates"),
            alpha = 0.2
        ) +
        geom_line(
            data = dfMT_iter,
            aes(x = year_0 + 1951, y = estimate__, color = "Mortality Estimates", linetype = "Mortality Estimates"),
            size = 1
        ) +
        # SB line + ribbon (scaled) with dashed line
        geom_ribbon(
            data = dfSB_iter,
            aes(x = year_0 + 1951, ymin = lower__ * scale_factor, ymax = upper__ * scale_factor, fill = "Summit Bid Window"),
            alpha = 0.2
        ) +
        geom_line(
            data = dfSB_iter,
            aes(x = year_0 + 1951, y = estimate__ * scale_factor, color = "Summit Bid Window", linetype = "Summit Bid Window"),
            size = 1
        ) +
        # Facets per PEAKID
        facet_wrap(~PEAKID, scales = "free_y", nrow = 4) +
        # Dual axis
        scale_y_continuous(
            name = "Mortality Estimates",
            sec.axis = sec_axis(~ . / scale_factor, name = "Summit Bid Window")
        ) +
        # Manual legend for color and fill
        scale_color_manual(
            name = "Estimate Type",
            values = c("Mortality Estimates" = "#3f0027", "Summit Bid Window" = "#003f25")
        ) +
        scale_fill_manual(
            name = "Estimate Type",
            values = c("Mortality Estimates" = "#3f0027", "Summit Bid Window" = "#003f25")
        ) +
        scale_linetype_manual(
            name = "Estimate Type",
            values = c("Mortality Estimates" = "solid", "Summit Bid Window" = "dashed")
        ) +
        # Axis title colors
        theme_bw() +
        theme(
            axis.title.y = element_text(color = "#3f0027"),
            axis.title.y.right = element_text(color = "#003f25"),
            # Legend at top in single row
            legend.position = "top",
            legend.direction = "horizontal",
            legend.key.width = unit(2, "cm")
        ) +
        labs(x = "Year")
    ggsave(
        BidMortality_gg,
        file = FName,
        width = ifelse(PKNames == 1, 24, 48), height = ifelse(PKNames == 1, 20, 40)
    )
})


## Figure 2 - Climate Trends at Summits -----------------------------------
lapply(1:length(PeakGroups), FUN = function(PKNames) {
    FName <- file.path(Dir.Exports, paste0("Figure2_", names(PeakGroups)[PKNames], ".png"))
    if (file.exists(file.path(Dir.Exports, paste0("Figure2_", names(PeakGroups)[PKNames], ".png")))) {
        return(FName)
    }
    VarPlots <- lapply(1:nrow(MainVars), FUN = function(i) {
        # i = 1
        # PKNames <- 1
        ### Prepare plot data -----
        plot_df <- peaks_ts_df[peaks_ts_df$Variable == MainVars$VarName[i], ] ## limit to variable
        plot_df <- plot_df[plot_df$PeakID %in% PeakGroups[[PKNames]], ] ## limit to 8k peaks
        plot_df$Date <- as.POSIXct(plot_df$Date)
        plot_df$Year <- as.numeric(substr(plot_df$Date, 1, 4))
        plot_df <- plot_df %>%
            left_join(
                summits_df %>% select(PKNAME, ID),
                by = c("PeakID" = "PKNAME")
            )
        breaks <- seq(1051, 2021, by = 10)
        plot_df$YearBin <- cut(
            plot_df$Year,
            breaks = breaks,
            include.lowest = TRUE,
            right = FALSE,
            labels = paste0(breaks[-length(breaks)], ":", breaks[-1])
        )

        ### Line plots of change over years -----
        #### Prepare Labels +++++
        # Prepare label information for each ID and Season
        labelInfo <- lapply(split(plot_df, list(plot_df$ID, plot_df$Season)), function(dat) {
            lmmod <- lm(Value ~ Year, data = dat)
            label_year <- if (unique(dat$Season) == "Pre") min(dat$Year) else max(dat$Year)
            Label <- predict(lmmod, newdata = data.frame(Year = label_year))
            data.frame(
                ID = unique(dat$ID),
                Season = unique(dat$Season),
                LabelPos = Label,
                LabelYear = label_year,
                Pval = summary(lmmod)$coefficients["Year", "Pr(>|t|)"],
                Effect = summary(lmmod)$coefficients["Year", "Estimate"]
            )
        })
        labelInfo <- do.call(rbind, labelInfo)
        labelInfo$Fill <- as.character(labelInfo$Pval < 0.5)

        season_trends <- lapply(unique(plot_df$Season), function(season) {
            dat <- plot_df %>% filter(Season == season)
            lmmod <- lm(Value ~ Year, data = dat)
            # Label positions: start for pre, end for post
            label_year <- if (season == "Pre") min(dat$Year) else max(dat$Year)
            Label <- predict(lmmod, newdata = data.frame(Year = label_year))
            data.frame(
                ID = "Region",
                Season = season,
                LabelPos = Label,
                LabelYear = label_year,
                Pval = summary(lmmod)$coefficients["Year", "Pr(>|t|)"],
                Fill = "Region",
                Effect = summary(lmmod)$coefficients["Year", "Estimate"]
            )
        })
        season_trends <- do.call(rbind, season_trends)
        labelInfo <- rbind(labelInfo, season_trends)

        # Merge labels back to main dataframe for significance
        plot_df <- plot_df %>%
            left_join(labelInfo, by = c("ID", "Season")) %>%
            mutate(StatSig = Pval < 0.05)
        plot_df$Effect <- as.numeric(plot_df$Effect)
        labelInfo$Effect <- as.numeric(labelInfo$Effect)

        #### Actual Plot +++++
        Return_gg <- ggplot(plot_df, aes(x = Year, y = Value, group = interaction(ID, Season))) +
            # Trendlines for each peak-season
            geom_smooth(aes(linetype = StatSig, col = Effect), method = "lm", alpha = 0.3, linewidth = 0.8, show.legend = TRUE) +
            # Overall trendlines per season
            geom_smooth(aes(group = Season), method = "lm", linewidth = 2, col = "#1a1b1b", show.legend = FALSE) +
            # Labels for each ID
            geom_label_repel(
                data = labelInfo,
                aes(x = LabelYear, y = LabelPos, label = ID, fill = Effect),
                color = "#000000",
                size = 8,
                hjust = ifelse(labelInfo$Season == "Pre", 1.5, -1),
                direction = "y",
                nudge_x = 0,
                show.legend = FALSE
            ) +
            # Linetype legend for StatSig
            scale_linetype_manual(
                name = "Statistical Significance of Trend",
                values = c(`FALSE` = "dashed", `TRUE` = "solid"),
                guide = guide_legend(
                    title.position = "top",
                    label.position = "bottom",
                    nrow = 1,
                    override.aes = list(
                        size = 1.5, # thicker lines
                        color = "black" # force legend lines to black
                    ),
                    keywidth = 4,
                    keyheight = 1,
                    order = 2 # shows second in legend
                )
            ) +
            # Continuous color bar for Effect
            scale_color_viridis_c(
                direction = -1,
                name = "Mean Change Year over Year",
                guide = guide_colorbar(
                    title.position = "top",
                    barwidth = 15,
                    barheight = 1.5,
                    order = 1 # shows first in legend
                ),
            ) +
            # Fill scale for labels (optional if you want label fill colors)
            scale_fill_viridis_c(direction = -1) +
            # Axis labels
            labs(
                x = "Year",
                y = MainVars$ClearName[i]
            ) +
            # Plot limits
            lims(x = c(min(plot_df$Year) - 5, max(plot_df$Year) + 5)) +
            # Theme and legend layout
            theme_bw() +
            theme(
                legend.position = "bottom",
                legend.box = "horizontal",
                legend.direction = "horizontal",
                legend.justification = "center",
                legend.box.spacing = unit(5, "pt") # small spacing between legends
            ) +
            # Remove unwanted fill legend for labels
            guides(
                fill = "none"
            )

        ### Density plots of change over years -----
        #### Individual peaks +++++
        Var_density <- ggplot(plot_df, aes(
            x = Value,
            color = YearBin,
            fill = YearBin
        )) +
            geom_density(alpha = 0.05, adjust = 1) +
            facet_grid(PeakID ~ Season, scales = "free_x") +
            scale_color_viridis_d(option = "C", direction = -1) +
            scale_fill_viridis_d(option = "C", direction = -1) +
            labs(
                x = MainVars$ClearName[i],
                y = "Density",
                color = "10-year bins",
                fill = "10-year bins"
            ) +
            theme_bw() +
            theme(
                legend.position = "bottom",
                legend.box = "horizontal",
                legend.direction = "horizontal",
                legend.justification = "center",
                legend.box.spacing = unit(0, "pt")
            ) +
            guides(
                color = guide_legend(nrow = 1, byrow = TRUE),
                fill  = guide_legend(nrow = 1, byrow = TRUE)
            )

        #### Entire Region +++++
        plot_df$Season <- factor(plot_df$Season, levels = c("Pre", "Post"))
        plot_df$Season <- factor(
            plot_df$Season,
            levels = c("Pre", "Post"), # original values in your data
            labels = c("Pre-Monsoon", "Post-Monsoon") # desired facet labels
        )
        season_plots <- ggplot(plot_df, aes(x = Value, color = YearBin, fill = YearBin)) +
            geom_density(alpha = 0.05, adjust = 1) +
            # Continuous color scales for decades
            scale_color_viridis_d(option = "C", direction = -1, name = "10-year bins") +
            scale_fill_viridis_d(option = "C", direction = -1, name = "10-year bins") +
            # Facet by Season in one column
            facet_wrap(~Season, ncol = 1, scales = "free_y") +
            # Labels
            labs(
                x = MainVars$ClearName[i],
                y = "Density"
            ) +
            # Theme
            theme_bw() +
            theme(
                legend.position = "bottom",
                legend.box = "horizontal",
                legend.direction = "horizontal",
                legend.justification = "center",
                legend.box.spacing = unit(0, "pt"),
                strip.text = element_text(size = 12, face = "bold")
            ) +
            guides(
                color = guide_legend(nrow = 1, byrow = TRUE),
                fill  = guide_legend(nrow = 1, byrow = TRUE)
            )

        ### Return plots -----
        #### Saving individual peak plot +++++
        clean_x <- gsub("\\[.*?\\]", "", MainVars$ClearName[i])
        dir.create(file.path(Dir.Exports, clean_x))
        ggsave(
            Var_density,
            file = file.path(Dir.Exports, clean_x, paste0("Raw_Density_", clean_x, "_", names(PeakGroups)[PKNames], ".png")),
            width = 12, height = 20
        )

        ### Making panel plot +++++
        Up_gg <- plot_grid(
            Return_gg, season_plots,
            nrow = 1, rel_widths = c(1, 0.75)
        )

        Up_gg
    })
    ggsave(
        plot_grid(plotlist = VarPlots, ncol = 1),
        file = FName,
        width = ifelse(PKNames == 1, 24, 48), height = ifelse(PKNames == 1, 24, 48)
    )
})

## Figure 3 - Extreme Weather Events at Summits ---------------------------
stop("make this happen in one loop for all figures produced here and show (1) count of extremes, (2) length of consecutive extremes and (3) degree of ectreemness")

## Figure 4 - Shift in Seasonal Predictability ----------------------------
stop("show change in predictability change over time in season")

## Figure 5 - Shift in Seasonal Usage Patterns ----------------------------
stop("show change in expeditions in seasons over time as well as amount of people in the mountains")
# + Predictability

## More supplementary figures ---------------------------------------------
stop("should keep figures that produce maps!!!")


### Climate Change -------
message("Climate Change Visualizations for the Full Region...")

symmetric_range <- function(x) {
    max_abs <- max(abs(range(x, na.rm = TRUE)))
    c(-max_abs, max_abs)
}

GG_ClimChangeRegion_ls <- lapply(CDSData_ls, FUN = function(RasterIter) {
    # RasterIter <- CDSData_ls[[1]]
    LegendTitle <- climate_vars$display[climate_vars$name == varnames(RasterIter)]
    ClimateChangeTitle <- gsub("\\s*\\[.*?\\]", "", climate_vars$display[climate_vars$name == varnames(RasterIter)])

    print(paste("                                            ...", LegendTitle))

    StartMean <- mean(RasterIter[[format(time(RasterIter), "%Y") %in% head(unique(format(time(RasterIter), "%Y")), 20)]])
    StopMean <- mean(RasterIter[[format(time(RasterIter), "%Y") %in% tail(unique(format(time(RasterIter), "%Y")), 20)]])

    Currentday_gg <- ggplot() +
        geom_spatraster(data = StopMean) +
        geom_sf(data = sevenks_sf, shape = 2, colour = "white") +
        ggrepel::geom_text_repel(
            data = summits_df[summits_df$HEIGHTM >= 8000, ],
            aes(x = LON, y = LAT, label = PKNAME),
            max.overlaps = 30,
            colour = "white"
        ) +
        scale_fill_viridis_c(
            option = "C",
            name = LegendTitle,
            guide = guide_colourbar(title.vjust = 0.75)
        ) +
        labs(
            x = "Longitude",
            y = "Latitude",
            title = paste("Average", ClimateChangeTitle, "During the Last 20 Years on Record"),
        ) +
        theme(
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.width = unit(2, "cm"),
            legend.key.height = unit(1, "cm"),
            panel.background = element_rect(fill = "#2c2c2c", color = "#2c2c2c")
        )

    ClimateChange_gg <- ggplot() +
        geom_spatraster(data = StopMean - StartMean) +
        geom_sf(data = sevenks_sf, shape = 2, colour = "white") +
        ggrepel::geom_text_repel(
            data = summits_df[summits_df$HEIGHTM >= 8000, ],
            aes(x = LON, y = LAT, label = PKNAME),
            max.overlaps = 30,
            colour = "white"
        ) +
        scale_fill_gradient2(
            low = "cyan",
            mid = "grey30",
            high = "red",
            midpoint = 0,
            limits = symmetric_range(values(StopMean - StartMean)), # This ensures symmetric color scale
            name = paste("Δ", LegendTitle),
            guide = guide_colourbar(title.vjust = 0.75)
        ) +
        labs(
            x = "Longitude",
            y = "Latitude",
            title = paste(ClimateChangeTitle, "Change Between First and Last 20 Years on Record"),
        ) +
        theme(
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.width = unit(2, "cm"),
            legend.key.height = unit(1, "cm"),
            panel.background = element_rect(fill = "#2c2c2c", color = "#2c2c2c")
        )

    Boxplot_df <- rbind(
        data.frame(
            Values = values(StartMean)[, 1],
            Time = "Beginning"
        ),
        data.frame(
            Values = values(StopMean)[, 1],
            Time = "Ending"
        )
    )

    Boxplot_gg <- ggplot(Boxplot_df, aes(x = Time, y = Values)) +
        geom_violin() +
        geom_boxplot(width = 0.05) +
        stat_compare_means(comparisons = list(c("Beginning", "Ending")), paired = TRUE) +
        theme_bw() +
        scale_x_discrete(labels = c(
            "Beginning" = paste(head(unique(format(time(RasterIter), "%Y")), 20)[1], "-", tail(head(unique(format(time(RasterIter), "%Y")), 20), 1)),
            "Ending" = paste(tail(unique(format(time(RasterIter), "%Y")), 20)[1], "-", tail(tail(unique(format(time(RasterIter), "%Y")), 20), 1))
        )) +
        labs(x = "Time-Windows", y = LegendTitle)

    ret_gg <- plot_grid(
        plot_grid(
            Currentday_gg, ClimateChange_gg,
            ncol = 2
        ),
        Boxplot_gg,
        ncol = 1, rel_heights = c(1.7, 1)
    )

    var <- climate_vars$name[climate_vars$name == varnames(RasterIter)]
    ggsave(
        ret_gg,
        filename = file.path(Dir.Exports, var, paste0("ClimateChange_", var, ".png")),
        width = 24 * 1.7, height = 16 * 1.7, units = "cm"
    )

    ret_gg
})
names(GG_ClimChangeRegion_ls) <- names(CDSData_ls)

### Predictability -------
message("Predictability Change Visualizations for the Full Region...")
GG_PredictabilityRegion_ls <- lapply(CDSData_ls, FUN = function(RasterIter) {
    # RasterIter <- CDSData_ls[[1]]
    LegendTitle <- climate_vars$display[climate_vars$name == varnames(RasterIter)]
    ClimateChangeTitle <- gsub("\\s*\\[.*?\\]", "", climate_vars$display[climate_vars$name == varnames(RasterIter)])

    print(paste("                                                   ...", LegendTitle))

    Start_rast <- RasterIter[[format(time(RasterIter), "%Y") %in% head(unique(format(time(RasterIter), "%Y")), 20)]]
    Stop_rast <- RasterIter[[format(time(RasterIter), "%Y") %in% tail(unique(format(time(RasterIter), "%Y")), 20)]]


    Predictability_ls <- pblapply(c(1, 2, 3, 5, 10), FUN = function(k) {
        # k = 1
        ## ARs
        BeginAr <- app(Start_rast,
            fun = function(x) ar(x, lag = k)
        )

        EndAr <- app(Stop_rast,
            fun = function(x) ar(x, lag = k)
        )

        Currentday_gg <- ggplot() +
            geom_spatraster(data = EndAr) +
            geom_sf(data = sevenks_sf, shape = 2, colour = "white") +
            ggrepel::geom_text_repel(
                data = summits_df[summits_df$HEIGHTM >= 8000, ],
                aes(x = LON, y = LAT, label = PKNAME),
                max.overlaps = 30,
                colour = "white"
            ) +
            scale_fill_viridis_c(
                option = "C",
                name = paste(ClimateChangeTitle, "AR", k),
                guide = guide_colourbar(title.vjust = 0.75)
            ) +
            labs(
                x = "Longitude",
                y = "Latitude",
                title = paste("Average", ClimateChangeTitle, k, "Day Autocorrelation During the Last 20 Years on Record"),
            ) +
            theme(
                legend.position = "bottom",
                legend.direction = "horizontal",
                legend.key.width = unit(2, "cm"),
                legend.key.height = unit(1, "cm"),
                panel.background = element_rect(fill = "#2c2c2c", color = "#2c2c2c")
            )

        ClimateChange_gg <- ggplot() +
            geom_spatraster(data = EndAr - BeginAr) +
            geom_sf(data = sevenks_sf, shape = 2, colour = "white") +
            ggrepel::geom_text_repel(
                data = summits_df[summits_df$HEIGHTM >= 8000, ],
                aes(x = LON, y = LAT, label = PKNAME),
                max.overlaps = 30,
                colour = "white"
            ) +
            scale_fill_gradient2(
                low = "cyan",
                mid = "grey30",
                high = "red",
                midpoint = 0,
                limits = symmetric_range(values(EndAr - BeginAr)), # This ensures symmetric color scale
                name = paste("Δ", ClimateChangeTitle, "AR", k),
                guide = guide_colourbar(title.vjust = 0.75)
            ) +
            labs(
                x = "Longitude",
                y = "Latitude",
                title = paste("Change in", k, "day Autocorrelation of", ClimateChangeTitle, "Between First and Last 20 Years on Record"),
            ) +
            theme(
                legend.position = "bottom",
                legend.direction = "horizontal",
                legend.key.width = unit(2, "cm"),
                legend.key.height = unit(1, "cm"),
                panel.background = element_rect(fill = "#2c2c2c", color = "#2c2c2c")
            )

        Boxplot_df <- rbind(
            data.frame(
                Values = values(BeginAr)[, 1],
                Time = "Beginning"
            ),
            data.frame(
                Values = values(EndAr)[, 1],
                Time = "Ending"
            )
        )

        Boxplot_gg <- ggplot(Boxplot_df, aes(x = Time, y = Values)) +
            geom_violin() +
            geom_boxplot(width = 0.05) +
            stat_compare_means(comparisons = list(c("Beginning", "Ending")), paired = TRUE) +
            theme_bw() +
            scale_x_discrete(labels = c(
                "Beginning" = paste(head(unique(format(time(RasterIter), "%Y")), 20)[1], "-", tail(head(unique(format(time(RasterIter), "%Y")), 20), 1)),
                "Ending" = paste(tail(unique(format(time(RasterIter), "%Y")), 20)[1], "-", tail(tail(unique(format(time(RasterIter), "%Y")), 20), 1))
            )) +
            labs(x = "Time-Windows", y = paste(ClimateChangeTitle, "AR", k))

        plot_grid(Currentday_gg, Boxplot_gg, ClimateChange_gg, ncol = 3, rel_widths = c(1, 0.7, 1))
    })
    names(Predictability_ls) <- paste0("AR", c(1, 2, 3, 5, 10))

    var <- climate_vars$name[climate_vars$name == varnames(RasterIter)]
    ggsave(
        plot_grid(plotlist = Predictability_ls, ncol = 1),
        filename = file.path(Dir.Exports, var, paste0("PredictabilityChange_", var, ".png")),
        width = 30 * 2.1, height = 42 * 2.1, units = "cm"
    )

    Predictability_ls
})
names(GG_PredictabilityRegion_ls) <- names(CDSData_ls)

## Summits ---------------------------------------------------------
stop("cleaned to here - remake plots and models")


### Climate Change -------
peaks_ts_df
### Predictability -------
pred_df
### Extremes -------
peaks_ts_df

## Humans ----------------------------------------------------------
ModelData_df
### Expeditions over Time -------
# Calculate expedition duration and create year column
ModelData_df$Duration <- as.numeric(difftime(as.Date(ModelData_df$TERMDATE),
    as.Date(ModelData_df$BCDATE),
    units = "days"
))
ModelData_df$Year <- as.numeric(format(as.Date(ModelData_df$BCDATE), "%Y"))

# Create daily presence data with member counts
expedition_days <- do.call(rbind, lapply(1:nrow(ModelData_df), function(i) {
    dates <- seq(as.Date(ModelData_df$BCDATE[i]),
        as.Date(ModelData_df$TERMDATE[i]),
        by = "day"
    )
    data.frame(
        Date = dates,
        Year = as.numeric(format(dates, "%Y")),
        DayOfYear = as.numeric(format(dates, "%j")),
        Members = ModelData_df$TOTMEMBERS[i] # Add member count
    )
}))

# Create complete grid of years and days
all_years <- 1951:2021
all_days <- 1:365
complete_grid <- expand.grid(Year = all_years, DayOfYear = all_days)

# Count expeditions and total members per day
expedition_counts <- aggregate(
    cbind(Expeditions = 1, Members = Members) ~ Year + DayOfYear,
    data = expedition_days,
    FUN = sum
)
expedition_counts$Year <- as.numeric(as.character(expedition_counts$Year))
expedition_counts$DayOfYear <- as.numeric(as.character(expedition_counts$DayOfYear))

# Merge with complete grid to include all years and days
expedition_counts <- merge(
    complete_grid,
    expedition_counts,
    by = c("Year", "DayOfYear"),
    all.x = TRUE
)

# Replace NAs with zeros
expedition_counts$Expeditions[is.na(expedition_counts$Expeditions)] <- 0
expedition_counts$Members[is.na(expedition_counts$Members)] <- 0

# Reshape data to long format for faceting
expedition_counts <- tidyr::pivot_longer(
    expedition_counts,
    cols = c(Expeditions, Members),
    names_to = "type",
    values_to = "value"
)
expedition_counts <- expedition_counts[expedition_counts$value != 0, ]

# Create separate plots for expeditions and members
ExpeditionPlots_ls <- lapply(unique(expedition_counts$type), function(plot_type) {
    data_subset <- expedition_counts[expedition_counts$type == plot_type, ]

    # Create title based on type
    title <- if (plot_type == "Expeditions") {
        "Number of Active Expeditions Throughout the Year"
    } else {
        "Total Number of Expedition Members Throughout the Year"
    }

    ggplot(
        data_subset,
        aes(
            x = DayOfYear,
            y = as.factor(Year),
            fill = value
        )
    ) +
        # Add vertical lines for pre and post monsoon seasons
        geom_vline(xintercept = c(60, 151), color = "white", alpha = 0.3, linetype = "dashed") + # March (60) to May (151)
        geom_vline(xintercept = c(244, 304), color = "white", alpha = 0.3, linetype = "dashed") + # September (244) to November (304)
        geom_tile() +
        scale_fill_viridis_c(
            name = "", # plot_type,
            option = "C"
        ) +
        scale_y_discrete(
            name = "Year",
            limits = rev(as.character(1951:2021)), # Reverse to show earliest at bottom
            breaks = rev(as.character(seq(1951, 2021, by = 5))) # Show only every 5th year
        ) +
        scale_x_continuous(
            name = "Month",
            breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
            labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"),
            limits = c(1, 365), # Set limits to show full year
            expand = c(0, 0) # Remove spacing at edges
        ) +
        labs(
            title = title,
            subtitle = "Vertical lines indicate pre-monsoon and post-monsoon seasons"
        ) +
        theme_bw() +
        theme(
            legend.position = "right",
            legend.key.width = unit(0.5, "cm"),
            legend.key.height = unit(2, "cm"),
            panel.background = element_rect(fill = "#2c2c2c", color = "#2c2c2c"),
            panel.grid = element_line(color = "grey30"),
            text = element_text(color = "white"),
            axis.text = element_text(color = "white"),
            plot.background = element_rect(fill = "#2c2c2c", color = "#2c2c2c"),
            legend.background = element_rect(fill = "#2c2c2c"),
            legend.text = element_text(color = "white"),
            legend.title = element_text(color = "white"),
            panel.spacing = unit(0, "lines"), # Remove spacing between tiles
            axis.ticks = element_blank() # Remove axis ticks for cleaner look
        )
})

# Name the list elements
names(ExpeditionPlots_ls) <- unique(expedition_counts$type)

# Combine plots vertically
ExpeditionDurations_gg <- plot_grid(
    ExpeditionPlots_ls$Expeditions,
    ExpeditionPlots_ls$Members,
    ncol = 1,
    align = "v",
    rel_heights = c(1, 1)
)
ExpeditionDurations_gg

# Save the plot
ggsave(
    ExpeditionDurations_gg,
    filename = file.path(Dir.Exports, "ExpeditionDurations.png"),
    width = 24, height = 24, units = "cm"
)

### Summit Bid Window -------
### Mortality -------

## Kriging ----------------------------------------------------------
### Covariates -------
### Kriging -------
### Visualization -------






post_df <- aggregate(Windspeed_Season == "Post" ~ YEAR, ModelData_df, FUN = sum)
colnames(post_df) <- c("Year", "Count")
post_df$Total <- table(ModelData_df$YEAR)
post_df$Proportion <- post_df$Count / post_df$Total

BarTime <- ggplot(post_df, aes(x = factor(Year), y = Total)) +
    geom_bar(stat = "identity", fill = "#156082") +
    theme_bw() +
    labs(y = paste("Number of Expeditions on Mountains >8,000 masl"), x = "Year") +
    geom_label_repel(aes(label = Total),
        box.padding   = 0.35,
        point.padding = 0.5,
        segment.color = "grey50"
    )

ggsave(BarTime,
    filename = file.path(Dir.Exports, "BarTime.png"),
    width = 26, height = 8
)

SeasonsExpeds <- ggplot(post_df, aes(x = Year, y = Proportion)) +
    geom_point() +
    geom_label_repel(aes(label = Total),
        box.padding   = 0.35,
        point.padding = 0.5,
        segment.color = "grey50"
    ) +
    geom_smooth(method = "loess", fill = "#156082", col = "#156082") +
    theme_bw() +
    labs(y = paste("Proportion of Post-Monsoon Season Expeditions"), x = "Year")

ggsave(SeasonsExpeds,
    filename = file.path(Dir.Exports, "SeasonsExpeds .png"),
    width = 26, height = 8
)

post_df <- aggregate(STORM ~ YEAR + Windspeed_Season, ModelData_df, FUN = sum)
colnames(post_df) <- c("Year", "Season", "Count")
ModelData_df$Dummy <- 1
post_df$Total <- aggregate(Dummy ~ YEAR + Windspeed_Season, ModelData_df, FUN = sum)[, 3]
post_df$Proportion <- post_df$Count / post_df$Total

Storm_df <- post_df
Storm_df$Condition <- "STORM"


post_df <- aggregate(AVALANCHE ~ YEAR + Windspeed_Season, ModelData_df, FUN = sum)
colnames(post_df) <- c("Year", "Season", "Count")
ModelData_df$Dummy <- 1
post_df$Total <- aggregate(Dummy ~ YEAR + Windspeed_Season, ModelData_df, FUN = sum)[, 3]
post_df$Proportion <- post_df$Count / post_df$Total
post_df$Condition <- "AVALANCHE"

AVASTORMS <- ggplot(rbind(Storm_df, post_df), aes(x = Year, y = Proportion)) +
    geom_point() +
    geom_label_repel(aes(label = Total),
        box.padding = 0.1,
        max.overlaps = 30,
        segment.color = "grey50"
    ) +
    facet_grid(factor(Condition, levels = c("STORM", "AVALANCHE")) ~ factor(Season, levels = c("Pre", "Post"))) +
    geom_smooth(fill = "#156082", col = "#156082") +
    theme_bw() +
    labs(y = paste("Proportion of Expeditions per Season Reporting Hazard Events"))
# +
# lims(y = c(0, 1))

ggsave(AVASTORMS,
    filename = file.path(Dir.Exports, "AVALANCHES_STORMS .png"),
    width = 22 * 1.2, height = 8 * 1.2
)
