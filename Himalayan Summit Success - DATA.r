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
options(timeout = 600) # set timeout for downloading data and packages - my home internet connection is off-grid and so can drop at times, this fixes timeout issues

## Packages ---------------------------------------------------------------
packages <- list(
    core = c("readr", "dplyr", "tidyr", "terra", "sf", "sp", "devtools", "foreign"),
    viz = c("ggplot2", "viridis", "cowplot", "mapview", "ggrepel", "tidyterra", "ggpubr", "grid", "png", "ggh4x"),
    spatial = c("rnaturalearth", "rnaturalearthdata"),
    stats = c("brms", "tidybayes", "broom", "changepoint", "lme4", "lmerTest"), # , "purr"
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
    name = c("2m_temperature", "snow_cover", "windspeed"),
    display = c("Air Temperature [K]", "Snow Cover [%]", "Windspeed [m/s]"),
    subset = TRUE,
    stringsAsFactors = FALSE
)

# For backward compatibility
VNames_vec <- climate_vars$display
SubsetVariables_vec <- climate_vars$name[climate_vars$subset]

### Plotting Settings ------
PreColour <- "#AEC647"
PostColour <- "#3F6B8F"

### Death Types ------
DEATHTYPE_LABELS <- c(
    `4` = "Fall", `7` = "Avalanche", `6` = "Icefall/\nserac",
    `8` = "Rockfall", `5` = "Crevasse"
)
CAUSE_ORDER <- c("Fall", "Avalanche", "Icefall/\nserac", "Rockfall", "Crevasse")

# A. DATA LOADING =========================================================
message("#### Loading Data from Disk ##################################")

## Expedition DATA -------------------------------------------------
Expeditions_df <- read.csv(file.path(Dir.Data, "CleanedExpeditions.csv"))
df_members <- read.dbf(file.path(Dir.Data,"members.DBF"))

## Summits as Spatial Objects --------------------------------------
summits_df <- read_csv(file.path(Dir.Data, "selected_peaks_coordinates_counts.csv")) # load positions and names of summits
summits_df <- summits_df[!duplicated(summits_df$ID), ] # check for duplicates and eventually delete them
## making sp object of summits
summits_sp <- summits_df
coordinates(summits_sp) <- ~ LON + LAT
proj4string(summits_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
## making into sf object
summits_sf <- st_as_sf(summits_sp)
PeakswithEnoughExpeds <- names(table(Expeditions_df$PEAKID[Expeditions_df$YEAR > 1950]))[table(Expeditions_df$PEAKID[Expeditions_df$YEAR > 1950]) > 25]
eightks_sf <- summits_sf[summits_sf$HEIGHTM >= 8000, ]
eightks_sf <- eightks_sf[eightks_sf$ID %in% PeakswithEnoughExpeds, ]
TargetIDs <- eightks_sf$ID
Years_vec <- 1951:2021

## ERA5-Land Data Download -----------------------------------------
if (file.exists(file.path(Dir.Exports, "peaks_time_series_DAILY.rds"))) {
    message("#### Data already prepared ##################################")
} else {
    message("#### Download Data from CDS ##################################")
    Variables_vec <- c(
        "2m_temperature",
        # "skin_temperature",
        "10m_u_component_of_wind",
        "10m_v_component_of_wind",
        "snow_cover"
        # ,
        # "snow_density",
        # "snow_depth",
        # "snow_depth_water_equivalent",
        # "snow_evaporation",
        # "snowfall",
        # "snowmelt",
        # "temperature_of_snow_layer"
    )

    CDSData_ls <- lapply(Variables_vec, FUN = function(Var_Iter) {
        message(Var_Iter)
        if (file.exists(file.path(Dir.Data, paste0(Var_Iter, "_Raw.nc")))) {
            print("Already prepared")
            return(rast(file.path(Dir.Data, paste0(Var_Iter, "_Raw.nc"))))
        }

        while (!file.exists(file.path(Dir.Data, paste0(Var_Iter, "_Raw.nc")))) {
            tryCatch(
                {
                    Raw_rast <- KrigR::CDownloadS(
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
                },
                error = function(e) {
                    message("Error downloading data for ", Var_Iter, ": ", e$message)
                    message("Retrying in 10 seconds...")
                    Sys.sleep(10) # wait for a minute before retrying
                }
            )
        }

        # Var_data <- KrigR:::Temporal.Aggr(Raw_rast,
        #     BaseResolution = "hour",
        #     BaseStep = 1,
        #     TResolution = "day",
        #     TStep = 1,
        #     FUN = mean,
        #     Cores = numberOfCores,
        #     TZone = "UTC"
        # )
        # terra::writeCDF(Var_data, file = file.path(Dir.Data, paste0(Var_Iter, ".nc")))
        # Var_data
    })
    names(CDSData_ls) <- Variables_vec

    # ## Windspeed Calculation -------------------------------------------
    # ### Hourly -------
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
    # CSData_ls$windspeed <- rast(file.path(Dir.Data, "windspeed_RAW.nc"))

    ## Data Subsetting -------------------------------------------------
    # CDSData_ls <- CDSData_ls[SubsetVariables_vec]

    # B. DATA EXTRACTION ======================================================
    message("#### Extracting Data for Analyses ############################")

    ## Time-Series Extraction for Peaks --------------------------------
    print("Creating time series for peaks...")
    ### Extract daily time series -------
    print("                               ... extracting data")
    peak_ts_ls <- pblapply(names(CDSData_ls), FUN = function(VarName) {
        VarIter <- CDSData_ls[[VarName]]
        values_mat <- t(terra::extract(
            VarIter,
            eightks_sf
            # ,
            # method = "bilinear",
            # fun = "mean"
        ))
        values_df <- as.data.frame(values_mat[-1, ])
        colnames(values_df) <- eightks_sf$PKNAME
        values_df$Date <- time(VarIter)
        values_df$Variable <- VarName
        values_df$Month <- as.numeric(format(values_df$Date, "%m"))
        values_df$Season <- "Out Of Season"
        values_df$Season[values_df$Month %in% 3:5] <- "Pre"
        values_df$Season[values_df$Month %in% 9:11] <- "Post"
        values_df
    })
    names(peak_ts_ls) <- names(CDSData_ls)
    peak_ts_ls$windspeed <- peak_ts_ls$`10m_u_component_of_wind`
    peak_ts_ls$windspeed$Variable <- "windspeed"
    peak_ts_ls$windspeed[, 1:nrow(eightks_sf)] <- sqrt(peak_ts_ls$`10m_u_component_of_wind`[, 1:nrow(eightks_sf)]^2 + peak_ts_ls$`10m_v_component_of_wind`[, 1:nrow(eightks_sf)]^2)

    ### Calculate seasonal bounds -------
    ## loop here twice, once for hourly and once for daily, to get seasonal bounds for both and write final data frames to disk
    lapply(1:2, FUN = function(i) {
        if (i == 1) {
            print("                               ... calculating seasonal bounds for HOURLY data")
        } else {
            print("                               ... calculating seasonal bounds for DAILY data")
            peak_ts_ls <- lapply(peak_ts_ls, FUN = function(x) {
                ## build daily mean based on Date for the different peaks, assign to new data frame and attach Variable, Month and Season columns with first value for each day
                x <- x %>%
                    dplyr::mutate(Day = as.Date(Date)) %>%
                    dplyr::group_by(Day) %>%
                    dplyr::summarise(
                        Variable = dplyr::first(Variable),
                        Month = dplyr::first(Month),
                        Season = dplyr::first(Season),
                        dplyr::across(dplyr::all_of(eightks_sf$PKNAME), ~ mean(.x, na.rm = TRUE)),
                        .groups = "drop"
                    ) %>%
                    dplyr::ungroup() %>%
                    dplyr::rename(Date = Day) %>%
                    dplyr::select(dplyr::all_of(eightks_sf$PKNAME), Date, Variable, Month, Season)
                x
            })
        }

        seasonal_bounds <- do.call(rbind, pblapply(names(peak_ts_ls), FUN = function(VarName) {
            var_df <- peak_ts_ls[[VarName]]
            do.call(rbind, lapply(eightks_sf$PKNAME, FUN = function(peak) {
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
        # print(seasonal_bounds)

        ### Create final time series dataframe -------
        print("                              ... making final data frame")
        FinalVars <- names(peak_ts_ls) ## but remove 10m_u_component_of_wind and 10m_v_component_of_wind, as they are not needed anymore
        FinalVars <- FinalVars[!FinalVars %in% c("10m_u_component_of_wind", "10m_v_component_of_wind")]
        peaks_ts_df <- do.call(rbind, pblapply(FinalVars, FUN = function(VarName) {
            var_df <- peak_ts_ls[[VarName]]
            var_long <- tidyr::pivot_longer(
                var_df,
                cols = eightks_sf$PKNAME,
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
                    var_long$ExtremeRatioLow <- bounds$Value[bounds$Bound == "Lower"] / var_long$Value
                }
            }

            var_long
        }))

        ## sort by Variable, PeakID, and Date
        peaks_ts_df <- peaks_ts_df[order(peaks_ts_df$Variable, peaks_ts_df$PeakID, peaks_ts_df$Date), ]

        ### Save time series data ---------------------------------------
        saveRDS(
            peaks_ts_df,
            file = file.path(Dir.Exports, paste0("peaks_time_series_", ifelse(i == 1, "HOURLY", "DAILY"), ".rds"))
        )
    }) # temporal resolution lapply
}
Peaks_ts_HOURLY_df <- readRDS(file.path(Dir.Exports, "peaks_time_series_HOURLY.rds"))
Peaks_ts_DAILY_df <- readRDS(file.path(Dir.Exports, "peaks_time_series_DAILY.rds"))

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
    pred_df <- do.call(rbind, pblapply(unique(Peaks_ts_DAILY_df$PeakID), FUN = function(peak) {
        peak_data <- Peaks_ts_DAILY_df[Peaks_ts_DAILY_df$PeakID == peak, ]

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
    Expeditions_df <- Expeditions_df[Expeditions_df$PKNAME %in% eightks_sf$PKNAME, ] # reduce to only those summits for which we have coordinates, losing 31 rows of data in Expeditions out of a total of 1917
    Expeditions_df <- Expeditions_df[as.numeric(substr(Expeditions_df$TERMDATE, 1, 4)) > 1950, ] # climate data only available for 1951 onwards, losing a further 5 expeditions
    Expeditions_df <- Expeditions_df[as.numeric(substr(Expeditions_df$BCDATE, 1, 4)) < 2022, ] # climate data only available until end of 2021, losing a further 26 expeditions
    Expeditions_df <- Expeditions_df[format(as.Date(Expeditions_df$TERMDATE), "%m") %in% c("03", "04", "05", "09", "10", "11") & format(as.Date(Expeditions_df$BCDATE), "%m") %in% c("03", "04", "05", "09", "10", "11"), ] # make sure data falls into seasons

    ModelData_ls <- pblapply(1:nrow(Expeditions_df), FUN = function(expedition) {
        # expedition = 72
        # print(expedition)
        Iter_df <- Expeditions_df[expedition, ]
        weather_df <- Peaks_ts_DAILY_df[
            Peaks_ts_DAILY_df$PeakID == Iter_df$PKNAME &
                Peaks_ts_DAILY_df$Date >= as.Date(Iter_df$BCDATE) &
                Peaks_ts_DAILY_df$Date <= as.Date(Iter_df$TERMDATE),
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
