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
    utils = c("pbapply")
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
message("#### Extracting Data for Analyses ###########################")

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

## Entire Region ---------------------------------------------------
### Climate Change -------
### Predictability -------

## Summits ---------------------------------------------------------
### Climate Change -------
### Predictability -------
### Extremes -------

## Humans ----------------------------------------------------------
### Expeditions over Time -------
### Summit Bid Window -------
### Mortality -------

## Kriging ----------------------------------------------------------
### Covariates -------
### Kriging -------
### Visualization -------





stop("cleaned to here - remake plots and models")
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
