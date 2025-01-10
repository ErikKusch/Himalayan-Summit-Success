#' ####################################################################### #
#' PROJECT: [Himalayan Summit Success - Climate Data]
#' CONTENTS:
#'  - Download, aggregate, and analyse climate data in Himalayas
#'  DEPENDENCIES:
#'  - None
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PREAMBLE ================================================================
rm(list = ls()) # some may not like it, but it helps my workflow

## Packages ---------------------------------------------------------------
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}
### CRAN PACKAGES ----
package_vec <- c(
  "pbapply", # for parallel lapply
  "readr", # for reading csv
  "ggplot2", # for plotting
  "tidyr", # for turning rasters into ggplot-dataframes
  "viridis", # colour palettes
  "cowplot", # gridding multiple plots
  "rnaturalearth", # for shapefiles
  "rnaturalearthdata", # for high-resolution shapefiles
  "mapview", # for generating mapview outputs
  "sp", # for legacy spatial handling
  "sf", # for more efficient handling of sp data
  "terra", # for more efficient handling of raster data
  "ggplot2", # ggploting engine
  "ggrepel", # repelled labels
  "tidyterra", # ggploting of terra files
  "dplyr", # for reshaping time series extractions
  "ggpubr" # for statistical comparisons in plots
)
sapply(package_vec, install.load.package)

### NON-CRAN PACKAGES ----
#' needs rewrite for dev version of KrigR
if ("KrigR" %in% rownames(installed.packages()) == FALSE) { # KrigR check
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
  devtools::install_github("ErikKusch/KrigR")
}
library(KrigR)
package_vec <- c("KrigR", package_vec)

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
rm(Dirs)

# DATA ====================================================================
## Summit windows ---------------------------------------------------------
if (file.exists(file.path(Dir.Data, "CleanedExpeditions.csv"))) {
  Expeditions_df <- read.csv(file.path(Dir.Data, "CleanedExpeditions.csv"))
} else {
  Expeditions_df <- read_csv(file.path(Dir.Data, "individual_expeditions_higher_7000.csv"))
  ## throw out everything that lacks all date specifications or only have one specified date (due to date format ambiguity, this may me severly misleading)
  Expeditions_df <- Expeditions_df[rowSums(is.na(
    Expeditions_df[, c("BCDATE", "SMTDATE", "TERMDATE")]
  )) < 2, ]
  ## when summit days and or total days have been recorded, then fill date fields with them
  Expeditions_df$TOTDAYS <- ifelse(Expeditions_df$TOTDAYS == 0, Expeditions_df$SMTDAYS, Expeditions_df$TOTDAYS)
  DaysCheck <- which(rowSums(Expeditions_df[, c("SMTDAYS", "TOTDAYS")]) > 0)
  for (Days_i in DaysCheck) {
    # Days_i <- 491
    Days_iter <- Expeditions_df[Days_i, c("BCDATE", "SMTDATE", "TERMDATE", "SMTDAYS", "TOTDAYS")]
    DaysDiff <- unlist(c(Days_iter[4:5]))
    DaysDiff <- c(DaysDiff, diff(DaysDiff))
    ## try to grab and fix date issues here
    if (sum(!is.na(Days_iter[1:3])) > 1) {
      ActualDiff <- Days_iter[which(!is.na(Days_iter[1:3]))][2] - Days_iter[which(!is.na(Days_iter[1:3]))][1]
      RecordedDiff <- DaysDiff[max(which(!is.na(Days_iter[1:3])))]
      if (ActualDiff != RecordedDiff) {
        Days_iter$BCDATE <- as.Date(as.character(Days_iter$BCDATE), format = "%Y-%d-%m")
        Days_iter$SMTDATE <- as.Date(as.character(Days_iter$SMTDATE), format = "%Y-%d-%m")
        Days_iter$TERMDATE <- as.Date(as.character(Days_iter$TERMDATE), format = "%Y-%d-%m")
      }
    }
    if (!is.na(Days_iter$BCDATE)) {
      if (is.na(Days_iter$SMTDATE)) {
        Days_iter$SMTDATE <- Days_iter$BCDATE + Days_iter$SMTDAYS
      }
      if (is.na(Days_iter$TERMDATE)) {
        Days_iter$TERMDATE <- Days_iter$BCDATE + Days_iter$TOTDAYS
      }
    }
    if (!is.na(Days_iter$TERMDATE)) {
      if (is.na(Days_iter$BCDATE)) {
        Days_iter$BCDATE <- Days_iter$TERMDATE - Days_iter$TOTDAYS
      }
      if (is.na(Days_iter$SMTDATE)) {
        Days_iter$SMTDATE <- Days_iter$BCDATE + Days_iter$SMTDAYS
      }
    }
    if (!is.na(Days_iter$SMTDATE)) {
      if (is.na(Days_iter$BCDATE)) {
        Days_iter$BCDATE <- Days_iter$SMTDATE - Days_iter$SMTDAYS
      }
      if (is.na(Days_iter$TERMDATE)) {
        Days_iter$TERMDATE <- Days_iter$BCDATE + Days_iter$TOTDAYS
      }
    }
    Expeditions_df[Days_i, c("BCDATE", "SMTDATE", "TERMDATE", "SMTDAYS", "TOTDAYS")] <- Days_iter
  }
  ## check data formatting errors (we have a mix of american YYYY-DD-MM and sensible YYYY-MM-DD formatting); if TERMDate or SMTDate predate BCDate or each other, we know there is an issue
  DateCheck <- which(Expeditions_df$SMTDATE < Expeditions_df$BCDATE |
    Expeditions_df$TERMDATE < Expeditions_df$BCDATE |
    Expeditions_df$TERMDATE < Expeditions_df$SMTDATE)
  Expeditions_df[DateCheck, ]$BCDATE <- as.Date(as.character(Expeditions_df[DateCheck, ]$BCDATE), format = "%Y-%d-%m")
  Expeditions_df[DateCheck, ]$SMTDATE <- as.Date(as.character(Expeditions_df[DateCheck, ]$SMTDATE), format = "%Y-%d-%m")
  Expeditions_df[DateCheck, ]$TERMDATE <- as.Date(as.character(Expeditions_df[DateCheck, ]$TERMDATE), format = "%Y-%d-%m")
  ## calculate median times between Basecamp and summit and basecamp and termination by decade
  Expeditions_df$DIFF_BCtoSMT <- Expeditions_df$SMTDATE - Expeditions_df$BCDATE
  Expeditions_df$DIFF_BCtoTERM <- Expeditions_df$TERMDATE - Expeditions_df$BCDATE
  Diff_SMT <- aggregate(DIFF_BCtoSMT ~ DECADE, FUN = median, data = Expeditions_df)
  Diff_TERM <- aggregate(DIFF_BCtoTERM ~ DECADE, FUN = median, data = Expeditions_df)
  ## Impute BC and termination date for each expedition lacking this information
  Expeditions_df$IMPUTATION <- FALSE
  ImputeCheck <- which(is.na(Expeditions_df$BCDATE) |
    (is.na(Expeditions_df$TERMDATE) & is.na(Expeditions_df$SMTDATE)))
  for (Impute_i in ImputeCheck) {
    Expeditions_df$IMPUTATION[Impute_i] <- TRUE
    Impute_df <- Expeditions_df[Impute_i, c("BCDATE", "SMTDATE", "TERMDATE")]
    Impute_decade <- Expeditions_df$DECADE[Impute_i]
    ### Basecamp date imputation
    if (is.na(Impute_df$BCDATE)) {
      ImSMT <- Impute_df$SMTDATE - Diff_SMT$DIFF_BCtoSMT[Diff_SMT$DECADE == Impute_decade]
      ImTERM <- Impute_df$TERMDATE - Diff_TERM$DIFF_BCtoTERM[Diff_TERM$DECADE == Impute_decade]
      if (sum(is.na(c(ImSMT, ImTERM)) != 2)) { ## if only one could be computed
        BC_imp <- na.omit(c(ImSMT, ImTERM))
        ImSource <- c("SMT", "TERM")[which(!is.na(c(ImSMT, ImTERM)))]
      } else { ## otherwise take earlier step
        BC_imp <- c(ImSMT, ImTERM)[which.min(c(ImSMT, ImTERM))]
        ImSource <- c("SMT", "TERM")[which.min(c(ImSMT, ImTERM))]
      }
      Expeditions_df$IMPUTATION[Impute_i] <- paste0(
        Expeditions_df$IMPUTATION[Impute_i],
        "_BC-", ImSource
      )
      Expeditions_df$BCDATE[Impute_i] <- BC_imp
    }
    ### Termination date imputation
    if (is.na(Impute_df$TERMDATE) & is.na(Impute_df$SMTDATE)) {
      TERM_imp <- Expeditions_df$BCDATE[Impute_i] +
        Diff_TERM$DIFF_BCtoTERM[Diff_TERM$DECADE == Impute_decade]
      Expeditions_df$IMPUTATION[Impute_i] <- paste0(Expeditions_df$IMPUTATION[Impute_i], "_TERM")
      Expeditions_df$TERMDATE[Impute_i] <- TERM_imp
    }
  }
  Expeditions_df <- Expeditions_df[which(rowSums(is.na(Expeditions_df[, c("BCDATE", "SMTDATE", "TERMDATE")])) == 0), ]
  Expeditions_df <- Expeditions_df[!(Expeditions_df$SMTDATE > Expeditions_df$TERMDATE), ]
  write.csv(Expeditions_df, file.path(Dir.Data, "CleanedExpeditions.csv"))
}

## find start and stop for analysis
TimeWindows <- data.frame(
  Start = as.POSIXct(Expeditions_df$BCDATE, tz = "UTC") - 14,
  Stop = as.POSIXct(Expeditions_df$TERMDATE, tz = "UTC")
)

## Summits as Spatial Objects ---------------------------------------------
Nepal_shp <- ne_states(country = "Nepal") # shape data of Nepal from rnaturalearth package.

summits_df <- read_csv(file.path(Dir.Data, "selected_peaks_coordinates_counts.csv")) # load positions and names of summits; more peaks in: mountains_df2.csv
summits_df <- summits_df[!duplicated(summits_df$ID), ] # check for duplicates and eventually delete them
## making sp object of summits
summits_sp <- summits_df
coordinates(summits_sp) <- ~ LON + LAT
proj4string(summits_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
summits_sf <- st_as_sf(summits_sp)
eightks_sf <- summits_sf[summits_sf$HEIGHTM >= 8000, ]
narrow_buffer <- KrigR::Buffer.pts(summits_sf, 2e4) # equates to roughly three grid cells in either direction
wide_buffer <- KrigR::Buffer.pts(summits_sf, 1e5) # equates to roughly three grid cells in either direction

## ERA5-Land --------------------------------------------------------------
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

Data_ls <- lapply(Variables_vec, FUN = function(Var_Iter) {
  message(paste("###", Var_Iter))
  if (file.exists(file.path(Dir.Data, paste0(Var_Iter, ".nc")))) {
    print("Already prepared")
    Var_data <- rast(file.path(Dir.Data, paste0(Var_Iter, ".nc")))
  } else {
    Raw_rast <- CDownloadS(
      Variable = Var_Iter,
      CumulVar = ifelse(Var_Iter %in% "snowfall", TRUE, FALSE),
      DataSet = "reanalysis-era5-land",
      Type = NA,
      DateStart = paste0(Years_vec[1], "-01-01 00:00"),
      DateStop = paste0(tail(Years_vec, 1), "-12-31 23:00"),
      TResolution = "hour",
      TStep = 1,
      Extent = ext(c(78.83, 89.33, 25.143, 31.543)), # ext(Nepal_shp) + c(-1, 1, -1, 1)
      Dir = Dir.Data,
      FileName = paste0(Var_Iter, "_Raw"),
      API_User = API_User, # API User Number
      API_Key = API_Key, # API User Key
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
  }
  Var_data
})
names(Data_ls) <- Variables_vec

## calculate windspeed from u and v component
if (!file.exists(file.path(Dir.Data, "windspeed_RAW.nc"))) {
  RAW_ls <- lapply(c("10m_u_component_of_wind", "10m_v_component_of_wind"), FUN = function(i) {
    # print(i)
    # rast(
    #   list.files(file.path(Dir.Data, i), full.names = TRUE)
    # )
    rast(file.path(Dir.Data, paste0(i, "_RAW.nc")))
  })
  Indices <- ceiling((1:terra::nlyr(RAW_ls[[1]])) / 2e4)
  ru_ls <- terra::split(x = RAW_ls[[1]], f = Indices)
  rv_ls <- terra::split(x = RAW_ls[[2]], f = Indices)
  ret_ls <- pblapply(1:length(ru_ls), FUN = function(BASE_iter) {
    ret_rast <- sqrt(abs(ru_ls[[BASE_iter]])^2 + abs(rv_ls[[BASE_iter]])^2) # make m/s total windspeed
    ret_rast
  })
  Windspeed <- do.call(c, ret_ls)
  rm("RAW_ls", "ru_ls", "rv_ls", "ret_ls")
  terraOptions(memmax = 10)
  writeCDF(Windspeed[[1:nlyr(Windspeed) / 2]], filename = file.path(Dir.Data, "windspeed_RAW.nc"))
  terraOptions(memfrac = .9)
}

Data_ls$Windspeed <- sqrt(abs(Data_ls$`10m_u_component_of_wind`)^2 + abs(Data_ls$`10m_v_component_of_wind`)^2)

VNames_vec <- c(
  "Air Temperature [K]",
  "Skin Temperature [K]",
  "Eastward windspeed [m/s]",
  "Northward windspeed [m/s]",
  "Snow Cover [%]",
  "Snow Density [kg/m^3]",
  "Snow Depth [m]",
  "Snow Depth Water Equivalent [m]",
  "Snow Evaporation [m]",
  "Snowfall [m]",
  "Snowmelt [m]",
  "Temperature of Snow Layer [K]",
  "Windspeed [m/s]"
)

# LAYER-TIME Identification ===============================================
PreMonsoon_ls <- lapply(Data_ls, FUN = function(x) {
  x[[format(terra::time(Data_ls[[1]]), "%m") %in% c("03", "04", "05")]] # March - May
})

PostMonsoon_ls <- lapply(Data_ls, FUN = function(x) {
  x[[format(terra::time(Data_ls[[1]]), "%m") %in% c("09", "10", "11")]] # September - November
})

# TimeCheck <- c()
# for (Time_i in 1:nlyr(Data_ls[[1]])) {
#   Time_i <- terra::time(Data_ls[[1]])[Time_i]
#   print(Time_i)
#   TimeCheck <- c(
#     TimeCheck,
#     sum(rowSums((cbind(Time_i >= TimeWindows$Start, Time_i <= TimeWindows$Stop))) > 1)
#   )
# }
# TimeCheck <- as.logical(TimeCheck)
# sum(TimeCheck) / length(TimeCheck)

# CLIMATE CHANGE ==========================================================
ClimChange_ls <- pblapply(1:length(Data_ls), FUN = function(i) {
  # print(i)
  # print(Variables_vec[i])
  ## limit to region around peaks
  data_rast <- KrigR::Handle.Spatial(Data_ls[[i]], wide_buffer)

  ## Mean -----------------------------------------------------------------
  ## Make Yearly Aggregates
  MeanAnnual <- KrigR:::Temporal.Aggr(
    data_rast,
    "day", 1,
    "year", 1,
    FUN = mean,
    TZone = "UTC"
  )
  Var <- VNames_vec[[i]]

  ## Analyse change in region by 2-decade-intervals
  plot_df <- data.frame(
    Values = c(
      as.vector(terra::values(MeanAnnual[[1:20]])),
      as.vector(terra::values(MeanAnnual[[(nlyr(MeanAnnual) - 19):nlyr(MeanAnnual)]]))
    ),
    Time = rep(c("Beginning", "Ending"), each = length(as.vector(terra::values(MeanAnnual[[1:20]]))))
  )
  climatechange_timewindow <- ggplot(plot_df, aes(x = Time, y = Values)) +
    geom_violin() +
    geom_boxplot(width = 0.1) +
    stat_compare_means(comparisons = list(c("Beginning", "Ending")), paired = TRUE) +
    theme_bw() +
    scale_x_discrete(labels = c("Beginning" = "1951 - 1970", "Ending" = "2002 - 2021")) +
    labs(x = "Time-Windows", y = Var)

  ## Analyse change in region by year
  plot_df <- data.frame(
    Year = as.numeric(format(terra::time(MeanAnnual), "%Y")),
    mean = terra::global(MeanAnnual, mean, na.rm = TRUE),
    sd = terra::global(MeanAnnual, sd, na.rm = TRUE)
  )
  climatechange_region <- ggplot(plot_df, aes(x = Year, y = mean)) +
    geom_point() +
    stat_smooth(method = "lm") +
    theme_bw() +
    labs(y = paste("Annual", Var, "Mean"))

  # ## Analyse change at summits
  # means_summits <- terra::extract(MeanAnnual, eightks_sf, method = "bilinear")
  # rownames(means_summits) <- eightks_sf$ID
  # plot_df <- data.frame(
  #   Summit = rep(eightks_sf$ID, each = nlyr(MeanAnnual)),
  #   mean = unlist(as.vector(means_summits)[-1]),
  #   Year = as.numeric(format(terra::time(MeanAnnual), "%Y"))
  #   )
  # climchange_summits <- ggplot(plot_df, aes(x = Year, y = mean
  #                                           ,
  #                                           colour = as.factor(Summit)
  #                                           )) +
  #   geom_point() +
  #   stat_smooth(method = "lm", alpha = 0.1) +
  #   scale_colour_viridis_d() +
  #   theme_bw() +
  #   theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 1)) +
  #   labs(y = "Annual Mean", colour = "Summit ID")
  # climchange_summits

  save_gg <- cowplot::plot_grid(climatechange_timewindow,
    climatechange_region,
    nrow = 1
  )
  ggsave(save_gg,
    filename = file.path(Dir.Exports, paste0("ClimChange_", sub(" \\[.*", "", Var), ".png")),
    width = 32, height = 12, units = "cm"
  )

  ## Predictability -------------------------------------------------------
  ar <- function(x, lag = 1) {
    # Remove NA values
    ts <- x[!is.na(x)]

    # Check if there is sufficient data for the specified lag
    if (length(ts) <= lag) {
      return(NA) # Not enough data for the specified lag
    }

    # Create lagged time series
    ts_lagged <- ts[1:(length(ts) - lag)] # Exclude the last `lag` values
    ts_original <- ts[(lag + 1):length(ts)] # Exclude the first `lag` values

    # Compute the correlation
    cor(ts_lagged, ts_original, use = "complete.obs")
  }

  Predcitability_ls <- lapply(c(1, 2, 3, 5, 10), FUN = function(k) {
    ## ARs
    BeginAr <- app(data_rast[[format(terra::time(data_rast), "%Y") %in% as.character(1951:1970)]],
      fun = function(x) ar(x, lag = k)
    )

    EndAr <- app(data_rast[[format(terra::time(data_rast), "%Y") %in% as.character(2002:2021)]],
      fun = function(x) ar(x, lag = k)
    )

    ## map
    map_gg <- KrigR::Plot.SpatRast(EndAr - BeginAr,
      Dates = paste0("AR ", k, " Difference (", Var, ")"), SF = summits_sf, Shape = 2, Size = 3,
      Legend = paste("AR", k)
    ) +
      ggrepel::geom_text_repel(
        data = summits_df[summits_df$HEIGHTM >= 8000, ],
        aes(x = LON, y = LAT, label = PKNAME),
        max.overlaps = 30
      )

    ## boxplot
    plot_df <- data.frame(
      Values = c(terra::values(BeginAr), terra::values(EndAr)),
      Time = rep(c("Beginning", "Ending"), each = length(terra::values(BeginAr)))
    )
    bp_gg <- ggplot(plot_df, aes(x = Time, y = Values)) +
      geom_violin() +
      geom_boxplot(width = 0.1) +
      stat_compare_means(comparisons = list(c("Beginning", "Ending")), paired = TRUE) +
      theme_bw() +
      scale_x_discrete(labels = c("Beginning" = "1951 - 1970", "Ending" = "2002 - 2021")) +
      labs(x = "Time-Windows", y = paste0("AR ", k, " (", Var, ")"))

    ## combined plot
    save_gg <- cowplot::plot_grid(bp_gg,
      map_gg,
      nrow = 1
    )
    save_gg
  })

  ggsave(cowplot::plot_grid(plotlist = Predcitability_ls, ncol = 1),
    filename = file.path(Dir.Exports, paste0("Predictability_", sub(" \\[.*", "", Var), ".png")),
    width = 39.5, height = 46, units = "cm"
  )
})


# EXTREMES ================================================================
Seasons_ls <- list(
  `Pre-Monsoon` = PreMonsoon_ls,
  `Post-Monsoon` = PostMonsoon_ls
)

pblapply(1:length(Seasons_ls[[1]]), FUN = function(VarI) {
  VName <- VNames_vec[VarI]
  # message(VarI)
  # print(VName)

  Seasons_ggs <- lapply(1:length(Seasons_ls), FUN = function(SeasonI) {
    SeasonName <- names(Seasons_ls)[SeasonI]
    data_rast <- Seasons_ls[[SeasonName]][[VarI]]

    ## peaks
    data_df <- data.frame(t(terra::extract(data_rast, eightks_sf, method = "bilinear", df = TRUE)[, -1]))
    colnames(data_df) <- eightks_sf$PKNAME
    rownames(data_df) <- c() # terra::time(data_rast)

    quantiles_df <- apply(X = data_df, MARGIN = 2, FUN = quantile, na.rm = TRUE, probs = c(0.05, 0.95))


    Big_df <- tidyr::pivot_longer(data_df, cols = colnames(data_df))

    Big_df$time <- as.numeric(rep(format(terra::time(data_rast), "%Y"), each = ncol(data_df)))
    Big_df$QL <- as.numeric(quantiles_df[1, match(Big_df$name, colnames(quantiles_df))])
    Big_df$QU <- as.numeric(quantiles_df[2, match(Big_df$name, colnames(quantiles_df))])
    Big_df$Upper <- Big_df$value > Big_df$QU
    Big_df$Lower <- Big_df$value < Big_df$QL

    ## average exceeding lengths
    # Function to calculate average run length
    calculate_avg_run_length <- function(data, col_name) {
      rle_result <- rle(data[[col_name]])
      mean(rle_result$lengths[rle_result$values])
    }

    # Group by time and name, then calculate the average run length
    result <- Big_df %>%
      group_by(time, name) %>%
      summarise(
        avg_upper_run_length = calculate_avg_run_length(cur_data(), "Upper"),
        avg_lower_run_length = calculate_avg_run_length(cur_data(), "Lower"),
        .groups = "drop"
      )

    Lower <- aggregate(Big_df, Lower ~ time * name, FUN = sum)
    Lower$Bound <- "0.05"
    colnames(Lower) <- c("Year", "Peak", "Exceeded", "Quantile")
    Lower$Length <- result$avg_lower_run_length

    Upper <- aggregate(Big_df, Upper ~ time * name, FUN = sum)
    Upper$Bound <- "0.95"
    colnames(Upper) <- c("Year", "Peak", "Exceeded", "Quantile")
    Upper$Length <- result$avg_upper_run_length

    # combining data
    plot_df <- rbind(
      Lower,
      Upper
    )

    # number exceeding quantile bounds
    Legend <- get_legend(ggplot(plot_df, aes(
      x = as.numeric(Year), y = Exceeded,
      shape = Quantile, colour = Peak,
      lty = Quantile
    )) +
      geom_point() +
      stat_smooth(method = "lm", alpha = 0.05) +
      scale_colour_viridis_d() +
      theme_bw() +
      labs(x = "Year", y = "Days Exceeding Quantile Bounds", title = SeasonName) +
      theme(legend.position = "bottom"))


    save_gg <- plot_grid(
      ggplot(plot_df, aes(
        x = as.numeric(Year), y = Exceeded,
        shape = Quantile, colour = Peak,
        lty = Quantile
      )) +
        geom_point() +
        stat_smooth(method = "lm", alpha = 0.05) +
        scale_colour_viridis_d() +
        theme_bw() +
        labs(x = "Year", y = "Days Exceeding Quantile Bounds", title = SeasonName) +
        theme(legend.position = "none"),
      ggplot(plot_df, aes(
        x = as.numeric(Year), y = Length,
        shape = Quantile, colour = Peak,
        lty = Quantile
      )) +
        geom_point() +
        stat_smooth(method = "lm", alpha = 0.05) +
        scale_colour_viridis_d() +
        theme_bw() +
        labs(x = "Year", y = "Average Number of Consecutive Days Exceeding Quantile Bounds", title = SeasonName) +
        theme(legend.position = "none") +
        lims(y = c(0, NA)),
      ncol = 1
    )
    list(
      save_gg,
      Legend
    )
  })

  save_gg <- plot_grid(
    plot_grid(plotlist = lapply(Seasons_ggs, "[[", 1), nrow = 1),
    Seasons_ggs[[1]][[2]],
    ncol = 1, rel_heights = c(1, 0.1)
  )


  ggsave(save_gg,
    filename = file.path(Dir.Exports, paste0("Extremes_", sub(" \\[.*", "", VName), ".png")),
    width = 32 * 1.5, height = 24 * 1.5, units = "cm"
  )
})

# COMPOUND EVENTS =========================================================
## thresholds and ideas modelled after  ISBN-13 ‏ : ‎ 978-0071370264




# FUSING WITH EXPEDITION DATA =============================================




# KRIGING =================================================================
## Raw Data ---------------------
Data <- Data_ls$"2m_temperature"[[which(terra::time(Data_ls[["2m_temperature"]]) == "1996-05-10")]] # krakauer storm

## Covariates ---------------------
COP_DEM <- rast(file.path(Dir.Covariates, "COP30.tif"))
# COP_DEM <- terra::aggregate(COP_DEM, fact = 2) # have to upscale to circumvent "coordinate intervals are not constant" when creating spatialpixels data frame in krigR

## Kriging Buffer ---------------------
Peaks_buffer <- st_buffer(eightks_sf[, "ID"], dist = 30 * 1e3)
st_crs(Peaks_buffer) <- terra::crs(COP_DEM)
st_crs(eightks_sf) <- terra::crs(COP_DEM)

## Preparing Rasters -------------------
# Data <- crop(Data, Peaks_buffer)
# COP_DEM <- crop(COP_DEM, Peaks_buffer)
COP_coarse <- resample(COP_DEM, Data)
COP_fine <- KrigR:::Handle.Spatial(COP_DEM, Peaks_buffer)
rm(COP_DEM)

## Actual Kriging ---------------------
COP_krig <- Kriging(
  Data = Data,
  Covariates_training = COP_coarse,
  Covariates_target = COP_fine,
  Keep_Temporary = FALSE,
  Cores = 1, # parallel::detectCores(),
  Equation = "COP30",
  nmax = 40,
  FileName = paste0("DEMO", "_KRIG"),
  Dir = Dir.Exports
)

## Plotting ---------------------
COP_kriged <- COP_krig$Prediction - 273.15
PlotData <- crop(Data, ext(COP_kriged)) - 273.15

Map_gg <- ggplot() +
  geom_spatraster(data = PlotData) +
  geom_spatraster(data = COP_kriged) +
  geom_sf(data = st_union(Peaks_buffer), color = "black", fill = "transparent") +
  geom_sf(data = eightks_sf, shape = 2, colour = "white") +
  ggrepel::geom_text_repel(
    data = summits_df[summits_df$HEIGHTM >= 8000, ],
    aes(x = LON, y = LAT, label = PKNAME),
    max.overlaps = 30,
    colour = "white"
  ) +
  scale_fill_viridis_c(option = "A", na.value = "transparent", name = "[°C]") +
  labs(title = "Air Temperature on 1996-05-10", y = "Latitude", x = "Longitude") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "darkgrey", color = "darkgrey")) +
  theme(legend.position = "bottom", legend.key.width = unit(5, "cm"), legend.key.height = unit(1, "cm"))
Map_gg

ggsave(Map_gg,
  filename = file.path(Dir.Exports, "DEMO_Interpolation.png"),
  width = 16, height = 10
)
