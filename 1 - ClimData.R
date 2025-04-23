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
  "ggpubr", # for statistical comparisons in plots
  "brms", # for some bayesian models
  "tidybayes" # fore viz of bayesian models
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
eightks_sf <- summits_sf[summits_sf$HEIGHTM >= 7000, ]
sevenks_sf <- summits_sf[summits_sf$HEIGHTM >= 7000, ]
narrow_buffer <- KrigR::Buffer.pts(eightks_sf, 2e4) # equates to roughly three grid cells in either direction
wide_buffer <- KrigR::Buffer.pts(eightks_sf, 1e5) # equates to roughly three grid cells in either direction

## ERA5-Land --------------------------------------------------------------
message("#### Raw Data ############################################")
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
  message(Var_Iter)
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
Variables_vec <- c(Variables_vec, "windspeed")

VNames_vec <- VNames_vec[c(1, 5, 7, 10, 13)]
Variables_vec <- Variables_vec[c(1, 5, 7, 10, 13)]
Data_ls <- Data_ls[c(1, 5, 7, 10, 13)]

# LAYER-TIME Identification ===============================================
PreMonsoon_ls <- lapply(Data_ls, FUN = function(x) {
  x[[format(terra::time(Data_ls[[1]]), "%m") %in% c("03", "04", "05")]] # March - May
})

PostMonsoon_ls <- lapply(Data_ls, FUN = function(x) {
  x[[format(terra::time(Data_ls[[1]]), "%m") %in% c("09", "10", "11")]] # September - November
})

# CLIMATE CHANGE ==========================================================
message("#### Climate Change ############################################")
ClimChange_ls <- for (i in 1:length(Data_ls)) {
  # i = 1
  # print(i)
  message(Variables_vec[i])

  Dir.Var <- file.path(Dir.Exports, Variables_vec[i])
  if (!dir.exists(Dir.Var)) {
    dir.create(Dir.Var)
  }
  FName <- file.path(Dir.Var, paste0("ClimChange_", sub(" \\[.*", "", Variables_vec[i]), ".png"))

  if (file.exists(paste0(tools::file_path_sans_ext(FName), "_Predictability.png"))) {
    next()
  }
  ## limit to region around peaks
  data_rast <- KrigR::Handle.Spatial(Data_ls[[i]], wide_buffer)

  ## Mean -----------------------------------------------------------------
  print("Mean Change")
  ## Make Yearly Aggregates
  MeanAnnual <- KrigR:::Temporal.Aggr(
    data_rast,
    "day", 1,
    "year", 1,
    FUN = mean,
    TZone = "UTC"
  )
  SDAnnual <- KrigR:::Temporal.Aggr(
    data_rast,
    "day", 1,
    "year", 1,
    FUN = sd,
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

  ## Build models of change linearly
  ModelMean <- brms::brm(formula = mean ~ Year, data = plot_df)
  ModelSD <- brms::brm(formula = sd ~ Year, data = plot_df)

  modelplot_df <- data.frame(
    Value = c(
      unlist(lapply(brms::as_draws(ModelMean), "[[", "b_Year")),
      unlist(lapply(brms::as_draws(ModelSD), "[[", "b_Year"))
    ),
    Outcome = rep(c("Mean", "SD"), each = length(as_draws(ModelSD)[[1]][[1]]) * length(as_draws(ModelSD)))
  )

  BRM_region_gg <- ggplot(modelplot_df, aes(y = Outcome, x = Value)) +
    stat_halfeye() +
    theme_bw() +
    geom_vline(xintercept = 0)

  save_gg <- plot_grid(
    plot_grid(climatechange_timewindow, climatechange_region, nrow = 1),
    BRM_region_gg,
    nrow = 2
  )

  ggsave(save_gg,
    filename = paste0(tools::file_path_sans_ext(FName), "_Region.png"),
    width = 32, height = 24, units = "cm"
  )

  ## Analyse change at summits
  message("Change at Summits")
  means_summits <- t(terra::extract(MeanAnnual, eightks_sf, method = "bilinear", fun = "mean"))[-1, ]
  sd_summits <- t(terra::extract(SDAnnual, eightks_sf, method = "bilinear", fun = "mean"))[-1, ]
  colnames(means_summits) <- eightks_sf$PKNAME
  plot_df <- data.frame(
    Summit = rep(colnames(means_summits), each = nlyr(MeanAnnual)),
    mean = unlist(as.vector(means_summits)),
    sd = unlist(as.vector(sd_summits)),
    Year = as.numeric(format(terra::time(MeanAnnual), "%Y"))
  )

  labelInfo <- split(plot_df, plot_df$Summit)
  labelInfo <- lapply(labelInfo, function(dat) {
    meanLabel <- predict(lm(mean ~ Year, data = dat), newdata = data.frame(Year = max(dat$Year)))
    sdLabel <- predict(lm(sd ~ Year, data = dat), newdata = data.frame(Year = max(dat$Year)))
    summit <- unique(dat$Summit)
    data.frame(meanLabel = meanLabel, sdLabel = sdLabel, Summit = summit)
  })
  labelInfo <- do.call(rbind, labelInfo)

  ### Line plots of trends at summits
  Mean_gg <- ggplot(plot_df, aes(x = Year, y = mean)) +
    geom_smooth(method = "lm", fill = "#480472", col = "#480472") +
    geom_smooth(aes(group = Summit), col = "#535353", method = "lm", alpha = 0.2) +
    geom_label_repel(
      data = labelInfo,
      aes(
        x = max(plot_df$Year), y = meanLabel,
        label = Summit
      ),
      color = "#535353",
      nudge_x = 7
    ) +
    theme_bw() +
    labs(y = "Mean")

  SD_gg <- ggplot(plot_df, aes(x = Year, y = sd)) +
    geom_smooth(method = "lm", fill = "#480472", col = "#480472") +
    geom_smooth(aes(group = Summit), col = "#535353", method = "lm", alpha = 0.2) +
    geom_label_repel(
      data = labelInfo,
      aes(
        x = max(plot_df$Year), y = sdLabel,
        label = Summit
      ),
      color = "#535353",
      nudge_x = 7
    ) +
    theme_bw() +
    labs(y = "Standard Deviation")

  ### BRMS models of trends at all summits and indvidual summits
  ModelMean <- brms::brm(formula = mean ~ Year, data = plot_df)
  ModelSD <- brms::brm(formula = sd ~ Year, data = plot_df)

  modelplot_df <- data.frame(
    Value = c(
      unlist(lapply(brms::as_draws(ModelMean), "[[", "b_Year")),
      unlist(lapply(brms::as_draws(ModelSD), "[[", "b_Year"))
    ),
    Outcome = rep(c("Mean", "SD"), each = length(as_draws(ModelSD)[[1]][[1]]) * length(as_draws(ModelSD))),
    Summit = "ALL"
  )

  Peaks_gg <- lapply(unique(plot_df$Summit), function(PeakIter) {
    Iter_df <- plot_df[plot_df$Summit == PeakIter, ]

    ModelMean <- brms::brm(formula = mean ~ Year, data = Iter_df)
    ModelSD <- brms::brm(formula = sd ~ Year, data = Iter_df)
    modelplot_df <- data.frame(
      Mean = unlist(lapply(brms::as_draws(ModelMean), "[[", "b_Year")),
      SD = unlist(lapply(brms::as_draws(ModelSD), "[[", "b_Year"))
    )

    modelplot_df <- data.frame(
      Value = c(
        unlist(lapply(brms::as_draws(ModelMean), "[[", "b_Year")),
        unlist(lapply(brms::as_draws(ModelSD), "[[", "b_Year"))
      ),
      Outcome = rep(c("Mean", "SD"), each = length(as_draws(ModelSD)[[1]][[1]]) * length(as_draws(ModelSD))),
      Summit = PeakIter
    )
  })
  Peaks_gg <- do.call(rbind, Peaks_gg)
  modelplot_df <- rbind(modelplot_df, Peaks_gg)
  StatSig <- aggregate(Value ~ Summit + Outcome, modelplot_df, FUN = quantile, c(0.05, 0.95))
  StatSig$Direction <- sign(StatSig$Value[, 1]) + sign(StatSig$Value[, 2])
  StatSig$Sig <- abs(StatSig$Direction) == 2
  modelplot_df$StatSig <- NA
  for (i in 1:nrow(modelplot_df)) {
    modelplot_df$StatSig[i] <- StatSig$Sig[intersect(which(StatSig$Outcome == modelplot_df$Outcome[i]), which(StatSig$Summit == modelplot_df$Summit[i]))]
  }
  modelplot_df$col <- ifelse(modelplot_df$StatSig, "green", "red")
  modelplot_df$fill <- ifelse(modelplot_df$Summit == "ALL", "#480472", "#535353")

  BMean_gg <-
    ggplot(
      modelplot_df[modelplot_df$Outcome == "Mean", ],
      aes(
        y = factor(Summit, levels = rev(c("ALL", rev(labelInfo$Summit[order(labelInfo$meanLabel)])))),
        x = Value
      )
    ) +
    stat_halfeye(aes(fill = fill)) +
    scale_fill_manual(values = c("#480472", "#535353")) +
    geom_boxplot(aes(col = col), width = 0.3, lwd = 1.1) +
    scale_color_manual(values = c("#003b05", "#5c0000"), breaks = c("green", "red")) +
    geom_vline(xintercept = 0) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "BRMS Model Coefficient Posterior Samples", y = "")

  BSD_gg <-
    ggplot(
      modelplot_df[modelplot_df$Outcome == "SD", ],
      aes(
        y = factor(Summit, levels = rev(c("ALL", rev(labelInfo$Summit[order(labelInfo$sdLabel)])))),
        x = Value
      )
    ) +
    stat_halfeye(aes(fill = fill)) +
    scale_fill_manual(values = c("#480472", "#535353")) +
    geom_boxplot(aes(col = col), width = 0.3, lwd = 1.1) +
    scale_color_manual(values = c("#003b05", "#5c0000"), breaks = c("green", "red")) +
    geom_vline(xintercept = 0) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "BRMS Model Coefficient Posterior Samples", y = "")

  ### Fusing and final export
  ggsave(plot_grid(Mean_gg, BMean_gg, SD_gg, BSD_gg, ncol = 2),
    filename = paste0(tools::file_path_sans_ext(FName), "_Summits.png"),
    width = 32 * 1.2, height = 22 * 1.2, units = "cm"
  )

  ## Predictability -------------------------------------------------------
  message("Predictability Change")
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

  Predcitability_ls <- pblapply(c(1, 2, 3, 5, 10), FUN = function(k) {
    # k = 1
    ## ARs
    BeginAr <- app(data_rast[[format(terra::time(data_rast), "%Y") %in% as.character(1951:1970)]],
      fun = function(x) ar(x, lag = k)
    )

    EndAr <- app(data_rast[[format(terra::time(data_rast), "%Y") %in% as.character(2002:2021)]],
      fun = function(x) ar(x, lag = k)
    )

    ## map
    map_gg <- KrigR::Plot.SpatRast(EndAr - BeginAr,
      Dates = paste0("AR ", k, " Difference (", Var, ")"), SF = eightks_sf, Shape = 2, Size = 3,
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
    list(Box = bp_gg, Map = EndAr - BeginAr)
    # save_gg <- cowplot::plot_grid(bp_gg,
    #   map_gg,
    #   nrow = 1
    # )
    # save_gg
  })

  map_gg <- Plot.SpatRast(do.call(c, lapply(Predcitability_ls, "[[", "Map")),
    Dates = paste0("AR ", c(1, 2, 3, 5, 10), " Difference (", Var, ")"), SF = eightks_sf, Shape = 2, Size = 3,
    Legend = paste("AR", k),
    ncol = 1
  ) +
    ggrepel::geom_text_repel(
      data = summits_df[summits_df$HEIGHTM >= 8000, ],
      aes(x = LON, y = LAT, label = PKNAME),
      max.overlaps = 30
    )

  ggsave(plot_grid(plot_grid(plotlist = lapply(Predcitability_ls, "[[", "Box"), ncol = 1), map_gg, ncol = 2),
    filename = paste0(tools::file_path_sans_ext(FName), "_Predictability.png"),
    width = 45, height = 46, units = "cm"
  )
}

# SEASONS & EXTREMES ======================================================
message("#### Seasons & Climate Extremes ################################")
Seasons_ls <- list(
  `Pre-Monsoon` = PreMonsoon_ls,
  `Post-Monsoon` = PostMonsoon_ls
)

for (VarI in 1:length(Seasons_ls[[1]])) {
  # VarI <- 1
  VName <- VNames_vec[VarI]
  message(VName)
  Dir.Var <- file.path(Dir.Exports, Variables_vec[VarI])
  if (!dir.exists(Dir.Var)) {
    dir.create(Dir.Var)
  }
  FName <- file.path(Dir.Var, paste0("Extremes_", sub(" \\[.*", "", Variables_vec[VarI]), ".png"))
  print(FName)

  if (file.exists(paste0(tools::file_path_sans_ext(FName), "_Length.png"))) {
    next()
  }

  Seasons_ggs <- lapply(1:length(Seasons_ls), FUN = function(SeasonI) {
    # SeasonI = 1
    SeasonName <- names(Seasons_ls)[SeasonI]
    print(SeasonName)
    data_rast <- Seasons_ls[[SeasonName]][[VarI]]

    ## climate change in seasons
    data_rast <- KrigR::Handle.Spatial(data_rast, wide_buffer)

    ## Mean -----------------------------------------------------------------
    print("Mean Change")
    ## Make Yearly Aggregates
    MeanAnnual <- terra::tapp(data_rast, index = as.numeric(substr(time(data_rast), 1, 4)), fun = mean, na.rm = TRUE)
    SDAnnual <- terra::tapp(data_rast, index = as.numeric(substr(time(data_rast), 1, 4)), fun = sd, na.RM = TRUE)
    time(MeanAnnual) <- time(SDAnnual) <- as.POSIXct(paste0(unique(as.numeric(substr(time(data_rast), 1, 4))), "-01-01"))
    Var <- VNames_vec[[i]]

    ## Analyse change in region by 2-decade-intervals
    Climchange_df <- data.frame(
      Values = c(
        as.vector(terra::values(MeanAnnual[[1:20]])),
        as.vector(terra::values(MeanAnnual[[(nlyr(MeanAnnual) - 19):nlyr(MeanAnnual)]]))
      ),
      Time = rep(c("Beginning", "Ending"), each = length(as.vector(terra::values(MeanAnnual[[1:20]])))),
      Season = SeasonName
    )
    climatechange_timewindow <- ggplot(plot_df, aes(x = Time, y = Values)) +
      geom_violin() +
      geom_boxplot(width = 0.1) +
      stat_compare_means(comparisons = list(c("Beginning", "Ending")), paired = TRUE) +
      theme_bw() +
      scale_x_discrete(labels = c("Beginning" = "1951 - 1970", "Ending" = "2002 - 2021")) +
      labs(x = "Time-Windows", y = Var)

    ## Analyse change in region by year
    ClimchangeRegion_df <- plot_df <- data.frame(
      Year = as.numeric(format(terra::time(MeanAnnual), "%Y")),
      mean = terra::global(MeanAnnual, mean, na.rm = TRUE),
      sd = terra::global(MeanAnnual, sd, na.rm = TRUE),
      Season = SeasonName
    )
    climatechange_region <- ggplot(plot_df, aes(x = Year, y = mean)) +
      geom_point() +
      stat_smooth(method = "lm") +
      theme_bw() +
      labs(y = paste("Annual", Var, "Mean"))

    ## Build models of change linearly
    ModelMean <- brms::brm(formula = mean ~ Year, data = plot_df)
    ModelSD <- brms::brm(formula = sd ~ Year, data = plot_df)

    ClimchangeModel_df <- modelplot_df <- data.frame(
      Value = c(
        unlist(lapply(brms::as_draws(ModelMean), "[[", "b_Year")),
        unlist(lapply(brms::as_draws(ModelSD), "[[", "b_Year"))
      ),
      Outcome = rep(c("Mean", "SD"), each = length(as_draws(ModelSD)[[1]][[1]]) * length(as_draws(ModelSD))),
      Season = SeasonName
    )

    BRM_region_gg <- ggplot(modelplot_df, aes(y = Outcome, x = Value)) +
      stat_halfeye() +
      theme_bw() +
      geom_vline(xintercept = 0)

    climchange_gg <- plot_grid(climatechange_timewindow, climatechange_region, BRM_region_gg, nrow = 1) ## needs returning to be fused with

    ## Analyse change at summits
    message("Change at Summits")
    means_summits <- t(terra::extract(MeanAnnual, eightks_sf, method = "bilinear", fun = "mean"))[-1, ]
    sd_summits <- t(terra::extract(SDAnnual, eightks_sf, method = "bilinear", fun = "mean"))[-1, ]
    colnames(means_summits) <- eightks_sf$PKNAME
    SummitsChange_df <- plot_df <- data.frame(
      Summit = rep(colnames(means_summits), each = nlyr(MeanAnnual)),
      mean = unlist(as.vector(means_summits)),
      sd = unlist(as.vector(sd_summits)),
      Year = as.numeric(format(terra::time(MeanAnnual), "%Y")),
      Season = SeasonName
    )

    labelInfo <- split(plot_df, plot_df$Summit)
    labelInfo <- lapply(labelInfo, function(dat) {
      meanLabel <- predict(lm(mean ~ Year, data = dat), newdata = data.frame(Year = max(dat$Year)))
      sdLabel <- predict(lm(sd ~ Year, data = dat), newdata = data.frame(Year = max(dat$Year)))
      summit <- unique(dat$Summit)
      data.frame(meanLabel = meanLabel, sdLabel = sdLabel, Summit = summit)
    })
    labelInfo <- do.call(rbind, labelInfo)
    labelInfo$Season <- SeasonName
    SummitLabels <- labelInfo

    ### Line plots of trends at summits
    Mean_gg <- ggplot(plot_df, aes(x = Year, y = mean)) +
      geom_smooth(method = "lm", fill = "#480472", col = "#480472") +
      geom_smooth(aes(group = Summit), col = "#535353", method = "lm", alpha = 0.2) +
      geom_label_repel(
        data = labelInfo,
        aes(
          x = max(plot_df$Year), y = meanLabel,
          label = Summit
        ),
        color = "#535353",
        nudge_x = 7
      ) +
      theme_bw() +
      labs(y = "Mean")

    SD_gg <- ggplot(plot_df, aes(x = Year, y = sd)) +
      geom_smooth(method = "lm", fill = "#480472", col = "#480472") +
      geom_smooth(aes(group = Summit), col = "#535353", method = "lm", alpha = 0.2) +
      geom_label_repel(
        data = labelInfo,
        aes(
          x = max(plot_df$Year), y = sdLabel,
          label = Summit
        ),
        color = "#535353",
        nudge_x = 7
      ) +
      theme_bw() +
      labs(y = "Standard Deviation")

    ### BRMS models of trends at all summits and indvidual summits
    ModelMean <- brms::brm(formula = mean ~ Year, data = plot_df)
    ModelSD <- brms::brm(formula = sd ~ Year, data = plot_df)

    modelplot_df <- data.frame(
      Value = c(
        unlist(lapply(brms::as_draws(ModelMean), "[[", "b_Year")),
        unlist(lapply(brms::as_draws(ModelSD), "[[", "b_Year"))
      ),
      Outcome = rep(c("Mean", "SD"), each = length(as_draws(ModelSD)[[1]][[1]]) * length(as_draws(ModelSD))),
      Summit = "ALL",
      Season = SeasonName
    )

    Peaks_gg <- lapply(unique(plot_df$Summit), function(PeakIter) {
      Iter_df <- plot_df[plot_df$Summit == PeakIter, ]

      ModelMean <- brms::brm(formula = mean ~ Year, data = Iter_df)
      ModelSD <- brms::brm(formula = sd ~ Year, data = Iter_df)
      modelplot_df <- data.frame(
        Mean = unlist(lapply(brms::as_draws(ModelMean), "[[", "b_Year")),
        SD = unlist(lapply(brms::as_draws(ModelSD), "[[", "b_Year"))
      )

      modelplot_df <- data.frame(
        Value = c(
          unlist(lapply(brms::as_draws(ModelMean), "[[", "b_Year")),
          unlist(lapply(brms::as_draws(ModelSD), "[[", "b_Year"))
        ),
        Outcome = rep(c("Mean", "SD"), each = length(as_draws(ModelSD)[[1]][[1]]) * length(as_draws(ModelSD))),
        Summit = PeakIter,
        Season = SeasonName
      )
    })
    Peaks_gg <- do.call(rbind, Peaks_gg)
    modelplot_df <- rbind(modelplot_df, Peaks_gg)
    StatSig <- aggregate(Value ~ Summit + Outcome, modelplot_df, FUN = quantile, c(0.05, 0.95))
    StatSig$Direction <- sign(StatSig$Value[, 1]) + sign(StatSig$Value[, 2])
    StatSig$Sig <- abs(StatSig$Direction) == 2
    modelplot_df$StatSig <- NA
    for (i in 1:nrow(modelplot_df)) {
      modelplot_df$StatSig[i] <- StatSig$Sig[intersect(which(StatSig$Outcome == modelplot_df$Outcome[i]), which(StatSig$Summit == modelplot_df$Summit[i]))]
    }
    modelplot_df$col <- ifelse(modelplot_df$StatSig, "green", "red")
    modelplot_df$fill <- ifelse(modelplot_df$Summit == "ALL", "#480472", "#535353")
    SummitsChangemodelplot_df <- modelplot_df

    BMean_gg <-
      ggplot(
        modelplot_df[modelplot_df$Outcome == "Mean", ],
        aes(
          y = factor(Summit, levels = rev(c("ALL", rev(labelInfo$Summit[order(labelInfo$meanLabel)])))),
          x = Value
        )
      ) +
      stat_halfeye(aes(fill = fill)) +
      scale_fill_manual(values = c("#480472", "#535353")) +
      geom_boxplot(aes(col = col), width = 0.3, lwd = 1.1) +
      scale_color_manual(values = c("#003b05", "#5c0000"), breaks = c("green", "red")) +
      geom_vline(xintercept = 0) +
      theme_bw() +
      theme(legend.position = "none") +
      labs(x = "BRMS Model Coefficient Posterior Samples", y = "")

    BSD_gg <-
      ggplot(
        modelplot_df[modelplot_df$Outcome == "SD", ],
        aes(
          y = factor(Summit, levels = rev(c("ALL", rev(labelInfo$Summit[order(labelInfo$sdLabel)])))),
          x = Value
        )
      ) +
      stat_halfeye(aes(fill = fill)) +
      scale_fill_manual(values = c("#480472", "#535353")) +
      geom_boxplot(aes(col = col), width = 0.3, lwd = 1.1) +
      scale_color_manual(values = c("#003b05", "#5c0000"), breaks = c("green", "red")) +
      geom_vline(xintercept = 0) +
      theme_bw() +
      theme(legend.position = "none") +
      labs(x = "BRMS Model Coefficient Posterior Samples", y = "")

    SummitsMean_gg <- plot_grid(Mean_gg, BMean_gg, ncol = 2)
    SummitsSD_gg <- plot_grid(SD_gg, BSD_gg, ncol = 2)

    ## Predictability -------------------------------------------------------
    message("Predictability Change")
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

    Predcitability_ls <- pblapply(c(1, 2, 3, 5, 10), FUN = function(k) {
      # k = 1
      ## ARs
      BeginAr <- app(data_rast[[format(terra::time(data_rast), "%Y") %in% as.character(1951:1970)]],
        fun = function(x) ar(x, lag = k)
      )

      EndAr <- app(data_rast[[format(terra::time(data_rast), "%Y") %in% as.character(2002:2021)]],
        fun = function(x) ar(x, lag = k)
      )

      ## map
      map_gg <- KrigR::Plot.SpatRast(EndAr - BeginAr,
        Dates = paste0("AR ", k, " Difference (", Var, ")"), SF = eightks_sf, Shape = 2, Size = 3,
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
      list(Box = bp_gg, Map = EndAr - BeginAr)
      # save_gg <- cowplot::plot_grid(bp_gg,
      #   map_gg,
      #   nrow = 1
      # )
      # save_gg
    })

    map_gg <- Plot.SpatRast(do.call(c, lapply(Predcitability_ls, "[[", "Map")),
      Dates = paste0("AR ", c(1, 2, 3, 5, 10), " Difference (", Var, ")"), SF = eightks_sf, Shape = 2, Size = 3,
      Legend = paste("AR", k),
      ncol = 1
    ) +
      ggrepel::geom_text_repel(
        data = summits_df[summits_df$HEIGHTM >= 8000, ],
        aes(x = LON, y = LAT, label = PKNAME),
        max.overlaps = 30
      )
    Predictability_gg <- plot_grid(plot_grid(plotlist = lapply(Predcitability_ls, "[[", "Box"), ncol = 1), map_gg, ncol = 2)

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
    ### Line plots of trends at summits
    Quant_ls <- lapply(unique(plot_df$Quantile), FUN = function(QuanIter) {
      # print(QuanIter)
      Iter_df <- Iter_df <- plot_df[plot_df$Quantile == QuanIter, ]
      labelInfo <- split(Iter_df, Iter_df$Peak)
labelInfo <- lapply(labelInfo, function(dat) {
    ELabel <- predict(lm(Exceeded ~ Year, data = dat), newdata = data.frame(Year = max(dat$Year)))
    LLabel <- predict(lm(Length ~ Year, data = dat), newdata = data.frame(Year = max(dat$Year)))
    summit <- unique(dat$Peak)
    data.frame(ELabel = ELabel, LLabel = LLabel, Summit = summit)
})
labelInfo <- do.call(rbind, labelInfo)

ELine_gg <- ggplot(Iter_df, aes(x = Year, y = Exceeded)) +
    geom_smooth(aes(group = Peak), col = "#535353", method = "lm", alpha = 0.2) +
    geom_smooth(method = "lm", fill = "#480472", col = "#480472") +
    geom_label_repel(
        data = labelInfo,
        aes(
            x = max(Iter_df$Year), y = ELabel,
            label = Summit
        ),
        color = "#535353",
        nudge_x = 7
    ) +
    theme_bw() +
    labs(x = "Year", y = "# Days Exceeding Quantile Bounds", title = paste(SeasonName, QuanIter, "Quantile")) +
    expand_limits(y = 0)

LLine_gg <- ggplot(Iter_df, aes(x = Year, y = Length)) +
    geom_smooth(aes(group = Peak), col = "#535353", method = "lm", alpha = 0.2) +
    geom_smooth(method = "lm", fill = "#480472", col = "#480472") +
    geom_label_repel(
        data = labelInfo,
        aes(
            x = max(Iter_df$Year), y = LLabel,
            label = Summit
        ),
        color = "#535353",
        nudge_x = 7
    ) +
    theme_bw() +
    labs(x = "Year", y = "Continuous Days Exceeding Quantile Bounds", title = paste(SeasonName, QuanIter, "Quantile")) +
    expand_limits(y = 0)

      ### BRMS models of trends at all summits and indvidual summits
      Big_df$Exceeded <- as.numeric(data.frame(Big_df)[, ifelse(QuanIter == "0.05", "Lower", "Upper")])
      ModelExceeded <- brms::brm(formula = Exceeded ~ time, data = Big_df, family = bernoulli(link = "logit")) # this really should be each day 0 or 1 for extreme or not
      run_lengths_df <- Big_df %>% # calculate length of each exceeding run (multiple per year)
        group_by(time) %>% 
  arrange(time) %>%
  summarise(runs = list({
    r <- rle(Exceeded)
    lengths <- r$lengths[r$values == TRUE]
    if (length(lengths) == 0) NA else lengths
  })) %>%
  unnest_longer(runs) %>%
  filter(!is.na(runs))
ModelLength <- brms::brm(formula = runs ~ time, data = run_lengths_df, family = "poisson")

modelplot_df <- data.frame(
    Value = c(
        unlist(lapply(brms::as_draws(ModelExceeded), "[[", "b_time")),
        unlist(lapply(brms::as_draws(ModelLength), "[[", "b_time"))
    ),
    Outcome = rep(c("Exceeded", "Length"), each = length(as_draws(ModelExceeded)[[1]][[1]]) * length(as_draws(ModelExceeded))),
    Summit = "ALL"
)

Peaks_gg <- lapply(unique(Big_df$name), function(PeakIter) {
    message(PeakIter)
    message(QuanIter)
    message(SeasonName)

    Peak_df <- Big_df[Big_df$name == PeakIter, ]
    run_lengths_df <- Peak_df %>% # calculate length of each exceeding run (multiple per year)
  group_by(time) %>% 
  arrange(time) %>%
  summarise(runs = list({
    r <- rle(Exceeded)
    lengths <- r$lengths[r$values == TRUE]
    if (length(lengths) == 0) NA else lengths
  })) %>%
  unnest_longer(runs) %>%
  filter(!is.na(runs))

    if (sum(Peak_df$Exceeded != 0) > 3) {
          E_BRM <- brms::brm(formula = Exceeded ~ time, data = Peak_df, family = bernoulli(link = "logit")) # this really should be each day 0 or 1 for extreme or not
          Exceeded <- exp(unlist(lapply(brms::as_draws(E_BRM), "[[", "b_time")))
        } else {
          Exceeded <- rep(NA, 4000)
        }
        if (nrow(run_lengths_df)>3) {
          L_BRM <- exp(brms::brm(formula = runs ~ time, data = run_lengths_df, family = "poisson"))
          Length <- unlist(lapply(brms::as_draws(L_BRM), "[[", "b_time"))
        } else {
          Length <- rep(NA, 4000)
        }

      data.frame(
      Value = c(
        Exceeded,
        Length
      ),
      Outcome = c(rep("Exceeded", length(Exceeded)), rep("Length", length(Length))),
      Summit = PeakIter
        )
})
names(Peaks_gg) <- unique(Iter_df$Peak)

Peaks_gg <- do.call(rbind, Peaks_gg)
modelplot_df <- rbind(modelplot_df, Peaks_gg)
StatSig <- aggregate(Value ~ Summit + Outcome, modelplot_df, FUN = quantile, c(0.05, 0.95))
StatSig$Direction <- sign(StatSig$Value[, 1]) + sign(StatSig$Value[, 2])
StatSig$Sig <- abs(StatSig$Direction) == 2
modelplot_df$StatSig <- NA
for (i in 1:nrow(modelplot_df)) {
    modelplot_df$StatSig[i] <- StatSig$Sig[intersect(which(StatSig$Outcome == modelplot_df$Outcome[i]), which(StatSig$Summit == modelplot_df$Summit[i]))]
}
modelplot_df$col <- ifelse(modelplot_df$StatSig, "green", "red")
modelplot_df$fill <- ifelse(modelplot_df$Summit == "ALL", "#480472", "#535353")

BExceeded_gg <-
    ggplot(
        modelplot_df[modelplot_df$Outcome == "Exceeded", ],
        aes(
            y = factor(Summit, levels = rev(c("ALL", rev(labelInfo$Summit[order(labelInfo$ELabel)])))),
            x = Value
        )
    ) +
    stat_halfeye(aes(fill = fill)) +
    scale_fill_manual(values = c("#480472", "#535353")) +
    geom_boxplot(aes(col = col), width = 0.3, lwd = 1.1) +
    scale_color_manual(values = c("#003b05", "#5c0000"), breaks = c("green", "red")) +
    geom_vline(xintercept = 0) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "BRMS Model Coefficient Posterior Samples", y = "")

BLength_gg <-
    ggplot(
        modelplot_df[modelplot_df$Outcome == "Length", ],
        aes(
            y = factor(Summit, levels = rev(c("ALL", rev(labelInfo$Summit[order(labelInfo$LLabel)])))),
            x = Value
        )
    ) +
    stat_halfeye(aes(fill = fill)) +
    scale_fill_manual(values = c("#480472", "#535353")) +
    geom_boxplot(aes(col = col), width = 0.3, lwd = 1.1) +
    scale_color_manual(values = c("#003b05", "#5c0000"), breaks = c("green", "red")) +
    geom_vline(xintercept = 0) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "BRMS Model Coefficient Posterior Samples", y = "")

SummitsExceeded_gg <- plot_grid(ELine_gg, BExceeded_gg, ncol = 2)
SummitsLength_gg <- plot_grid(LLine_gg, BLength_gg, ncol = 2)

## return plots
list(
  Predictability = Predictability_gg,
  Exceeded = SummitsExceeded_gg,
  Length = SummitsLength_gg,
)
      }) # quantile loop
      names(Quant_ls) <- unique(plot_df$Quantile)

      list(ClimChange = climchange_gg,
      Mean = SummitsMean_gg,
      SD = SummitsSD_gg,
      Quants = Quant_ls
      )
  }) # season loop
  stop("Fuse and export plots")
  FName
} # variable loop


# COMPOUND EVENTS =========================================================
# abandoned as data does not show compound events well
# message("#### Compound Events #########################################")
# ## thresholds and ideas modelled after  ISBN-13 ‏ : ‎ 978-0071370264
# if (file.exists(file.path(Dir.Data, "blizzard_binary.nc")) & file.exists(file.path(Dir.Data, "blizzard_continuous.nc"))) {
#   blizzard_Continuous <- rast(file.path(Dir.Data, "blizzard_continuous.nc"))
#   blizzard_Binary <- rast(file.path(Dir.Data, "blizzard_binary.nc"))
# } else {
#   ## Snowfall threshold ---------------------
#   snow_raw <- rast(file.path(Dir.Data, "snowfall_RAW.nc"))
#   snow_cumsum <- roll(snow_raw, 24, "sum", type = "to")
#   snow_thresh <- snow_cumsum / 0.1 # m; anything greater than 1 is a storm: https://www.weather.gov/media/meg/WinterStormCriteriaMEG.pdf

#   ## Wind threshold ---------------------
#   wind_raw <- rast(file.path(Dir.Data, "windspeed_RAW.nc"))
#   wind_raw <- wind_raw[[time(wind_raw) %in% time(snow_raw)]]
#   wind_thresh <- wind_raw / 7.7 # m/s; anything greater than 1 is a storm: DOI: 10.1175/1520-0450(1997)036<0205:EOTWSF>2.0.CO;2

#   ## Blizzard conditions ---------------------
#   if (file.exists(file.path(Dir.Data, "blizzard_binary.nc"))) {
#     blizzard_Binary <- rast(file.path(Dir.Data, "blizzard_binary.nc"))
#   } else {
#     blizzard_Binary <- ((snow_thresh > 1) + (wind_thresh > 1)) == 2
#     time(blizzard_Binary) <- time(snow_thresh)
#     blizzard_Binary <- ClimHub:::WriteRead.NC(
#       NC = blizzard_Binary,
#       FName = file.path(Dir.Data, "blizzard_binary.nc"),
#       Variable = "BlizzardConditions",
#       LongVar = "BlizzardConditions",
#       Unit = "",
#       Attrs = NULL,
#       Compression = NA,
#       Write = TRUE
#     )
#   }

#   if (file.exists(file.path(Dir.Data, "blizzard_continuous.nc"))) {
#     blizzard_Continuous <- rast(file.path(Dir.Data, "blizzard_continuous.nc"))
#   } else {
#     blizzard_Continuous <- (snow_thresh + wind_thresh) / 2
#     time(blizzard_Continuous) <- time(snow_thresh)
#     blizzard <- Continuous <- ClimHub:::WriteRead.NC(
#       NC = blizzard_Continuous,
#       FName = file.path(Dir.Data, "blizzard_continuous.nc"),
#       Variable = "BlizzardConditions",
#       LongVar = "BlizzardConditions",
#       Unit = "",
#       Attrs = NULL,
#       Compression = NA,
#       Write = TRUE
#     )
#   }
# }

# ## Analysis for peaks ---------------------
# binary_summits <- t(blizzard_Binary[
#   terra::cellFromXY(blizzard_Binary[[1]], sf::st_coordinates(eightks_sf))
# ])
# # binary_summits <- t(terra::extract(blizzard_Binary, eightks_sf, method = "bilinear", fun = "mean"))[-1, ]

# continuous_summits <- t(blizzard_Continuous[
#   terra::cellFromXY(blizzard_Continuous[[1]], sf::st_coordinates(eightks_sf))
# ])
# # continuous_summits <- t(terra::extract(blizzard_Continuous), eightks_sf, method = "bilinear", fun = "mean"))[-1, ]

# colnames(binary_summits) <- colnames(continuous_summits) <- eightks_sf$PKNAME

# plot_df <- data.frame(
#   Summit = rep(colnames(binary_summits), each = nlyr(blizzard_Binary)),
#   binary = unlist(as.vector(binary_summits)),
#   continuous = unlist(as.vector(continuous_summits)),
#   Year = as.numeric(format(terra::time(blizzard_Binary), "%Y")),
#   Date = terra::time(blizzard_Binary)
# )

# save(plot_df, file = file.path(Dir.Data, "SnowStormPeaks.RData"))


# ## count storm events per year
# count_gg <- aggregate(binary ~ Year + Summit, data = plot_df, FUN = sum, na.rm = TRUE)
# ggplot(count_gg, aes(x = Year, y = binary, color = Summit)) +
#   geom_point() +
#   stat_smooth(method = "lm") +
#   theme_bw()

# cont_gg <- aggregate(continuous ~ Year + Summit, data = plot_df, FUN = mean, na.rm = TRUE)
# ggplot(cont_gg, aes(x = Year, y = continuous, color = Summit)) +
#   geom_point() +
#   stat_smooth(method = "lm") +
#   theme_bw()




# Peak_ls <- pblapply(unique(cont_gg$Summit), FUN = function(PeakIter) {
#   print(PeakIter)
#   Peak_df <- cont_gg[cont_gg$Summit == PeakIter, ]
#   C_BRM <- brms::brm(formula = continuous ~ Year, data = Peak_df)
#   unlist(lapply(brms::as_draws(C_BRM), "[[", "b_Year"))
# })
# names(Peak_ls) <- unique(cont_gg$Summit)

# modelplot_df <- data.frame(
#   Exceeded = unlist(lapply(Peak_ls, "[[", "Exceeded")),
#   Length = unlist(lapply(Peak_ls, "[[", "Length")),
#   Peak = rep(names(Peak_ls), each = length(Peak_ls[[1]][[1]]))
# )

# BRM_E_gg <- ggplot(modelplot_df, aes(y = Peak, x = Exceeded, fill = Peak)) +
#   stat_halfeye() +
#   scale_fill_viridis_d() +
#   theme_bw() +
#   geom_vline(xintercept = 0) +
#   theme(legend.position = "bottom")




# Peak_ls <- pblapply(unique(plot_df$Summit), FUN = function(PeakIter) {
#   print(PeakIter)
#   Peak_df <- na.omit(plot_df[plot_df$Summit == PeakIter, ])
#   C_BRM <- brms::brm(formula = continuous ~ Year, data = Peak_df)
#   # summary(C_BRM)
#   unlist(lapply(brms::as_draws(C_BRM), "[[", "b_Year"))
# })
# names(Peak_ls) <- unique(plot_df$Summit)

# ggplot(plot_df, aes(x = Year, y = continuous, color = Summit)) +
#   geom_point() +
#   stat_smooth(method = "lm") +
#   theme_bw()


# modelplot_df <- data.frame(
#   Exceeded = unlist(lapply(Peak_ls, "[[", "Exceeded")),
#   Length = unlist(lapply(Peak_ls, "[[", "Length")),
#   Peak = rep(names(Peak_ls), each = length(Peak_ls[[1]][[1]]))
# )

# BRM_E_gg <- ggplot(modelplot_df, aes(y = Peak, x = Exceeded, fill = Peak)) +
#   stat_halfeye() +
#   scale_fill_viridis_d() +
#   theme_bw() +
#   geom_vline(xintercept = 0) +
#   theme(legend.position = "bottom")



## count storm events per season

## strength of storms (only where storm conditions are met) per year and season
## length of storms (only where storm conditions are met) per year and season

## count storm conditions met per day to create data for model data frame

# FUSING WITH EXPEDITION DATA =============================================
message("#### Model Data Frame ############################################")
Expeditions_df <- Expeditions_df[Expeditions_df$PKNAME %in% summits_sf$PKNAME, ] # reduce to only those summits for which we have coordinates, losing 31 rows of data in Expeditions out of a total of 1917
Expeditions_df <- Expeditions_df[as.numeric(substr(Expeditions_df$TERMDATE, 1, 4)) > 1950, ] # climate data only available for 1951 onwards, losing a further 5 expeditions
Expeditions_df <- Expeditions_df[as.numeric(substr(Expeditions_df$BCDATE, 1, 4)) < 2022, ] # climate data only available until end of 2021, losing a further 26 expeditions

Datadf_ls <- pblapply(names(Data_ls), FUN = function(VarName) {
  # VarName = names(Data_ls)[1]
  VarIter <- Data_ls[[VarName]]
  VarTime <- t(terra::extract(
    VarIter, # variable during time between BC and TERM
    summits_sf[summits_sf$PKNAME %in% Expeditions_df$PKNAME, ],
    method = "bilinear", fun = "mean"
  ))
  VarTime <- data.frame(VarTime[-1, ])
  colnames(VarTime) <- summits_sf$PKNAME[summits_sf$PKNAME %in% Expeditions_df$PKNAME]
  rownames(VarTime) <- as.character(time(VarIter))
  VarTime
})
names(Datadf_ls) <- names(Data_ls)

quant_df <- do.call(rbind, lapply(names(Datadf_ls), FUN = function(VarIter) {
  quant_df <- do.call(rbind, lapply(colnames(Datadf_ls[[VarIter]]), FUN = function(PeakIter) {
    data_df <- data.frame(var = Datadf_ls[[VarIter]][, PeakIter], time = rownames(Datadf_ls[[VarIter]]))
    data.frame(
      Var = rep(VarIter, 4),
      Season = c("Pre", "Pre", "Post", "Post"),
      Bounds = c("L", "U", "L", "U"),
      Values = c(
        quantile(data_df[as.numeric(substr(data_df$time, 6, 7)) %in% 3:5, "var"], probs = c(0.05, 0.95), na.rm = TRUE),
        quantile(data_df[as.numeric(substr(data_df$time, 6, 7)) %in% 9:11, "var"], probs = c(0.05, 0.95), na.rm = TRUE)
      ),
      Peak = rep(PeakIter, 4)
    )
  }))
}))

Extract_df <- pbapply(Expeditions_df, MARGIN = 1, FUN = function(ExIter) {
  # ExIter <- Expeditions_df[1, ]
  # print(ExIter[1])
  Peak <- ExIter["PKNAME"]
  BC <- ExIter["BCDATE"]
  TERM <- ExIter["TERMDATE"]
  # print(BC)
  # print(TERM)

  Season <- c("Pre", "Post")[which.min(abs(c(
    4 - as.numeric(substr(as.character(BC), 6, 7)),
    8 - as.numeric(substr(as.character(BC), 6, 7))
  )))]
  # print(Season)

  # if (as.numeric(substr(as.character(BC), 6, 7)) %in% 3:5) {
  #   Season <- "Pre"
  # }
  # if (as.numeric(substr(as.character(BC), 6, 7)) %in% 9:11) {
  #   Season <- "Post"
  # }
  # if (!exists("Season")) {
  #   stop("not in seasons")
  # }

  df_ls <- lapply(names(Datadf_ls), FUN = function(VarName) {
    # VarName = names(Datadf_ls)[1]
    VarIter <- Datadf_ls[[VarName]]

    VarTime <- VarIter[which(rownames(VarIter) >= as.character(BC))[1]:tail(which(rownames(VarIter) <= as.character(TERM)), 1), which(colnames(VarIter) == Peak)]

    Bounds_df <- quant_df[quant_df$Peak == as.character(Peak) & quant_df$Season == Season & quant_df$Var == VarName, ]

    df1 <- data.frame(
      "mean" = mean(as.numeric(VarTime)),
      "sd" = sd(as.numeric(VarTime)),
      "ExtremeLow" = sum(VarTime < Bounds_df$Values[Bounds_df$Bounds == "L"]),
      "ExtremeHigh" = sum(VarTime > Bounds_df$Values[Bounds_df$Bounds == "U"])
    )
    colnames(df1) <- paste(VarName, colnames(df1), sep = "_")
    df1
  })
  do.call(cbind, df_ls)
})

Expeditions_Export <- cbind(Expeditions_df, do.call(rbind, Extract_df))
write.csv(Expeditions_Export, file = file.path(Dir.Exports, "ModelData.csv"))
head(Expeditions_Export)

stop("Ready for Chris")

# KRIGING =================================================================
message("#### Kriging Showcase ############################################")

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
    data = summits_df[summits_df$HEIGHTM >= 7000, ],
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
