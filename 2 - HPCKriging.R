#' ####################################################################### #
#' PROJECT: [Himalayan Summit Success - Kriging on HPC] 
#' CONTENTS: 
#'  - Kriging of raw data products
#'  DEPENDENCIES:
#'  - NETCDF files in Dta directory
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PREAMBLE ================================================================
rm(list=ls()) # some may not like it, but it helps my workflow

## Packages ---------------------------------------------------------------
install.load.package <- function(x) {
	if (!require(x, character.only = TRUE))
		install.packages(x, repos='http://cran.us.r-project.org')
	require(x, character.only = TRUE)
}
### CRAN PACKAGES ----
package_vec <- c(
	"readr", # for reading csv
	"rnaturalearth", # for shapefiles
	"rnaturalearthdata", # for high-resolution shapefiles
	"sf", # for more efficient handling of sp data
	"terra", # for more efficient handling of raster data
	## for kriging
	"ncdf4", 
	"stringr", 
	"raster",
	"rgdal",
	"doParallel",
	"foreach",
	"doSNOW", 
	"automap",
	"lubridate",
	"sp",
	"sf",
	"fasterize",
	"stars",
	"httr",
	"terra"
)
sapply(package_vec, install.load.package)

# ### NON-CRAN PACKAGES ----
# if("KrigR" %in% rownames(installed.packages()) == FALSE){ # KrigR check
# 	Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
# 	devtools::install_github("ErikKusch/KrigR")
# }
# library(KrigR)
# 
# package_vec <- c("KrigR", package_vec)

## Directories ------------------------------------------------------------
### Define directories in relation to project directory
Dir.Base <- getwd() # identifying the current directory
Dir.Data <- file.path(Dir.Base, "Data") # folder path for data
Dir.Covariates <- file.path(Dir.Base, "Covariates") # folder path for covariates
Dir.Exports <- file.path(Dir.Base, "Exports") # folder path for exports
### create directories, if they don't exist yet
Dirs <- sapply(c(Dir.Data, Dir.Covariates, Dir.Exports), 
							 function(x) if(!dir.exists(x)) dir.create(x))
rm(Dirs)

## Kriging Function -------------------------------------------------------
#' Just in case KrigR cannot be loaded on HPC, fill in kriging functions here
krigR <- function(Data = NULL, Covariates_coarse = NULL, Covariates_fine = NULL, KrigingEquation = "ERA ~ DEM", Cores = detectCores(), Dir = getwd(), FileName, Keep_Temporary = TRUE, SingularTry = 10, Variable, PrecipFix = FALSE, Type = "reanalysis", DataSet = "era5-land", DateStart, DateStop, TResolution = "month", TStep = 1, FUN = 'mean', Extent, Buffer = 0.5, ID = "ID", API_Key, API_User, Target_res, Source = "USGS", nmax = Inf,  TryDown = 10, verbose = TRUE, TimeOut = 36000, SingularDL = FALSE, ...){
	## CALL LIST (for storing how the function as called in the output) ----
	if(is.null(Data)){
		Data_Retrieval <- list(Variable = Variable,
													 Type = Type,
													 PrecipFix = PrecipFix,
													 DataSet = DataSet,
													 DateStart = DateStart,
													 DateStop = DateStop,
													 TResolution = TResolution,
													 TStep = TStep,
													 Extent = Extent)
	}else{
		Data_Retrieval <- "None needed. Data was not queried via krigR function, but supplied by user."
	}
	## CLIMATE DATA (call to download_ERA function if no Data set is specified) ----
	if(is.null(Data)){ # data check: if no data has been specified
		Data <- download_ERA(Variable = Variable, PrecipFix = PrecipFix, Type = Type, DataSet = DataSet, DateStart = DateStart, DateStop = DateStop, TResolution = TResolution, TStep = TStep, FUN = FUN, Extent = Extent, API_User = API_User, API_Key = API_Key, Dir = Dir, TryDown = TryDown, verbose = verbose, ID = ID, Cores = Cores, TimeOut = TimeOut, SingularDL = SingularDL)
	} # end of data check
	
	## COVARIATE DATA (call to download_DEM function when no covariates are specified) ----
	if(is.null(Covariates_coarse) & is.null(Covariates_fine)){ # covariate check: if no covariates have been specified
		if(class(Extent) == "SpatialPolygonsDataFrame" | class(Extent) == "data.frame"){ # Extent check: if Extent input is a shapefile
			Shape <- Extent # save shapefile for use as Shape in masking covariate data
		}else{ # if Extent is not a shape, then extent specification is already baked into Data
			Shape <- NULL # set Shape to NULL so it is ignored in download_DEM function when masking is applied
		} # end of Extent check
		Covs_ls <- download_DEM(Train_ras = Data, Target_res = Target_res, Shape = Shape, Buffer = Buffer, ID = ID, Keep_Temporary = Keep_Temporary, Dir = Dir, Source = Source)
		Covariates_coarse <- Covs_ls[[1]] # extract coarse covariates from download_DEM output
		Covariates_fine <- Covs_ls[[2]] # extract fine covariates from download_DEM output
	} # end of covariate check
	
	## KRIGING FORMULA (assure that KrigingEquation is a formula object) ----
	KrigingEquation <- as.formula(KrigingEquation)
	
	## CALL LIST (for storing how the function as called in the output) ----
	Call_ls <- list(Data = SummarizeRaster(Data),
									Covariates_coarse = SummarizeRaster(Covariates_coarse),
									Covariates_fine = SummarizeRaster(Covariates_fine),
									KrigingEquation = KrigingEquation,
									Cores = Cores,
									FileName = FileName,
									Keep_Temporary = Keep_Temporary,
									nmax = nmax,
									Data_Retrieval = Data_Retrieval)
	
	## SANITY CHECKS (step into check_Krig function to catch most common error messages) ----
	Check_Product <- check_Krig(Data = Data, CovariatesCoarse = Covariates_coarse, CovariatesFine = Covariates_fine, KrigingEquation = KrigingEquation)
	KrigingEquation <- Check_Product[[1]] # extract KrigingEquation (this may have changed in check_Krig)
	DataSkips <- Check_Product[[2]] # extract which layers to skip due to missing data (this is unlikely to ever come into action)
	Terms <- unique(unlist(strsplit(labels(terms(KrigingEquation)), split = ":"))) # identify which layers of data are needed
	
	## DATA REFORMATTING (Kriging requires spatially referenced data frames, reformatting from rasters happens here) ---
	Origin <- raster::as.data.frame(Covariates_coarse, xy = TRUE) # extract covariate layers
	Origin <- Origin[, c(1:2, which(colnames(Origin) %in% Terms))] # retain only columns containing terms
	
	Target <- raster::as.data.frame(Covariates_fine, xy = TRUE) # extract covariate layers
	Target <- Target[, c(1:2, which(colnames(Target) %in% Terms))] # retain only columns containing terms
	Target <- na.omit(Target)
	suppressWarnings(gridded(Target) <- ~x+y) # establish a gridded data product ready for use in kriging
	Target@grid@cellsize[1] <- Target@grid@cellsize[2] # ensure that grid cells are square
	
	## SET-UP TEMPORARY DIRECTORY (this is where kriged products of each layer will be saved) ----
	Dir.Temp <- file.path(Dir, paste("Kriging", FileName, sep="_"))
	if(!dir.exists(Dir.Temp)){dir.create(Dir.Temp)}
	
	## KRIGING SPECIFICATION (this will be parsed and evaluated in parallel and non-parallel evaluations further down) ----
	looptext <- "
  OriginK <- cbind(Origin, raster::extract(x = Data[[Iter_Krige]], y = Origin[,1:2], df=TRUE)[, 2]) # combine data of current data layer with training covariate data
  OriginK <- na.omit(OriginK) # get rid of NA cells
  colnames(OriginK)[length(Terms)+3] <- c(terms(KrigingEquation)[[2]]) # assign column names
  suppressWarnings(gridded(OriginK) <-  ~x+y) # generate gridded product
  OriginK@grid@cellsize[1] <- OriginK@grid@cellsize[2] # ensure that grid cells are square

  Iter_Try = 0 # number of tries set to 0
  kriging_result <- NULL
  while(class(kriging_result)[1] != 'autoKrige' & Iter_Try < SingularTry){ # try kriging SingularTry times, this is because of a random process of variogram identification within the automap package that can fail on smaller datasets randomly when it isn't supposed to
    try(invisible(capture.output(kriging_result <- autoKrige(formula = KrigingEquation, input_data = OriginK, new_data = Target, nmax = nmax))), silent = TRUE)
    Iter_Try <- Iter_Try +1
  }
  if(class(kriging_result)[1] != 'autoKrige'){ # give error if kriging fails
    message(paste0('Kriging failed for layer ', Iter_Krige, '. Error message produced by autoKrige function: ', geterrmessage()))
  }

  ## retransform to raster
  try( # try fastest way - this fails with certain edge artefacts in meractor projection and is fixed by using rasterize
    Krig_ras <- raster(x = kriging_result$krige_output, layer = 1), # extract raster from kriging product
    silent = TRUE
  )
  try(
    Var_ras <- raster(x = kriging_result$krige_output, layer = 3), # extract raster from kriging product
    silent = TRUE
  )
  if(!exists('Krig_ras') & !exists('Var_ras')){
    Krig_ras <- rasterize(x = kriging_result$krige_output, y = Covariates_fine[[1]])[[2]] # extract raster from kriging product
    Var_ras <- rasterize(x = kriging_result$krige_output, y = Covariates_fine)[[4]] # extract raster from kriging product
  }
  crs(Krig_ras) <- crs(Data) # setting the crs according to the data
  crs(Var_ras) <- crs(Data) # setting the crs according to the data

  if(Cores == 1){
  Ras_Krig[[Iter_Krige]] <- Krig_ras
  Ras_Var[[Iter_Krige]] <- Var_ras
  } # stack kriged raster into raster list if non-parallel computing

  terra::writeCDF(x = as(brick(Krig_ras), 'SpatRaster'), filename = file.path(Dir.Temp, paste0(str_pad(Iter_Krige,4,'left','0'), '_data.nc')), overwrite = TRUE)
  # writeRaster(x = Krig_ras, filename = file.path(Dir.Temp, paste0(str_pad(Iter_Krige,4,'left','0'), '_data.nc')), overwrite = TRUE, format='CDF') # save kriged raster to temporary directory
  terra::writeCDF(x = as(brick(Var_ras), 'SpatRaster'), filename = file.path(Dir.Temp, paste0(str_pad(Iter_Krige,4,'left','0'), '_SE.nc')), overwrite = TRUE)
 # writeRaster(x = Var_ras, filename = file.path(Dir.Temp, paste0(str_pad(Iter_Krige,4,'left','0'), '_SE.nc')), overwrite = TRUE, format='CDF') # save kriged raster to temporary directory

  if(Cores == 1){ # core check: if processing non-parallel
    if(Count_Krige == 1){ # count check: if this was the first actual computation
      T_End <- Sys.time() # record time at which kriging was done for current layer
      Duration <- as.numeric(T_End)-as.numeric(T_Begin) # calculate how long it took to krig on layer
      message(paste('Kriging of remaining ', nlayers(Data)-Iter_Krige, ' data layers should finish around: ', as.POSIXlt(T_Begin + Duration*nlayers(Data), tz = Sys.timezone(location=TRUE)), sep='')) # console output with estimate of when the kriging should be done
      ProgBar <- txtProgressBar(min = 0, max = nlayers(Data), style = 3) # create progress bar when non-parallel processing
      Count_Krige <- Count_Krige + 1 # raise count by one so the stimator isn't called again
    } # end of count check
    setTxtProgressBar(ProgBar, Iter_Krige) # update progress bar with number of current layer
  } # end of core check
  "
	
	## KRIGING PREPARATION (establishing objects which the kriging refers to) ----
	Ras_Krig <- as.list(rep(NA, nlayers(Data))) # establish an empty list which will be filled with kriged layers
	Ras_Var <- as.list(rep(NA, nlayers(Data))) # establish an empty list which will be filled with kriged layers
	
	if(verbose){message("Commencing Kriging")}
	## DATA SKIPS (if certain layers in the data are empty and need to be skipped, this is handled here) ---
	if(!is.null(DataSkips)){ # Skip check: if layers need to be skipped
		for(Iter_Skip in DataSkips){ # Skip loop: loop over all layers that need to be skipped
			Ras_Krig[[Iter_Skip]] <- Data[[Iter_Skip]] # add raw data (which should be empty) to list
			terra::writeCDF(x = as(brick(Ras_Krig[[Iter_Skip]]), 'SpatRaster'), filename = file.path(Dir.Temp, str_pad(Iter_Skip,4,'left','0')), overwrite = TRUE)
			# writeRaster(x = Ras_Krig[[Iter_Skip]], filename = file.path(Dir.Temp, str_pad(Iter_Skip,4,'left','0')), overwrite = TRUE, format = 'CDF') # save raw layer to temporary directory, needed for loading back in when parallel processing
		} # end of Skip loop
		Layers_vec <- 1:nlayers(Data) # identify vector of all layers in data
		Compute_Layers <- Layers_vec[which(!Layers_vec %in% DataSkips)] # identify which layers can actually be computed on
	}else{ # if we don't need to skip any layers
		Compute_Layers <- 1:nlayers(Data) # set computing layers to all layers in data
	} # end of Skip check
	
	
	## ACTUAL KRIGING (carry out kriging according to user specifications either in parallel or on a single core) ----
	if(Cores > 1){ # Cores check: if parallel processing has been specified
		### PARALLEL KRIGING ---
		ForeachObjects <- c("Dir.Temp", "Cores", "Data", "KrigingEquation", "Origin", "Target", "Covariates_coarse", "Covariates_fine", "Terms", "SingularTry", "nmax") # objects which are needed for each kriging run and are thus handed to each cluster unit
		pb <- txtProgressBar(max = length(Compute_Layers), style = 3)
		progress <- function(n){setTxtProgressBar(pb, n)}
		opts <- list(progress = progress)
		cl <- makeCluster(Cores) # Assuming Cores node cluster
		registerDoSNOW(cl) # registering cores
		foreach(Iter_Krige = Compute_Layers, # kriging loop over all layers in Data, with condition (%:% when(...)) to only run if current layer is not present in Dir.Temp yet
						.packages = c("raster", "stringr", "automap", "ncdf4", "rgdal", "terra"), # import packages necessary to each itteration
						.export = ForeachObjects,
						.options.snow = opts) %:% when(!paste0(str_pad(Iter_Krige,4,"left","0"), '_data.nc') %in% list.files(Dir.Temp)) %dopar% { # parallel kriging loop
							Ras_Krig <- eval(parse(text=looptext)) # evaluate the kriging specification per cluster unit per layer
						} # end of parallel kriging loop
		close(pb)
		stopCluster(cl) # close down cluster
		Files_krig <- list.files(Dir.Temp)[grep(pattern = "_data.nc", x = list.files(Dir.Temp))]
		Files_var <- list.files(Dir.Temp)[grep(pattern = "_SE.nc", x = list.files(Dir.Temp))]
		for(Iter_Load in 1:length(Files_krig)){ # load loop: load data from temporary files in Dir.Temp
			Ras_Krig[[Iter_Load]] <- raster(file.path(Dir.Temp, Files_krig[Iter_Load])) # load current temporary file and write contents to list of rasters
			Ras_Var[[Iter_Load]] <- raster(file.path(Dir.Temp, Files_var[Iter_Load])) # load current temporary file and write contents to list of rasters
		} # end of load loop
	}else{ # if non-parallel processing has been specified
		### NON-PARALLEL KRIGING ---
		Count_Krige <- 1 # Establish count variable which is targeted in kriging specification text for producing an estimator
		for(Iter_Krige in Compute_Layers){ # non-parallel kriging loop over all layers in Data
			if(paste0(str_pad(Iter_Krige,4,'left','0'), '_data.nc') %in% list.files(Dir.Temp)){ # file check: if this file has already been produced
				Ras_Krig[[Iter_Krige]] <- raster(file.path(Dir.Temp, paste0(str_pad(Iter_Krige,4,'left','0'), '_data.nc'))) # load already produced kriged file and save it to list of rasters
				Ras_Var[[Iter_Krige]] <- raster(file.path(Dir.Temp, paste0(str_pad(Iter_Krige,4,'left','0'), '_SE.nc')))
				if(!exists("ProgBar")){ProgBar <- txtProgressBar(min = 0, max = nlayers(Data), style = 3)} # create progress bar when non-parallel processing}
				setTxtProgressBar(ProgBar, Iter_Krige) # update progress bar
				next() # jump to next layer
			} # end of file check
			T_Begin <- Sys.time() # record system time when layer kriging starts
			eval(parse(text=looptext)) # evaluate the kriging specification per layer
		} # end of non-parallel kriging loop
	} # end of Cores check
	
	## SAVING FINAL PRODUCT ----
	if(is.null(DataSkips)){ # Skip check: if no layers needed to be skipped
		# convert list of kriged layers in actual rasterbrick of kriged layers
		names(Ras_Krig) <- names(Data)
		if(class(Ras_Krig) != "RasterBrick"){Ras_Krig <- brick(Ras_Krig)}
		Krig_terra <- as(Ras_Krig, "SpatRaster")
		names(Krig_terra) <- names(Data)
		terra::writeCDF(x = Krig_terra, filename = file.path(Dir, paste0(FileName, ".nc")), overwrite = TRUE)
		# writeRaster(x = Ras_Krig, filename = file.path(Dir, FileName), overwrite = TRUE, format="CDF") # save final product as raster
		# convert list of kriged layers in actual rasterbrick of kriged layers
		names(Ras_Var) <- names(Data)
		if(class(Ras_Var) != "RasterBrick"){Ras_Var <- brick(Ras_Var)}
		Var_terra <- as(Ras_Var, "SpatRaster")
		names(Var_terra) <- names(Data)
		
		terra::writeCDF(x = Var_terra, filename = file.path(Dir, paste0("SE_", paste0(FileName, ".nc"))), overwrite = TRUE)
		# writeRaster(x = Ras_Var, filename = file.path(Dir, paste0("SE_",FileName)), overwrite = TRUE, format="CDF") # save final product as raster
	}else{ # if some layers needed to be skipped
		warning(paste0("Some of the layers in your raster could not be kriged. You will find all the individual layers (kriged and not kriged) in ", Dir, "."))
		Keep_Temporary <- TRUE # keep temporary files so kriged products are not deleted
	} # end of Skip check
	
	### REMOVE FILES FROM HARD DRIVE ---
	if(Keep_Temporary == FALSE){ # cleanup check
		unlink(Dir.Temp, recursive = TRUE)
	}  # end of cleanup check
	
	Krig_ls <- list(Ras_Krig, Ras_Var, Call_ls)
	names(Krig_ls) <- c("Kriging_Output", "Kriging_SE", "Call")
	return(Krig_ls) # return raster or list of layers
}

check_Krig <- function(Data, CovariatesCoarse, CovariatesFine, KrigingEquation){
	### RESOLUTIONS ----
	if(res(CovariatesFine)[1] < res(Data)[1]/10){
		warning("It is not recommended to use kriging for statistical downscaling of more than one order of magnitude. You are currently attempting this. Kriging will proceed.")
	}
	if(all.equal(res(CovariatesCoarse)[1], res(Data)[1]) != TRUE){
		stop(paste0("The resolution of your data (", res(Data)[1], ") does not match the resolution of your covariate data (", res(CovariatesCoarse)[1], ") used for training the kriging model. Kriging can't be performed!" ))
	}
	### EXTENTS ----
	# if(extent(Data) == extent(-180, 180, -90, 90)){
	#   stop("You are attempting to use kriging at a global extent. For reasons of computational expense and identity of relationships between covariates and variables not being homogenous across the globe, this is not recommended. Instead, try kriging of latitude bands if global kriging is really your goal.")
	# }
	if(!all.equal(extent(CovariatesCoarse), extent(Data))){
		stop("The extents of your data and training covariates don't match. Kriging can't be performed!")
	}
	
	### DATA AVAILABILITY ----
	DataSkips <- NULL # data layers without enough data to be skipped in kriging
	Data_vals <- base::colSums(matrix(!is.na(values(Data)), ncol = nlayers(Data))) # a value of 0 indicates a layer only made of NAs
	if(length(which(Data_vals < 2)) > 0){
		if(length(which(Data_vals < 2)) != nlayers(Data)){
			warning(paste0("Layer(s) ", paste(which(Data_vals == 0), collapse=", "), " of your data do(es) not contain enough data. Kriging will result in a raster identical do the input for this layer."))
			DataSkips <- which(Data_vals < 2)
		}else{
			stop("Your Data does not contain enough values. Kriging can't be performed!")
		}
	}
	CovCo_vals <- base::colSums(matrix(!is.na(values(CovariatesCoarse)), ncol = nlayers(CovariatesCoarse))) # a value of 0 indicates a layer only made of NAs
	if(length(which(CovCo_vals < 2)) > 0){
		if(length(which(CovCo_vals < 2)) != nlayers(CovariatesCoarse)){
			warning(paste0("Layer(s) ", paste(which(CovCo_vals < 2), collapse=", "), " of your covariates at training resolution do(es) not contain enough data. This/these layer(s) is/are dropped. The Kriging equation might get altered."))
			CovariatesCoarse <- CovariatesCoarse[[-which(CovCo_vals < 2)]]
		}else{
			stop("Your covariate data at training resolution does not contain enough values. Kriging can't be performed!")
		}
	}
	CovFin_vals <- base::colSums(matrix(!is.na(values(CovariatesFine)), ncol = nlayers(CovariatesFine))) # a value of 0 indicates a layer only made of NAs
	if(length(which(CovFin_vals < 2)) > 0){
		if(length(which(CovFin_vals < 2)) != nlayers(CovariatesFine)){
			warning(paste0("Layer(s) ", paste(which(CovFin_vals == 0), collapse=", "), " of your covariates at target resolution do(es) not contain enough data. This/these layer(s) is/are dropped."))
			CovariatesFine <- CovariatesFine[[-which(CovFin_vals < 2)]]
		}else{
			stop("Your covariate data at target resolution does not contain enough values. Kriging can't be performed!")
		}
	}
	### EQUATION ----
	Terms <- unlist(strsplit(labels(terms(KrigingEquation)), split = ":")) # identify parameters called to in formula
	Terms_Required <- unique(Terms) # isolate double-references (e.g. due to ":" indexing for interactions)
	Terms_Present <- Reduce(intersect, list(Terms_Required, names(CovariatesCoarse), names(CovariatesFine))) # identify the terms that are available and required
	if(sum(Terms_Required %in% Terms_Present) != length(Terms_Required)){
		if(length(Terms_Present) == 0){ # if none of the specified terms were found
			KrigingEquation <- paste0("Data ~ ", paste(names(CovariatesCoarse), collapse = "+"))
			warn <- paste("None of the terms specified in your KrigingEquation are present in the covariate data sets. The KrigingEquation has been altered to include all available terms in a linear model:", KrigingEquation)
		}else{ # at least some of the specified terms were found
			KrigingEquation <- paste0("Data ~ ", paste(Terms_Present, collapse = "+"))
			warn <- paste("Not all of the terms specified in your KrigingEquation are present in the covariate data sets. The KrigingEquation has been altered to include all available and specified terms in a linear model:", KrigingEquation)
		}
		Cotinue <- menu(c("Yes", "No"), title=paste(warn, "Do you wish to continue using the new formula?"))
		if(Cotinue == 2){ # break operation if user doesn't want this
			stop("Kriging terminated by user due to formula issues.")
		}
	}
	### NA DATA IN LAYERS ----
	# CovariatesFine <- CovariatesFine[[which(names(CovariatesFine) %in% Terms_Present)]] # only look at layers that the krigignequation targets
	# if(nlayers(CovariatesFine) > 1){
	#   MaskedPix <- length(which(values(sum(CovariatesFine, na.rm = TRUE)) != 0)) # number of non-masked pixels in which data is present in at least one layer
	#   MissingPix <- length(which(!is.na(values(sum(CovariatesFine, na.rm = FALSE))))) # number of pixels in which all layers have data
	#   if(MissingPix < MaskedPix){ # when there are any pixels for which data is absent for at least one layer
	#     stop("One or more more of your target covariate layers is missing data in locations where data is present for other layers. Please either fill these pixels with data or omit terms targeting these layers from your Kriging equation.")
	#   }
	# }
	return(list(as.formula(KrigingEquation), DataSkips))
}

SummarizeRaster <- function(Object_ras = NULL){
	Summary_ls <- list(Class = class(Object_ras),
										 Dimensions = list(nrow = nrow(Object_ras),
										 									ncol = ncol(Object_ras),
										 									ncell = ncell(Object_ras)),
										 Extent = Object_ras@extent,
										 CRS = crs(Object_ras),
										 layers = names(Object_ras))
	return(Summary_ls)
}


# DATA ====================================================================
## Loading from Disk ------------------------------------------------------
Nepal_shp <- ne_states(country = "Nepal") # shape data of Nepal from rnaturalearth package.
Nepal_shp <- as(Nepal_shp, "Spatial")

Data <- rast(list.files(Dir.Data, pattern = "2m_temperature.nc", full.names = TRUE))

COP_DEM <- rast(file.path(Dir.Covariates, "COP30.tif"))
COP_DEM <- terra::aggregate(COP_DEM, fact = 2) # have to upscale to circumvent "coordinate intervals are not constant" when creating spatialpixels data frame in krigR

Peaks_df <- read.csv(file.path(Dir.Data, "selected_peaks_coordinates_counts.csv"))
colnames(Peaks_df)[colnames(Peaks_df) %in% c("LON", "LAT")] <- c("Lon", "Lat")
Peaks_sf <- st_as_sf(Peaks_df, coords = c("Lon", "Lat"))
Peaks_buffer <- st_buffer(Peaks_sf[,"ID"], dist = 0.4, endCapStyle = "SQUARE")
st_crs(Peaks_buffer) <- terra::crs(COP_DEM)
st_crs(Peaks_sf) <- terra::crs(COP_DEM)

Data <- crop(Data, Peaks_buffer)
COP_DEM <- crop(COP_DEM, Peaks_buffer)
COP_coarse <- resample(COP_DEM, Data)
COP_fine <- mask(COP_DEM, Peaks_buffer)
rm(COP_DEM)

# Peaks_buffer <- st_union(Peaks_buffer)

# KRIGING =================================================================
## Time-Windows -----------------------------------------------------------
expeditions_df <- read_csv(file.path(Dir.Data, "individual_expeditions_higher_7000.csv"))
## throw out everything that lacks all date specifications or only have one specified date (due to date format ambiguity, this may me severly misleading)
expeditions_df <- expeditions_df[rowSums(is.na(
	expeditions_df[, c("BCDATE", "SMTDATE", "TERMDATE")]
	)) < 2, ]
## when summit days and or total days have been recorded, then fill date fields with them
expeditions_df$TOTDAYS <- ifelse(expeditions_df$TOTDAYS == 0, expeditions_df$SMTDAYS, expeditions_df$TOTDAYS)
DaysCheck <- which(rowSums(expeditions_df[, c("SMTDAYS", "TOTDAYS")]) > 0)
for(Days_i in DaysCheck){
	# Days_i <- 491
	Days_iter <- expeditions_df[Days_i, c("BCDATE", "SMTDATE", "TERMDATE", "SMTDAYS", "TOTDAYS")]
	DaysDiff <- unlist(c(Days_iter[4:5]))
	DaysDiff <- c(DaysDiff, diff(DaysDiff))
	## try to grab and fix date issues here
	if(sum(!is.na(Days_iter[1:3]))>1){
		ActualDiff <- Days_iter[which(!is.na(Days_iter[1:3]))][2]-Days_iter[which(!is.na(Days_iter[1:3]))][1]
		RecordedDiff <- DaysDiff[max(which(!is.na(Days_iter[1:3])))]
		if(ActualDiff != RecordedDiff){
			Days_iter$BCDATE <- as.Date(as.character(Days_iter$BCDATE), format = "%Y-%d-%m")
			Days_iter$SMTDATE <- as.Date(as.character(Days_iter$SMTDATE), format = "%Y-%d-%m")
			Days_iter$TERMDATE <- as.Date(as.character(Days_iter$TERMDATE), format = "%Y-%d-%m")
		}
	}
	if(!is.na(Days_iter$BCDATE)){
		if(is.na(Days_iter$SMTDATE)){Days_iter$SMTDATE <- Days_iter$BCDATE+Days_iter$SMTDAYS}
		if(is.na(Days_iter$TERMDATE)){Days_iter$TERMDATE <- Days_iter$BCDATE+Days_iter$TOTDAYS}
	}
	if(!is.na(Days_iter$TERMDATE)){
		if(is.na(Days_iter$BCDATE)){Days_iter$BCDATE <- Days_iter$TERMDATE-Days_iter$TOTDAYS}
		if(is.na(Days_iter$SMTDATE)){Days_iter$SMTDATE <- Days_iter$BCDATE+Days_iter$SMTDAYS}
	}
	if(!is.na(Days_iter$SMTDATE)){
		if(is.na(Days_iter$BCDATE)){Days_iter$BCDATE <- Days_iter$SMTDATE-Days_iter$SMTDAYS}
		if(is.na(Days_iter$TERMDATE)){Days_iter$TERMDATE <- Days_iter$BCDATE+Days_iter$TOTDAYS}
	}
	expeditions_df[Days_i, c("BCDATE", "SMTDATE", "TERMDATE", "SMTDAYS", "TOTDAYS")] <- Days_iter
}
## check data formatting errors (we have a mix of american YYYY-DD-MM and sensible YYYY-MM-DD formatting); if TERMDate or SMTDate predate BCDate or each other, we know there is an issue
DateCheck <- which(expeditions_df$SMTDATE<expeditions_df$BCDATE | 
	expeditions_df$TERMDATE<expeditions_df$BCDATE |
	expeditions_df$TERMDATE<expeditions_df$SMTDATE)
expeditions_df[DateCheck,]$BCDATE <- as.Date(as.character(expeditions_df[DateCheck,]$BCDATE), format = "%Y-%d-%m")
expeditions_df[DateCheck,]$SMTDATE <- as.Date(as.character(expeditions_df[DateCheck,]$SMTDATE), format = "%Y-%d-%m")
expeditions_df[DateCheck,]$TERMDATE <- as.Date(as.character(expeditions_df[DateCheck,]$TERMDATE), format = "%Y-%d-%m")
## calculate median times between Basecamp and summit and basecamp and termination by decade
expeditions_df$DIFF_BCtoSMT <- expeditions_df$SMTDATE-expeditions_df$BCDATE
expeditions_df$DIFF_BCtoTERM <- expeditions_df$TERMDATE-expeditions_df$BCDATE
Diff_SMT <- aggregate(DIFF_BCtoSMT ~ DECADE, FUN = median, data = expeditions_df)
Diff_TERM <- aggregate(DIFF_BCtoTERM ~ DECADE, FUN = median, data = expeditions_df)
## Impute BC and termination date for each expedition lacking this information
expeditions_df$IMPUTATION <- FALSE
ImputeCheck <- which(is.na(expeditions_df$BCDATE) | 
										 	(is.na(expeditions_df$TERMDATE) & is.na(expeditions_df$SMTDATE)))
for(Impute_i in ImputeCheck){
	expeditions_df$IMPUTATION[Impute_i] <- TRUE
	Impute_df <- expeditions_df[Impute_i, c("BCDATE", "SMTDATE", "TERMDATE")]
	Impute_decade <- expeditions_df$DECADE[Impute_i]
	### Basecamp date imputation
	if(is.na(Impute_df$BCDATE)){
		ImSMT <- Impute_df$SMTDATE - Diff_SMT$DIFF_BCtoSMT[Diff_SMT$DECADE == Impute_decade]
		ImTERM <- Impute_df$TERMDATE - Diff_TERM$DIFF_BCtoTERM[Diff_TERM$DECADE == Impute_decade]
		if(sum(is.na(c(ImSMT, ImTERM)) != 2)){ ## if only one could be computed
			BC_imp <- na.omit(c(ImSMT, ImTERM))
			ImSource <- c("SMT", "TERM")[which(!is.na(c(ImSMT, ImTERM)))]
		}else{ ## otherwise take earlier step
			BC_imp <- c(ImSMT, ImTERM)[which.min(c(ImSMT, ImTERM))]
			ImSource <- c("SMT", "TERM")[which.min(c(ImSMT, ImTERM))]
		}
		expeditions_df$IMPUTATION[Impute_i] <- paste0(expeditions_df$IMPUTATION[Impute_i], 
																									"_BC-", ImSource)
		expeditions_df$BCDATE[Impute_i] <- BC_imp
	}
	### Termination date imputation
	if(is.na(Impute_df$TERMDATE) & is.na(Impute_df$SMTDATE)){
		TERM_imp <- expeditions_df$BCDATE[Impute_i] + 
			Diff_TERM$DIFF_BCtoTERM[Diff_TERM$DECADE == Impute_decade]
		expeditions_df$IMPUTATION[Impute_i] <- paste0(expeditions_df$IMPUTATION[Impute_i], "_TERM")
		expeditions_df$TERMDATE[Impute_i] <- TERM_imp
	}
}
## find start and stop for kriging
Compare_df <- cbind(
	as.character(expeditions_df$SMTDATE),
	as.character(expeditions_df$TERMDATE)
)
Compare_df[is.na(Compare_df)] <- "0000-00-00"
KrigingTimeWindows <- data.frame(Start = expeditions_df$BCDATE-14,
																 Stop = ifelse(Compare_df[,1] > Compare_df[,2], Compare_df[,1], Compare_df[,2]))
KrigingTimeWindows <- na.omit(KrigingTimeWindows) # imputation failes for a few decades because there is only one expedition: 1900 and 1940

KrigCheck <- c()
for(Time_i in 1:nlyr(Data)){
	Time_i <- terra::time(Data)[Time_i]
	# print(Time_i)
	KrigCheck <- c(KrigCheck, 
								 sum(rowSums((cbind(Time_i >= KrigingTimeWindows$Start, Time_i <= KrigingTimeWindows$Stop)))>1)
								 )
}
KrigCheck <- as.logical(KrigCheck)
sum(KrigCheck > 1)/length(KrigCheck)


stop("DoubleCheck the above with Christian at a new date")

## Computations -----------------------------------------------------------
rm("Peaks_df", "Nepal_shp", "COP_DEM")
gc()
# closeAllConnections()

COP_krig <- krigR(Data = stack(Data[[as.logical(KrigCheck)]]),
													Covariates_coarse = stack(COP_coarse),
													Covariates_fine = stack(COP_fine),
													Keep_Temporary = TRUE,
													Cores = 1, # parallel::detectCores(),
													KrigingEquation = "mean ~ COP30",
													nmax = 40,
													FileName = paste0(terra::varnames(Data), "_KRIG"),
													Dir = Dir.Exports
)

