#' ####################################################################### #
#' PROJECT: [Himalayan Summit Success - Climate Data] 
#' CONTENTS: 
#'  - Supporting functionality for analysis script
#'  DEPENDENCIES:
#'  - None
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# ERA5-Land download on yearly basis ======================================
FUN.RawDown <- function(Variable, 
												Extent, 
												Years, 
												Dir, 
												parallel,
												API_User = API_User, # API User Number
												API_Key = API_Key # API User Key
												){
	if(parallel == 1){parallel <- NULL} # no parallelisation
	if(!is.null(parallel)){ # parallelisation
		message("Registering cluster for parallel processing")
		print("Registering cluster")
		parallel <- parallel::makeCluster(parallel)
		on.exit(stopCluster(parallel))
		print("R Objects loading to cluster")
		parallel::clusterExport(parallel, varlist = c(
			"package_vec", "install.load.package",
			"Variable", "Years", "Dir", "Extent",
			"API_User", "API_Key"
		), envir = environment())
		print("R Packages loading on cluster")
		clusterpacks <- clusterCall(parallel, function() sapply(package_vec, install.load.package))
		message("Starting parallel download")
	}
	
	Var_ls <- pblapply(Years, 
										 cl = parallel,
										 function(Year_Iter){
										 	
										 	print(Year_Iter)
										 	
										 	if(file.exists(file.path(Dir, paste0(paste(Variable, Year_Iter, sep="-"), ".nc")))){
										 		message("Already downloaded")
										 		Var_Year_ras <- stack(file.path(Dir, paste0(paste(Variable, Year_Iter, sep="-"), ".nc")))
										 	}else{
										 		
										 		StartYear <- paste0(Year_Iter, "-01-01")
										 		if(Year_Iter == 1950){StartYear <- paste0(Year_Iter, "-02-01")}
										 		
										 		Var_Year_ras <- download_ERA(
										 			Variable = Variable, # target variable
										 			DataSet = "era5-land", # data set
										 			DateStart = StartYear, # starting date of time-window
										 			DateStop = paste0(Year_Iter, "-12-31"), # final date of time-window
										 			TResolution = "hour",
										 			TStep = 1,
										 			Extent = Extent, # the spatial preference
										 			Dir = Dir, # where to store data
										 			FileName = paste(Variable, Year_Iter, sep="-"), # a name downloaded file
										 			API_User = API_User, # API User Number
										 			API_Key = API_Key, # API User Key
										 			SingularDL = TRUE
										 		)	
										 	}
										 })
	stack(Var_ls)
}

