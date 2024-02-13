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
                        parallel = NULL,
                        API_User = API_User, # API User Number
                        API_Key = API_Key # API User Key
){
  in_parallel <- parallel
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
      "API_User", "API_Key", "in_parallel"
      #, "substrRight"
    ), envir = environment())
    print("R Packages loading on cluster")
    clusterpacks <- clusterCall(parallel, function() sapply(package_vec, install.load.package))
    message("Starting parallel download")
  }
  
  Var_ls <- pblapply(Years, 
                     cl = parallel,
                     function(Year_Iter){
                       
                       sink(file.path(Dir, paste0(Year_Iter, ".txt")))
                       print("Working on it")
                       sink()
                       
                       if(file.exists(file.path(Dir, paste0(paste(Variable, Year_Iter, sep="-"), ".nc")))){
                         message("Already downloaded")
                         Var_Year_ras <- stack(file.path(Dir, paste0(paste(Variable, Year_Iter, sep="-"), ".nc")))
                       }else{
                         Y_diff <- (Year_Iter - Years[1])
                         Sys.sleep(Year_Iter-Years[1] - in_parallel*floor(Y_diff/in_parallel))
                         ## error message handling to keep reiterating downloads until login has been validated
                         regex.escape <- function(string) {   gsub("([][{}()+*^${|\\\\?])", "\\\\\\1", string) }
                         errmsg <- "Default"
                         ## loop for as long as error message is no actual error (first iteration) or as long as it relates to validation errors on login
                         while(grepl(regex.escape(errmsg), pattern = "validate") | errmsg == "Default"){
                         	try(invisible(capture.output(Var_Year_ras <- download_ERA(
                         		Variable = Variable, # target variable
                         		DataSet = "era5-land", # data set
                         		DateStart = paste0(Year_Iter, "-01-01"), # starting date of time-window
                         		DateStop = paste0(Year_Iter, "-12-31"), # final date of time-window
                         		TResolution = "hour",
                         		TStep = 1,
                         		Extent = Extent, # the spatial preference
                         		Dir = Dir, # where to store data
                         		FileName = paste(Variable, Year_Iter, sep="-"), # a name downloaded file
                         		API_User = API_User, # API User Number
                         		API_Key = API_Key, # API User Key
                         		SingularDL = TRUE,
                         		TryDown = 1
                         	), silent = TRUE)))
                         	
                         	if(!exists("Var_Year_ras")){ # if download fails
                         		message(paste0('Download failing for ', 
                         									 Year_Iter, '. Error: ', 
                         									 geterrmessage()
                         		)
                         		)
                         		errmsg <- geterrmessage()
                         	}
                         } # while loop
                       } # else statement on file check
                       unlink(file.path(Dir, paste0(Year_Iter, ".txt")))
                       if(exists("Var_Year_ras")){ # if download doesn't fail
                       Var_Year_ras
                       }
                     }) # year pbapply
  
  names(Var_ls) <- Years
  
  Processed_ls <- pblapply(Var_ls, 
                           cl = ifelse(in_parallel > 4, 4, in_parallel),
                           function(Ras_Iter){
                             
                             Year_Iter <- substr(names(Ras_Iter)[1], start = 2, stop = 5)
                             message(Year_Iter)
                             
                             Processed_data <- stackApply(Ras_Iter, 
                                                          indices = rep(1:(nlayers(Ras_Iter)/24), each = 24),
                                                          fun = mean,
                                                          progress = "text",
                             )
                             
                             writeRaster(readAll(Processed_data), 
                                         filename = file.path(Dir,  paste("DAILY", Variable, paste0(Year_Iter, ".nc"), sep="-")),
                                         format = "CDF", overwrite = TRUE)
                             
                             Processed_data
                           })
  
  
  Processed_data <- rast(list.files(Dir, pattern = "DAILY", full.names = TRUE))
  terra::time(Processed_data) <- seq(as.Date(paste0(Years[1], "-01-01")),
                                     as.Date(paste0(Years[length(Years)], "-12-31")),
                                     'days')
  terra::varnames(Processed_data) <- Variable
  
  terra::writeCDF(Processed_data, 
                  filename = file.path(Dir.Data, paste0(Variable, ".nc")), 
                  overwrite = TRUE)
  
  unlink(list.files(Dir, pattern = "DAILY", full.names = TRUE))
  
  Processed_data
}

