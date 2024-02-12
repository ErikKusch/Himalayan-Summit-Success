#' ####################################################################### #
#' PROJECT: [Himalayan Summit Success - Climate Data] 
#' CONTENTS: 
#'  - Download, downscale, and prepare climate data for study of success rates in Hiamalayan mountaineering expeditions
#'  DEPENDENCIES:
#'  - X - Functions.R
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
  "pbapply", #for parallel lapply
  "readr", # for reading csv
  "ggplot2", # for plotting
  "tidyr", # for turning rasters into ggplot-dataframes
  "viridis", # colour palettes
  "cowplot", # gridding multiple plots
  "rnaturalearth", # for shapefiles
  "rnaturalearthdata", # for high-resolution shapefiles
  "mapview", # for generating mapview outputs
  "sf", # for more efficient handling of sp data
  "terra" # for more efficient handling of raster data
)
sapply(package_vec, install.load.package)

### NON-CRAN PACKAGES ----
if("KrigR" %in% rownames(installed.packages()) == FALSE){ # KrigR check
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
  devtools::install_github("ErikKusch/KrigR")
}
library(KrigR)

package_vec <- c("KrigR", package_vec)

source("X - Functions.R")

## API Credentials --------------------------------------------------------
try(source("X - PersonalSettings.R")) 
if(!exists("API_Key") | !exists("API_User")){ # CS API check: if CDS API credentials have not been specified elsewhere
  API_User <- readline(prompt = "Please enter your Climate Data Store API user number and hit ENTER.")
  API_Key <- readline(prompt = "Please enter your Climate Data Store API key number and hit ENTER.")
} # end of CDS API check
# NUMBER OF CORES
if(!exists("numberOfCores")){ # Core check: if number of cores for parallel processing has not been set yet
  numberOfCores <- as.numeric(readline(prompt = paste("How many cores do you want to allocate to these processes? Your machine has", parallel::detectCores())))
} # end of Core check

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

# DATA ====================================================================
## Loading from Disk ------------------------------------------------------
Nepal_shp <- ne_states(country = "Nepal") # shape data of Nepal from rnaturalearth package.
Nepal_shp <- as(Nepal_shp, "Spatial")
#' summits_df <- read_csv(file.path(Dir.Data, "mountains_df2.csv")) # load positions and names of summits
#' summits_df <- summits_df[!duplicated(summits_df$ID), ] # check for duplicates and eventually delete them
#' #' "Palung Ri" in row 302 has been originally a duplicated but changed to Palung Ri 1 and 2, respectively.
#' ## making sp object of summits
#' summits_sp <- summits_df
#' coordinates(summits_sp) <- ~Lon+Lat
#' proj4string(summits_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

## ERA5-Land --------------------------------------------------------------
Variables_vec <- c("10m_u_component_of_wind", "10m_v_component_of_wind", 
                   "2m_temperature", "skin_temperature", 
                   "snow_cover", "snow_density", "snow_depth", "snow_depth_water_equivalent",
                   "snow_evaporation", "snowfall", "snowmelt", "temperature_of_snow_layer")
Years_vec <- 1951:2021

#substrRight <- function(x, n){
#	substr(x, nchar(x)-n+1, nchar(x))
#}

Data_ls <- lapply(rev(Variables_vec), FUN = function(Var_Iter){
  message(paste("###", Var_Iter))
  
  if(file.exists(file.path(Dir.Data, paste0(Var_Iter, ".nc")))){
    print("Already prepared")
    Var_data <- rast(file.path(Dir.Data, paste0(Var_Iter, ".nc")))
  }else{
    Dir.Var <- file.path(Dir.Data, Var_Iter)
    if(!dir.exists(Dir.Var)){dir.create(Dir.Var)}
    
    Var_data <- FUN.RawDown(Variable = Var_Iter, 
                            Extent = extent(Nepal_shp) + c(-1, 1, -1, 1), 
                            Years = Years_vec, 
                            Dir = Dir.Var, 
                            parallel = 12,
                            API_User = API_User, # API User Number
                            API_Key = API_Key # API User Key
    )
  }
  Var_data
})

# COVARIATES ==============================================================
# Data <- mean(rast(file.path(Dir.Data, list.files(Dir.Data, pattern = ".nc")[1])))
Data <- mean(rast(list.files(file.path(Dir.Data, list.files(Dir.Data, pattern = "temperature")[1]), pattern = ".nc", full.names = TRUE)[1]))

COP_DEM <- rast(file.path(Dir.Covariates, "COP30.tif"))
COP_DEM <- terra::aggregate(COP_DEM, fact = 3) # have to upscale to circumvent "coordinate intervals are not constant" when creating spatialpixels data frame in krigR

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

Peaks_buffer <- st_union(Peaks_buffer)

# KRIGING =================================================================
COP_krig <- krigR(Data = raster(Data),
									Covariates_coarse = raster(COP_coarse),
									Covariates_fine = raster(COP_fine),
									Keep_Temporary = TRUE,
									Cores = parallel::detectCores(),
									KrigingEquation = "mean ~ COP30",
									nmax = 40,
									FileName = "COP_krig",
									Dir = Dir.Exports
)

COP_kriged <- rast(file.path(Dir.Exports, "COP_krig.nc"))

library(ggplot2)
library(ggrepel)
library(tidyterra)

ggplot() +
	geom_spatraster(data = Data, aes(fill = mean)) +
	geom_spatraster(data = COP_kriged, aes(fill = COP_krig)) +
	geom_sf(data = Peaks_buffer, color = "black", fill = "transparent") + 
	geom_sf(data = Peaks_sf, shape = 2) + 
	ggrepel::geom_text_repel(data = Peaks_df, 
													 aes(x = Lon, y = Lat, label = PKNAME),
													 max.overlaps = 20) + 
	scale_fill_viridis_c(na.value = "transparent", name = "[K]") + 
	labs(title = "Temperature of Snow Layer") + 
	theme_bw() + 
	theme(legend.position = "bottom", legend.key.width = unit(2.5, "cm"))
	






