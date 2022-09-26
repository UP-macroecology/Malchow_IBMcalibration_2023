
###-----------------------------------------------------------------
##   READ AND PREPARE DATA FOR SIMULATION AND CALIBRATION
###-----------------------------------------------------------------


# The function readandpreparedata() reads the data necessary for the simulation and 
# calibration and assigns them to the parent (mostly global) environment.

# First, the prepared data for the creation of the initialisation scenarios is read.
# Then, the habitat maps are read from ASCII files and turned into matrices for use in RangeShiftR.
# Next, the observations are read and processed according to optional spatial and temporal aggregation.
# Finally, the indices of the observed cells and years are recorded for the comparison of 
# observed and simulated data in the calibration step.


readandpreparedata <- function() {

	#--- 1.) read initialisation data (stochastic & deterministic part)
	
	mhb_init  <<- readRDS("data/init/initcell_mhb99_2km.rds")
	pois_pred <<- readRDS("data/init/init_Pois_pred_2km.rds")
	
	
	#--- 2.) read habitat map(s):
	

	habitat_names <- paste0("data/habitatmaps/RK_hsi_2km_",sprintf("%02d",1998+1:nrYrSDM),".asc")
	habitatmap <- raster(habitat_names[20])
	habitat_matrices <- lapply(habitat_names, 
							   FUN = function(filename){
							   	matrix(raster(filename), 
							   		   ncol = ncol(habitatmap), 
							   		   nrow = nrow(habitatmap), byrow = TRUE)})

	Resolution <<- res(habitatmap)[1]
	LLC_coords <<- c(xmin(habitatmap),ymin(habitatmap))
	assign("habitat_matrices", habitat_matrices, envir = .GlobalEnv)
	

	#--- 3.) read observed data set:

	load(file = "data/abund/mhb_counts.Rdata", envir = .GlobalEnv)
	
	#sample_cells_habitat <- rowSums(MHBabunds[-1],na.rm = TRUE)>0 | values(habitatmap)[sample_cells_abd]>0
	sample_years_abd <<- (0:(nrYrMHB-1))   # sampled years to be compared with simulation
	simul_years <<- seq(from = OutStartPop-nrYrBurnin, by = OutIntPop, to = Years-nrYrBurnin-1)   # simulated years
	ix_years_abd <<- which(simul_years %in% sample_years_abd)   # indices of sampled years within vector of simulated years
}
