
###-----------------------------------------------------------------
##   RUN RANGESHIFTER SIMULATION
###-----------------------------------------------------------------


# The function runRSmodel() runs the RangeShifter simulation for a given set of parameters.

# It modifies the simulation parameters to be calibrated, validates the new parameter settings,
# draws a new seed, creates a new initialisation scenario, runs the simulation, and returns
# the output as a raster stack.


runRSmodel <- function(params){
	
	p <- ref_par
	p[par_sel] <- params
	
	s@land@K_or_DensDep <- p["DensDep"]
	#### reparametrise survival and development
	d1 <- p["Deve1"] * p["Surv1"]
	d2 <- p["Deve2"] * p["Surv2"]
	s1 <- p["Surv1"] - d1
	s2 <- p["Surv2"] - d2
	s@demog@StageStruct@TransMatrix <- matrix(c( 0,  0, 0, p["Fecund"],
												 1, s1,  0,           0,
												 0, d1, s2,           0,
												 0,  0, d2, p["Surv3"]),
											  nrow = 4, byrow = T)
	####
	# s@demog@StageStruct@TransMatrix <- matrix(c( 0,          0,          0, p["Fecund"],
	# 											 1, p["Surv1"],          0,           0,
	# 											 0, p["Deve1"], p["Surv2"],           0,
	# 											 0,          0, p["Deve2"], p["Surv3"]),
	# 										  nrow = 4, byrow = T)
	s@dispersal@Emigration@EmigProb[2,2] <- p["EmigProb"]
	if(class(s@dispersal@Transfer)[1]=="CorrRW") {
		s@dispersal@Transfer@StepLength <- p["DispDist"]
		s@dispersal@Settlement@Settle <- matrix(c(p["SettProb"],
												  par_SettSlop,
												  p["SettBeta"]/(p["DensDep"]*400)),
												nrow = 1)
	}
	if(class(s@dispersal@Transfer)[1]=="DispersalKernel") {
		s@dispersal@Transfer@Distances <- matrix(p["DispDist"])
	}
	
	#s@init@InitIndsList <- replicate(s@simul@Replicates, expr = createInitIndsFile(pop_size = p["InitPop"], stage_struct = s@demog@StageStruct, breedingpairs = TRUE), simplify = FALSE )
	s@init@InitIndsList <- replicate(s@simul@Replicates, expr = createInitIndsFile_Pois(stage_struct = s@demog@StageStruct, breedingpairs = TRUE), simplify = FALSE )
	
	s@control@seed <- runif(1, 1000000, 9999999)
	
	validateRS <- FALSE
	try(validateRS <- validateRSparams(s), silent = TRUE)
	if(validateRS)
		capture.output(
			out <- RunRS(s, dirpath = RSdir) #,type = c("output", "message")
		)
	else out <- NA
	
	return(out)
}
