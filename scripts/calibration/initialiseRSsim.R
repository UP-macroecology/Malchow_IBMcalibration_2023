
###-----------------------------------------------------------------
##   INITIALISE RANGESHIFTER SIMULATION
###-----------------------------------------------------------------


# The function initialiseRSsim() initialises the RangeShifter simulation.

# The model structure for the red kite is specified here. The function uses variables from
# the parent (global) environment to set the parameters. All calibration parameters are set
# to their reference values.

initialiseRSsim <- function() {
	
	simul <- Simulation(Years = Years, 
						Replicates = nr_replicates,
						Simulation = 1,
						OutStartPop = OutStartPop, 
						OutIntPop = OutIntPop,
						OutIntRange = 0,
						ReturnPopRaster = T,
						CreatePopFile = F )
	
	land <- ImportedLandscape(LandscapeFile = habitat_matrices,
							  OriginCoords = LLC_coords,
							  DynamicLandYears = c(0:(nrYrSDM-1)), # year 0=1999 (first year of MHB obs), 19=2018 (last year SDM)
							  HabPercent = TRUE,
							  Resolution = Resolution,
							  #SpDistFile = "RK_initdist.asc", SpDistResolution = Resolution,
							  K_or_DensDep = ref_par["DensDep"] )
	
	#### reparametrise survival and development
	d1 <- ref_par["Deve1"] * ref_par["Surv1"]
	d2 <- ref_par["Deve2"] * ref_par["Surv2"]
	s1 <- ref_par["Surv1"] - d1
	s2 <- ref_par["Surv2"] - d2
	tramat <- matrix(c( 0,  0, 0, ref_par["Fecund"],
						1, s1,  0,           0,
						0, d1, s2,           0,
						0,  0, d2, ref_par["Surv3"]),
					 nrow = 4, byrow = T)
	
	demog <- Demography(ReproductionType = 0,
						StageStruct = StageStructure(Stages = 4,
													 TransMatrix = tramat,
													 MaxAge = 12,
													 FecDensDep = TRUE,
													 FecStageWtsMatrix = matrix(c( 0,   0,   0,   0,
													 							   0,   0,   0,   0,
													 							   0,   0,   0,   0,
													 							   0,   0, 1.0, 1.0),
													 						    nrow = 4, byrow = T),
													 PostDestructn = TRUE
						)
	)
	
	#transfer <- DispersalKernel(Distances = ref_par["DispDist"])
	transfer <- CorrRW(StepLength = ref_par["DispDist"],
					   Rho = par_Rho)
	
	if(class(transfer)[1]=="CorrRW") settle <- Settlement(DensDep = TRUE,
														  Settle = matrix(c(ref_par["SettProb"],
														  				  par_SettSlop,
														  				  ref_par["SettBeta"]/(ref_par["DensDep"]*400)),
														  				nrow = 1),
														  MaxStepsYear = 5, MaxSteps = 10)

	if(class(transfer)[1]=="DispersalKernel") settle <- Settlement(Settle = 3)
	
	disp <- Dispersal(Emigration = Emigration(StageDep = TRUE,
											  EmigProb = matrix(c(0:3,0,ref_par["EmigProb"],0,0), nrow = 4 )),
					  Transfer = transfer,
					  Settlement = settle)
	
	init <- Initialise(InitType = 2, 
					   InitIndsFile = "NULL", # "initinds.txt", # List init
					   InitIndsList = replicate(nr_replicates,data.frame(0))
	)
	
	s <- RSsim(batchnum = batchnr, land = land, demog = demog, dispersal = disp, init = init, simul = simul)
	
	return(s)
}
