
###-----------------------------------------------------------------
##   RangeShiftR: calibration with simulated data
###-----------------------------------------------------------------

# This is the top-level script to run the sensitivity analysis and calibration 
# for the red kite model.
# 
# It starts by including the required libraries, functions and data, sets various 
# parameters that control the type of simulation, runs a sensitivity analysis, and
# the calibration with optimizer or MCMC.


### libraries ------------------------------------------------------------------

library(raster)
library(RangeShiftR)

### directories ----------------------------------------------------------------

RSdir <- "model/"

# load helper functions ------------------------------------

# these helper functions are mostly used to organise the code sections,
# as they heavily rely on parameters from the calling environment
source("scripts/calibration/readandpreparedata.R")
source("scripts/calibration/initialiseRSsim.R")
source("scripts/calibration/createInitIndsFile.R")
source("scripts/calibration/runRSmodel.R")
source("scripts/calibration/calculateFolds.R")


###-- 1.) Set simulation & calibration parameters -----------------------------------------------------------
{
#--- calibration parameters:
set_calib <- list(likdist = "GamPois", sampler = "DEzs", iter = 3e5, parallel = 3, InfVal = -1e6) # "Pois","GamPois","Gaus","TrGaus","LogNorm","Skell"
# set switches to choose which analyses to run
switch_analyses <- c("SA_oat" = FALSE, "SA_morris" = FALSE, "SA_bt" = FALSE,
                     "CA_mcmc" = TRUE, "CA_opt" = FALSE)

#---- spatial fold nr
spatial_fold <- 2

#---- data type:
simulatedData <- FALSE
aggr_blocks <- 14

#--- simulation parameters:
batchnr <- 25
nr_replicates <- 20
nrYrMHB <- 22 # nr of years in MHB data set (1999-2020)
nrYrSDM <- 20  # nr of years in SDM time series (1999-2018)
nrYrBurnin <- 2
Years <- nrYrMHB+nrYrBurnin  # years in output stack: 0...(Years-1), in Pop output: 0...Years
OutStartPop <- nrYrBurnin
OutIntPop <- 1
par_Rho <- 0.85
par_SettSlop <- -1.0
}

#--- model parameters:

{# set parameter ranges & reference params
  CA_params  <- data.frame(name = c( "DensDep","Fecund", "Surv1", "Surv2", "Surv3", "Deve1", "Deve2", "EmigProb","DispDist","SettBeta","SettProb","InitPop","sigma","GPsize"),
                           min  = c(   0.1e-2 ,   0.50 ,   0.01 ,   0.01 ,   0.01 ,   0.01 ,   0.01 ,      0.01 ,      400 ,     -15  ,     0.01 ,     200 ,  0.01,    1.00 ),
                           def  = c(   0.6e-2 ,   1.65 ,   0.42 ,   0.68 ,   0.80 ,   0.80 ,   0.55 ,      0.64 ,     2200 ,       4  ,     0.85 ,    1500 ,  0.20,   50.00 ),
                           max  = c(   2.0e-2 ,   5.00 ,   0.99 ,   0.99 ,   0.99 ,   0.99 ,   0.99 ,      0.99 ,     2400 ,      15  ,     0.99 ,    2500 ,  0.90,  500.00 ))
  
# select parameters for calibration
par_sel <- c(1:8,10,14)                           # ...InitPop is now replaced with Poisson model

# index of variance parameter(s)
sigma_parix <- which(par_sel==which(CA_params$name == "sigma"))
GPsize_parix <- which(par_sel==which(CA_params$name == "GPsize"))
#if(length(sigma_parix)!=1 || length(GPsize_parix)!=1) warning("variance parameter(s) must be selected for calibration") # sigma <- CA_params$def[CA_params$name == "sigma"]; GPsize <- CA_params$def[CA_params$name == "GPsize"]

# reference parameters (to be used when not chosen for calibration)
ref_par <- CA_params[["def"]]
names(ref_par) <- CA_params$name
if(length(sigma_parix)!=1)  sigma  <- ref_par["sigma"]
if(length(GPsize_parix)!=1) GPsize <- ref_par["GPsize"]
}


###-- 2.) Read and prepare data -------------------------------------------------------------

readandpreparedata()

###-- 3.) Create RS simulation -----------------------------------------------------

s <- initialiseRSsim()

# test runRSmodel():
#out <- runRSmodel(ref_par[par_sel])


###-- 4.) define likelihood / target function for calibration or optimisation

CAtarget <- function(params) 
{
	gof <- -Inf
	out = runRSmodel(params)
	if (class(out) == "RasterStack") {
		out <- stackApply(out, indices = as.integer(sapply(strsplit(names(out),	"year"), "[[", 2)) + 1 - nrYrBurnin, fun = mean)
		sim_data <- raster::extract(out[[ix_years_abd[ix_years_abd %in% seq(nlayers(out))]]], lookup[["cells"]], df = TRUE)
		sim_data$ID <- lookup$blockID
		names(sim_data)[names(sim_data) == "ID"] <- "blockID"
		notsim_years <- sum(!ix_years_abd %in% seq(nlayers(out)))
		if (notsim_years) 
			sim_data <- cbind(sim_data, data.frame(matrix(0, ncol = notsim_years, nrow = nrow(lookup))))
		sim_data <- sumBlocks(sim_data)
		sim_data[sim_data == 0] <- sigma
		GPsize <- params[GPsize_parix]
		gof <- sum(dnbinom(x = round(as.matrix(MHBabunds[, -1])), 
						   mu = as.matrix(sim_data[, -1]), size = GPsize, log = TRUE), 
				   na.rm = TRUE)
	}
	else return(set_calib$InfVal)
	return(gof)
}


# test CAtarget():
#CAtarget(CA_params$def[par_sel]*0.3)
#CAtarget(CA_params$def[par_sel])
#CAtarget(CA_params$max[par_sel]*0.8)

#replicate(20,CAtarget(CA_params$def[par_sel]))

calib_name_core <- paste0(ifelse(simulatedData,"sim","mhb"),"_",set_calib$likdist,ifelse(aggr_blocks,paste0("_agBlc",aggr_blocks),""),ifelse(spatial_fold,paste0("_spF",spatial_fold),""),"_rep",nr_replicates)


###-- 5.) Sensitivity analysis

##----- a) Local sensitivity analysis 

if(switch_analyses["SA_oat"]){
  source("scripts/calibration/sensitivityOAT.R")
  oneFactorSA <- sensitivityOAT(ll_points = 9, ll_replics = 9, savePlot = TRUE, minVal = set_calib$InfVal)
}


##----- b) Global sensitivity analysis with Morris screening

if(switch_analyses["SA_morris"]){
  require(sensitivity)
  require(BayesianTools)
  
  CAtarget_Pll <- generateParallelExecuter(CAtarget, parallel = ifelse(set_calib$parallel>1,set_calib$parallel,FALSE) ) # try to parallelize with BT
  
  morrisResult <- morris(model = CAtarget_Pll$parallelFun,
                         factors = CA_params$name[par_sel],
                         r = 15,
                         design = list(type = "oat", levels = 5, grid.jump = 3),
                         binf = CA_params$min[par_sel],
                         bsup = CA_params$max[par_sel],
                         scale = T) # scales calculated sensitivity in units of given interval, so that sensitivity is given w.r.t. its uncertainty
  #plot(morrisResult)
}

if(!any(switch_analyses[c("CA_mcmc","CA_opt")])) save.image(paste0("results/SA_out/SA_",calib_name_core,"_",batchnr,".RData"))

###-- 6.) Calibration

##----- a) Calibration with BayesianTools MCMC

if(switch_analyses["CA_mcmc"]){
  require(BayesianTools)

	# set the prior
	prior <- createTruncatedNormalPrior(mean = CA_params$def[par_sel], #(CA_params$min[par_sel]+(CA_params$max[par_sel]-CA_params$min[par_sel])/2),
                                      sd   = ifelse(rep(simulatedData,length(par_sel)),(CA_params$max[par_sel]-CA_params$min[par_sel])/c(2,4,4,4,4,3,3,3,2,2),c(0.0025,0.51,0.08,0.09,0.05,0.10,0.10,0.10,4,250)),
                                      lower = CA_params$min[par_sel],
                                      upper = CA_params$max[par_sel])

	setup <- createBayesianSetup(CAtarget,
                               prior = prior,
                               names = CA_params$name[par_sel],
                               parallel = ifelse(set_calib$parallel>1,set_calib$parallel,FALSE),
                               plotBest = CA_params$def[par_sel],
                               plotLower = CA_params$min[par_sel],
                               plotUpper = CA_params$max[par_sel])
  
  if(switch_analyses["SA_bt"]){
    pdf(file = paste0("results/SA_out/SA_bt_",calib_name_core,"_",batchnr,".pdf"), width = 12, height = 10)
      sen <- plotSensitivity(setup)
    dev.off()
  }
	
  # define MCMC settings
  if(TRUE){ 
  settings = list(iterations = set_calib$iter,
                  #f = ifelse(set_calib$sampler=="DEzs",1.7,2.38/sqrt(2*length(par_sel))),
                  #eps.add = 0.01,
                  startValue = ifelse(set_calib$parallel>1, set_calib$parallel, 3),
                  parallel = (set_calib$parallel>1)
                  )
  }else{ # run this instead to continue previous chain
    print("restart MCMC:\n")
    # restart DE sampler
    newZ_length <- 1e4
    restart_burnin <- set_calib$iter/2
    nrChains <- settings$startValue
    samples_x = getSample(MCMCout, start = restart_burnin)
    nrPar <- ncol(samples_x)
    #meansPost = apply(samples_x, 2, mean)
    #sdPost = apply(samples_x, 2, sd)
    rangePost = apply(samples_x, 2, range)
    newZ = matrix(runif(newZ_length*nrChains*nrPar, rangePost[1,], rangePost[2,]), ncol = nrPar, byrow = T)
    settings = list(iterations = set_calib$iter, Z = newZ, startValue = samples_x[(nrow(samples_x)-(nrChains-1)):nrow(samples_x), ])
  } 
  
  # run the chain
  MCMCout <- runMCMC(bayesianSetup = setup, sampler = set_calib$sampler, settings = settings)
  
  # produce output
  summary(MCMCout)
  getVolume(MCMCout, prior = T)
  
  save.image(paste0("results/CA_out/CA_",set_calib$sampler,"_",calib_name_core,"_it",set_calib$iter,"_",batchnr,".RData"))
  
  pdf(file = paste0("results/CA_out/CA_",set_calib$sampler,"_",calib_name_core,"_it",set_calib$iter,"_",batchnr,".pdf"), width = 12, height = 10)
    tracePlot(MCMCout, start=1)
    marginalPlot(MCMCout, start=1)
  dev.off()
}

