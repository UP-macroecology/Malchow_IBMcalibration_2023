
###-----------------------------------------------------------------
##   GENERATE OUTPUTS AND PLOTS FROM CALIBRATION RESULTS
###-----------------------------------------------------------------

# Main script for analysing the posterior sample and for plotting results; loads the functions stored in the other files in the same folder
# Note: The sensitivity analysis and its plotting are not included here, but can be found in its own script plot_sensitivity.R.

# Contents:
# 1.) Create and store combined MCMCs
# 2.) Plot posterior marginal distributions & its pairwise correlations
# 3.) MCMC diagnostic plots
# 4.) Spatial-block cross-validation  (Simulations and plotting)
# 5.) Projections  (Simulations and plotting)


# include packages

library(raster)
library(xtable)
library(ggplot2)
library(RangeShiftR)
library(coda)
library(bayesplot)
library(BayesianTools)
library(pROC)
library(Hmisc)
library(lattice)
library(latticeExtra)
library(grid)
library(RColorBrewer)


# read functions from other files in this folder

source("scripts/analysis/docu_colors.R")
source("scripts/analysis/combine_chains.R")
source("scripts/analysis/plot_posterior.R")

# choose type of graphics device: latex (TRUE) or pdf (FALSE)?
latex_device <- FALSE

# parameter names used in plot
if(latex_device){ par_names <- c("Density-dependence $b^{-1}$","Fecundity $\\phi_0$","Survival prob. $\\sigma_1$","Survival prob. $\\sigma_2$","Survival prob. $\\sigma_3$","Developm. prob. $\\gamma_{1\\rightarrow 2}$","Developm. prob. $\\gamma_{2\\rightarrow 3}$","Emigration prob. $e_1$","Settlem. par. $\\beta_s$","Dispersion $\\nu$")
}else{ par_names <- NULL } # if not latex, names are taken from MCMC object



###--- 1.) Create and store combined MCMCs

# Load and combine single MCMCs, generate some summary statistics, and store everything to new objects

## a) get info about the MCMC to identify the correct file names to read:

samplr <- "GamPois" # likelihood type
aggr_blocks <- 14   # amount of spatial aggregation
reps <- 20          # number of replicates in RS
iter <- 3e+05       # number of iterations in chain

# is the MCMC run on simulated or MHB data?
simulatedData <- FALSE # MHB data

# give batch IDs and fold numbers of single MCMCs 
if(simulatedData){
	batches <- matrix(101:112,nrow=3)
	fold <- NULL
}else{
	batches <- matrix(c(401:403,411:413,421:423,431:433,441:443,451:453),nrow=3) #+ 100
	fold <- 0:5
}

# parameters for chain treatment
it_burnin <- 5e4 # number of iterations treated as burn-in 
thinningrate <- 100 # rate of thinning

# are there two chains to be concatenated? (in case one chain was stopped and restarted by the second)
# (in this case, the respective second chains must have a batch number that is larger by 100 than their first chains )
concat <- TRUE # if TRUE, combineMCMCs() returns a Coda object; if FALSE, a BayesianTools::mcmcSampler object


## b) Create MCMC files: 
#     loop over columns in batches; for each batchID read single MCMC from file and combine MCMCs of all batchIDs in a column; then save to file

if(FALSE){ # deactivate the loop once combined MCMCs are generated
	for(batch in 1:ncol(batches) ) {
		
		simnames <- batches[,batch]
		foldID <- fold[batch]
		
		# read the calibration outputs for all runs and create combined MCMC sampler list
		sampler_list <- combineMCMCs(simulatedData,samplr,aggr_blocks,reps,iter,simnames,foldID,concatenateChains = concat)
		
		# retrieve info about MCMC
		if(concat) sampler_info <- posterior_df <- NULL else list2env(evalMCMC(sampler_list),globalenv())
		
		# write summary & diagnostics output to file
		if(concat){
			capture.output(
				effectiveSize(sampler_list),
				gelman.diag(sampler_list),
				(correlations <- coda::crosscorr(sampler_list)),
				autocorr.diag(sampler_list),
				file = paste0("results/CA_out/CA_DEzs_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_rep",reps,"_it",iter,"_",simnames[1],"-",simnames[3],"_concat_summary.txt"))
		}else{
			capture.output(
					correlations <- summary(sampler_list),
					gelmanDiagnostics(sampler_list, plot = F),
					file = paste0("results/CA_out/CA_DEzs_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_rep",reps,"_it",iter,"_",simnames[1],"-",simnames[3],"_summary.txt"))
		}

		# generate summarising plots
		if(TRUE){
			# open graphics device
			if(latex_device){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",ifelse(is.null(fold),paste0("bch_",batch),paste0("spF_",foldID)),ifelse(concat,"_concat",""),".tex"), width = 5.76, height = 3.84, standAlone = FALSE)
			}else{ pdf(file = paste0("results/CA_out/CA_DEzs_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_rep",reps,"_it",iter,"_",simnames[1],"-",simnames[3],ifelse(concat,"_concat",""),".pdf"),width = 12, height = 8) }
			if(concat){
				# coda plot
				plot(sampler_list, trace = !latex_device, density = TRUE, smooth = FALSE)
			}else{
				# plot marginals 
				myMarginalPlot(sampler_list, start = it_burnin, par_names = par_names)
				# plot chain traces
				if(!latex_device){ tracePlot(sampler_list) }
			}
			dev.off()
		}

		# save combined chain and info in R object, then remove
		save(sampler_list, sampler_info, posterior_df, correlations, 
			 file = paste0("results/CA_out/CA_DEzs_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_rep",reps,"_it",iter,"_",simnames[1],"-",simnames[3],ifelse(concat,"_concat",""),"_chain.RData"))
		
		rm(sampler_list, sampler_info, posterior_df, simnames, foldID, correlations)
	}
}
	 # batchcol = 1
	 # load(paste0("results/CA_out/CA_DEzs_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_rep",reps,"_it",iter,"_",batches[1,batchcol],"-",batches[3,batchcol],ifelse(concat,"_concat",""),"_chain.RData"))
	 # getVolume(sampler_list, prior = F)
	 # $priorVol  3.842592e-14
	 # $postVol   1.525253e-17




###--- 2.) Plot posterior marginal distributions & its pairwise correlations


if(FALSE){ # comment out the plotting

	# get parameter ranges
	CA_params  <- data.frame(name = c( "DensDep","Fecund", "Surv1", "Surv2", "Surv3", "Deve1", "Deve2", "EmigProb","DispDist","SettBeta","SettProb","InitPop","sigma","GPsize"),
							 min  = c(   0.1e-2 ,   0.50 ,   0.01 ,   0.01 ,   0.01 ,   0.01 ,   0.01 ,      0.01 ,      400 ,     -15  ,     0.01 ,     200 ,  0.01,    1.0  ),
							 def  = c(   0.6e-2 ,   1.65 ,   0.42 ,   0.68 ,   0.80 ,   0.80 ,   0.55 ,      0.64 ,     2200 ,       4  ,     0.85 ,    1500 ,  0.20,   50.0  ),
							 max  = c(   2.0e-2 ,   5.00 ,   0.99 ,   0.99 ,   0.99 ,   0.99 ,   0.99 ,      0.99 ,     2400 ,      15  ,     0.99 ,    2500 ,  0.90,  500.0 ))
	par_sel <- c(1:8,10,14) #

#- 2.1) Plot posterior marginal box plots

	##- 2.1.a) for all folds
		combinded_smpl <- sapply(0:max(fold), function(foldID, nPriorDraws = 20000, thin = "auto"){
			simnames <- batches[,fold==foldID]
			filename = paste0("results/CA_out/CA_DEzs_mhb_",samplr,"_agBlc",aggr_blocks,"_rep",reps,"_it",iter,"_",simnames[1],"-",simnames[3],ifelse(concat&foldID!=0,"_concat",""),"_chain.RData")
			attach(filename, name = "currChain")
			message(paste("Loading",filename,"..."))
			samplist <- get("sampler_list",pos="currChain")
			if(foldID==0) {samplMat = samplist[[1]]$setup$prior$sampler(nPriorDraws)  # draw prior from bayesianSetup, only for one fold as they all have the same prior
			}else{samplMat <- getSample(samplist, parametersOnly = TRUE, start = it_burnin, thin = thin)}
			detach("currChain", character.only = TRUE)
			samplMat
		})
	
		# file name extension
		filename_ext <- "marg_folds_boxplot"
		# open graphics device
		if(latex_device){ tikz(file = paste0(path_doc,"CA_mhb_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_",filename_ext,".tex"), width = 5.7, height = 6.44, standAlone = FALSE) 
		}else{ pdf(file = paste0(path_doc,"pdf/CA_mhb_",samplr,"_agBlc",aggr_blocks,min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_",filename_ext,".pdf"),width = 12, height = 8) }
		# plot
		multiMarginalBoxplot(combinded_smpl, par_names = par_names, xrange = t(as.matrix(CA_params[par_sel,c(2,4)])))
		# close graphics device
		dev.off()
		
	##- 2.1.b) only for the total fold (fold 0)

	# load the all-folds-combined MCMC (MHB data)
	foldID <- 0
	if(simulatedData){ simnames <- batches[,1]
	}else{ simnames <- batches[,fold==foldID] }
	
	filename = paste0("results/CA_out/CA_DEzs_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_rep",reps,"_it",iter,"_",simnames[1],"-",simnames[3],ifelse(concat,"_concat",""),"_chain.RData")
	attach(filename, name = "currChain")
	sampler_list <- get("sampler_list",pos="currChain")
	detach("currChain", character.only = TRUE)

	
		# get maximum-a-posteriori (MAP) sections of the joint posterior
		  #map_sect <- get_MAPsection(sampler_list, CIpercent = 99)
		
		filename_ext <- "marg_boxplot"
		if(latex_device){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",ifelse(is.null(fold),paste0("bch_",batch),paste0("spF_",foldID)),"_",filename_ext,".tex"), width = 5.7, height = 5.84, standAlone = FALSE)
		}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_rep",reps,"_it",iter,"_",simnames[1],"-",simnames[3],"_",filename_ext,".pdf"),width = 12, height = 8) }
		myMarginalPlot(sampler_list, coda_chain = coda_list, start = it_burnin, thin = "auto", par_names = par_names, xrange = t(as.matrix(CA_params[par_sel,c(2,4)])), singlePanel = TRUE, boxplots = TRUE )
		dev.off()
	
	
#- 2.2) Plot all parameters marginal distributions
	
	# all parameters in one panel
	filename_ext <- "marg_onepanel"
	if(latex_device){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",ifelse(is.null(fold),paste0("bch_",batch),paste0("spF_",foldID)),"_",filename_ext,".tex"), width = 5.7, height = 5.84, standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_rep",reps,"_it",iter,"_",simnames[1],"-",simnames[3],"_",filename_ext,".pdf"),width = 12, height = 8) }
	myMarginalPlot(sampler_list, coda_chain = coda_list, start = it_burnin, thin = thinningrate, par_names = par_names, xrange = t(as.matrix(CA_params[par_sel,c(2,4)])), singlePanel = TRUE, plot_ref = FALSE, plot_MAP = FALSE, MAP_section = NULL ) # MAP_section = map_sect
	dev.off()

	# one panel for each parameter
	filename_ext <- "marg_singlepanel"
	if(latex_device){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",ifelse(is.null(fold),paste0("bch_",batch),paste0("spF_",foldID)),"_",filename_ext,".tex"), width = 5.76, height = 3.84, standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_rep",reps,"_it",iter,"_",simnames[1],"-",simnames[3],"_",filename_ext,".pdf"),width = 12, height = 8) }
	myMarginalPlot(sampler_list, coda_chain = coda_list, start = it_burnin, thin = thinningrate, par_names = par_names, xrange = t(as.matrix(CA_params[par_sel,c(2,4)])), singlePanel = FALSE, plot_ref = FALSE, plot_MAP = FALSE, MAP_section = NULL)
	dev.off()
	
	
#- 2.3) Plot pairwise correlations
	
	# with BayesianTools
	filename_ext <- "post_correlations"
	if(FALSE){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",ifelse(is.null(fold),paste0("bch_",batch),paste0("spF_",foldID)),"_",filename_ext,".tex"), width = 5.76, height = 4.84, standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_rep",reps,"_it",iter,"_",simnames[1],"-",simnames[3],"_",filename_ext,".pdf"),width = 12, height = 12) }
	BayesianTools::correlationPlot(coda_list, start = it_burnin, thin = "auto")
	dev.off()
	
	
	# with Coda - simpler layout
	#post_sample <- getSample(sampler_list, coda = TRUE, parametersOnly = FALSE, start = it_burnin, thin = "auto") # turn sampler list into a Coda-chain-class object
	
	filename_ext <- "post_correlations_simple"
	if(FALSE){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",ifelse(is.null(fold),paste0("bch_",batch),paste0("spF_",foldID)),"_",filename_ext,".tex"), width = 5.76, height = 4.84, standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_rep",reps,"_it",iter,"_",simnames[1],"-",simnames[3],"_",filename_ext,".pdf"),width = 5, height = 5) }
	coda::crosscorr.plot(sampler_list)
	dev.off()
	
	
#- 2.3) Compare width of posterior vs prior
	
	plot_HPDIranges <- function(samples){
		# calculate HDPI ranges of posterior for different probability mass fractions 
		HPDI_probs <- c(0.01,seq(0.1,0.9,0.1),0.99)
		HPDI_ranges_post <- sapply(HPDI_probs, function(x, chain = samples){
			HPDIperChain <- data.frame(
				lapply(
					coda::HPDinterval(chain, prob=x ), 
					function(y){data.frame(range=y[,2]-y[,1])}))
			rowMeans(HPDIperChain)
		})
		# get rid of unneeded rows 
		HPDI_ranges_post <- HPDI_ranges_post[-which(rownames(HPDI_ranges_post)%in%c("Lposterior","Llikelihood","Lprior")),]
		 
		# get HPDIs for the Prior, too
		prior <- createTruncatedNormalPrior(mean = CA_params$def[par_sel],
											sd   = ifelse(rep(simulatedData,length(par_sel)),(CA_params$max[par_sel]-CA_params$min[par_sel])/c(2,4,4,4,4,3,3,3,2,2),c(0.0025,0.51,0.08,0.09,0.05,0.10,0.10,0.10,4,250)),
											lower = CA_params$min[par_sel],
											upper = CA_params$max[par_sel])
		prio_sample <- as.mcmc(prior$sampler(coda::niter(samples)))
		HPDI_ranges_prio <- as.data.frame(
			sapply(HPDI_probs, function(x, chain = prio_sample){
				tmp <- coda::HPDinterval( chain , prob = x )
				return( data.frame(range=tmp[,2]-tmp[,1]) )
			}))
		
		# normalise parameter ranges
		HPDI_ranges_post <- HPDI_ranges_post / (CA_params[par_sel,4]-CA_params[par_sel,2])
		HPDI_ranges_prio <- HPDI_ranges_prio / (CA_params[par_sel,4]-CA_params[par_sel,2])
		
		# plot
		if(FALSE){
			# all in one plot
			plot(NULL, xlim = c(0,1), ylim = c(0,1), main = "HPDI ranges", xlab = "Probability mass fraction", ylab = "Parameter range fraction")
			for(i in 1:length(par_sel)) points(HPDI_probs,HPDI_ranges_post[i,], type = "b")
			for(i in 1:length(par_sel)) points(HPDI_probs,HPDI_ranges_prio[i,], type = "b", col = "red", lty = 2)
		}else{
			# multi-plot for each parameter
			if(is.null(par_names)) par_mains = rownames(HPDI_ranges_post) else par_mains = par_names
			par(mfrow = c(4,3), oma = c(2,3,2.5,1), mar = c(2, 2, 2.3, 1))
			for(i in 1:length(par_sel)){
				plot(NULL, xlim = c(0,1), xaxt = 'n', ylim = c(0,1), yaxt = 'n', main = par_mains[i])
				axis(1, at = seq(0,1,.5))
				axis(2, at = seq(0,1,.5))
				points(HPDI_probs,HPDI_ranges_post[i,], type = "b", col = "black", lty = 1)
				points(HPDI_probs,HPDI_ranges_prio[i,], type = "b", col = "red", lty = 2)
			}
			par(fig = c(0, 1, 0, 1), oma = c(2,2,2,2), mar = c(0, 0, 0, 0), new = TRUE)
			#title(main = "HPDI ranges", xlab = "Probability mass fraction", ylab = "Parameter range fraction")
			mtext("HPDI ranges", 3, cex = 1.3)
			mtext("Probability mass fraction", 1)
			mtext("Parameter range fraction", 2)
			legend('bottomright', c('posterior', 'prior'), lty = c(1,2), col = c(1,2), bg = "white", cex = 1.3, pt.cex = 2.5, inset = c(0.0,-0.14))
		}
	}
	
	# make plots
	filename_ext <- "HPDI_ranges"
	if(latex_device){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",ifelse(is.null(fold),paste0("bch_",batch),paste0("spF_",foldID)),"_",filename_ext,".tex"), width = 5.76, height = 5.5, standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_rep",reps,"_it",iter,"_",simnames[1],"-",simnames[3],"_",filename_ext,".pdf"),width = 12, height = 12) }
	plot_HPDIranges(samples = post_sample)
	dev.off()
	
	
	
###--- 3.) MCMC diagnostic plots 

#- 3.1) Plot Gelman diagnostic for each parameter
	filename_ext <- "diag_gelman"
	if(latex_device){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",ifelse(is.null(fold),paste0("bch_",batch),paste0("spF_",foldID)),"_",filename_ext,".tex"), width = 5.76, height = 4.84, standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_rep",reps,"_it",iter,"_",simnames[1],"-",simnames[3],"_",filename_ext,".pdf"),width = 5.76, height = 6) }
	par(mfrow=c(4,3),mar=c(2,2,2.5,0.5))
	coda::gelman.plot(coda_list, autoburnin = FALSE, auto.layout = FALSE, confidence = 0.99)
	dev.off()


#- 3.2) Plot Geweke diagnostic for each parameter
	filename_ext <- "diag_geweke"
	if(FALSE){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",ifelse(is.null(fold),paste0("bch_",batch),paste0("spF_",foldID)),"_",filename_ext,".tex"), width = 5.76, height = 4.84, standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_rep",reps,"_it",iter,"_",simnames[1],"-",simnames[3],"_",filename_ext,".pdf"),width = 12, height = 12) }
	coda::geweke.plot(coda_list)
	dev.off()
	
	
#- 3.3) evolution of the sample quantiles over number of iterations
	filename_ext <- "quantile_evolution"
	if(FALSE){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",ifelse(is.null(fold),paste0("bch_",batch),paste0("spF_",foldID)),"_",filename_ext,".tex"), width = 5.76, height = 4.84, standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_rep",reps,"_it",iter,"_",simnames[1],"-",simnames[3],"_",filename_ext,".pdf"),width = 12, height = 12) }
	coda::cumuplot(coda_list)
	dev.off()
	
	
#- 3.4) posterior marginals, coda style
	filename_ext <- "marg_coda"
	if(FALSE){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",ifelse(is.null(fold),paste0("bch_",batch),paste0("spF_",foldID)),"_",filename_ext,".tex"), width = 5.76, height = 4.84, standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_rep",reps,"_it",iter,"_",simnames[1],"-",simnames[3],"_",filename_ext,".pdf"),width = 12, height = 12) }
	par(mfrow = c(4,4))
	coda::densplot(coda_list)
	dev.off()
	
	
#- 3.5) Trace plots
	
	# get window from chains for plotting
	post_sample <- window(coda_list, start = it_burnin, thin = thinningrate)
	
	# with Coda
	filename_ext <- "chain_traces_coda"
	if(FALSE){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",ifelse(is.null(fold),paste0("bch_",batch),paste0("spF_",foldID)),"_",filename_ext,".tex"), width = 5.76, height = 4.84, standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_rep",reps,"_it",iter,"_",simnames[1],"-",simnames[3],"_",filename_ext,".pdf"),width = 5.76, height = 8) }
	par(mfrow = c(5,2),mar=c(2,2,2.5,0.8))
	coda::traceplot(post_sample)
	dev.off()

	# with BayesianTools
	filename_ext <- "chain_traces_bt"
	if(FALSE){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",ifelse(is.null(fold),paste0("bch_",batch),paste0("spF_",foldID)),"_",filename_ext,".tex"), width = 5.76, height = 4.84, standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_rep",reps,"_it",iter,"_",simnames[1],"-",simnames[3],"_",filename_ext,".pdf"),width = 12, height = 12) }
	BayesianTools::tracePlot(sampler_list, start = it_burnin, thin = thinningrate)
	dev.off()
	
	# with Bayesplot
	filename_ext <- "chain_traces_bp"
	if(FALSE){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",ifelse(is.null(fold),paste0("bch_",batch),paste0("spF_",foldID)),"_",filename_ext,".tex"), width = 5.76, height = 4.84, standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_rep",reps,"_it",iter,"_",simnames[1],"-",simnames[3],"_",filename_ext,".pdf"),width = 5.76, height = 9.54, pointsize = 9) }
	bayesplot_theme_update(strip.text = element_text(size = 9), axis.text = element_text( size = 8 ))
	bayesplot::mcmc_trace(coda_list, facet_args = list(ncol=2))
	#bayesplot::mcmc_trace_highlight(post_sample, highlight = 5)
	dev.off()
	
	
#- 3.6) Trace-rank plots
	
	filename_ext <- "chain_trank_overlay"
	if(latex_device){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",ifelse(is.null(fold),paste0("bch_",batch),paste0("spF_",foldID)),"_",filename_ext,".tex"), width = 5.76, height = 7.76, pointsize = 9, standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_rep",reps,"_it",iter,"_",simnames[1],"-",simnames[3],"_",filename_ext,".pdf"),width = 5.76, height = 7.76, pointsize = 9) }
	bayesplot_theme_update(strip.text = element_text(size = 9), axis.text = element_text( size = 8 ))
	bayesplot::mcmc_rank_overlay(coda_list, iter1 = it_burnin, facet_args = list(ncol=2))
	dev.off()
	
	filename_ext <- "chain_trank_perChain"
	if(FALSE){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",ifelse(is.null(fold),paste0("bch_",batch),paste0("spF_",foldID)),"_",filename_ext,".tex"), width = 5.76, height = 4.84, standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_rep",reps,"_it",iter,"_",simnames[1],"-",simnames[3],"_",filename_ext,".pdf"),width = 12, height = 12) }
	bayesplot::mcmc_rank_hist(coda_list)
	dev.off()


rm(CA_params, par_sel)

}




###--- 4.) Spatial-block cross-validation


#- 4.1) Run Simulations for Validation


	#- 4.1.1) Get observed data for validation
	
	# get MHB data
	MHBdata <- readRDS("data/abund/validation_counts.rds")
	
	# get lookup table of which landscape cell belongs to which spatial fold
	lookup <- readRDS(paste0("data/spatial_blocking/spatial_blocks_lookup_",aggr_blocks,".rds"))
	# retain only those cells that are also in MHB set...
	lookup <- lookup[lookup$mhb_cell,]

	# reshape to long format
	MHBdata_l <- reshape(MHBdata, direction = "long", varying = list(2:23), timevar = "Year", idvar = "cells", v.names = "Count")
	
	
	#--- Simulation of posterior predictions or load results
	
	if(FALSE){ # comment out simulation loop
		
		
		# 4.1.2) Prepare RS setup
		
		# number of draws from the posterior
		nr_postdraws <- 999
		
		# load one MCMC chain from the calibration results to get simulation setup that was used for calibration
		# filename <- paste0("results/CA_out/CA_DEzs_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,ifelse(is.null(foldID)|foldID==0,"",paste0("_spF",foldID)),"_rep",reps,"_it",iter,"_",simnames[1],".RData")
		# alternatively, save objects created in 'RangeShiftR_calibration.R', lines (1-213) and load that:
		filename <- "results/RSparam_Master.Rdata"
		
		attach(filename, name = "CAsetup")
		#s <- s
		if( !identical(aggr_blocks,get("aggr_blocks",pos="CAsetup")) ||
			!identical(simulatedData,get("simulatedData",pos="CAsetup"))||
			!identical(CA_params,get("CA_params",pos="CAsetup"))||
			!identical(par_sel,get("par_sel",pos="CAsetup")) ) {
			warning("Loaded calibration results don't match info in header.")
			detach("CAsetup", character.only = TRUE)
		}
		
		# set number of RS replicates to 1
		s@simul@Replicates <- 1
		
		# time frame for prediction and number of output years
		s@simul@Years <- 24        # 1997 - 2020
		#s@simul@OutStartPop <-  2 # = Year 1999
		
		
		
		#- 4.1.3) Simulation loop
		
		# loop over spatial folds
		for (foldID in 1:max(fold)) {
			
			# find cells of current fold
			fold_cells <- lookup[lookup$foldID == foldID,"cells"]
			
			# load chain into own namespace
			simnames <- batches[,fold==foldID]
			filename <- paste0("results/CA_out/CA_DEzs_mhb_",samplr,"_agBlc",aggr_blocks,"_rep",reps,"_it",iter,"_",simnames[1],"-",simnames[3],ifelse(concat,"_concat",""),"_chain.RData")
			attach(filename, name = "currChain")
			message(paste("Loading",filename,"..."))
			
			# get sample from joint posterior
			post_sample <- getSample(get("sampler_list",pos="currChain"), start = it_burnin, numSamples = nr_postdraws )
			
			# detach namespace again
			detach("currChain", character.only = TRUE)
			
			# loop over samples
			for (smpl in 1:nr_postdraws){
				
				if(smpl %% 45 == 1) print(paste("fold",foldID,", sample #",smpl))
				
				# generate prediction -> returns raster stack
				smpl_pred <- runRSmodel(post_sample[smpl,])
				
				# extract MHB cells within current fold and for surveyed years
				smpl_pred <- raster::extract(smpl_pred[[ix_years_abd[ix_years_abd %in% seq(nlayers(smpl_pred))]]], fold_cells, df = TRUE)
				smpl_pred$ID <- fold_cells
					
				# add all predicted counts and their squares together to calculate mean & variance at the end
				if(smpl==1) {
					fold_pred_mean <- fold_pred_var  <- smpl_pred
					fold_pred_var[,-1]  <- (smpl_pred[,-1])^2
				}else{
					fold_pred_mean[,-1] <- fold_pred_mean[,-1] + smpl_pred[,-1]
					fold_pred_var[,-1] <- fold_pred_var[,-1] + (smpl_pred[,-1])^2
				}
			}
			
			# calculate final values of mean & variance in MHB cells for this fold
			fold_pred_mean[,-1] <- fold_pred_mean[,-1]/nr_postdraws
			fold_pred_var[,-1] <- (fold_pred_var[,-1]/nr_postdraws) - (fold_pred_mean[,-1])^2
			
			# store in data frame that holds results for all folds
			if(foldID==1) {
				post_pred_mean <- fold_pred_mean
				post_pred_var  <- fold_pred_var
			}else{
				post_pred_mean <- rbind( post_pred_mean, fold_pred_mean )
				post_pred_var  <- rbind( post_pred_var,  fold_pred_var  )
			}
		}
		
		
		#- 4.1.4) Process results
		
		# Validation results: 
		# reshape dataframes to long format
		post_pred_mean <- reshape(post_pred_mean, direction = "long", varying = list(2:23), timevar = "Year", idvar = "ID", v.names = "Abund")
		post_pred_var  <- reshape(post_pred_var , direction = "long", varying = list(2:23), timevar = "Year", idvar = "ID", v.names = "varAbund")
		
		if(identical(post_pred_mean$ID,post_pred_var$ID)){post_pred <- cbind(post_pred_mean,post_pred_var[3])
		}else{warning("Format error in prediction output")}
		
		pred_mhbcells <- merge.data.frame(MHBdata_l,post_pred,by.x = c("cells","Year"), by.y = c("ID","Year"))
		
		
		# detach namespace CA_setup again
		detach("CAsetup", character.only = TRUE)
		
		
		#- 4.1.5) Save all results to file
		save(pred_mhbcells, 
			 file = paste0("results/predict/CA_DEzs_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_valid.RData"))
	}else{
		# load predictions from file
		load(paste0("results/predict/CA_DEzs_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_valid.RData"))
	}
	
	
	pred_mhbcells$sdAbund <- sqrt(pred_mhbcells$varAbund)
	pred_mhbcells$uppAbund <- pred_mhbcells$Abund + pred_mhbcells$sdAbund
	pred_mhbcells$lowAbund <- pred_mhbcells$Abund - pred_mhbcells$sdAbund
	pred_mhbcells$lowAbund[pred_mhbcells$lowAbund<0] <- 0

	
	
#- 4.2) Validate all MHB samples (all site-year)

#- 4.2.a) with ROC & AUC
	
	require(pROC)
	pred_occ <- data.frame(pred_mhbcells, 
						   mhb_occ = as.numeric(pred_mhbcells$Count>0), 
						   pred_occ = (1-dpois(0, lambda = pred_mhbcells$Abund)))
	(occ_roc <- roc(mhb_occ~pred_occ, data = pred_occ, na.rm = TRUE, smooth = FALSE, ci = TRUE, ret = c("roc", "coords", "all_coords"), plot = TRUE ))
	
	occ_roc_df <- data.frame(spec = occ_roc$specificities, sens = occ_roc$sensitivities, thres = occ_roc$thresholds)
	
	# make plot
	filename_ext <- "AUCallFold"
	if(latex_device){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_",filename_ext,".tex"), width = 2.52, height = 1.9, standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_",filename_ext,".pdf"),width = 6, height = 4) }
	plot(occ_roc)
	dev.off()
	
	
	# repeat this but with only those cells that change occupancy
	MHBocc <- cbind(MHBdata[[1]], ifelse(MHBdata[,-1]>0,1,0) )
	cell_chg <- MHBocc[apply(MHBocc[,-1],1,var,na.rm=T)>0,1]
	pred_occ_chg <- pred_occ[pred_occ$cells%in%cell_chg,]
	(occ_roc_chg <- roc(mhb_occ~pred_occ, data = pred_occ_chg, na.rm = TRUE, smooth = FALSE, ci = TRUE, ret = c("roc", "coords", "all_coords"), plot = TRUE ))
	
	
	# AUC for individual folds
	auc_perFold <- sapply(1:max(fold), function(foldID){
		fold_cells <- lookup[lookup$foldID == foldID,"cells"]
		roc(mhb_occ~pred_occ, data = pred_occ[pred_occ$cells%in%fold_cells,], na.rm = TRUE, smooth = FALSE, ci = TRUE, plot = FALSE, ret = c("roc", "coords"))
	})
	# add a column for all folds combined
	auc_perFold <- cbind(roc(mhb_occ~pred_occ, data = pred_occ, na.rm = TRUE, smooth = FALSE, ci = TRUE, plot = FALSE, ret = c("roc", "coords")), auc_perFold)
	
	# keep only AUC and 95% CIs
	auc_perFold_stats <- sapply(1:ncol(auc_perFold), function(fold){ as.numeric(auc_perFold['ci',fold]$ci) })
	
	# make dataframe
	auc_perFold_df <- data.frame(factor(c("all",paste("Fold",1:5)), levels=rev(c("all",paste("Fold",1:5))), ordered = TRUE), t(auc_perFold_stats) )
	names(auc_perFold_df) <- c("Fold ID","lower bound","AUC","upper bound")
	

#- 4.2.b) with C-index
	
# C-index for individual folds
	cix_perFold <- sapply(1:max(fold), function(foldID){
		fold_cells <- lookup[lookup$foldID == foldID,"cells"]
		with(pred_mhbcells[pred_mhbcells$cells%in%fold_cells,], rcorr.cens(Abund, Count))
	})
	cix_perFold_stats <- cix_perFold[c("C Index","S.D."),]
	
	# add a column for all folds combined
	cix_perFold_stats <- cbind(with(pred_mhbcells,rcorr.cens(Abund, Count))[c("C Index","S.D.")],cix_perFold_stats)
	
	# make dataframe
	cix_perFold_df <- data.frame(factor(c("all",paste("Fold",1:5)), levels=rev(c("all",paste("Fold",1:5))), ordered = TRUE), t(cix_perFold_stats) )
	cix_perFold_df$'95% CI' <- 1.96 * cix_perFold_df$S.D. / 2  # c.f. example of rcorr.cens()
	names(cix_perFold_df)[1:2] <- c("Fold ID","C index")
	

#- 4.3. make tables and plots
	
	# export latex table
	options(xtable.floating = FALSE)
	options(xtable.include.rownames = FALSE)
	options(xtable.booktabs = TRUE)
	print(xtable(auc_perFold_df, 
				 type = "latex",
				 include.rownames = FALSE,
				 align = "ccccc",
				 floating = FALSE),
		  file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_AUCtable.tex"))
	
	print(xtable(cix_perFold_df[!names(cix_perFold_df) %in% c('S.D.')],
				 type = "latex",
				 include.rownames = FALSE,
				 align = "cccc",
				 floating = FALSE),
		  file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_CIndextable.tex"))
	
	# ggplot
	auc_ggplot <- 
		ggplot(auc_perFold_df,
			   aes(
			   	x = AUC,
			   	y = `Fold ID`
			   )) +
		geom_point() +
		geom_errorbarh(aes(xmax = `upper bound`, xmin = `lower bound`, height = .2)) +
		xlim(NA, 1) +
		theme_bw() +
		theme(plot.title = element_text(hjust = 0.5),
			  panel.grid.major.x = element_blank(),
			  panel.grid.minor.x = element_blank() ) +
		#coord_flip() + 
		labs(x = "AUC", y="")
			 #title="AUC for all spatial folds")
	
	cix_ggplot <- 
		ggplot(cix_perFold_df,
			   aes(
			   	x = `C index`,
			   	y = `Fold ID`
			   )) +
		geom_point() +
		geom_errorbarh(aes(xmax = `C index`+`95% CI`, xmin = `C index`-`95% CI`, height = .2)) +
		xlim(NA, 1) +
		theme_bw() +
		theme(plot.title = element_text(hjust = 0.5),
			  panel.grid.major.x = element_blank(),
			  panel.grid.minor.x = element_blank() ) +
		#coord_flip() + 
		labs(x = "C Index", y="")
	#title="AUC for all spatial folds")
	
	# make plots
	filename_ext <- "AUCperFold"
	if(latex_device){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_",filename_ext,".tex"), width = 2.52, height = 1.9, standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_",filename_ext,".pdf"),width = 2.52, height = 1.9) }
	plot(auc_ggplot)
	dev.off()
	filename_ext <- "CIxperFold"
	if(latex_device){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_",filename_ext,".tex"), width = 2.52, height = 1.9, standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_",filename_ext,".pdf"),width = 2.52, height = 1.9) }
	plot(cix_ggplot)
	dev.off()
	
	
#- 4.3) Validate total abundances
	
	# calculate total abundances, MHB and posterior predictions, over all folds
	total_abunds <- lapply(1:max(fold), function(foldID){
		fold_cells <- lookup[lookup$foldID == foldID,"cells"]
		fold_data <- pred_mhbcells[pred_mhbcells$cells %in% fold_cells, c("Year","Count","Abund")]
		fold_sums <- aggregate(cbind(Count,Abund) ~ Year, dat = fold_data, sum)
		cbind(foldID = foldID , fold_sums[c("Year","Count","Abund")])
	})
	total_abunds <- do.call(rbind,total_abunds)
	
	# add rows for combined folds (labelled as foldID = 0)
	total_abunds <- rbind(total_abunds, cbind(foldID = 0, aggregate(cbind(Count,Abund) ~ Year, dat = pred_mhbcells[c("Year","Count","Abund")], sum)))
	
	
	# C-index for times series of total abundance within folds
	cindex_perFold <- sapply(fold, function(foldID){
		fold_data <- total_abunds[total_abunds$foldID==foldID,]
		rcorr.cens(fold_data$Abund, fold_data$Count)
	})
	
	# make dataframe
	cix_perFold_df <- data.frame(factor(c("all",paste("Fold",1:5)), levels=rev(c("all",paste("Fold",1:5))), ordered = TRUE), t(cindex_perFold[c("C Index","S.D."),]) )
	cix_perFold_df$'95% CI' <- 1.96 * cix_perFold_df$S.D. / 2  # c.f. example of rcorr.cens()
	names(cix_perFold_df)[1:2] <- c("Fold ID","C index")

	
	# export latex table
	options(xtable.floating = FALSE)
	options(xtable.include.rownames = FALSE)
	options(xtable.booktabs = TRUE)
	print(xtable(cix_perFold_df[!names(cix_perFold_df) %in% c('S.D.')],
				 type = "latex",
				 include.rownames = FALSE,
				 align = "cccc",
				 floating = FALSE),
		  file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_CItotalsTable.tex"))
	
	# ggplot
	cix_tots_ggplot <- 
		ggplot(cix_perFold_df,
			   aes(
			   	x = `C index`,
			   	y = `Fold ID`
			   )) +
		geom_point() +
		geom_errorbarh(aes(xmax = `C index`+`95% CI`, xmin = `C index`-`95% CI`, height = .2)) +
		xlim(layer_scales(cix_ggplot)$x$get_limits()) + # get limits from the other c-index plot for comparability
		theme_bw() +
		theme(plot.title = element_text(hjust = 0.5),
			  panel.grid.major.x = element_blank(),
			  panel.grid.minor.x = element_blank() ) +
		#coord_flip() + 
		labs(x = "C Index", y="")
	
	# make plots
	filename_ext <- "CItotals"
	if(latex_device){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_",filename_ext,".tex"), width = 2.52, height = 1.9, standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_",filename_ext,".pdf"),width = 2.52, height = 1.9) }
	plot(cix_tots_ggplot)
	dev.off()
	
	
	
#- 4.4) Validate cells with highest variance in abundance over time
	
	# calculate count variance in each cell over time and sort
	cell_var <- cbind(MHBdata[,1],apply(MHBdata[,-1], 1, var, na.rm = T))
	cell_var <- cell_var[order(cell_var[,2], decreasing = TRUE),]
	
	# select the 15% highest variance cells
	#high_var_cell <- MHBdata$cells[cell_var>0.3]
	high_var_cell <- cell_var[1:as.integer(nrow(cell_var)*0.15),1]
	
	# calculate AUC & C-index
	(roc(mhb_occ~pred_occ, data = pred_occ[pred_occ$cells%in%high_var_cell,], na.rm = TRUE, smooth = FALSE, ci = TRUE, ret = c("roc", "coords", "all_coords"), plot = TRUE))
	with(pred_mhbcells[pred_mhbcells$cells%in%high_var_cell,], rcorr.cens(Abund, Count))[c("C Index","S.D.")]
	
	# C-index for individual folds
	
	cix_perFold <- sapply(0:max(fold), function(foldID){
		if(foldID) fold_cells <- lookup[lookup$foldID == foldID,"cells"]
		else       fold_cells <- lookup[,"cells"]
		MHBsubset <- MHBdata[MHBdata[,1] %in% fold_cells, ]
		cell_var <- cbind(MHBsubset[,1],apply(MHBsubset[,-1], 1, var, na.rm = T))
		cell_var <- cell_var[order(cell_var[,2], decreasing = TRUE),]
		high_var_cell <- cell_var[1:as.integer(nrow(cell_var)*0.15),1]
		with(pred_mhbcells[pred_mhbcells$cells%in%high_var_cell,], rcorr.cens(Abund, Count))
	})
	cix_perFold_stats <- cix_perFold[c("C Index","S.D."),]
	
	# make dataframe
	cix_perFold_df <- data.frame(factor(c("all",paste("Fold",1:5)), levels=rev(c("all",paste("Fold",1:5))), ordered = TRUE), t(cix_perFold_stats) )
	cix_perFold_df$'95% CI' <- 1.96 * cix_perFold_df$S.D. / 2  # c.f. example of rcorr.cens()
	names(cix_perFold_df)[1:2] <- c("Fold ID","C index")
	
	# export latex table
	options(xtable.floating = FALSE)
	options(xtable.include.rownames = FALSE)
	options(xtable.booktabs = TRUE)
	print(xtable(cix_perFold_df[!names(cix_perFold_df) %in% c('S.D.')],
				 type = "latex",
				 include.rownames = FALSE,
				 align = "cccc",
				 floating = FALSE),
		  file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_CIhighvarTable.tex"))
	
	# ggplot
	cix_highvar_ggplot <- 
		ggplot(cix_perFold_df,
			   aes(
			   	x = `C index`,
			   	y = `Fold ID`
			   )) +
		geom_point() +
		geom_errorbarh(aes(xmax = `C index`+`95% CI`, xmin = `C index`-`95% CI`, height = .2)) +
		#xlim(layer_scales(cix_ggplot)$x$get_limits()) + # get limits from the other c-index plot for comparability
		xlim(NA, 1) +
		theme_bw() +
		theme(plot.title = element_text(hjust = 0.5),
			  panel.grid.major.x = element_blank(),
			  panel.grid.minor.x = element_blank() ) +
		#coord_flip() + 
		labs(x = "C Index", y="")
	
	# make plots
	filename_ext <- "CIhighvar"
	if(latex_device){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_",filename_ext,".tex"), width = 2.52, height = 1.9, standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_",filename_ext,".pdf"),width = 2.52, height = 1.9) }
	plot(cix_highvar_ggplot)
	dev.off()










###--- 5.) Projections

	
	#- 5.1) Run simulations for Projections
	
	# select calibration to full dataset for making projections
	foldID = 0
	simnames <- batches[,fold==foldID]
	
	# Prior predictions?
	prior_predictions <- FALSE
	
	# Use a habitat suitability scenario?
	HSI_scenario <- FALSE
	
	if(FALSE){ # comment out simulation loop
		
		# 5.1.1) Prepare RS setup
		
		# number of draws from the posterior
		nr_postdraws <- 999
		
		# load one MCMC chain from the calibration results to get simulation setup that was used for calibration
		  #filename <- paste0("results/CA_out/CA_DEzs_mhb_GamPois_agBlc14_spF2_rep20_it3e+05_521.RData")
		# alternatively, save objects created in 'RangeShiftR_calibration.R', lines (1-213) and load that:
		filename <- "results/RSparam_Master.Rdata"
		
		attach(filename, name = "CAsetup")
		#s <- s
		if( !identical(aggr_blocks,get("aggr_blocks",pos="CAsetup")) ||
			!identical(simulatedData,get("simulatedData",pos="CAsetup"))||
			!identical(CA_params,get("CA_params",pos="CAsetup"))||
			!identical(par_sel,get("par_sel",pos="CAsetup")) ) {
			warning("Loaded calibration results don't match info in header.")
			detach("CAsetup", character.only = TRUE)
		}

		
		# set number of RS replicates to 1
		s@simul@Replicates <- 1
		
		# Extend simulation time by 30 years into future
		s@simul@Years <- 54        # 1997 - 2050
		#s@simul@OutStartPop <-  2 # = Year 1999
		
		
		# Prepare objects for storage
		
		# prepare empty raster stack to store mean predicted abundances over all posterior draws
		
		# load one of the habitatmaps to use as a template
		habitatmap <- raster(paste0("data/habitatmaps/RK_hsi_2km_2000.asc"))
		pred_raster_empty <- habitatmap # copy map
		pred_raster_empty[!is.na(pred_raster_empty)] <- 0 # reset all values
		pred_raster_mean <- stack(replicate((s@simul@Years-s@simul@OutStartPop),pred_raster_empty)) # replicate for each output year
		names(pred_raster_mean) <- paste0("year_",as.character(1997+(s@simul@OutStartPop:(s@simul@Years-1))))
		pred_raster_var <- pred_raster_mean # same template to store variances
		
		# mask of Swiss borders
		bg <- raster('/vsicurl/https://damariszurell.github.io/SDM-Intro/CH_mask.tif')
		bg <- raster::crop(bg,habitatmap)
		bg <- aggregate(bg,2,sum)
		
		# empty array to store time series of total abundance for each posterior draw
		pred_total <- array(NA, dim = c(nr_postdraws, nlayers(pred_raster_mean)))
		dimnames(pred_total)[[2]] <- names(pred_raster_mean)
		
		# Scenario?
		 # For counter-factual scenario with altered land-use change:
		#HSI_scenario <- "plus" # "mnus"
		if(is.character(HSI_scenario)){
			# Set dynamic landscapes to false:
			#s@land@DynamicLandYears <- 0
			# Set landscape to first landscape only:
			
			habitat_names <- paste0("data/habitatmaps/RK_hsi_",HSI_scenario,"_2km_",sprintf("%02d",1998+1:nrYrSDM),".asc")
			habitatmap <- raster(habitat_names[20])
			habitat_matrices <- lapply(habitat_names,
									   FUN = function(filename){
									   	matrix(raster(filename),
									   		   ncol = ncol(habitatmap),
									   		   nrow = nrow(habitatmap), byrow = TRUE)})
			s@land@LandscapeFile <- habitat_matrices
		}
		
		#- 5.1.a) Prior predictions
		if(prior_predictions){
			
			# draw from prior from bayesianSetup
			nr_priordraws <- nr_postdraws
			prior_sample <- get("setup","CAsetup")$prior$sampler(nr_priordraws)
			
			# loop over samples
			for (smpl in 1:nr_priordraws){
				
				if(smpl %% 50 == 1) print(paste("sample #",smpl))
				
				# generate prediction -> returns raster stack
				smpl_pred <- runRSmodel(prior_sample[smpl,])
				
				# if simulation ends early, fill up raster stack with empty maps
				if( nlayers(smpl_pred) < nlayers(pred_raster_mean) ){
					smpl_pred <- stack(smpl_pred, stack(replicate(nlayers(pred_raster_mean)-nlayers(smpl_pred),pred_raster_empty)) )
				}
				
				# store mean & var of predicted abundances for each cell in this fold
				pred_raster_mean <- pred_raster_mean + smpl_pred
				pred_raster_var <- pred_raster_var + (smpl_pred)^2
				
				# get total abundance
				smpl_pred <- raster::mask(smpl_pred, bg) # crop to Switzerland to exclude population in buffer area
				pred_total[smpl,] <- sapply(1:nlayers(smpl_pred), FUN = function(lyr_ix){ sum(values(smpl_pred[[lyr_ix]]), na.rm=T )})
			}

		}else{
		
		#- 5.1.b) Posterior predictions
			
			# load chain into own namespace
			filename <- paste0("results/CA_out/CA_DEzs_mhb_",samplr,"_agBlc",aggr_blocks,"_rep",reps,"_it",iter,"_",simnames[1],"-",simnames[3],ifelse(concat,"_concat",""),"_chain.RData")
			attach(filename, name = "currChain")
			
			# get sample from joint posterior
			post_sample <- getSample(get("sampler_list",pos="currChain"), start = it_burnin, numSamples = nr_postdraws )
			
			# detach namespace again
			detach("currChain", character.only = TRUE)
			
			# loop over samples
			for (smpl in 1:nr_postdraws){
				
				if(smpl %% 50 == 1) print(paste("fold",foldID,", sample #",smpl))
				
				# generate prediction -> returns raster stack
				smpl_pred <- runRSmodel(post_sample[smpl,])
				
				# store mean & var of predicted abundances for each cell in this fold
				pred_raster_mean <- pred_raster_mean + smpl_pred
				pred_raster_var <- pred_raster_var + (smpl_pred)^2
				
				# get total abundance
				smpl_pred <- raster::mask(smpl_pred, bg) # crop to Switzerland to exclude population in buffer area from the total
				pred_total[smpl,] <- sapply(1:nlayers(smpl_pred), FUN = function(lyr_ix){ sum(values(smpl_pred[[lyr_ix]]), na.rm=T )})
				
			}
		}
		
		#- 5.2.) Process results
		#- (same for prior and posterior predictions)
		
		# Calculate mean & var of total abundance over samples
		pred_total_means <- apply(pred_total, 2, mean, na.rm = TRUE)
		pred_total_vars  <- apply(pred_total, 2, var,  na.rm = TRUE)
		
		# Calculate final values of mean & variance for raster stacks
		pred_raster_mean <- pred_raster_mean/nr_postdraws
		pred_raster_var <- pred_raster_var/nr_postdraws - (pred_raster_mean)^2
		
		
		# save all results to file
		save(pred_total, pred_total_means, pred_total_vars, pred_raster_mean, pred_raster_var,
			 file = paste0("results/predict/CA_DEzs_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(simnames),"-",max(simnames),ifelse(prior_predictions,"_prior","_post"),"_predict",ifelse(is.character(HSI_scenario),paste0("_",HSI_scenario),""),".RData"))
		
		# detach namespace CA_setup again
		detach("CAsetup", character.only = TRUE)
		
	}else{
		
		# load predictions from file
		load(paste0("results/predict/CA_DEzs_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(simnames),"-",max(simnames),ifelse(prior_predictions,"_prior","_post"),"_predict",ifelse(is.character(HSI_scenario),paste0("_",HSI_scenario),""),".RData"))
	}






#- 5.2) Maps

	# put together the maps to plot
	plot_years <- c(1,21,41)        					 # Years 1999, 2019, 2019+20
	plot_stack <- pred_raster_mean[[plot_years]]
	plot_stack[[4]] <- plot_stack[[2]]-plot_stack[[1]]   # Difference 1999 - 2019
	plot_stack[[5]] <- plot_stack[[3]]-plot_stack[[2]]   # Difference 2019 - 2019+20
	names(plot_stack)[c(4:5)] <- c("diff_1999_2019","diff_2019_2039")
	spplot(plot_stack)
	
	# load Swiss borders and hillshade
	bg <- raster('/vsicurl/https://damariszurell.github.io/SDM-Intro/CH_mask.tif')
	bg <- raster::crop(bg,plot_stack)
	load("data/hillshade/hill_ch.Rdata")
	hill_ch <- projectRaster(hill_ch, bg, method = 'ngb')
	values(hill_ch$layer)[is.na(values(hill_ch$layer))] <- 1
	hill_ch <- raster::mask(hill_ch, bg)
	#hill_ch <- raster::crop(hill_ch,plot_stack)
	
	# crop to Swiss borders
	bg <- aggregate(bg,2,sum)
	plot_stack <- raster::mask(plot_stack, bg)
	
	# multiply by 25 to get units of breeding pairs per 100km²
	plot_stack <- 25*plot_stack
	
	# settings
	plot_names <- c("Year 1999","Year 2019","Year 2019+20","Difference 1999 - 2019","Difference 2019 - 2019+20")
	names(plot_stack) <- plot_names
	my.theme <-
	    list(strip.background = list(col="gray90"), #,strip.border=list(col="transparent")
	    	 par.strip.text = list(cex=ifelse(latex_device,.9,.7)),
	    	 layout.heights =
		         list(top.padding = 0.2,
		              main.key.padding = 0,
		              key.axis.padding = 0,
		              axis.xlab.padding = 0,
		              xlab.key.padding = 0,
		              key.sub.padding = 0,
		              bottom.padding = 0.72),
	         layout.widths =
		         list(left.padding = 0.2,
		              key.ylab.padding = 0,
		              ylab.axis.padding = 0,
		              axis.key.padding = 1,
		         	  key.left <- 1.2,
		              right.padding = 0.2))
	col_palette_grain <- 7
	
	
	# plot absolute abundances
	
	north_arrow <- list("SpatialPolygonsRescale", layout.north.arrow(type=2), offset = c(496000, 236000), scale = 22*2000, which = 1)
	scale_bar   <- list("SpatialPolygonsRescale", layout.scale.bar(), offset = c(742000, 98000), scale = 36*2000, fill=c("white","black"), which = 1)
	scale_bar_text0 <- list("sp.text", c(744000, 87000), "0", cex = .85, which = 1)
	scale_bar_text1 <- list("sp.text", c(812000, 87000), "2000 m", cex = .85, which = 1)
	
	max_abs <- max(abs(values(plot_stack[[c(1:3)]])),na.rm=T)
	brk <- do.breaks(c(0, max_abs), col_palette_grain)
	# Colors
	cols <- colorRampPalette(c("#EEEEEE","#000E76"))(col_palette_grain)
	my.theme$layout.widths$left.padding <- 0.2
	my.theme$layout.heights$top.padding <- 0.2
	my.theme$layout.heights$bottom.padding <- 0.72
	my.theme$layout.widths$key.left <- 1.2
	filename_ext <- "abund_maps_abs"
	if(latex_device){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(simnames),"-",max(simnames),"_",filename_ext,ifelse(prior_predictions,"_prior","_post"),".tex"), width = 2.9, height = 5.4, standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(simnames),"-",max(simnames),"_",filename_ext,ifelse(prior_predictions,"_prior","_post"),".pdf"), width = 3.9, height = 6.2) }
		spplot(plot_stack[[c(1:3)]],
		   col.regions = sapply(cols, t_col, 30),
		   at = brk,
		   names.attr = plot_names[c(1:3)],
		   colorkey = list(space = "left", width = .75,
		   				title = list(
		   					label = ifelse(latex_device,"Population density [BP / 100m$^2$]","Population density [BP / 100m^2]"),
		   					cex = 1.0,
		   					rot = 90),
		   				title.control = list(side = "left", padding = 3),
		   				axis.text = list(cex = 1.0)),
		   layout = c(1,3),
		   par.settings = my.theme,
		   sp.layout = list(north_arrow, scale_bar, scale_bar_text0, scale_bar_text1)) +
		as.layer(spplot(hill_ch$layer, col.regions = hcl.colors(16, palette = "Light Grays")), under = TRUE)
	dev.off()
	
	
	# plot abundance differences
	
	max_abs <- max(abs(values(plot_stack[[c(4:5)]])),na.rm=T)
	brk <- do.breaks(c(-max_abs, max_abs), 9)
	cols <- brewer.pal(10, "PiYG")
	cols <- sapply(cols,makeColDarker,factor = 0.25)
	cols[6] <- "#EEEEEE"
	my.theme$layout.widths$right.padding <- 1.4
	my.theme$layout.widths$left.padding <- 0.2
	filename_ext <- "abund_maps_diff"
	if(latex_device){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(simnames),"-",max(simnames),"_",filename_ext,ifelse(prior_predictions,"_prior","_post"),".tex"), width = 2.9, height = 3.6, standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(simnames),"-",max(simnames),"_",filename_ext,ifelse(prior_predictions,"_prior","_post"),".pdf"), width = 3.9, height = 4.2) }
	spplot(plot_stack[[c(4:5)]],
		   layout = c(1,2),
		   col.regions = sapply(cols, t_col, 30),
		   at = brk,
		   names.attr = plot_names[c(4:5)],
		   colorkey = list(space = "right", width = .75,
		   				   title = list(
		   				   		label = ifelse(latex_device,"Population density [BP / 100m$^2$]","Population density [BP / 100m^2]"),
		   				   		cex = 1.0,
		   				   		rot = -90),
		   				   title.control = list(side = "right", padding = 2),
		   				   axis.text = list(cex = 1.0)),
		   par.settings = my.theme,
		   sp.layout = list(north_arrow, scale_bar, scale_bar_text0, scale_bar_text1)) +
		as.layer(spplot(hill_ch$layer, col.regions = hcl.colors(16, palette = "Light Grays")), under = TRUE)
	dev.off()
	
	

#- 5.3) Time series
	
	# get quantiles to plot
	total_abund <- data.frame(years = 1999:2050,
							  t(apply(pred_total, 2, quantile, probs = c(0.025, 0.500, 0.975), na.rm = FALSE))
	)
	names(total_abund)[2:4] <- c("pred_lower", "pred_mean", "pred_upper")
	
		# earlier version: plotting mean and SD
		#total_abund <- data.frame(years = 1999:2050,
		#						  pred_mean  = pred_total_means["all_fold",], 
		#						  pred_upper = pred_total_means["all_fold",] + 2*sqrt(pred_total_vars["all_fold",]),
		#						  pred_lower = pred_total_means["all_fold",] - 2*sqrt(pred_total_vars["all_fold",]))
		
	
	# add BBI data
	load(file = "../Observed_data/breeding_bird_index.rds")
	total_abund$mhb_rel <- c(bbi_RK$Model, rep(NA,31))
	
	# calculate relative abundances
	ref_pred_abund <- total_abund$pred_lower[1] * 0.95
	total_abund$pred_mean_rel <- total_abund$pred_mean / ref_pred_abund
	total_abund$pred_upper_rel <- total_abund$pred_upper / ref_pred_abund
	total_abund$pred_lower_rel <- total_abund$pred_lower / ref_pred_abund
	
	# Colors
	pred_col <- "#62708EFF"
	mhb_col <- "#C9704EFF"
	index_col <- "#AFA800FF"
	
	# plot
	filename_ext <- "PredTotalTS"
	if(latex_device){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(simnames),"-",max(simnames),"_",filename_ext,ifelse(prior_predictions,"_prior","_post"),".tex"), width = 3.8, height = ifelse(prior_predictions,1.7,2.2), standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(simnames),"-",max(simnames),"_",filename_ext,ifelse(prior_predictions,"_prior","_post"),".pdf"),width = 4.8, height = ifelse(prior_predictions,2.2,2.8)) }
		par(mar=c(3,3,0.2,3), mgp = c(2.0, 0.66, 0) )
		ylim = c(0,9.2)
		if(latex_device){ par(lwd = 1.7) }
		map_years <-  c(1999, 2019, 2019+20) # year index that marks end of time window for aggregation
		if(!prior_predictions){ 
			total_plot <- with(total_abund, {
				plot(NULL, xlim = range(years), ylim = ylim, # ylim = range(c(pred_lower_rel,pred_upper_rel,mhb_rel),na.rm = T),
					 xlab = "Year", ylab="Relative abundance", axes = FALSE,
					 main = "" ) # "Total Swiss red kite abundance"
				axis(1, at = c(1999,2009,2019,2029,2039,2049), labels = c('1999','','2019','','2019+20',''))
				axis(2)
				box()
				polygon(c(years,rev(years)),c(pred_lower_rel,rev(pred_upper_rel)), col = t_col(pred_col,70), border = NA )
				for(year in map_years) abline(v=year, lty = 2, lwd = 1.7, col = index_col)
				lines(years,pred_lower_rel, type = "l", col = t_col(pred_col,30), lwd = 1)
				lines(years,pred_upper_rel, type = "l", col = t_col(pred_col,30), lwd = 1)
				lines(years,pred_mean_rel,  type = "p", col = pred_col, pch = 16, cex = 0.5, lwd = 2)
				lines(years,mhb_rel, type = "b", col = mhb_col, lwd = 2)
				par(new = TRUE) # Add new plot
				#plot(years,mhb_rel, type = "b", col = "red", lwd = 2, axes = FALSE, xlab = "", ylab = "") # Create second plot without axes
				plot(NULL, axes = FALSE, xlab = "", ylab = "", 
					 xlim = range(years), ylim = ylim *ref_pred_abund ) # ylim = range(c(pred_lower_rel,pred_upper_rel,mhb_rel),na.rm = T)*ref_pred_abund) # Create second plot without axes
				axis(side = 4, at = c(0,2000,4500,7000)) # second y-axis
				mtext("Absolute abundance", side = 4, line = 2)
				legend(x=ifelse(latex_device,2020,2021), y=ifelse(latex_device,2300,1800), legend=c("Observed","Simulated"), col=c( mhb_col, pred_col), pch=c(16,16), bty="n")
			})
		}else{
			par(mar=c(0.1,3,0.2,3), mgp = c(2.0, 0.66, 0) )
			total_plot <- with(total_abund, {
				plot(NULL, xlim = range(years), ylim = ylim, # ylim = range(c(pred_lower_rel,pred_upper_rel)),
					 xaxt = "n", xlab = "", #"Year", 
					 ylab="Relative abundance",
					 main = "" ) # "Total Swiss red kite abundance"
				polygon(c(years,rev(years)),c(pred_lower_rel,rev(pred_upper_rel)), col = t_col(pred_col,70), border = NA )
				for(year in map_years) abline(v=year, lty = 2, lwd = 1.7, col = index_col)
				lines(years,pred_lower_rel, type = "l", col = t_col(pred_col), lwd = 1)
				lines(years,pred_upper_rel, type = "l", col = t_col(pred_col,50), lwd = 1)
				lines(years,pred_mean_rel,  type = "p", col = pred_col, pch = 16, cex = 0.5, lwd = 2)
				par(new = TRUE) # Add new plot
				plot(NULL, axes = FALSE, xlab = "", ylab = "", 
					 xlim = range(years), ylim = ylim *ref_pred_abund ) #  ylim = range(c(pred_lower_rel,pred_upper_rel,mhb_rel),na.rm = T)*ref_pred_abund) # Create second plot without axes
				axis(side = 4, at = c(0,2000,4500,7000)) # second y-axis
				mtext("Absolute abundance", side = 4, line = 2)
			})
		}
	dev.off()

	
	
#- 5.4) Validation
	
	# For comparison with cross-validation, calculate AUC and C-index with calibration from full dataset
	
	# get MHB data
	MHBdata <- readRDS("data/abund/validation_counts.rds")
	
	# reshape to long format
	MHBdata_l <- reshape(MHBdata, direction = "long", varying = list(2:23), timevar = "Year", idvar = "cells", v.names = "Count")
	
	# extract mean predictions from raster stack
	pred_mhb <- raster::extract(pred_raster_mean[[1:22]], MHBdata$cells, df = TRUE)
	pred_mhb$ID <- MHBdata$cells
	
	# reshape predictions to long format
	pred_mhb_l <- reshape(pred_mhb, direction = "long", varying = list(2:23), timevar = "Year", idvar = "ID", v.names = "Pred")
	
	# merge observations and predictions
	pred_mhbcells <- merge.data.frame(MHBdata_l, pred_mhb_l, by.x = c("cells","Year"), by.y = c("ID","Year"))
	
	# add observed and predicted occupancy
	pred_mhbcells <- data.frame(pred_mhbcells,
								mhb_occ = as.numeric(pred_mhbcells$Count>0),
								pred_occ = (1-dpois(0, lambda = pred_mhbcells$Pred)))
	
	
	## AUC for calibration with full data set
	
	(occ_roc <- roc(mhb_occ~pred_occ, data = pred_mhbcells, na.rm = TRUE, smooth = FALSE, ci = TRUE, ret = c("roc", "coords", "all_coords"), plot = TRUE ))
	#occ_roc_df <- data.frame(spec = occ_roc$specificities, sens = occ_roc$sensitivities, thres = occ_roc$thresholds)
	# AUC: 0.8906  (0.8812-0.9000)
	
	## C-index for calibration with full data set (spatio-temporal)
	with(pred_mhbcells, rcorr.cens(Pred, Count))
	# c index: 0.8817  (SD: 0.0090 )
	
	## C-index for total abundance
	total_abunds <- aggregate(cbind(Count,Pred) ~ Year, dat = pred_mhbcells[c("Year","Count","Pred")], sum)
	with(total_abunds, rcorr.cens(Pred, Count))
	# c index: 0.9437  (SD: 0.0452 )
	
	## C-index for highly-variable cells
	cell_var <- cbind(MHBdata[,1],apply(MHBdata[,-1], 1, var, na.rm = T))  # calculate count variance in each cell over time and sort
	cell_var <- cell_var[order(cell_var[,2], decreasing = TRUE),]          # sort by highest variance
	high_var_cell <- cell_var[1:as.integer(nrow(cell_var)*0.15),1]         # select the 15% highest variance cells
	# calculate AUC & C-index
	(roc(mhb_occ~pred_occ, data = subset(pred_mhbcells, cells %in% high_var_cell), na.rm = TRUE, smooth = FALSE, ci = TRUE, ret = c("roc", "coords", "all_coords"), plot = TRUE))
	with(subset(pred_mhbcells, cells %in% high_var_cell), rcorr.cens(Pred, Count))
	# c index: 0.6549  (SD: 0.0299 )
	
	

	
#- 5.5) Plot comparative time series of scenarios
#
	# 5.5.1) load predictions
	
	# a) plus scenario
	attach(paste0("results/predict/CA_DEzs_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(simnames),"-",max(simnames),ifelse(prior_predictions,"_prior","_post"),"_predict_plus.RData"),
		   name = "readPred")
	# save predictions in global env.
	pred_total_plus <- get("pred_total", pos = "readPred")
	# detach namespace again
	detach("readPred", character.only = TRUE)
	
	# b) minus scenario
	attach(paste0("results/predict/CA_DEzs_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(simnames),"-",max(simnames),ifelse(prior_predictions,"_prior","_post"),"_predict_mnus.RData"),
		   name = "readPred")
	# save predictions in global env.
	pred_total_mnus <- get("pred_total", pos = "readPred")
	# detach namespace again
	detach("readPred", character.only = TRUE)

	# c) factual scenario
	attach(paste0("results/predict/CA_DEzs_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(simnames),"-",max(simnames),ifelse(prior_predictions,"_prior","_post"),"_predict.RData"),
		   name = "readPred")
	# save predictions in global env.
	pred_total_real <- get("pred_total", pos = "readPred")
	# detach namespace again
	detach("readPred", character.only = TRUE)
	
	
	
	# 5.5.2) plot time series
	 

	# extract time series to plot
	total_abund <- data.frame(years = 1999:2050,
							  plus = as.numeric(apply(pred_total_plus, 2, mean, na.rm = FALSE)),
							  mnus = as.numeric(apply(pred_total_mnus, 2, mean, na.rm = FALSE)),
							  real = as.numeric(apply(pred_total_real, 2, mean, na.rm = FALSE))
	)
	
	# plot
	scen_col <- "#769C1D"
	filename_ext <- "ScenariosTS"
	if(latex_device){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(simnames),"-",max(simnames),"_",filename_ext,".tex"), width = 4.5, height = 3.2, standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(simnames),"-",max(simnames),"_",filename_ext,".pdf"),width = 6, height = 4) }
		with(total_abund, {
			plot(NULL, xlim = range(years), ylim = c(0,max(total_abund[-1],na.rm = T)), axes = FALSE,
				 xlab = "Year", ylab="Absolute abundance",
				 main = "" ) #ifelse(latex_device,"","Total Swiss red kite abundance"))
			axis(1, at = c(1999,2009,2019,2029,2039,2049), labels = c('1999','','2019','','2019+20',''))
			axis(2)
			box()
			lines(years,plus, type = "l", col = scen_col, lwd = 1.5)
			lines(years,mnus, type = "l", col = scen_col, lwd = 1.5)
			lines(years,real, type = "l", col = pred_col, lwd = 1.5)
			legend('bottomright', c('factual', 'counter-factual'), lty = 1, lwd = 2, col = c(pred_col,scen_col), bg = "white", inset = c(0.06,0.06))
		})
	dev.off()
	