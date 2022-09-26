
###-----------------------------------------------------------------
##   GENERATE OUTPUTS AND PLOTS FROM CALIBRATION RESULTS
###-----------------------------------------------------------------

# Main script for analysing the posterior sample and for plotting results; loads the functions stored in the other files in the same folder
# Note: The sensitivity analysis and its plotting are not included here, but can be found in its own script plot_sensitivity.R.

# Contents:
# 1.) Create and store combined MCMCs
# 2.) Plot posterior marginal distributions & its pairwise correlations
# 3.) MCMC diagnostic plots
# 4.) Run Simulations for Validation & Prediction
# 5.) Validation
# 6.) Predictions


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


# read functions from other files in this folder

source("scripts/analysis/docu_colors.R")
source("scripts/analysis/combine_chains.R")
source("scripts/analysis/get_MAPsection.R")
source("scripts/analysis/plot_posterior.R")

# choose type of graphics device: latex (TRUE) or pdf (FALSE)?
latex_device <- TRUE

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

concat <- TRUE # if TRUE, combineMCMCs() returns a Cosa object; if FALSE, a BayesianTools::mcmcSampler object 


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



###--- 2.) Plot posterior marginal distributions & its pairwise correlations

# get parameter ranges
CA_params  <- data.frame(name = c( "DensDep","Fecund", "Surv1", "Surv2", "Surv3", "Deve1", "Deve2", "EmigProb","DispDist","SettBeta","SettProb","InitPop","sigma","GPsize"),
						 min  = c(   0.1e-2 ,   0.50 ,   0.01 ,   0.01 ,   0.01 ,   0.01 ,   0.01 ,      0.01 ,      400 ,     -15  ,     0.01 ,     200 ,  0.01,    1.0  ),
						 def  = c(   0.6e-2 ,   1.65 ,   0.42 ,   0.68 ,   0.80 ,   0.80 ,   0.55 ,      0.64 ,     2200 ,       4  ,     0.85 ,    1500 ,  0.20,   50.0  ),
						 max  = c(   2.0e-2 ,   5.00 ,   0.99 ,   0.99 ,   0.99 ,   0.99 ,   0.99 ,      0.99 ,     2400 ,      15  ,     0.99 ,    2500 ,  0.90,  500.0 ))
par_sel <- c(1:8,10,14) #


if(FALSE){ # comment out the plotting

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
	
	filename = paste0("results/CA_out/CA_DEzs_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_rep",reps,"_it",iter,"_",simnames[1],"-",simnames[3],"_chain.RData")
	attach(filename, name = "currChain")
	sampler_list <- get("sampler_list",pos="currChain")
	detach("currChain", character.only = TRUE)
	if(concat) {
		filename = paste0("results/CA_out/CA_DEzs_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_rep",reps,"_it",iter,"_",simnames[1],"-",simnames[3],"_concat_chain.RData")
		attach(filename, name = "currChain")
		coda_list <- get("sampler_list",pos="currChain")
		detach("currChain", character.only = TRUE)
	}
	
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
}




###--- 4.) Run Simulations (for Validation & Prediction)


#- 4.1) Prepare RS setup

	# number of draws from the posterior
	nr_postdraws <- 999
	
	# load one of the calibration results files to get simulation setup that was used for calibration
	foldID = 0; simnames = 401:403
	filename <- paste0("results/CA_out/CA_DEzs_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,ifelse(is.null(foldID)|foldID==0,"",paste0("_spF",foldID)),"_rep",reps,"_it",iter,"_",simnames[1],".RData")
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
	s@simul@Years <- 54        # 1997 - 2050
	#s@simul@OutStartPop <-  2 # = Year 1999
	#s@simul@OutStartPop <- 23 # = Year 2020

#- 4.2) Prepare objects for Validation
	
	# get MHB data
	MHBdata <- readRDS("../Observed_data/MHB_data/abunds_wide_1999_2020.rds")
	# reshape to long format
	MHBdata_l <- reshape(MHBdata, direction = "long", varying = list(2:23), timevar = "Year", idvar = "cells", v.names = "Count")
	
	# prepare empty data frames to store mean and variance of simulated abundance at MHB cells
	post_pred_mean <- post_pred_var <- data.frame()
	
	# get lookup table of which landscape cell belongs to which spatial fold
	lookup <- readRDS(paste0("scripts/calibration/spatial_blocks_lookup_",aggr_blocks,".rds"))

#- 4.3) Prepare objects for Prediction
	
	# load one of the habitatmaps to use as a template
	habitatmap <- raster(paste0("../Observed_data/Atlas_data/SDM_results/asciimaps/RK_hsi_2km_2000.asc"))
	
	# prepare empty raster stack to store predicted abundances
	pred_raster_empty <- habitatmap # copy map
	pred_raster_empty[!is.na(pred_raster_empty)] <- 0 # reset all values
	pred_raster_mean <- stack(replicate((s@simul@Years-s@simul@OutStartPop),pred_raster_empty)) # replicate for each output year
	names(pred_raster_mean) <- paste0("year_",as.character(1997+(s@simul@OutStartPop:(s@simul@Years-1)))) #paste0("Year ",1999:2050)
	pred_raster_var <- pred_raster_mean # same template to store variances
	
	# empty array to store time series of total abundance in each fold for all posterior draws (this is needed to later get the variance right for time series of total Swiss abundance)
	post_pred_total <- array(NA, dim = c(nr_postdraws, ncol = max(fold)+1, nlayers(pred_raster_mean)))
	dimnames(post_pred_total)[[2]] <- c(paste0("fold_", 1:max(fold)),"all_fold")
	dimnames(post_pred_total)[[3]] <- names(pred_raster_mean)
	
	# mask of Swiss borders
	bg <- raster('/vsicurl/https://damariszurell.github.io/SDM-Intro/CH_mask.tif')
	bg <- raster::crop(bg,habitatmap)
	bg <- aggregate(bg,2,sum)
	
	# switch to produce prior predictions
	prior_predictions <- FALSE
	
#- 4.4) Run simulations

if(FALSE){ # comment out simulation loop
	# loop over spatial folds
	for (foldID in 1:max(fold)) {
		# find cells of current fold
		fold_cells <- lookup[lookup$foldID == foldID,"cells"]
		# load chain into own namespace
		simnames <- batches[,fold==foldID]
		filename <- paste0("results/CA_out/CA_DEzs_mhb_",samplr,"_agBlc",aggr_blocks,"_rep",reps,"_it",iter,"_",simnames[1],"-",simnames[3],"_chain.RData")
		attach(filename, name = "currChain")
		message(paste("Loading",filename,"..."))
		# get sample from joint posterior
		post_sample <- getSample(get("sampler_list",pos="currChain"), start = it_burnin, numSamples = nr_postdraws )
		# loop over samples
		for (smpl in 1:nr_postdraws){
			if(smpl %% 45 == 1) print(paste("fold",foldID,", sample #",smpl))
			# generate prediction -> returns raster stack
			smpl_pred <- runRSmodel(post_sample[smpl,])
			# take mean over replicates -- not needed here as only 1 rep
			#smpl_pred <- stackApply(smpl_pred, indices = as.integer(sapply(strsplit(names(smpl_pred),"year"), "[[", 2)) + 1 - nrYrBurnin, fun = mean)
			
			# store mean & var of predicted abundances for each cell in this fold
			pred_raster_mean[fold_cells] <- pred_raster_mean[fold_cells] + smpl_pred[fold_cells]
			pred_raster_var[fold_cells] <- pred_raster_var[fold_cells] + (smpl_pred[fold_cells])^2
			
			# get total abundance for this fold (this gets stored for each posterior sample)
			smpl_pred <- raster::mask(smpl_pred, bg) # crop to Switzerland to exclude population in buffer area from the total
			post_pred_total[smpl,foldID,] <- sapply(1:nlayers(smpl_pred), FUN = function(lyr_ix){ sum(values(smpl_pred[[lyr_ix]])[fold_cells], na.rm=T )})
			
			# extract MHB cells within current fold and for surveyed years
			smpl_pred <- raster::extract(smpl_pred[[ix_years_abd[ix_years_abd %in% seq(nlayers(smpl_pred))]]], fold_cells[fold_cells %in% MHBdata$cells], df = TRUE)
			smpl_pred$ID <- fold_cells[fold_cells %in% MHBdata$cells]
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
		# transfer these results to full results data frames
		post_pred_mean <- rbind(post_pred_mean, fold_pred_mean )
		post_pred_var  <- rbind(post_pred_var,  fold_pred_var  )
		
		# detach namespace again
		detach("currChain", character.only = TRUE)
	}
	if(prior_predictions){ # run this instead of the preceding block to make prior predictions
		# draw from prior from bayesianSetup
		nr_priordraws <- nr_postdraws
		prior_sample <- get("setup","CAsetup")$prior$sampler(nr_priordraws)
		prior_pred_total <- post_pred_total[,1,]
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
			# get total abundance for this fold (this gets stored for each posterior sample)
			smpl_pred <- raster::mask(smpl_pred, bg) # crop to Switzerland to exclude population in buffer area
			prior_pred_total[smpl,] <- sapply(1:nlayers(smpl_pred), FUN = function(lyr_ix){ sum(values(smpl_pred[[lyr_ix]]), na.rm=T )})
		}
		# calculate mean & var of total abundance over samples
		prior_total_means <- apply(prior_pred_total, 2, mean, na.rm = TRUE)
		prior_total_vars <- apply(prior_pred_total, 2, var, na.rm = TRUE)
		
		# save all results to file
		save(prior_pred_total, prior_total_means, prior_total_vars, pred_raster_mean, pred_raster_var,
			 file = paste0("results/predict/CA_DEzs_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_prior_predict.RData"))
	}
	
	detach("CAsetup", character.only = TRUE)
	
	
	#- 4.4) Process results
	
	# Validation results: 
	# reshape dataframes to long format
	post_pred_mean <- reshape(post_pred_mean, direction = "long", varying = list(2:23), timevar = "Year", idvar = "ID", v.names = "Abund")
	post_pred_var  <- reshape(post_pred_var , direction = "long", varying = list(2:23), timevar = "Year", idvar = "ID", v.names = "varAbund")
	
	if(identical(post_pred_mean$ID,post_pred_var$ID)){post_pred <- cbind(post_pred_mean,post_pred_var[3])
	}else{warning("Format error in prediction output")}
	
	pred_mhbcells <- merge.data.frame(MHBdata_l,post_pred,by.x = c("cells","Year"), by.y = c("ID","Year"))
	
	# Prediction results:
	# calculate sum of total abundance over all folds, then mean & var over posterior draws for each fold
	post_pred_total[,"all_fold",] <- apply(post_pred_total[,1:max(fold),], c(1,3), sum)
	pred_total_means <- apply(post_pred_total, c(2,3), mean)
	pred_total_vars <- apply(post_pred_total, c(2,3), var)
	# calculate final values of mean & variance for raster stacks
	pred_raster_mean <- pred_raster_mean/nr_postdraws
	pred_raster_var <- pred_raster_var/nr_postdraws - (pred_raster_mean)^2
	
	
	#- 4.5) Save all results to file
	save(pred_mhbcells, post_pred_total, pred_total_means, pred_total_vars, pred_raster_mean, pred_raster_var,
		 file = paste0("results/predict/CA_DEzs_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_predict.RData"))
}else{
	# load predictions from file
	if(prior_predictions){
		load(paste0("results/predict/CA_DEzs_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_prior_predict.RData"))
	}else{
		load(paste0("results/predict/CA_DEzs_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_predict.RData"))
	}
}


###--- 5.) Validation

	pred_mhbcells$sdAbund <- sqrt(pred_mhbcells$varAbund)
	pred_mhbcells$uppAbund <- pred_mhbcells$Abund + pred_mhbcells$sdAbund
	pred_mhbcells$lowAbund <- pred_mhbcells$Abund - pred_mhbcells$sdAbund
	pred_mhbcells$lowAbund[pred_mhbcells$lowAbund<0] <- 0

#- 5.1) Validate all MHB samples (all site-year)

#- 5.1.a) with ROC & AUC
	
	require(pROC)
	pred_occ <- data.frame(pred_mhbcells, 
						   mhb_occ = as.numeric(pred_mhbcells$Count>0), 
						   pred_occ = (1-dpois(0, lambda = pred_mhbcells$Abund)))
	(occ_roc <- roc(mhb_occ~pred_occ, data = pred_occ, na.rm = TRUE, smooth = FALSE, ci = TRUE, ret = c("roc", "coords", "all_coords"), plot = TRUE ))
	
	occ_roc_df <- data.frame(spec = occ_roc$specificities, sens = occ_roc$sensitivities, thres = occ_roc$thresholds)
	
	# repeat this but with only those sellc that change occupancy
	MHBocc <- cbind(MHBdata[,1], ifelse(MHBdata[,-1]>0,1,0) )
	cell_chg <- MHBocc[apply(MHBocc[,-1],1,var,na.rm=T)>0,1]
	pred_occ_chg <- pred_occ[pred_occ$cells%in%cell_chg,]
	(occ_roc <- roc(mhb_occ~pred_occ, data = pred_occ, na.rm = TRUE, smooth = FALSE, ci = TRUE, ret = c("roc", "coords", "all_coords"), plot = TRUE ))
	
	# make plot
	filename_ext <- "AUCallFold"
	if(latex_device){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_",filename_ext,".tex"), width = 2.52, height = 1.9, standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_",filename_ext,".pdf"),width = 6, height = 4) }
	plot(occ_roc)
	dev.off()
	
	# AUC for individual folds
	auc_perFold <- sapply(1:max(fold), function(foldID){
		fold_cells <- lookup[lookup$foldID == foldID,"cells"]
		roc(mhb_occ~pred_occ, data = pred_occ[pred_occ$cells%in%fold_cells,], na.rm = TRUE, smooth = FALSE, ci = TRUE, plot = FALSE, ret = c("roc", "coords"))
	})
	
	auc_perFold_stats <- sapply(1:max(fold), function(foldID){ as.numeric(auc_perFold['ci',foldID]$ci) })
	
	auc_perFold_df <- data.frame(factor(paste("Fold",1:5)), t(auc_perFold_stats) )
	names(auc_perFold_df) <- c("Fold ID","lower bound","AUC","upper bound")
	

#- 5.1.b) with C-index
	
# C-index for individual folds
	cix_perFold <- sapply(1:max(fold), function(foldID){
		fold_cells <- lookup[lookup$foldID == foldID,"cells"]
		with(pred_mhbcells[pred_mhbcells$cells%in%fold_cells,], rcorr.cens(Abund, Count))
	})
	cix_perFold_stats <- cix_perFold[c("C Index","S.D."),]
	
	# add a column for all folds combined
	cix_perFold_stats <- cbind(with(pred_mhbcells,rcorr.cens(Abund, Count))[c("C Index","S.D.")],cix_perFold_stats)
	
	# make dataframe
	cix_perFold_df <- data.frame(factor(c("all",paste("Fold",1:5))), t(cix_perFold_stats) )
	names(cix_perFold_df) <- c("Fold ID","C index","2 SD")
	

#- 5.2. make tables and plots
	
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
	
	print(xtable(cix_perFold_df,
				 type = "latex",
				 include.rownames = FALSE,
				 align = "cccc",
				 floating = FALSE),
		  file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_CIndextable.tex"))
	
	# ggplot
	names(auc_perFold_df) <- c("Fold_ID","lower_bound","AUC","upper_bound")
	auc_perFold_df$Fold_ID <- factor(auc_perFold_df$Fold_ID,levels = rev(levels(auc_perFold_df$Fold_ID)),ordered = TRUE)
	auc_ggplot <- 
		ggplot(auc_perFold_df,
			   aes(
			   	x = AUC,
			   	y = Fold_ID # reorder(Fold_ID, 1:5),
			   )) +
		geom_point() +
		geom_errorbarh(aes(xmax = upper_bound, xmin = lower_bound, height = .2)) +
		xlim(NA, 1) +
		theme_bw() +
		theme(plot.title = element_text(hjust = 0.5),
			  panel.grid.major.x = element_blank(),
			  panel.grid.minor.x = element_blank() ) +
		#coord_flip() + 
		labs(x = "AUC", y="")
			 #title="AUC for all spatial folds")
	
	names(cix_perFold_df) <- c("Fold_ID","C_Index","SD")
	cix_perFold_df$Fold_ID <- factor(cix_perFold_df$Fold_ID,levels = rev(levels(cix_perFold_df$Fold_ID)),ordered = TRUE)
	cix_ggplot <- 
		ggplot(cix_perFold_df,
			   aes(
			   	x = C_Index,
			   	y = Fold_ID # reorder(Fold_ID, 1:5),
			   )) +
		geom_point() +
		geom_errorbarh(aes(xmax = C_Index+SD, xmin = C_Index-SD, height = .2)) +
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
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_",filename_ext,".pdf"),width = 6, height = 4) }
	plot(auc_ggplot)
	dev.off()
	filename_ext <- "CIxperFold"
	if(latex_device){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_",filename_ext,".tex"), width = 2.52, height = 1.9, standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_",filename_ext,".pdf"),width = 6, height = 4) }
	plot(cix_ggplot)
	dev.off()
	
	
#- 5.2) Validate total abundances
	
	# calculate total MHB abundance over all folds
	mhb_totals <- sapply(1:max(fold), function(foldID){
		fold_cells <- lookup[lookup$foldID == foldID,"cells"]
		colSums(MHBdata[MHBdata$cells%in%fold_cells,-1], na.rm=T)
	})
	# add column for combined folds and transpose
	mhb_totals <- t(cbind(mhb_totals,colSums(MHBdata[,-1], na.rm=T)))
	
	# C-index for times series of total abundance within folds
	cindex_perFold <- sapply(1:nrow(pred_total_means), function(foldID){
		rcorr.cens(pred_total_means[foldID,1:22], mhb_totals[foldID,])
	})
	
	# make dataframe
	#cindex_perFold_df <- data.frame(factor(paste("Fold",1:5)), cindex_perFold["C Index",]-cindex_perFold["S.D.",] , cindex_perFold["C Index",] , cindex_perFold["C Index",]+cindex_perFold["S.D.",] )
	#names(cindex_perFold_df) <- c("Fold ID","lower bound","C-index","upper bound")
	cix_perFold_df <- data.frame(factor(c("all",paste("Fold",1:5))), t(cindex_perFold[c("C Index","S.D."),]) )
	names(cix_perFold_df) <- c("Fold ID","C index","2 SD")
	
	# export latex table
	options(xtable.floating = FALSE)
	options(xtable.include.rownames = FALSE)
	options(xtable.booktabs = TRUE)
	print(xtable(cix_perFold_df,
				 type = "latex",
				 include.rownames = FALSE,
				 align = "cccc",
				 floating = FALSE),
		  file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_CItotalsTable.tex"))
	
	# ggplot 
	names(cix_perFold_df) <- c("Fold_ID","C_Index","SD")
	cix_perFold_df$Fold_ID <- factor(cix_perFold_df$Fold_ID,levels = rev(levels(cix_perFold_df$Fold_ID)),ordered = TRUE)
	cix_tots_ggplot <- 
		ggplot(cix_perFold_df,
			   aes(
			   	x = C_Index,
			   	y = Fold_ID # reorder(Fold_ID, 1:5),
			   )) +
		geom_point() +
		geom_errorbarh(aes(xmax = C_Index+SD, xmin = C_Index-SD, height = .2)) +
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
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_",filename_ext,".pdf"),width = 6, height = 4) }
	plot(cix_tots_ggplot)
	dev.off()
	
	
#- 5.3) Validate cell with highest variance in abundance over time
	
	# AUC for calibration with full data set 
	
	# calculate count variance in each cell over time and sort
	cell_var <- cbind(MHBdata[,1],apply(MHBdata[,-1], 1, var, na.rm = T))
	cell_var <- cell_var[order(cell_var[,2], decreasing = TRUE),]
	
	# select the 15% highest variance cells
	#high_var_cell <- MHBdata$cells[cell_var>0.3]
	high_var_cell <- cell_var[1:as.integer(length(cell_var)*0.15),1]
	
	# calculate AUC & C-index
	(roc(mhb_occ~pred_occ, data = pred_occ[pred_occ$cells%in%high_var_cell,], na.rm = TRUE, smooth = FALSE, ci = TRUE, ret = c("roc", "coords", "all_coords"), plot = TRUE))

	
	# C-index for individual folds
	
	cix_perFold <- sapply(0:max(fold), function(foldID){
		if(foldID) fold_cells <- lookup[lookup$foldID == foldID,"cells"]
		else       fold_cells <- lookup[,"cells"]
		MHBsubset <- MHBdata[MHBdata[,1] %in% fold_cells, ]
		cell_var <- cbind(MHBsubset[,1],apply(MHBsubset[,-1], 1, var, na.rm = T))
		cell_var <- cell_var[order(cell_var[,2], decreasing = TRUE),]
		high_var_cell <- cell_var[1:as.integer(length(cell_var)*0.15),1]
		with(pred_mhbcells[pred_mhbcells$cells%in%high_var_cell,], rcorr.cens(Abund, Count))
	})
	cix_perFold_stats <- cix_perFold[c("C Index","S.D."),]
	
	# make dataframe
	cix_perFold_df <- data.frame(factor(c("all",paste("Fold",1:5))), t(cix_perFold_stats) )
	names(cix_perFold_df) <- c("Fold ID","C index","2 SD")
	
	# export latex table
	options(xtable.floating = FALSE)
	options(xtable.include.rownames = FALSE)
	options(xtable.booktabs = TRUE)
	print(xtable(cix_perFold_df,
				 type = "latex",
				 include.rownames = FALSE,
				 align = "cccc",
				 floating = FALSE),
		  file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_CIhighvarTable.tex"))
	
	# ggplot
	names(cix_perFold_df) <- c("Fold_ID","C_Index","SD")
	cix_perFold_df$Fold_ID <- factor(cix_perFold_df$Fold_ID,levels = rev(levels(cix_perFold_df$Fold_ID)),ordered = TRUE)
	cix_highvar_ggplot <- 
		ggplot(cix_perFold_df,
			   aes(
			   	x = C_Index,
			   	y = Fold_ID # reorder(Fold_ID, 1:5),
			   )) +
		geom_point() +
		geom_errorbarh(aes(xmax = C_Index+SD, xmin = C_Index-SD, height = .2)) +
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
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_",filename_ext,".pdf"),width = 6, height = 4) }
	plot(cix_highvar_ggplot)
	dev.off()


###--- 6.) Predictions

#- 6.1) Maps

	# sel_years = c(2,12,22,32,42,52)
	# spplot(pred_raster_mean[[sel_years]])
	# spplot(pred_raster_var[[sel_years]])
	# 
	# # get standard deviation from variance
	# pred_raster_sd <- sqrt(pred_raster_var)
	# pred_raster_ratio <- pred_raster_sd/pred_raster_mean
	# values(pred_raster_ratio)[values(pred_raster_sd)==0] <- 0
	# spplot(pred_raster_ratio[[sel_years]], zlim = c(0,2.0))
	# 
	# hist(pred_raster_mean, c(1,31), xlim = c(0,2.5), nclass = 20)
	# hist(pred_raster_sd, c(1,31), xlim = c(0,1.5), nclass = 20)
	# hist(pred_raster_ratio, c(1,31), xlim = c(0,2.5), nclass = 20)
	# 
	# diff_r <- pred_raster_mean[[31]]-pred_raster_mean[[1]]
	# values(diff_r)[abs(values(diff_r))<0.1] <- 0
	# spplot(diff_r)
	# 
	# diff_r2 <- pred_raster[[1]]
	# values(diff_r2)[values(diff_r)<0] <- -1
	# values(diff_r2)[values(diff_r)>0.1] <- 1
	# spplot(diff_r2)
	# 
	# thresh <- 1.2
	# diff_r <- pred_raster[[1]]
	# values(diff_r)[values(pred_raster_mean[[31]])>thresh] <- 2
	# values(diff_r)[values(pred_raster_mean[[1]]) >thresh] <- 1
	# spplot(diff_r)
	
	# put together the maps to plot
	plot_years <- c(2,22,42)
	plot_stack <- pred_raster_mean[[plot_years]]         # Years 2000, 2010, 2020,
	plot_stack[[4]] <- plot_stack[[2]]-plot_stack[[1]]   # Difference 2000-2020
	plot_stack[[5]] <- plot_stack[[3]]-plot_stack[[2]]   # Difference 2020-2040
	names(plot_stack)[c(4:5)] <- c("diff_2000_2020","diff_2020_2040")
	spplot(plot_stack)
	
	# crop to Swiss borders
	bg <- raster('/vsicurl/https://damariszurell.github.io/SDM-Intro/CH_mask.tif')
	bg <- raster::crop(bg,plot_stack)
	bg <- aggregate(bg,2,sum)
	plot_stack <- raster::mask(plot_stack, bg)
	
	# multiply by 25 to get units of breeding pairs per 100km²
	plot_stack <- 25*plot_stack
	
	# settings
	names(plot_stack) <- c("Year 2000","Year 2020","Year 2040","Difference 2000-2020","Difference 2020-2040")
	my.settings <- list(
		strip.background = list(col="gray90")#,strip.border=list(col="transparent")
	)
	col_palette_grain <- 7
	
	# plot absolute abundances
	max_abs <- max(abs(values(plot_stack[[c(1:3)]])),na.rm=T)
	brk <- do.breaks(c(0, max_abs), col_palette_grain)
	# Colors
	cols <- colorRampPalette(c("#E2E2E2",makeColDarker(cols_qualitative[3],factor = 0.25)))(col_palette_grain)
	#cols <- rev(hcl.colors(col_palette_grain, "Blues 2")) # "BrwnYl"
	filename_ext <- "abund_maps_abs"
	if(latex_device){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_",filename_ext,".tex"), width = 2.8, height = 4.8, standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_",filename_ext,".pdf"),
		width = 3.9, height = 6.7) }
	spplot(plot_stack[[c(1:3)]], col.regions = cols, at = brk, colorkey = list(space = "left", width = .75), 
		   layout = c(1,3), par.settings = my.settings, par.strip.text = list(cex=.8, lines=.8))
	dev.off()
	
	# plot abundance differences
	max_abs <- max(abs(values(plot_stack[[c(4:5)]])),na.rm=T)
	brk <- do.breaks(c(-max_abs, max_abs), (col_palette_grain*2-1))
	first_true <- which.max(brk > min(values(plot_stack[[c(2,4)]]),na.rm=T))
	brk <- brk[(first_true -1):length(brk)]
	brk <- brk[1:col_palette_grain+1]
	cols <- colorRampPalette(c("#E2E2E2",makeColDarker(cols_qualitative[8],factor = 0.25)))(col_palette_grain+1)
	#cols <- hcl.colors((col_palette_grain*2-1), "Blue-red")
	#cols <- cols[(first_true -1):length(cols)]
	#cols <- rev(cols)
	filename_ext <- "abund_maps_diff"
	if(latex_device){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_",filename_ext,".tex"), width = 2.8, height = 3.45, standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_",filename_ext,".pdf"),
		width = 3.9, height = 4.8) }
	spplot(plot_stack[[c(4:5)]], col.regions = cols, at = brk, colorkey = list(space = "right", width = .75),
		   layout = c(1,2), par.settings = my.settings, par.strip.text = list(cex=.8, lines=.8))
	dev.off()
	
	

#- 6.2) Time series
	
	# get quantiles to plot
	if(!prior_predictions){
		total_abund <- data.frame(years = 1999:2050,
								  t( sapply(seq(dim(post_pred_total)[3]), function(yr) {
								  	quantile(post_pred_total[,"all_fold",yr], probs = c(0.025, 0.500, 0.975), na.rm = FALSE)}))
		)	
		# earlier version: plotting mean and SD
		#total_abund <- data.frame(years = 1999:2050,
		#						  pred_mean  = pred_total_means["all_fold",], 
		#						  pred_upper = pred_total_means["all_fold",] + 2*sqrt(pred_total_vars["all_fold",]),
		#						  pred_lower = pred_total_means["all_fold",] - 2*sqrt(pred_total_vars["all_fold",]))
		
	}else{
		total_abund <- data.frame(years = 1999:2050,
								  t(apply(prior_pred_total, 2, quantile, probs = c(0.025, 0.500, 0.975), na.rm = FALSE))
		)
	}
	names(total_abund)[2:4] <- c("pred_lower", "pred_mean", "pred_upper")

	# add BBI data
	load(file = "../Observed_data/MHB_data/breeding_bird_index.rds")
	total_abund$mhb_rel <- c(bbi_RK$Model, rep(NA,31))
	
	# calculate relative abundances
	total_abund$pred_mean_rel <- total_abund$pred_mean / total_abund$pred_mean[1]
	total_abund$pred_upper_rel <- total_abund$pred_upper / total_abund$pred_mean[1]
	total_abund$pred_lower_rel <- total_abund$pred_lower / total_abund$pred_mean[1]
	
	# Colors
	colors = sapply(cols_qualitative,makeColDarker,factor = 0.25)
	pred_col <- colors[3]
	mhb_col <- colors[8]
	
	# plot
	filename_ext <- ifelse(prior_predictions,"PriorPredTotalTS","PredTotalTS")
	if(latex_device){ tikz(file = paste0(path_doc,"CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_",filename_ext,".tex"), width = 3.8, height =  ifelse(prior_predictions,1.7,2.2), standAlone = FALSE)
	}else{ pdf(file = paste0(path_doc,"pdf/CA_",ifelse(simulatedData,"sim","mhb"),"_",samplr,"_agBlc",aggr_blocks,"_",min(batches[,fold%in%(1:max(fold))]),"-",max(batches[,fold%in%(1:max(fold))]),"_",filename_ext,".pdf"),width = 6, height = 4) }
		par(mar=c(3,3,0.2,3), mgp = c(2.0, 0.66, 0) )
		if(latex_device){ par(lwd = 1.7) }
		map_years <-  c(2000, 2020, 2040) # year index that marks end of time window for aggregation
		if(!prior_predictions){ 
			total_plot <- with(total_abund, {
				plot(NULL, xlim = range(years), ylim = range(c(pred_lower_rel,pred_upper_rel)),
					 xlab = "Year", ylab="Posterior rel. pred. ab.", main = ifelse(latex_device,"","Total relative abundance Switzerland"))
				polygon(c(years,rev(years)),c(pred_lower_rel,rev(pred_upper_rel)), col = t_col(pred_col,70), border = NA )
				for(year in map_years) abline(v=year ,lty = 3, lwd = 3 , col = colors[10])
				lines(years,pred_lower_rel, type = "l", col = t_col(pred_col,50), lwd = 1)
				lines(years,pred_upper_rel, type = "l", col = t_col(pred_col,50), lwd = 1)
				lines(years,pred_mean_rel,  type = "p", col = pred_col, pch = 16, cex = 0.5, lwd = 2)
				lines(years,mhb_rel, type = "b", col = mhb_col, lwd = 2)
				par(new = TRUE) # Add new plot
				#plot(years,mhb_rel, type = "b", col = "red", lwd = 2, axes = FALSE, xlab = "", ylab = "") # Create second plot without axes
				plot(NULL, axes = FALSE, xlab = "", ylab = "", 
					 xlim = range(years), ylim = range(c(total_abund$pred_lower,total_abund$pred_upper))) # Create second plot without axes
				axis(side = 4, at = c(2000,4500,7000)) # second y-axis
				mtext("Absolute pred. ab.", side = 4, line = 2)
				legend(x=2020, y=3300, legend=c("Observed","Simulated"), col=c( mhb_col, pred_col), pch=c(16,16), bty="n")
			})
		}else{
			par(mar=c(0.1,3,0.2,3), mgp = c(2.0, 0.66, 0) )
			total_plot <- with(total_abund, {
				plot(NULL, xlim = range(years), ylim = range(c(pred_lower_rel,pred_upper_rel)),
					 xaxt = "n", xlab = "", #"Year", 
					 ylab="Prior rel. pred. ab.", main = ifelse(latex_device,"","Total relative abundance Switzerland"))
				polygon(c(years,rev(years)),c(pred_lower_rel,rev(pred_upper_rel)), col = t_col(pred_col,70), border = NA )
				for(year in map_years) abline(v=year ,lty = 3, lwd = 3 , col = colors[10])
				lines(years,pred_lower_rel, type = "l", col = t_col(pred_col,50), lwd = 1)
				lines(years,pred_upper_rel, type = "l", col = t_col(pred_col,50), lwd = 1)
				lines(years,pred_mean_rel,  type = "p", col = pred_col, pch = 16, cex = 0.5, lwd = 2)
				par(new = TRUE) # Add new plot
				plot(NULL, axes = FALSE, xlab = "", ylab = "", 
					 xlim = range(years), ylim = range(c(total_abund$pred_lower,total_abund$pred_upper))) # Create second plot without axes
				axis(side = 4, at = c(0,2000,4500)) # second y-axis
				mtext("Absolute pred. ab.", side = 4, line = 2)
			})
		}
	dev.off()
