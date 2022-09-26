
###-----------------------------------------------------------------
##   SUMMARISE AND COMBINE INDEPENDENT MCMCS
###-----------------------------------------------------------------




### Function of retrieve info about MCMC and posterior summary
######################################################################
#
# copied and slightly modified from BayesianTools::summary.mcmcSamplerList()
#
evalMCMC <- function(samplist, quants = c(0.025, 0.15, 0.50, 0.85, 0.975) ){
	require(BayesianTools)
	MAPvals <- round(MAP(samplist)$parametersMAP,3)
	gelDiag <- gelmanDiagnostics(samplist)
	psf <- round(gelDiag$psrf[,1], 3)
	
	if(length(quants)!=5){
		warning("quants must have 5 elements; using default")
		quants = c(0.025, 0.15, 0.05, 0.85, 0.975)
	}
	info_sample <- getSample(samplist, parametersOnly = T, coda = T)
	npar <- ncol(info_sample[[1]])
	lowerext <- lowerq <- medi <- upperq <- upperext <- numeric(npar)
	for(i in 1:npar){
		tmp <- unlist(info_sample[,i])
		tmp <- quantile(tmp, probs = quants )
		lowerext[i] <- round(tmp[1],3)
		lowerq[i] <- round(tmp[2],3)
		medi[i] <- round(tmp[3],3)
		upperq[i] <- round(tmp[4],3)
		upperext[i] <- round(tmp[5],3)
	}
	refpars <- round(samplist[[1]]$setup$info$plotBest,3)
	
	posterior_df <- cbind(psf, MAPvals, lowerext, lowerq, medi, upperq,  upperext, refpars)
	colnames(posterior_df) <- c("psf", "MAP", paste0(as.character(quants),"%"), "prior_best")
	if(quants[3]==0.5){colnames(posterior_df)[5] <- "median"}
	row.names(posterior_df) <- colnames(info_sample[[1]])	
	
	sampler_info <- list(nrChain = length(info_sample),
						 nrIter = nrow(info_sample[[1]]),
						 conv = round(gelDiag$mpsrf,3),
						 npar = npar,
						 parnames = colnames(info_sample[[1]]))
					 
	return(list(posterior_df = posterior_df,
				sampler_info = sampler_info))
}




### Function to combine MCMCs
#######################################

# Read the calibration outputs for repeat runs and create combined MCMC sampler list

combineMCMCs <- function(sim_mhb,sampler,spat_agg,RSrep,chain_iter,batchnrs,fold=NULL, concatenateChains = FALSE){
	require(BayesianTools)
	
	# construct file names from info given via parameters
	filename <- paste0("results/CA_out/CA_DEzs_",ifelse(sim_mhb,"sim","mhb"),"_",sampler,"_agBlc",spat_agg,ifelse(is.null(fold)|fold==0,"",paste0("_spF",fold)),"_rep",RSrep,"_it",chain_iter,"_",batchnrs,".RData")
	if(concatenateChains) filename2 <- paste0("results/CA_out/CA_DEzs_",ifelse(sim_mhb,"sim","mhb"),"_",sampler,"_agBlc",spat_agg,ifelse(is.null(fold)|fold==0,"",paste0("_spF",fold)),"_rep",RSrep,"_it",chain_iter,"_",(batchnrs+100),".RData")

	# get object 'MCMCout' (that contains the separate chain) from each file and store in a list
	MCMC_chains <- list()
	for (sim in 1:length(batchnrs)) {
		name <- paste0("run_",batchnrs[sim])
		attach(filename[sim], name = name)
		mcmcChain <- get("MCMCout", name)
		detach(name, character.only = TRUE)
		if(concatenateChains){
			name <- paste0("run_",batchnrs[sim]+100)
			attach(filename2[sim], name = name)
			mcmcChain2 <- get("MCMCout", name)
			detach(name, character.only = TRUE)
			# concatenate chains
			mcmcChain1 <- getSample(mcmcChain, coda = T)
			mcmcChain2 <- getSample(mcmcChain2, coda = T)
			mcmcChain <-  lapply(1:length(mcmcChain1), function(c){mcmc(rbind(mcmcChain1[[c]],mcmcChain2[[c]]))})
		}
		MCMC_chains <- c(MCMC_chains,mcmcChain)
	}
	
	# combine individual chains and return combined MCMC
	if(concatenateChains){return(mcmc.list(MCMC_chains))
	}else{return(createMcmcSamplerList(MCMC_chains)) }
}
