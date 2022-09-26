
###-----------------------------------------------------------------
##   RUN ONE-FACTOR-AT-A-TIME SENSITIVITY ANALYSIS
###-----------------------------------------------------------------


# The function sensitivityOAT() runs a one-factor-at-a-time sensitivity analysis with the
# RangeShiftR simulation specified in the parent environment.

# Each model parameter selected for calibration will be varied within its interval while
# all other model parameters are fixed at their respective reference values. The number of
# equi-distant points with the parameter interval (ll_pts) as well as the number of replicates
# at each point (ll_reps) can be set. Optionally, the results are plotted to a pdf file.


sensitivityOAT <- function(ll_points = 9, ll_replics = 9, savePlot = TRUE, minVal = -Inf){
	
	# get names and reference value to label plot
	par_names <- CA_params$name[par_sel]
	p_ref <- CA_params$def[par_sel]
	# set multiplot
	par(mfrow = c(3,4), mar = c(4,3,1,1))
	
	oneFactorSA <- lapply(seq(length(par_names)), function(factor_ix, ll_pts = ll_points, ll_reps = ll_replics ){ # ll_pts...nr of points from prior range, ll_reps...nr of replicates

		#reset reference parameters
		p_ref <- CA_params$def[par_sel]
		cat("parameter: ",as.character(par_names[factor_ix]),"\n")
		
		# make 1D grid for each factor
		par_ix <- par_sel[factor_ix]
		ll_out <- rep(0,ll_reps)
		ll_par <- seq(CA_params$min[par_ix],CA_params$max[par_ix],length.out = ll_pts)
		
		# calculate gof
		ll <- sapply(ll_par,
					 function(p){
					 	p_ref[factor_ix] <- p
					 	cat("values: ",signif(p_ref,digits=3),",\t rep: ")
					 	for(i in 1:ll_reps) {
					 		ll_out[i] <- CAtarget(p_ref)
					 		cat(i," ")
					 	}
					 	cat("\n")
					 	ll_out}
		)
		
		# plot
		if(TRUE){
		max_ll <- max(ll[ll>minVal])
		min_ll <- min(ll[ll>minVal])
		if(is.infinite(max_ll)) max_ll <- 0
		if(is.infinite(min_ll)) min_ll <- -1e3
		if(ll_reps>1){
			plot(colMeans(ll) ~ ll_par, type = "l", col = "red", lwd = 2, main = par_names[factor_ix], xlab = "", ylim = c(min_ll, max_ll))
			for(i in 1:ll_reps) points(ll[i,] ~ ll_par, col="grey50")
		}
		else{
			plot(ll ~ ll_par, type = "l", col = "red", lwd = 2, xlab = par_names[factor_ix], ylim = c(min_ll, max_ll))
			points(ll ~ ll_par, col="grey50")
		} 
		abline(v = p_ref[factor_ix])
		}
		
		# return gofs
		return(list(ll_par = ll_par, ll = ll))
	})
	
	if(savePlot){
		pdf(file = paste0("results/SA_out/SA_oat_",calib_name_core,"_",batchnr,".pdf"), width = 10, height = 7)
		# plot on equal scale
		temp_ll <- unlist(ifelse(length(oneFactorSA[-sigma_parix])==0,oneFactorSA,oneFactorSA[-sigma_parix]),recursive = F)
		temp_ll <- temp_ll[names(temp_ll)=="ll"]
		temp_ll <- unlist(temp_ll)
		temp_ll <- temp_ll[temp_ll>minVal]
		max_ll <- max(temp_ll)
		min_ll <- min(temp_ll)
		if(is.infinite(max_ll)) max_ll <- 0
		if(is.infinite(min_ll)) min_ll <- -1e3
		par(mfrow = c(3,4), mar = c(2,2,3,1))
		for (i in 1:length(oneFactorSA) ) {
			plot(colMeans(oneFactorSA[[i]]$ll) ~ oneFactorSA[[i]]$ll_par,
				 main = par_names[i], 
				 type = "l", col = "red", lwd = 2, 
				 ylim = c(min_ll, max_ll),#ifelse(rep(i==sigma_parix,2),c(min(oneFactorSA[[i]]$ll),max(oneFactorSA[[i]]$ll)),c(min_ll, max_ll)),
				 xlab = "",ylab = "")
			#plot(colMeans(oneFactorSA[[i]]$ll) ~ oneFactorSA[[i]]$ll_par, xlab = par_names[i], type = "l", col = "red", lwd = 2, ylim = c(min(oneFactorSA[[i]]$ll), max(oneFactorSA[[i]]$ll)))
			for(j in 1:ll_replics) points(oneFactorSA[[i]]$ll[j,] ~ oneFactorSA[[i]]$ll_par, col="grey50")
			abline(v = p_ref[i])
		}
		dev.off()
	}
	return(oneFactorSA)
}

# with output from BayesianTools SA

# par(mfrow = c(3,4), mar = c(4,3,1,1))
# for (i in 1:12 ) {
# 	plot(sen[[i]]$resp ~ sen[[i]]$par, xlab = par_names[i], type = "l", col = "red", lwd = 2, ylim = c(min(sen[[i]]$resp), max(sen[[i]]$resp)))
# 	abline(v = sen$reference[i])
# }

