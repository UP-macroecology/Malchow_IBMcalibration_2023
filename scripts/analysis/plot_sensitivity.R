library(tikzDevice)

source("scripts/analysis/docu_colors.R")


ll_replics <- 9
# get names and reference value to label plot
par_names <- c("Density-dependence $b^{-1}$","Fecundity $\\phi_0$","Survival prob. $\\sigma_1$","Survival prob. $\\sigma_2$","Survival prob. $\\sigma_3$","Developm. prob. $\\gamma_{1\\rightarrow 2}$","Developm. prob. $\\gamma_{2\\rightarrow 3}$","Emigration prob. $e_1$","Settlem. par. $\\beta_s$","Dispersion $\\nu$")
p_ref <- CA_params$def[par_sel]
# set multiplot


# plot on equal scale
temp_ll <- unlist(ifelse(length(oneFactorSA[-sigma_parix])==0,oneFactorSA,oneFactorSA[-sigma_parix]),recursive = F)
temp_ll <- temp_ll[names(temp_ll)=="ll"]
temp_ll <- unlist(temp_ll)
temp_ll <- temp_ll[temp_ll>minVal]
max_ll <- max(temp_ll)
min_ll <- min(temp_ll)
if(is.infinite(max_ll)) max_ll <- 0
if(is.infinite(min_ll)) min_ll <- -1e3

{
tikz(file = paste0(path_doc,"SA_oat_mhb.tex"), width = 5.76, height = 4.03, standAlone = FALSE)
	#par(mfrow = c(3,4), mar = c(2,2,3,1), oma = c(0,0,0,0))
	par(mfrow = c(3,4), mar = c(2.5,2.5,2,0.3), mgp=c(2.4,1,0), oma=c(0,1,0,0) , xpd=FALSE)
	for (i in 1:length(oneFactorSA) ) {
		#if(i %% 4 == 1) par(mar = c(2.5,2.7,2,0.3)) else par(mar = c(2.5,2.5,2,0.3))
		plot(NULL,
			 main = par_names[i], font.main = 1,
			 xlim = c(min(oneFactorSA[[i]]$ll_par),max(oneFactorSA[[i]]$ll_par)),
			 ylim = c(min_ll, max_ll),#ifelse(rep(i==sigma_parix,2),c(min(oneFactorSA[[i]]$ll),max(oneFactorSA[[i]]$ll)),c(min_ll, max_ll)),
			 xlab = "", ylab = '') # ifelse( i %% 4 == 1, 'likelihod', '') )
		#plot(colMeans(oneFactorSA[[i]]$ll) ~ oneFactorSA[[i]]$ll_par, xlab = par_names[i], type = "l", col = "red", lwd = 2, ylim = c(min(oneFactorSA[[i]]$ll), max(oneFactorSA[[i]]$ll)))
		lines(rep(p_ref[i],2), c(min_ll, max_ll), lwd = 3, col = cols_qualitative[1], lty = 2)
		for(j in 1:ll_replics) points(oneFactorSA[[i]]$ll[j,] ~ oneFactorSA[[i]]$ll_par, col=t_col(cols_qualitative[8],percent = 50))
		lines(colMeans(oneFactorSA[[i]]$ll) ~ oneFactorSA[[i]]$ll_par,type = "l", col = cols_qualitative[5], lwd = 2)
		#abline(v = p_ref[i], lwd = 3, col = cols_qualitative[1], lty = 2)
		if(i %% 4 == 1) mtext(text="likelihood",side=2,line=2.2,outer=FALSE,cex=0.75)
	}
	# overlay plot with empty plot to be able to place the legends freely
	par(fig = c(0, 1, 0, 1), oma = c(2,2,2,2), mar = c(0, 0, 0, 0), new = TRUE)
	plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
	leg.text <- c('mean', 'replicates', 'reference value')
	legend('bottomright', leg.text, xpd = TRUE, horiz = FALSE, inset = c(0, 0), lty = c(1,NA,2), lwd = c(2,NA,3),
		   bty = 'n', pch = c(NA,1,NA), col = c(cols_qualitative[5],cols_qualitative[8],cols_qualitative[1]), cex = 1.3, pt.cex = 1.5)
dev.off()
}



#######################################
### plot result of Morris screening 
#######################################
#
# copied and modified from BayesianTools::summary.mcmcSamplerList()
#
plot.morris <- function(x, identify = FALSE, atpen = FALSE,	y_col = NULL, y_dim3 = NULL, par_names = NULL, position = NULL, ...) {
	if (!is.null(x$ee)) {
		if(inherits(x$y, "numeric")){
			mu.star <- apply(x$ee, 2, function(x) mean(abs(x)))
			sigma <- apply(x$ee, 2, sd)
		} else if(inherits(x$y, "matrix")){
			if(is.null(y_col)) y_col <- 1
			if(!is.null(y_dim3)){
				warning("Argument \"y_dim3\" is ignored since the model output is ",
						"a matrix")
			}
			mu.star <- apply(x$ee[, , y_col, drop = FALSE], 2, function(x) mean(abs(x)))
			sigma <- apply(x$ee[, , y_col, drop = FALSE], 2, sd)
		} else if(inherits(x$y, "array")){
			if(is.null(y_col)) y_col <- 1
			if(is.null(y_dim3)) y_dim3 <- 1
			mu.star <- apply(x$ee[, , y_col, y_dim3, drop = FALSE], 2, function(x) mean(abs(x)))
			sigma <- apply(x$ee[, , y_col, y_dim3, drop = FALSE], 2, sd)
		}
		
		par(mar=c(4,4,0.4,0.4))
		plot(mu.star, sigma, pch = 20, xlab = "$\\mu^{\\star}$", ylab = "$\\sigma$", ...)
		
		if(is.null(par_names)) {labels = colnames(x$ee)} else {labels = par_names}
		if(is.null(position)) {position = 4} else {position = position}
		if (identify) {
			identify(mu.star, sigma, labels = labels, atpen = atpen)
		} else {
			text(mu.star, sigma, labels = labels, pos = position)
		}
	}
}

{
	tikz(file = paste0(path_doc,"SA_morris_mhb.tex"), width = 4.6, height = 3.2, standAlone = FALSE)
	plot.morris(morrisResult, par_names = par_names, position = c(4,1,4,2,2,4,4,4,4,4))
	dev.off()
}
