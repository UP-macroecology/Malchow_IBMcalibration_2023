
###-----------------------------------------------------------------
##   PLOT PRIOR AND POSTERIOR MARGINALS
###-----------------------------------------------------------------

###
### copied and modified from BayesianTools::marginalPlot() and BayesianTools::marginalPlotDensity()
###

# function to plot posterior and optionally prior marginal distributions from a sample list in different ways

myMarginalPlot <- function(samplist, coda_chain = NULL, prior = TRUE, xrange = NULL, nPriorDraws = 10000, start = 1, thin = "auto", col = NULL, singlePanel = FALSE, boxplots = FALSE, MAP_section = NULL,
						   mfrow = c(3,4), plot_MAP = TRUE, plot_ref = TRUE, par_names = NULL, ...){ # ...passed to getSample, e.g. numSamples (from posterior) or start (to discard burn-in)
	
	# sample from posterior
	if(is.null(coda_chain)) posteriorMat <- getSample(samplist, parametersOnly = TRUE, start = start, thin = thin, ...)
	else posteriorMat <- getSample(coda_chain, parametersOnly = TRUE, start = start, thin = thin, ...)
	
	# checking for parameter 'which'
	args <- list(...)   
	if("which" %in% names(args))
		which = args$which
	else
		which = 1:ncol(posteriorMat)
	nPar <- ncol(posteriorMat)
	
	if(is.null(col)) col = c(t_col(cols_qualitative[3],percent = 30),t_col(cols_qualitative[10],percent = 30))
	
	if(prior){
		priorMat = sampler_list[[1]]$setup$prior$sampler(nPriorDraws) # draw prior from bayesianSetup
	}else priorMat = NULL
	if (!is.null(priorMat)) {
		priorMat = priorMat[,which]
		if (ncol(posteriorMat) != ncol(priorMat)) stop("wrong dimensions of prior")
		colnames(priorMat) <- colnames(posteriorMat)    
	}
	
	# check xrange
	if (!is.null(xrange)) {
		if (!any(c('numeric', 'matrix') %in% class(xrange))) stop('xrange must be numeric or matrix, or NULL')
		if ('numeric' %in% class(xrange)) xrange <- matrix(rep(xrange), nPar, nrow = 2)
		else if ('matrix' %in% class(xrange)) {
			if (ncol(xrange) != ncol(posteriorMat)) stop('xrange must have as many colums as there are parameters')
			else if (nrow(xrange) != 2) stop('xrange must have two rows (min, max)')
		}
	} else {
		posteriorRanges <- apply(posteriorMat, 2, range)
		priorRanges <- if(!is.null(priorMat)) apply(priorMat, 2, range) else NULL
		xrange <- if (is.null(priorRanges)) posteriorRanges else apply(rbind(priorRanges, posteriorRanges), 2, range)
	}
	
	# check parameter names
	if (is.null(colnames(posteriorMat))) colnames(posteriorMat) <- paste('par', 1:nPar, sep = '')
	if (!is.null(priorMat)) colnames(priorMat) <- colnames(posteriorMat)
	if (is.null(par_names)) parNames <- colnames(posteriorMat)
	else parNames <- par_names
	
	# create densities for plotting
	posteriorDensities <- lapply(1:ncol(posteriorMat),
								 function(i) density(posteriorMat[,i], from = xrange[1,i], to = xrange[2,i], ...))
	priorDensities <- if(!is.null(priorMat)){
		lapply(1:ncol(priorMat),function(i) density(priorMat[,i], from = xrange[1,i], to = xrange[2,i], ...))
	}else NULL
	postXY <- lapply(posteriorDensities, function(d){
		xy <- cbind(c(d$x[1], d$x, d$x[length(d$x)]),
					c(     0, d$y,                0))
		colnames(xy) <- c('x', 'y')
		xy
	})												 
	priorXY <- if (!is.null(priorDensities)) lapply(priorDensities, function(d) {
		xy <- cbind(c(d$x[1], d$x, d$x[length(d$x)]),
					c(0, d$y, 0))
		colnames(xy) <- c('x', 'y')
		xy
	}) else NULL
	
	if(boxplots) plot_MAP <- plot_ref <- FALSE
	
	# get MAP parameter values
	if (plot_MAP) {
		MAPpars <- round(MAP(samplist)$parametersMAP,3)
		MAPcol  <- cols_qualitative[7]
	}
	else MAPcol <- "white"
	
	# get reference parameter values
	if (plot_ref) {
		REFpars <-  round(sampler_list[[1]]$setup$info$plotBest,3)
		REFcol <- cols_qualitative[1] #"grey80"
	}
	else REFcol <- "white"
	
	# PLOTTING
	if (singlePanel) {
		
		op <- par(mfrow = c(nPar,1), mar = c(1, 5, 1, 1), oma = c(4, 1, 0.5, 1))
		on.exit(par(op))
		
		if (boxplots) {
			#post_quants <- apply(posteriorMat, 2, quantile, probs =  c(0.025, 0.25, 0.5, 0.75, 0.975))
			#prior_quants <- apply(priorMat   , 2, quantile, probs =  c(0.025, 0.25, 0.5, 0.75, 0.975))
			
			for (i in 1:nPar) {
				
				postY <- posteriorMat[,i]
				if(prior) {
					priorY <- priorMat[,i]
					priorX <- 1 # to add legend
				}else priorY <- NULL
				
				plot(NULL, NULL, xlim = xrange[,i], ylim = c(0,4), main = NA,
					 xlab = NA, ylab = NA, bty = 'n', yaxt = 'n', xaxt = 'n')
				
				boxplot(list(priorY,postY), horizontal = TRUE, ann = TRUE, add = TRUE, 
						#names = c("Prior","Posterior"),
						frame.plot = FALSE, # draw outer box?
						yaxt = 'n',         # draw y-axis?
						outline = FALSE,    # draw outliers?
						col = rev(col))
				
				mtext(sprintf('%20s', parNames[i]), 2, las = 1, adj = ifelse(latex_device,0,1), line = ifelse(latex_device,5,0), cex = ifelse(latex_device,0.7,1.0))
			}
			
		}else{ # !boxplots
			
			for (i in 1:length(posteriorDensities)) {
				postX <- postXY[[i]][,1]
				postY <- postXY[[i]][,2]
				
				priorX <- if (!is.null(priorXY[[i]])) priorXY[[i]][,1] else NULL
				priorY <- if (!is.null(priorXY[[i]])) priorXY[[i]][,2] else NULL
				
				yrange <- if (is.null(priorX)) range(postY) else range(c(postY, priorY))
				
				plot(NULL, NULL, xlim = xrange[,i], ylim = yrange, main = NA,
					 xlab = NA, ylab = NA, bty = 'n', yaxt = 'n', xaxt = 'n')
				axis(1, at = xrange[,i], labels = NA, lwd.ticks=0)
				xticks <- axTicks(1)
				xticks <- xticks[xticks >= xrange[1,i] & xticks <= xrange[2,i]]
				
				axis(1, at = xticks)
				
				mtext(sprintf('%20s', parNames[i]), 2, las = 1, adj = 1, cex = ifelse(latex_device,0.7,1.0))
				
				# plot MAP section
				if(!is.null(MAP_section)){
					mapsectX <- rep(c(MAP_section[i,'low'],MAP_section[i,'upp']),2)[c(1,2,4,3)]
					mapsectY <- rep(yrange,each=2)
					polygon(mapsectX, mapsectY, col = t_col(MAPcol,percent = 30), border = 1)
				}
				
				# plot posterior & prior
				if (!is.null(priorX)) polygon(priorX, priorY, col = col[2], border = 1)
				polygon(postX, postY, col = col[1], border = 1)
				
				if (plot_ref) lines(rep(REFpars[i],2), yrange, lwd = 3, col = REFcol, lty = 2)
				if (plot_MAP) lines(rep(MAPpars[i],2), yrange, lwd = 2, col = MAPcol)
				
			}
			#mtext('Marginal parameter uncertainty', outer = TRUE, cex = 1.5)
		}
		
		# overlay plot with empty plot to be able to place the legends freely
		par(fig = c(0, 1, 0, 1), oma = c(2,2,2,2), mar = c(0, 0, 0, 0), new = TRUE)
		plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
		mtext('Parameter values', 1, las = 1)
		leg.ix <- c(1)
		if (!is.null(priorX)) leg.ix <- c(leg.ix,2)
		if (plot_ref) leg.ix <- c(leg.ix,3)
		if (plot_MAP) leg.ix <- c(leg.ix,4)
		leg.text <- c('posterior', 'prior', 'reference value', 'MAP value')
		legend('topright', leg.text[leg.ix], xpd = TRUE, horiz = FALSE, inset = c(0.0, ifelse(boxplots,0.01,0.0)), lty = c(NA,NA,2,1)[leg.ix], lwd = c(NA,NA,3,3)[leg.ix],
			   bty = 'o', bg = "white", pch = c(15,15,NA,NA)[leg.ix], col = c(col,REFcol,MAPcol)[leg.ix], cex = 1.3, pt.cex = 2.5) #pch = 15
		
	} else { # !singlePanel
		
		op <- par(mfrow = mfrow, mar = c(2.5,2.5,2,0.3), mgp=c(2.4,1,0), oma=c(0,1,0,0), xpd=NA)
		on.exit(par(op))
		
		for (i in 1:length(posteriorDensities)) {
			postX <- postXY[[i]][,1]
			postY <- postXY[[i]][,2]
			
			priorX <- if (!is.null(priorXY[[i]])) priorXY[[i]][,1] else NULL
			priorY <- if (!is.null(priorXY[[i]])) priorXY[[i]][,2] else NULL
			
			yrange <- if (is.null(priorX)) range(postY) else range(c(postY, priorY))
			
			if( i %% mfrow[2] == 1 ) ylab = 'density' else  ylab = ''
			plot(NULL, NULL, xlim = xrange[,i], ylim = yrange, main = parNames[i], font.main = 1, xlab = NA, ylab = ylab)
			
			# plot MAP section
			if(!is.null(MAP_section)){
				mapsectX <- rep(c(MAP_section[i,'low'],MAP_section[i,'upp']),2)[c(1,2,4,3)]
				mapsectY <- rep(yrange,each=2)
				polygon(mapsectX, mapsectY, col = t_col(cols_qualitative[7],percent = 30), border = 1)
			}
			
			# plot posterior & prior
			if (!is.null(priorX)) polygon(priorX, priorY, col = col[2], border = 1)
			polygon(postX, postY, col = col[1], border = 1)
			
			if (plot_ref) lines(rep(REFpars[i],2), yrange, lwd = 3, col = REFcol, lty = 2)
			if (plot_MAP) lines(rep(MAPpars[i],2), yrange, lwd = 2, col = MAPcol)
			
			#if (i %% 16 == 1) mtext('Marginal parameter uncertainty', outer = TRUE, cex = 1.5)
		}
		
		# overlay plot with empty plot to be able to place the legends freely
		par(fig = c(0, 1, 0, 1), oma = c(2,2,2,2), mar = c(0, 0, 0, 0), new = TRUE)
		plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
		leg.ix <- c(1)
		if (!is.null(priorX)) leg.ix <- c(leg.ix,2)
		if (plot_ref) leg.ix <- c(leg.ix,3)
		if (plot_MAP) leg.ix <- c(leg.ix,4)
		leg.text <- c('posterior', 'prior', 'reference value', 'MAP value')
		legend('bottomright', leg.text[leg.ix], xpd = TRUE, horiz = FALSE, inset = c(0, 0), lty = c(NA,NA,2,1)[leg.ix], lwd = c(NA,NA,3,3)[leg.ix],
			   bty = 'n', pch = c(15,15,NA,NA)[leg.ix], col = c(col,REFcol,MAPcol)[leg.ix], cex = 1.3, pt.cex = 2.5) #pch = 15
		
	} # end singlePanel
}


###-----------------------------------------------------------------
##   PLOT PRIOR AND MULTIPLE POSTERIOR MARGINALS AS BOXPLOTS
###-----------------------------------------------------------------

# function to plot multiple posterior and the prior marginal distributions from a sample list as boxplots

multiMarginalBoxplot <- function(samples, xrange = NULL, col = NULL, mfrow = c(3,4), par_names = NULL, ...){ # ...passed to boxplot()
	
	if(is.null(col)) col = c(t_col(cols_qualitative[3],percent = 30),t_col(cols_qualitative[10],percent = 30))
	
	nr_cols <- table(unlist(lapply(samples, ncol)))
	if(length(nr_cols)!=1) stop('all samples must have the same number of parameters / columns') 
	else nPar <- as.integer(names(nr_cols[1]))
	
	# check xrange
	if (!is.null(xrange)) {
		if (!any(c('numeric', 'matrix') %in% class(xrange))) stop('xrange must be numeric or matrix, or NULL')
		if ('numeric' %in% class(xrange)) xrange <- matrix(rep(xrange), nPar, nrow = 2)
		else if ('matrix' %in% class(xrange)) {
			if (ncol(xrange) != nPar) stop('xrange must have as many colums as there are parameters')
			else if (nrow(xrange) != 2) stop('xrange must have two rows (min, max)')
		}
	} else {
		xrange <- apply (matrix( apply( sapply(samples,apply,2,range), 1,range), ncol=nPar) ,2,range)
	}
	
	# check parameter names
	if(is.null(par_names)) {
		colnames <- lapply(samples, colnames)
		colnames_ix <- which(sapply(colnames,function(x){!is.null(x)}))
		if(length(colnames_ix)) parNames <- colnames[[colnames_ix[1]]] else parNames <- paste('par', 1:nPar, sep = '')
	}else{ parNames <- par_names }
	
	
	# PLOTTING
	
	op <- par(mfrow = c(nPar,1), mgp=c(2,0.5,0), mar = c(1, 5, .8, 0.1), oma = c(2, 1, 0.1, 0.2))
	on.exit(par(op))

	#post_quants <- apply(posteriorMat, 2, quantile, probs =  c(0.025, 0.25, 0.5, 0.75, 0.975))
	#prior_quants <- apply(priorMat   , 2, quantile, probs =  c(0.025, 0.25, 0.5, 0.75, 0.975))
	
	for (p in 1:nPar) {
		
		par_samplist <- lapply(samples, "[", i=TRUE, j=p)
		nBoxes <- length(par_samplist)
		
		plot(NULL, NULL, xlim = c(ifelse(p%in%c(1,2),0,xrange[1,p]),xrange[2,p]), ylim = c(0,nBoxes*1.1), main = NA,
			 xlab = NA, ylab = NA, bty = 'n', yaxt = 'n', xaxt = 'n')
		
		boxplot(par_samplist, horizontal = TRUE, ann = TRUE, add = TRUE, 
				#names = c("Prior","Posterior"),
				frame.plot = FALSE, # draw outer box?
				yaxt = 'n',         # draw y-axis?
				outline = FALSE,    # draw outliers?
				col = c(col[2],rep(col[1],(nBoxes-1))))
		
		mtext(sprintf('%20s', parNames[p]), 2, las = 1, adj = ifelse(latex_device,0,1), line = ifelse(latex_device,5,0), cex = ifelse(latex_device,1.0,1.0))
	}

	# overlay plot with empty plot to be able to place the legends freely
	par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
	plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
	mtext('Parameter values', 1, line = -1)
	leg.text <- c('posterior', 'prior')
	legend(0.7,ifelse(latex_device,1.078,1.07), leg.text, xpd = TRUE, horiz = FALSE, inset = c(0.0,0.0), pch = c(15,15), col = col,
		   bty = 'o', bg = "white", cex = 1.3, pt.cex = 2.5,y.intersp=1.0,x.intersp=1.7) #pch = 15
}
