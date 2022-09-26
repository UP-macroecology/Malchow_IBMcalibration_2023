
###-----------------------------------------------------------------
##   SETTINGS FOR DOCUMENTATION PLOTS
###-----------------------------------------------------------------

# load tikz device
library(tikzDevice)

# path for storing tikz figures
path_doc <- "doc/tex/figures/"


###-----------------------------------------------------------------
##   COLOR PALETTE
###-----------------------------------------------------------------

cols_sequential <- c('#00145e', '#00216f', '#002d80', '#0c3a8f', '#26459d', '#3851ab', '#485eb9', '#586ac7', '#6777d6', '#7684e5', '#8491f4', '#94a0fe', '#a5afff', '#b6bfff', '#c6cfff', '#d7dfff', '#e7efff')
cols_diverging  <- c('#00297b', '#1e4095', '#3c56b0', '#5b6dcb', '#7887dc', '#94a1ec', '#b2bcfc', '#d3d8fa', '#f5f5f5', '#f9d0d6', '#fea9b5', '#f08997', '#dc6a7a', '#c84b5d', '#ab3348', '#8d1a31', '#6e001a')
cols_qualitative <- c('#a0b869','#03a276','#00859a','#005faf','#072b8f','#003199','#9d1684','#d62d59','#e26f2c','#c9ad1e')


pal_sequential <- colorRampPalette(rev(cols_sequential), interpolate = "linear", space = "Lab")
pal_diverging  <- colorRampPalette(rev(cols_diverging), interpolate = "linear", space = "Lab")
pal_qualitative  <- colorRampPalette(rev(cols_qualitative), interpolate = "linear", space = "Lab")


filled.contour(volcano,
			   color = pal_qualitative,
			   asp = 1)




###-----------------------------------------------------------------
##   ADJUST COLOR
###-----------------------------------------------------------------

rgb_2_hex <- function(r,g,b){
	rgb(r, g, b, maxColorValue = 255, alpha = 255)
}

# col2rgb() # from package ’grDevices’. takes colorname or hexcode 


# make color darker 
makeColDarker <- function(color,factor=0.5,f_alpha=NULL){
	maxColorValue <- 255
	if(is.null(f_alpha)) f_alpha <- factor
	rgb <- col2rgb(color,alpha=TRUE)*c(rep(1-factor,3),1)
	rgb['alpha',] <- rgb['alpha',]+((maxColorValue-rgb['alpha',])*factor)
	return(rgb(rgb[1],rgb[2],rgb[3],rgb[4], maxColorValue = maxColorValue))
}


###-----------------------------------------------------------------
##   DEFINE COLOR TRANSPARENCY
###-----------------------------------------------------------------


## Transparent colors
## by Mark Gardener 2015
## www.dataanalytics.org.uk

# The function t_col() takes a R color name and a transparency value given in percent 
# and returns the RGB code of the corresponding transparent color

t_col <- function(color, percent = 50, name = NULL) {
	#      color = color name
	#    percent = % transparency
	#       name = an optional name for the color
	
	## Get RGB values for named color
	rgb.val <- col2rgb(color)
	
	## Make new color using input color as base and alpha set by transparency
	t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
				 max = 255,
				 alpha = (100 - percent) * 255 / 100,
				 names = name)
	
	## Save the color
	invisible(t.col)
}
## END