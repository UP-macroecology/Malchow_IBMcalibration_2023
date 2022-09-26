
###-----------------------------------------------------------------
##   GET TEMPORAL AGGREGATION AND RELATIVE ABUNDANCE
###-----------------------------------------------------------------

#
# The function sumBlocks() sums up the abundances in each spatial block (spatial aggregation) 
# and then calculates the mean within each temporal block (spatial aggregation).
#


sumBlocks <- function(input_df, timefolds = list(tf1 = 4:7, tf2 = 8:11, tf3 = 12:15, tf4 = 16:19, tf5 = 20:23) ){
	# 4 folds:  timefolds = list(tf1 = 4:8, tf2 = 9:13, tf3 = 14:18, tf4 = 19:23)
	# 5 folds:  timefolds = list(tf1 = 4:7, tf2 = 8:11, tf3 = 12:15, tf4 = 16:19, tf5 = 20:23)
	# 7 folds:  timefolds = list(tf1 = 3:5, tf2 = 6:8,  tf3 = 9:11,  tf4 = 12:14, tf5 = 15:17, tf6 = 18:20, tf7 = 21:23)
	
	# sum abundances in each spatial block
	input_df <- aggregate(x = input_df[-1], by = list(blockID = input_df$blockID), FUN = sum, na.rm = TRUE)
	
	# now, take means over years of each timefold, within each spatial block
	return( 
		data.frame( input_df[1],  # block IDs
					lapply(timefolds, FUN = function(tf){rowMeans(input_df[tf],na.rm = FALSE)})
	               )
	)
}
