
###-----------------------------------------------------------------
##   CREATE INITIAL INDIVIDUALS FILE
###-----------------------------------------------------------------


# The function createInitIndsFile() creates a list of individuals for initialisation and writes it to file.
# This file is needed for the initialisation option 'From initial individuals list file' (InitType = 2).

# The first argument of createInitIndsFile() is the desired size of the initial population given in #individuals.
# The distribution of ages and stages is based on the transition matrix and the maximum age. To get these, 
# createInitIndsFile() takes a StagesParams object as second parameter.
# Currently, the function is only defined for a 3-stage model with maximum age up to 50 yrs.

# The distribution of the individuals in space is dependent on 2 pre-defined data frames of x/y-coordinates,
# which are read from RDS.


createInitIndsFile <- function(pop_size, stage_struct, breedingpairs = TRUE, filename = "initinds.txt") {
	
	## check if stage_struct is of class "StagesParams"
	if(class(stage_struct)[1] != "StagesParams") {return(-1L)}
	
	## !! currently only defined for a 3-stage (4 RS stages) model!
	if(stage_struct@Stages != 4) {return(-2L)}
	
	## ! need a finite maximum age:
	if(stage_struct@MaxAge > 50) {return(-3L)}
	
	age_max <- stage_struct@MaxAge
	trans_mat <- stage_struct@TransMatrix
	
	### 1.) age- and stage- distributions
	
	age_probs <- sapply(0:(age_max-1), function(year){ # year = no. of transitions
		res <- rep(0,stage_struct@Stages-1)
		res[1] <- trans_mat[2,2]^year
		if(year>0){
			year1 <- 0:(year-1)
			year2 <- ((year-1)-year1)
			res[2] <- trans_mat[3,2]*sum(trans_mat[2,2]^year1*trans_mat[3,3]^year2) # year1 + year2 = year - 1
		}
		if(year>1){
			year1 <- rep(0:(year-2),((year-1):1))
			year2 <- NULL
			for (y in (year-2):0) {year2 <- c(year2, seq(0,y))}
			year3 <- ((year-2)-(year1+year2))
			res[3] <- trans_mat[3,2]*trans_mat[4,3]*sum(trans_mat[2,2]^year1*trans_mat[3,3]^year2*trans_mat[4,4]^year3) # year1 + year2 + year3 = year - 1
		}
		res
	})
	
	age_dist <- colSums(age_probs)/sum(colSums(age_probs)) # sum over stages and normalize
	stg_dist <- rowSums(age_probs)/sum(rowSums(age_probs)) # sum over ages and normalize
	
	if(breedingpairs){ # if pop_size is given in number of breeding pairs, adjust total pop. size according to stage dist.
		pop_size <- pop_size/stg_dist[3]
	}
	
	## apply distribution to total population:
	
	# define truncation fct
	trunc_prec = floor(log(pop_size, base=10))
	trunc <- function(x, ..., prec = trunc_prec) base::trunc(x * 10^prec, ...) / 10^prec;
	
	pop_ages <- as.integer( trunc(age_dist)*pop_size )
	# distribute the rest (add) randomly:
	add_ages <- table(sample(1:age_max, size = pop_size-sum(pop_ages), replace = TRUE, prob = age_dist ))
	pop_ages[as.integer(names(add_ages))] <- pop_ages[as.integer(names(add_ages))] + add_ages
	
	pop_stages <- as.integer( trunc(stg_dist)*pop_size )
	# distribute the rest (add) randomly:
	add_stages <- table(sample(1:(stage_struct@Stages-1), size = pop_size-sum(pop_stages), replace = TRUE, prob = stg_dist ))
	pop_stages[as.integer(names(add_stages))] <- pop_stages[as.integer(names(add_stages))] + add_stages
	
	
	age_vec <- rep(1:length(pop_ages), pop_ages)
	stg_vec <- rep(1:length(pop_stages), pop_stages)
	
	# check and correct for incompatible age-stage combinations:
	age_vec[age_vec-stg_vec<0] <- stg_vec[age_vec-stg_vec<0]
	
	# make groups per cell
	families <- table(stg_vec)
	# first version for family distribution:
	#nr_cells <- max(families) - nrow(mhb_init)
	#family_vec <- c(seq(families[1]),seq(families[2]),seq(families[3])) # needs adjustment if more stages allowed
	# new version that gives a cell to each individual of stage 2 and 3 and distributes the juveniles among the stg 3s
	family_vec <- c(seq(families[1]),seq((families[3]+1),(families[3]+families[2])),seq(families[3])) # needs adjustment if more stages allowed
	if(families[1]>families[3]){family_vec[(families[3]+1):families[1]] <- sample(1:families[3], size = (families[1]-families[3]), replace = TRUE) } # distribute 2nd juvenile per adult
	
	
	### 2.) initial location
	
	# load table of initial cells and their probabilities
	# now done in main script
	#hsi_vals <- readRDS("../Observed_data/MHB_data/initcell_probs_2km.rds")
	#mhb_init <- readRDS("../Observed_data/MHB_data/initcell_mhb99_2km.rds")
	
	# put first families into the constant cells from MHB of year 1999 and all further families draw a cell
	nr_cells <- families[["3"]] + families[["2"]] - nrow(mhb_init)
	hsivals_ix <- c(1:nrow(mhb_init), sample(seq(nrow(mhb_init)+1,nrow(mhb_init)+nrow(hsi_vals)), size = nr_cells, replace = TRUE, prob = hsi_vals$prob) )
	fam_cell_ix <- hsivals_ix[family_vec]
	
	# join the tables of potential cells, order them after family index, discard prob column
	hsi_vals <- rbind(mhb_init[,c('x','y')], hsi_vals[,c('x','y')])[fam_cell_ix,]
	
	
	### 3.) Create data frame and write to file
	
	# InitIndsFile column headers: 'Year', 'Species', 'X', 'Y', 'Ninds', 'Age', 'Stage'
	# note that x,y coordinates start counting at 0, and y-coordinate is reversed:
	ind_table <- data.frame(Year = 0, Species = 0, X = hsi_vals[,'x']-1, Y = 114-hsi_vals[,'y'], Ninds = 1, Age = age_vec, Stage = stg_vec)
	
	# write table ro file
	#write.table(ind_table, paste0(dirWorking,"Inputs/",filename) , row.names = FALSE, quote = FALSE, sep = "\t")
	
	# finish successfully
	#return(NULL)
	return(ind_table)
}



createInitIndsFile_Pois <- function(stage_struct, breedingpairs = TRUE) {
	
	## check if stage_struct is of class "StagesParams"
	if(class(stage_struct)[1] != "StagesParams") {return(-1L)}
	
	## !! currently only defined for a 3-stage (4 RS stages) model!
	if(stage_struct@Stages != 4) {return(-2L)}
	
	## ! need a finite maximum age:
	if(stage_struct@MaxAge > 50) {return(-3L)}
	
	age_max <- stage_struct@MaxAge
	trans_mat <- stage_struct@TransMatrix
	
	### 1.) age- and stage- distributions
	
	age_probs <- sapply(0:(age_max-1), function(year){ # year = no. of transitions
		res <- rep(0,stage_struct@Stages-1)
		res[1] <- trans_mat[2,2]^year
		if(year>0){
			year1 <- 0:(year-1)
			year2 <- ((year-1)-year1)
			res[2] <- trans_mat[3,2]*sum(trans_mat[2,2]^year1*trans_mat[3,3]^year2) # year1 + year2 = year - 1
		}
		if(year>1){
			year1 <- rep(0:(year-2),((year-1):1))
			year2 <- NULL
			for (y in (year-2):0) {year2 <- c(year2, seq(0,y))}
			year3 <- ((year-2)-(year1+year2))
			res[3] <- trans_mat[3,2]*trans_mat[4,3]*sum(trans_mat[2,2]^year1*trans_mat[3,3]^year2*trans_mat[4,4]^year3) # year1 + year2 + year3 = year - 1
		}
		res
	})
	
	age_dist <- colSums(age_probs)/sum(colSums(age_probs)) # sum over stages and normalize
	stg_dist <- rowSums(age_probs)/sum(rowSums(age_probs)) # sum over ages and normalize
	
	
	### 2.) draw number of initial individuals in each cell from Poisson distribution
	
	# duplicate rows so that one row for each (adult) individual
	adults <- pois_pred[ rep(1:nrow(pois_pred),rpois(nrow(pois_pred),lambda = pois_pred$Pois_pred)) , c("X","Y")]
	
	# total pop
	pop_size <- nrow(adults)
	
	if(breedingpairs){ # if pop_size is given in number of breeding pairs, adjust total pop. size according to stage dist.
		pop_size <- pop_size/stg_dist[3]
	}
	
	# if survival is very low, pop_size can get very large for breedingpairs=TRUE, as the expected proportion of adults becomes very small
	if(pop_size>1e5) warning("Warning: Large initial population size!")
	if(pop_size>1e6) {pop_size <- 1e6; warning("...setting init pop size to 1e6.",call. = FALSE)}
	
	### 3.) apply stage and age distributions to total population:
	
	# define truncation fct
	trunc_prec = floor(log(pop_size, base=10))
	trunc <- function(x, ..., prec = trunc_prec) base::trunc(x * 10^prec, ...) / 10^prec;
	
	pop_ages <- as.integer( trunc(age_dist)*pop_size )
	# distribute the rest (add) randomly:
	add_ages <- table(sample(1:age_max, size = pop_size-sum(pop_ages), replace = TRUE, prob = age_dist ))
	pop_ages[as.integer(names(add_ages))] <- pop_ages[as.integer(names(add_ages))] + add_ages
	
	pop_stages <- as.integer( trunc(stg_dist)*pop_size )
	# distribute the rest (add) randomly:
	add_stages <- table(sample(1:(stage_struct@Stages-1), size = pop_size-sum(pop_stages), replace = TRUE, prob = stg_dist ))
	pop_stages[as.integer(names(add_stages))] <- pop_stages[as.integer(names(add_stages))] + add_stages
	
	age_vec <- rep(1:length(pop_ages), pop_ages)
	stg_vec <- rep(1:length(pop_stages), pop_stages)
	
	# check and correct for incompatible age-stage combinations:
	age_vec[age_vec-stg_vec<0] <- stg_vec[age_vec-stg_vec<0]
	
	
	### 4.) Match individuals with cells
	
	# make sure condition from MHB year 1999 is represented in adults
	adults$Stage <- 3
	adults_mhb <- merge(adults, mhb_init, by.x = c("X","Y"), by.y = c("x","y"), all = TRUE)
	adultsmatch <- sum(!is.na(adults_mhb$Stage+adults_mhb$count)) # how many MHB cells are already occupied by Poisson draw?
	adultsrdm <- which(is.na(adults_mhb$count)) # indices of adults drawn from Poisson model that are non-MHB records
	adultsrdm <- adultsrdm[-sample(1:length(adultsrdm),(nrow(mhb_init)-adultsmatch))] # remove as many randomly drawn adults as are missing from MHB
	adultsmhb <- which(!is.na(adults_mhb$count)) # indices of MHB cells that are not in the set of adults drawn from Poisson model
	adultsmhb <- adultsmhb[which(!duplicated(adults_mhb[adultsmhb,1:2]))] # check for double-accounting
	adults <- adults_mhb[sample(c(adultsmhb,adultsrdm)),1:3] # merge adults from MHB-1999 records with those from Poisson model
	
	# now assign age to adults
	adults$Age <- age_vec[stg_vec==3][1:nrow(adults)]
	if(pop_stages[3]<nrow(adults)) adults$Age[is.na(adults$Age)] <- 6
	adults$Stage <- 3
	
	# then distribute juveniles among adults and assign age
	juveniles <- adults[sample(1:nrow(adults), pop_stages[1], replace = TRUE),c("X","Y")]
	juveniles$Age <- age_vec[stg_vec==1]
	juveniles$Stage <- 1
	
	# now distribute sub-adults among all cells and assign age
	subadults <- pois_pred[sample(1:nrow(pois_pred), pop_stages[2], replace = TRUE, prob = pois_pred$Pois_pred),c("X","Y")]
	subadults$Age <- age_vec[stg_vec==2]
	subadults$Stage <- 2
	
	# join together
	ind_table <- rbind(adults,subadults,juveniles)

	
	### 5.) Create data frame and write to file
	
	# InitIndsFile column headers: 'Year', 'Species', 'X', 'Y', 'Ninds', 'Age', 'Stage'
	ind_table <- data.frame(Year = 0, Species = 0, X = ind_table['X'], Y = ind_table['Y'], Ninds = 1, Age = ind_table$Age, Stage = ind_table$Stage)
	
	# write table ro file
	#write.table(ind_table, paste0(dirWorking,"Inputs/",filename) , row.names = FALSE, quote = FALSE, sep = "\t")
	
	# finish successfully
	#return(NULL)
	return(ind_table)
}
