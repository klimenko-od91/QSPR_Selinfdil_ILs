### get unique arepresentations of ILs

load(file = "datafolder/unique_ILs.R")
load(file = "datafolder/external_test_set/unique_ILs_external.R")

### merge them andget unique cation and anion representations
all_unique <- unique(c(unique_ILs, unique_df_full))

### commas and dots in phosphorous anions


cation <- 0
anion <- 0

for (i in 1:length(all_unique)) {
cation[i] <- strsplit(all_unique[i], split = "_")[[1]][1]
anion[i] <- strsplit(all_unique[i], split = "_")[[1]][2]

}

cation_unique <- unique(cation)
anion_unique <- unique(anion)

### remove a duplicate anion representation [FSI] vs "bis(fluorosulfonyl)imide" 
anion_unique <- anion_unique[-which(anion_unique == "bis(fluorosulfonyl)imide" ) ]

### Get the confusion matrix for presence or absence of selinfdil data for particular ILs.
presence_matrix <- as.data.frame(matrix(FALSE, nrow = length(cation_unique), ncol = length(anion_unique)  ))
rownames(presence_matrix) <- cation_unique
colnames(presence_matrix) <- anion_unique

for (i in 1:length(all_unique)) {

cation_pos <- which(rownames(presence_matrix) %in% strsplit(all_unique[i], split = "_")[[1]][1] )
anion_pos <- which(colnames(presence_matrix) %in% strsplit(all_unique[i], split = "_")[[1]][2] )

presence_matrix[cation_pos,anion_pos] <- TRUE

	} # good for visualization in a publication

save(presence_matrix, file = "datafolder/combinatorial_library_ILs/presence_matrix_ILs.R") 

### ONLY 249 OUT OF 5200 (4.8%) COMBINATIONS WERE EXPERIMENTALLY TESTED!

### get the names of the combinations that were not tried
IL_grid <- expand.grid(cation_unique, anion_unique, stringsAsFactors = FALSE)
ILs_new <- paste(IL_grid [,1], IL_grid [,2], sep = "_")
ILs_new <- ILs_new[!ILs_new %in% all_unique] #removing previously tried combinations

### generate descriptors
## load colnames
load(file = "datafolder/colnames_Selectivity_df.R")

### load the descriptor data

## download external desriptor data 
descr <- list()
component <- c("cation", "anion", "Solute")

for (i in 1:length(component)) {


descr[[i]] <- read.csv( paste("datafolder", component[i],"_small_descr_processed.txt",  sep = ""),
 header=TRUE, check.names=FALSE, sep= "\t",
                 as.is=TRUE,  row.names=NULL) 



	}

### download Paduz. descriptor data
df_full <- read.csv("datafolder/Selectivity_df_QSPR.csv", 
 header=TRUE, check.names=FALSE,
                 as.is=TRUE,  row.names=2)

Solute_descr <- read.csv('datafolder/aniline_descr.txt' ,
 header=TRUE, check.names=FALSE, sep= "\t",
                 as.is=TRUE,  row.names=NULL)

##############################################################################START OF THE DESCRIPTOR LOOP################################################################################################

## empty matrix
descr_df <- as.data.frame(matrix(NA, nrow = length(ILs_new), ncol = length(x) )) # do it at fixed temperature
colnames(descr_df) <- x
rownames(descr_df) <- ILs_new

descr_df[,"T"] <- 298

### fill the matrix. 
## Solute

descr_df[,44:60] <- Solute_descr[1, 2:ncol(Solute_descr)]

## raffinate
## descr file upload

x <- which(df_full[,'solute'] == 'n-dodecane')[1]
descr_df[,61:77] <- df_full[x,44:60]


## remove unnecessary columns
df_full[,44:ncol(df_full)] <- NULL

## cation and anion from external test set
cation_df <- 0
anion_df <- 0

for (i in 1:length(ILs_new)) {
cation_df[i] <- strsplit(ILs_new[i], split = "_")[[1]][1]
anion_df[i] <- strsplit(ILs_new[i], split = "_")[[1]][2]

}

#cation
for (i in 1:length(descr[[1]][,1]) ) {

descr_df[which(cation_df %in% descr[[1]][i,1]) ,27:43] <- descr[[1]][i,3:19]

	}

#anion

for (i in 1:length(descr[[2]][,1]) ) {

descr_df[which(anion_df %in% descr[[2]][i,1]) ,10:26] <- descr[[2]][i,3:19]

	}

## cation and anion from the Paduz.
for (i in 1:length(cation_unique)) {


if (sum(df_full[,"cation"] %in% cation_unique[i]) != 0 ) {


descr_df[cation_df %in% cation_unique[i]  , 27:43] <- df_full[which(df_full[,"cation"] %in% cation_unique[i] )[1] , 27:43]

	}

		}


for (i in 1:length(anion_unique)) {


if (sum(df_full[,"anion"] %in% anion_unique[i]) != 0 ) {


descr_df[anion_df %in% anion_unique[i]  , 10:26] <- df_full[which(df_full[,"anion"] %in% anion_unique[i] )[1] , 10:26]

	}

		}


### checking for NA's

smthfun <- function(x) {

result <- sum(is.na(x))
return(result)

		}

na_vector <- apply(descr_df[,10:ncol(descr_df)], 1, smthfun) #17, why?


### remove unneeded columns
descr_df[,1:8] <- NULL

## test 2
print(sum(apply(descr_df, 1, is.na) )) #equals to zero, as it should

### send unprocessed descriptors to the server for AD calculation
save(descr_df, file = "datafolder/more_az_descr.R")
