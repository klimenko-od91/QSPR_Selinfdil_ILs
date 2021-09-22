### get descr of predicted azeotropes

load('datafolder/more_az_descr.R')
sss <- vector("numeric", length = nrow(descr_new_az_all))
descr_new_az_all <- cbind(sss, descr_new_az_all)
df_external <- as.matrix(descr_new_az_all)

### get SelinfDB set
df_full <- read.csv("datafolder/Selectivity_df_QSPR.csv", 
 header=TRUE, check.names=FALSE,
                 as.is=TRUE,  row.names=2)

### optimization set labeling
load("datafolder/ts_index.R")
ts_index <- ts_index[-1]

## determining the thresholds, non-normalised values
colMin <- apply(df_full[-ts_index,9:ncol(df_full)],2,min)
colMax <- apply(df_full[-ts_index,9:ncol(df_full)],2,max)

AD_matrix <- as.data.frame(matrix(0, nrow = nrow(df_external), ncol = (ncol(df_external)-1)))
## predicting structural AD for the test set

for (i in 1:nrow(df_external) ) {

AD_matrix[i,] <- df_external[i,2:ncol(df_external)] >= colMin & df_external[i,2:ncol(df_external)] <= colMax
	}

AD_outcome <- apply(AD_matrix,1,prod)
save(AD_outcome, file = "datafolder/Bounding_box_more_az.R")


