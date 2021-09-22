
df_full <- read.csv("datafolder/Selectivity_df_QSPR.csv", 
 header=TRUE, check.names=FALSE,
                 as.is=TRUE,  row.names=2)

### optimization set selection
load("datafolder/ts_index.R")
ts_index <- ts_index[-1]

## determining the thresholds, non-normalised values
colMin <- apply(df_full[-ts_index,9:ncol(df_full)],2,min)
colMax <- apply(df_full[-ts_index,9:ncol(df_full)],2,max)

AD_matrix <- as.data.frame(matrix(0, nrow = length(ts_index), ncol = (ncol(df_full)-8)))
## predicting structural AD for the optimization set
for (i in 1:length(ts_index) ) {
AD_matrix[i,] <- df_full[ts_index[i],9:ncol(df_full)] >= colMin & df_full[ts_index[i],9:ncol(df_full)] <= colMax

	}

AD_outcome <- apply(AD_matrix,1,prod)
save(AD_outcome, file = "datafolder/Bounding_box.R")



