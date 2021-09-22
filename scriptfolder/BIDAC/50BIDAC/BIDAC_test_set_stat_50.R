### download data

df_full <- read.csv("datafolder/Selectivity_df_QSPR.csv", 
 header=TRUE, check.names=FALSE,
                 as.is=TRUE,  row.names=2)

### test set selection
load("datafolder/ts_index.R")
ts_index <- ts_index[-1]

### give AD
load("datafolder/Bounding_box.R") #100% in AD - no need for additional calculations

### assessment with a decision function
#download stats
source("datafolder/stat_2020.R")

best_results_list <- as.data.frame(matrix(0, nrow = 72, ncol = 10))

for (k in 1:72) {

###########################################################################roster loop starts##############################################################################################################

df_pred <- read.csv(paste("scriptfolder/BIDAC/50BIDAC/", k, "_ts_pred_BIDAC.csv", sep = "" ) , 
 header=TRUE, check.names=FALSE,
                 as.is=TRUE,  row.names=NULL) # 

df_test <- df_full[ts_index, 1:7]

df_test <- cbind(df_test, df_pred)

### change the probabilities into classes
df_test[df_test[,"pred_test"] >= 0.5,  "pred_test"] <- 1
df_test[df_test[,"pred_test"] < 0.5,  "pred_test"] <- 0

### calculate BA
best_results_list[k,1:8] <- as.vector(getBinaryStat(df_test[,"big IDAC"], df_test[,"pred_test"]))

### calculate accuracy per IL, because some ILs only have all positive or all negative obs or predicted results

unique_ILs <- unique(df_test[,"ILs"])

IL_index <- 0
Acc <- list(0)

for (i in 1:length(unique_ILs)) {

IL_index <- df_test[df_test[,"ILs"] %in% unique_ILs[i] , ]

Acc[i] <- (sum(IL_index[,"big IDAC"] == IL_index[,"pred_test"]))/nrow(IL_index)

	}

names(Acc) <- unique_ILs

best_results_list[k,9] <- mean(as.numeric(Acc))

### geometric mean for BA and Acc per compound
best_results_list[k,10] <- prod(best_results_list[k,c(7,9)])**(1/2)

### output
colnames(best_results_list) <- c("TN","TP","FN","FP", "sensitivity", "specificity", "BA", "Acc", "Acc (per compound)", "decision")

write.csv(best_results_list, file = "scriptfolder/BIDAC/50BIDAC/best_results_list.csv")

		} 

which(best_results_list[, 'decision'] == max(best_results_list[, 'decision'] )) #29 is the best

