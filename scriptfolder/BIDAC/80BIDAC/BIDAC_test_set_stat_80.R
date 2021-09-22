### download data

df_full <- read.csv("datafolder/Selectivity_df_QSPR.csv", 
 header=TRUE, check.names=FALSE,
                 as.is=TRUE,  row.names=2)

### optimization set selection
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

df_pred <- read.csv(paste("scriptfolder/BIDAC/80BIDAC/", k, "_ts_pred_BIDAC.csv", sep = "" ) , 
 header=TRUE, check.names=FALSE,
                 as.is=TRUE,  row.names=NULL) # 

df_test <- df_full[ts_index, 1:7]

df_test <- cbind(df_test, df_pred)

### 'perfect prediction probability' does not allow to apply binary crossentropy

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

write.csv(best_results_list, file = "scriptfolder/BIDAC/80BIDAC/best_results_list.csv")

		}


### get CV results for the best model
k = 33


Internal_Val_stat <- as.data.frame(matrix(0, nrow = 7, ncol = 3)) 
rownames(Internal_Val_stat) <- c("fold1", "fold2", "fold3", "fold4", "fold5", "mean", "sd")
colnames(Internal_Val_stat) <- c("BinCross (TrS)", "BinCross (Val)". "Acc (Val)")

for (j in 1:5) {

Trs_BinCross <- read.csv(paste("scriptfolder/BIDAC/80BIDAC/", k, "_", j, "_BIGIDAC_model_fitted_min.csv", sep = "" ) , 
 header=TRUE, check.names=FALSE,
                 as.is=TRUE,  row.names=1)

Val_BinCross <- read.csv(paste("scriptfolder/BIDAC/80BIDAC/", k, "_", j, "_BIGIDAC_model_fitted_df.csv", sep = "" ) , 
 header=TRUE, check.names=FALSE,
                 as.is=TRUE,  row.names=1)

Chosen_Val_BinCross <- Val_BinCross[which(Val_BinCross[,1] %in% Trs_BinCross[1,1] ) , 3]
Chosen_Val_Acc <- Val_BinCross[which(Val_BinCross[,1] %in% Trs_BinCross[1,1] ) , 4] #the fold model for prediction was chosen based on the lowest TrS binary cross entropy 

Internal_Val_stat[j,1] <- Trs_BinCross[1,1]
Internal_Val_stat[j,2] <- Chosen_Val_BinCross
Internal_Val_stat[j,3] <- Chosen_Val_Acc

	}

Internal_Val_stat["mean",] <- round(apply(Internal_Val_stat[1:5,], 2, mean) ,5)
Internal_Val_stat["sd",] <- round(apply(Internal_Val_stat[1:5,], 2, sd) ,5)

write.csv(Internal_Val_stat, file = "scriptfolder/BIDAC/80BIDAC/BIDAC_33_Internal_Val_stat.csv")


