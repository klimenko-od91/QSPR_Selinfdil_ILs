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
source("datafolder/stat.R")

best_results_list <- as.data.frame(matrix(0, nrow = 72, ncol = 5)) 

for (k in 1:72) {

###########################################################################roster loop starts##############################################################################################################

df_pred <- read.csv(paste("scriptfolder/log10S/20logS/", k, "_ts_pred_Selectivity.csv", sep = "" ) , 
 header=TRUE, check.names=FALSE,
                 as.is=TRUE,  row.names=NULL)

df_test <- df_full[ts_index, 1:8]

df_test <- cbind(df_test, df_pred)

df_test[,11] <- round(df_test[,11], 2)

best_results_list[k,1] <- mae(df_test[,10], df_test[,11])

# mae per IL
unique_ILs <- unique(df_test[,"ILs"])

IL_index <- 0
mae_per_comp <- list(0)

for (i in 1:length(unique_ILs)) {

IL_index <- df_test[df_test[,"ILs"] %in% unique_ILs[i] , ]

mae_per_comp[i] <- mae(IL_index[,10], IL_index[,11])


	}

names(mae_per_comp) <- unique_ILs
save(mae_per_comp, file = "scriptfolder/log10S/20logS/mae_per_comp.R") 
best_results_list[k,2] <- mean(as.numeric(mae_per_comp))

### range
range_obs <- max(df_test[,10]) - min(df_test[,10])
range_pred <- max(df_test[,11]) - min(df_test[,11])
range_opt <- abs(range_pred - range_obs )
best_results_list[k,3] <- range_opt

### Calculating covariance. Temperature is normalised.
covariance_obs <- 0
covariance_pred <- 0

covariance_obs <- cov(df_test[,10], df_test[,12])
covariance_pred <- cov(df_test[,11],df_test[,12])

cov_difference_opt <- abs(covariance_obs- covariance_pred)
best_results_list[k,4] <- cov_difference_opt

### calculating the function
best_results_list[k,5] <- prod(best_results_list[k,1:4])**(1/4)

	} #end of the loop

### output
colnames(best_results_list) <- c("MAE","MAE per IL","range","covariance", "decision")

write.csv(best_results_list, file = "scriptfolder/log10S/20logS/best_results_list.csv")

# 20 is the best

### get CV results for the best model
k = 20

Internal_Val_stat <- as.data.frame(matrix(0, nrow = 7, ncol = 2)) 
rownames(Internal_Val_stat) <- c("fold1", "fold2", "fold3", "fold4", "fold5", "mean", "sd")
colnames(Internal_Val_stat) <- c("MAE_TrS", "MAE_CV")

for (j in 1:5) {

Trs_MAE <- read.csv(paste("scriptfolder/log10S/20logS/", k, "_", j, "_model_fitted_min.csv", sep = "" ) , 
 header=TRUE, check.names=FALSE,
                 as.is=TRUE,  row.names=1)

Val_MAE <- read.csv(paste("scriptfolder/log10S/20logS/", k, "_", j, "_model_fitted_df.csv", sep = "" ) , 
 header=TRUE, check.names=FALSE,
                 as.is=TRUE,  row.names=1)

Chosen_Val_MAE <- Val_MAE[which(Val_MAE[,1] %in% Trs_MAE[1,1] ) , 3] #the fold model for prediction was chosen based on the lowest TrS MAE 

Internal_Val_stat[j,1] <- Trs_MAE[1,1]
Internal_Val_stat[j,2] <- Chosen_Val_MAE

	}

Internal_Val_stat["mean",] <- round(apply(Internal_Val_stat[1:5,], 2, mean) ,2)
Internal_Val_stat["sd",] <- round(apply(Internal_Val_stat[1:5,], 2, sd) ,2)

write.csv(Internal_Val_stat, file = "scriptfolder/log10S/20logS/Selectivity_20_Internal_Val_stat.csv")

######################################################################Classification interpretation od the results##################################################################################
rm(list = ls())

k = 20

df_pred <- read.csv(paste("scriptfolder/log10S/20logS/", k, "_ts_pred_Selectivity.csv", sep = "" ) , 
 header=TRUE, check.names=FALSE,
                 as.is=TRUE,  row.names=NULL) # this will be called differently! 72_ts_pred_Selectivity.csv

df_test <- df_full[ts_index, 1:8]

df_test <- cbind(df_test, df_pred)

df_test[,11] <- round(df_test[,11], 2)

### calculate classification interpretation
df_test[, "lg [K]"] <- df_test[, "lg [K]"] >= 1
df_test[, "pred_test"] <- df_test[, "pred_test"] >= 1

### classification stats
Acc <- sum(df_test[, "lg [K]"] == df_test[, "pred_test"])/nrow(df_test)
Sens <-  sum(df_test[, "lg [K]"] == "TRUE" & df_test[, "pred_test"] == "TRUE")   /   sum(df_test[, "lg [K]"])
Spec <- sum(df_test[, "lg [K]"] == "FALSE" & df_test[, "pred_test"] == "FALSE")   /   sum(df_test[, "lg [K]"] == "FALSE") 
BA <- (Sens+Spec)/2
PPV <- sum(df_test[, "lg [K]"] == "TRUE" & df_test[, "pred_test"] == "TRUE")   /   sum(df_test[, "pred_test"])
NPV <- sum(df_test[, "lg [K]"] == "FALSE" & df_test[, "pred_test"] == "FALSE")   / sum(df_test[, "pred_test"] == "FALSE")
TP <- sum(df_test[, "lg [K]"] == "TRUE" & df_test[, "pred_test"] == "TRUE")
TN <- sum(df_test[, "lg [K]"] == "FALSE" & df_test[, "pred_test"] == "FALSE")
FP <- sum(df_test[, "lg [K]"] == "FALSE" & df_test[, "pred_test"] == "TRUE")
FN <- sum(df_test[, "lg [K]"] == "TRUE" & df_test[, "pred_test"] == "FALSE")

stat_class <- c(Acc, Sens, Spec, BA, PPV, NPV, TP, TN, FP, FN)
names(stat_class) <- c("Acc", "Sens", "Spec", "BA", "PPV", "NPV", "TP", "TN", "FP", "FN")
write.csv(stat_class, file = "scriptfolder/log10S/20logS/20_stat_class_test_no_AD.csv")

