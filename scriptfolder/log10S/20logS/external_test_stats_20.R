df_test <- read.csv("datafolder/Selectivity_external_for_SQL.csv", 
 header=TRUE, check.names=FALSE,
                 as.is=TRUE,  row.names=1)

k = 20 #best model

df_pred <- read.csv(paste("scriptfolder/log10S/20logS/", k, "_external_ts_pred_Selectivity.csv", sep = "" ) , 
 header=TRUE, check.names=FALSE,
                 as.is=TRUE,  row.names=NULL) 

df_test <- cbind(df_test, df_pred)

df_test[,"pred_test"] <- round(df_test[,"pred_test"], 2)

### checking for potential training or test set compounds data. It was done manually
## get names of the Selectivity_df ILs
df_full <- read.csv("datafolder/Selectivity_df_QSPR.csv", 
 header=TRUE, check.names=FALSE,
                 as.is=TRUE,  row.names=2)

unique_df_full <- unique(df_full[, "ILs"])

unique_ILs <- unique(df_test[,"Extracion phase"])

## manual
x <- grep("IM-6,1_CL", unique_df_full) # vary cation
unique_df_full[x]# found presence of : "PYR-4,1_DCA" "IM-6,1_CL"

### removing results for 2 IL present in SelinfDB

df_test <- df_test[which(df_test[,"Extracion phase"] != unique_ILs[5]) , ]
df_test <- df_test[which(df_test[,"Extracion phase"] != unique_ILs[27]) , ] 

unique_ILs <- unique_ILs[-c(5, 27)]

#download stats
source("datafolder/stat.R")

##calculate - no AD
stat_vector <- vector("numeric", length = 5)
names(stat_vector) <- c("MAE","MAE per IL","range","covariance", "decision")

stat_vector["MAE"] <- mae(df_test[,"Selectivity [lg]"], df_test[,"pred_test"])

IL_index <- 0
mae_per_comp <- list(0)

for (i in 1:length(unique_ILs)) {

IL_index <- df_test[df_test[,"Extracion phase"] %in% unique_ILs[i] , ]

mae_per_comp[i] <- mae(IL_index[,"Selectivity [lg]"], IL_index[,"pred_test"]) 

	}

stat_vector["MAE per IL"] <- mean(as.numeric(mae_per_comp))
names(mae_per_comp) <- unique_ILs

range_obs <- max(df_test[,"Selectivity [lg]"]) - min(df_test[,"Selectivity [lg]"])
range_pred <- max(df_test[,"pred_test"]) - min(df_test[,"pred_test"])
range_opt <- abs(range_pred - range_obs )
stat_vector["range"] <- range_opt 

### Calculating covariance.
covariance_obs <- 0
covariance_pred <- 0

covariance_obs <- cov(df_test[,"Selectivity [lg]"], df_test[,ncol(df_test)])
covariance_pred <- cov(df_test[,"pred_test"],df_test[,ncol(df_test)])

cov_difference_opt <- abs(covariance_obs- covariance_pred)
stat_vector["covariance"] <- cov_difference_opt # predicted covariance is bigger than observed one. In the future do not use covariance or atleast use covariance per compound

stat_vector["decision"] <- prod(stat_vector[1:4])**(1/4)

### calculate fraction of big IDAC for every compound
BIGIDAC_fraction <- 0

for (i in 1:length(unique_ILs)) {


IL_index <- which(df_test[,"Extracion phase"] %in% unique_ILs[i]) 
BIGIDAC_fraction[i]  <- sum(df_test[IL_index, "big IDAC flag"])/length(IL_index)

	} 

names(BIGIDAC_fraction) <- unique_ILs

BIGIDAC_fraction <- BIGIDAC_fraction[order(as.numeric(mae_per_comp))]

write.csv(stat_vector, file = "scriptfolder/log10S/20logS/stat_external_no_AD_20.csv")

### calculate classification interpretation
colnames(df_test)[8] <- "obs"
df_test[, "obs"] <- df_test[, "obs"] >= 1
df_test[, "pred_test"] <- df_test[, "pred_test"] >= 1

nrow(df_test[df_test[, "obs"] == "TRUE" & df_test[, "big IDAC flag"] == 0 , ])

### classification stats
Acc <- sum(df_test[, "obs"] == df_test[, "pred_test"])/nrow(df_test)
Sens <-  sum(df_test[, "obs"] == "TRUE" & df_test[, "pred_test"] == "TRUE")   /   sum(df_test[, "obs"])
Spec <- sum(df_test[, "obs"] == "FALSE" & df_test[, "pred_test"] == "FALSE")   /   sum(df_test[, "obs"] == "FALSE") 
BA <- (Sens+Spec)/2
PPV <- sum(df_test[, "obs"] == "TRUE" & df_test[, "pred_test"] == "TRUE")   /   sum(df_test[, "pred_test"])
NPV <- sum(df_test[, "obs"] == "FALSE" & df_test[, "pred_test"] == "FALSE")   / sum(df_test[, "pred_test"] == "FALSE")
TP <- sum(df_test[, "obs"] == "TRUE" & df_test[, "pred_test"] == "TRUE")
TN <- sum(df_test[, "obs"] == "FALSE" & df_test[, "pred_test"] == "FALSE")
FP <- sum(df_test[, "obs"] == "FALSE" & df_test[, "pred_test"] == "TRUE")
FN <- sum(df_test[, "obs"] == "TRUE" & df_test[, "pred_test"] == "FALSE")

stat_class <- c(Acc, Sens, Spec, BA, PPV, NPV, TP, TN, FP, FN)
names(stat_class) <- c("Acc", "Sens", "Spec", "BA", "PPV", "NPV", "TP", "TN", "FP", "FN")
write.csv(stat_class, file = "scriptfolder/log10S/20logS/stat_class_external_no_AD_20.csv")
