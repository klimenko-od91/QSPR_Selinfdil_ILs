df_test <- read.csv("datafolder/Selectivity_external_for_SQL.csv", 
 header=TRUE, check.names=FALSE,
                 as.is=TRUE,  row.names=1) 

k = 9 #best model

df_pred <- read.csv(paste("scriptfolder/BIDAC/20BIDAC/", k, "_external_ts_pred_BIDAC.csv", sep = "" ) , 
 header=TRUE, check.names=FALSE,
                 as.is=TRUE,  row.names=NULL)

df_test <- cbind(df_test, df_pred)

### change the probabilities into classes
df_test[df_test[,"pred_test"] >= 0.5,  "pred_test"] <- 1
df_test[df_test[,"pred_test"] < 0.5,  "pred_test"] <- 0

### checking for potential training or test set compounds data. had to inspect manually
## get names of the Selectivity_df ILs
df_full <- read.csv("scriptfolder/BIDAC/20BIDAC/Selectivity_df_QSPR.csv", 
 header=TRUE, check.names=FALSE,
                 as.is=TRUE,  row.names=2)

unique_df_full <- unique(df_full[, "ILs"])

unique_ILs <- unique(df_test[,"Extracion phase"])

x <- grep("IM-6,1_CL", unique_df_full) # vary cation
unique_df_full[x]# found presence of : "PYR-4,1_DCA" "IM-6,1_CL"

### removing results for 2 IL present in SelinfDB

df_test <- df_test[which(df_test[,"Extracion phase"] != unique_ILs[5]) , ]
df_test <- df_test[which(df_test[,"Extracion phase"] != unique_ILs[27]) , ]

unique_ILs <- unique_ILs[-c(5, 27)]

### download stats
source("datafolder/stat_2020.R")

### make stat vector
best_results_list <- 0

### calculate BA
best_results_list <- as.vector(getBinaryStat(df_test[,"big IDAC flag"], df_test[,"pred_test"]))

### calculate accuracy per IL, because some ILs only have all positive or all negative obs or predicted results

unique_ILs <- unique(df_test[,"Extracion phase"])

IL_index <- 0
Acc <- list(0)

for (i in 1:length(unique_ILs)) {

IL_index <- df_test[df_test[,"Extracion phase"] %in% unique_ILs[i] , ]

Acc[i] <- (sum(IL_index[,"big IDAC flag"] == IL_index[,"pred_test"]))/nrow(IL_index)

	}

names(Acc) <- unique_ILs

best_results_list[9] <- mean(as.numeric(Acc))

### geometric mean for BA and Acc per compound
best_results_list[10] <- prod(as.numeric(best_results_list[c(7,9)]))**(1/2)

### names 
names(best_results_list) <- c("TN","TP","FN","FP", "sensitivity", "specificity", "BA", "Acc", "Acc (per compound)", "decision")
write.csv(best_results_list, file = "scriptfolder/BIDAC/20BIDAC/9_stat_external.csv")
