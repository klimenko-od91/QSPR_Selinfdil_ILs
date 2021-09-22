# Activate the library
library('tensorflow')
library('keras')

### download data
df_full <- read.csv("datafolder/Selectivity_df_QSPR.csv", 
 header=TRUE, check.names=FALSE,
                 as.is=TRUE,  row.names=2)

### optimization set selection
load("datafolder/ts_index.R")

ts_index <- ts_index[-1]

# load external script containing neccessary functions
source("datafolder/stat.R") 

### removing non-numeric data columns
df_full <- df_full[,-(1:7)]

### setting aside the trainingset
df_test <- df_full[ts_index,] 

df_test <- as.matrix(df_test)

### normalizing optimization data

load(file ="datafolder/Selectivity_k.R")
load(file ="datafolder/Selectivity_c.R")

for (i in 2:ncol(df_test)) {

df_test[,i] <- round(	(k[i]*df_test[,i]+c[i]),	1)

}

for (j in 1:72) {

### load model of interest
pred_test <- matrix(0, nrow = nrow(df_test), ncol = 5)

for (i in 1:5) {

model <- load_model_tf(paste("scriptfolder/log10S/50logS/", j,"_", i,"_","Selectivity_NN_50.tf", sep = ""), c('PReLU' = layer_activation_parametric_relu), compile = TRUE)

pred_test_fold <-  predict(model,df_test[, 2:ncol(df_test)])

pred_test[,i] <- pred_test_fold

}

### averaging 5-fold results
pred_test <- apply(pred_test,1,mean)
pred_test_df <- data.frame(df_test[, 1], pred_test ,df_test[, 2])

### export
write.csv(pred_test_df, file = paste("scriptfolder/log10S/50logS/", j,"_","ts_pred_Selectivity.csv", sep = "") ) 

	} #end of the roster loop

