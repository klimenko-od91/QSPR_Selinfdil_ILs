# Activate the library
library('tensorflow')
library('keras')

df_test <- read.csv("datafolder/external_descr_full.csv", 
 header=TRUE, check.names=FALSE,
                 as.is=TRUE,  row.names=1)

df_test <- as.matrix(df_test)

### optimization test data

### importing

load(file ="datafolder/Selectivity_k.R")
load(file ="datafolder/Selectivity_c.R")

for (i in 2:ncol(df_test)) {

df_test[,i] <- round(	(k[i]*df_test[,i]+c[i]),	1)

}

#j is a artifact from trying all models on the best test set
j = 29

### load model of interest
pred_test <- matrix(0, nrow = nrow(df_test), ncol = 5)

for (i in 1:5) {

model <- load_model_tf(paste("scriptfolder/BIDAC/50BIDAC/", j,"_", i,"_","BIGIDAC_NN_50.tf", sep = ""), c('PReLU' = layer_activation_parametric_relu), compile = TRUE)

# test
pred_test_fold <-  predict(model,df_test[, 2:ncol(df_test)])

pred_test[,i] <- pred_test_fold

}

### averaging 5-fold results
pred_test <- apply(pred_test,1,mean)
pred_test_df <- data.frame(df_test[, 1], pred_test ,df_test[, 2])

### export
write.csv(pred_test_df, file = paste("scriptfolder/BIDAC/50BIDAC/", j,"_","external_ts_pred_Selectivity.csv", sep = "") ) 

	} #end of the old roster loop
