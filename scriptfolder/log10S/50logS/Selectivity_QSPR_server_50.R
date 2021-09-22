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
dict_df <- df_full[,1:7]

df_full <- df_full[,-(1:7)]

### checking if columns have constant values
df_full <- df_full[, apply(df_full[-ts_index,], 2, sd) !=0 ] #69 descriptors, including T

### setting aside the trainingset
df_train <- df_full[-ts_index,]

### randomise the data so that models would be less biased
df_train <- df_train[sample(nrow(df_train),nrow(df_train)) , ]
x <- rownames(df_train)

# remove extra
rm(dict_df,df_full)

## ALTERNATIVE: by minmax
k <- 0
c <- 0
###  k and c must be saved for normalisation of future compounds 
for (i in 2:ncol(df_train)) {
k[i] <- 1/(max(df_train[,i])-min(df_train[,i]))
c[i] <- 1-k[i]*max(df_train[,i])

df_train[,i] <- round(	(k[i]*df_train[,i]+c[i]),	1)

}

### exporting
save(k, file ="datafolder/Selectivity_k.R")
save(c, file ="datafolder/Selectivity_c.R")

## other preparations

colnames(df_train) <- NULL
df_train <- as.matrix(df_train)

roster <- expand.grid(c("rmsprop", "adamax","adadelta", "sgd"), c("lecun_uniform", "glorot_uniform", "he_uniform"), c(100,200), c(10, 20, 40)   )
colnames(roster) <- c("optimization", "initializer", "nodes1", "nodes2")
write.csv(roster,file = "datafolder/roster_Selectivity.csv")

################################################################################roster loop##################################################################################################################

for (j in 1:nrow(roster)) {


### global variable, just in case
model_fitted <- 0

## making function. input_shape ,creating framework

build_model <- function() {

  model <- keras_model_sequential() %>%  
layer_dense(units = roster[j, "nodes1"], kernel_initializer= roster[j, "initializer"],
                activation = layer_activation_parametric_relu(), input_shape = (ncol(df_train)-1) ) %>%  
layer_dense(units =roster[j, "nodes2"], activation = "relu") %>%
layer_dense(units = 1,activation = "linear") 
  
  model %>% compile(
    loss = 'mean_absolute_error',
    optimizer = roster[j, "optimization"],
    metrics = list('mean_absolute_error', 'mean_absolute_percentage_error')
  )
  
  model
}

model <- build_model()

### loop will start here
for (i in 1:5) {

### continue splitting data because of memory error
spl_val <- seq(from = i, to = nrow(df_train), by = 2	)

## Creating callbacks for saving and stopping
some_path <- paste("scriptfolder/log10S/50logS/", j,"_",i,"_","Selectivity_NN_50.tf", sep = "") 

checkpoint_NN <- callback_model_checkpoint(some_path, monitor = "loss", verbose = 1, save_best_only = TRUE, save_weights_only = FALSE, mode = "min", period = 1)
early_stop_NN <- callback_early_stopping(monitor = "loss", min_delta = 0,patience = 50, verbose = 1, mode = "min")

### Fitting
model_fitted <- model %>% fit(
  df_train[-spl_val,2:ncol(df_train)] , #will be var column
  df_train[-spl_val,1]	,#will be 1st column
  epochs = 500,
  
  validation_data = list(df_train[spl_val,2:ncol(df_train)] , df_train[spl_val,1]), 
  verbose = 1,
  callbacks = list(early_stop_NN,checkpoint_NN),
shuffle = TRUE,
batch_size = 2000 )

### saving model development data
model_fitted_df <- data.frame(model_fitted$metrics$mean_absolute_error,
model_fitted$metrics$mean_absolute_percentage_error,
model_fitted$metrics$val_mean_absolute_error,
model_fitted$metrics$val_mean_absolute_percentage_error)

model_fitted_min <- min(model_fitted$metrics$mean_absolute_error)

write.csv(model_fitted_df, file = paste("scriptfolder/log10S/50logS/", j,"_",i,"_", "model_fitted_df",".csv", sep = ""))
write.csv(model_fitted_min, file = paste("scriptfolder/log10S/50logS/", j,"_",i,"_", "model_fitted_min",".csv", sep = ""))

}

	} ### end of roster loop
