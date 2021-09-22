getBinaryStat <- function (obs, pred) {
  if (length(obs) != length(pred)) {
    return(NULL)
  }
  tbl <- table(obs, pred)
  lst <- list()
  lst$TN <- tbl[1,1]
  lst$TP <- tbl[2,2]
  lst$FN <- tbl[2,1]
  lst$FP <- tbl[1,2]
  lst$sensitivity <- lst$TP / (lst$TP + lst$FN)
  lst$specificity <- lst$TN / (lst$TN + lst$FP)
  lst$balanced.acc <- (lst$sensitivity + lst$specificity) / 2
  lst$acc <- mean(obs == pred)
  return(lst)
}

rmse <- function(obs, pred) {
  if (length(obs) != length(pred)) {
    return(NULL)
  }
  return(sqrt(sum((obs - pred) ** 2) / length(obs)))
}

min.dist <- function(ws.df, pred) {
  if (ncol(ws.df) != length(pred)) {
    return(NULL)
  }
  return(min(apply(ws.df, 1, function(x) {sqrt(sum((log1p(x)-log1p(pred))**2))})))
}

R2test <- function(ts.obs, ts.pred, ws.obs) {
  return(1 - sum((ts.pred - ts.obs)**2)/sum((ts.obs - mean(ws.obs))**2))
}

mae <- function(obs, pred) 
{return (sum(abs(obs-pred))/length(obs))
}

mape <- function(obs, pred) 
{return (sum(abs(obs-pred)/abs(obs))/length(obs)*100)
}

binary_crossenthropy <- function(obs,pred) {

H = round( (-1/length(obs) )* (sum( obs*log(pred) + (1-obs)*log(1-pred) ) ) ,9)
return(H)
	}


categorical_crossenthropy <- function(obs,pred) {

enthropy_vector <- vector("numeric", length = nrow(obs) ) 

for (j in 1:nrow(obs)) {

enthropy_vector[j] <- sum( obs[j,]*log(pred[j,]) + (1-obs[j,])*log(1-pred[j,]) )
	}
H = round( (-1/nrow(obs))*sum(enthropy_vector) ,3)
return(H)
		}# rows in input data frames must be observations. Columns are classes. every class in prediction values must have probability > 0.0

