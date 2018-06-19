#' Building predictive models
#'
#' This function give us the tools to build predictive models for time series.
#'
#' Returns an object \code{modl} which stores all the information related to
#' the final chosen model (errors, parameters, model).
#'
#' Currently this function covers two different methods: the widely know ARIMA
#' and the "not so used for prediction" data mining. For the data mining we make
#' use of the \code{caret} package.
#'
#' The \code{caret} package offers plenty of data mining algorithms.
#' For the data splitting here we use a rolling forecasting origin technique, wich
#' works better on time series.
#'
#' @param tserie A ts or prep object.
#' @param method A string. Current methods available are "arima" and "dataMining". Method "arima" is set as default.
#' @param algorithm A string. In case \code{method} is "dataMining", pick the algorithm you want to use. There is a complete list of available algorithms here (only regression type allowed): \url{http://topepo.github.io/caret/train-models-by-tag.html}.
#' @param formula An integer vector. Contains the indexes from the time series wich will indicate how to extract the features. The last value will be the class index. Default value: c(1:16)
#' @param initialWindow An integer. The initial number of consecutive values in each training set sample. Default value: 30.
#' @param horizon An integer. The number of consecutive values in test set sample. Default value: 15.
#' @param fixedWindow A logical: if FALSE, the training set always start at the first sample and the training set size will vary over data splits. Default value: TRUE.
#' @return A list is returned of class \code{modl} containing:
#' \item{tserie}{Original time serie.}
#' \item{tserieDF}{Time serie converted to data frame.}
#' \item{method}{Method used to build the model.}
#' \item{algorithm}{If method is data mining, indicates wich algorithm was used.}
#' \item{horizon}{Horizon for the splitting.}
#' \item{model}{Model result from \code{caret}. It is a list, result of the \code{caret::train} function.}
#' \item{errors}{Contains three different metrics to evaluate the model.}
#' @author Alberto Vico Moreno
#' @seealso{
#' \code{\link{prep}}
#' \code{\link{modl.arima}},
#' \code{\link{modl.tsToDataFrame}},
#' \code{\link{modl.trControl}},
#' \code{\link{modl.dataMining}}
#' }
#' @export
#' @references \url{http://topepo.github.io/caret/index.html}
#' @examples
#' p <- prep(AirPassengers)
#' modl(p,method='arima')
#' \donttest{modl(p,method='dataMining',algorithm='rpart')}
modl <- function(tserie, method='arima'
                  ,algorithm=NULL
                  ,formula=NULL
                  ,initialWindow=NULL,horizon=NULL,fixedWindow=NULL){
  if(methods::is(tserie,"prep")) tserie <- tserie$tserie
  else if(!stats::is.ts(tserie)) stop('Not a ts or prep object')

  model <- NULL
  errors <- NULL
  tserieDF <- NULL
  k <- 60 #minimum number of observations to make good predictions

  #arima model
  if(method=='arima'){
    model <- modl.arima(tserie)

    errors <- c(forecast::accuracy(model)[2],forecast::accuracy(model)[3],forecast::accuracy(model)[4]) #RMSE, MAE and MAPE errors
  }

  #data mining model
  else if(method=='dataMining'){
    if(is.null(algorithm)) stop('Data mining regression algorithm required. Consult names here:\nhttp://topepo.github.io/caret/train-models-by-tag.html')

    if(is.null(fixedWindow)) fixedWindow <- TRUE
    if(is.null(initialWindow)) initialWindow <- as.integer(k/2)
    if(is.null(horizon)) horizon <- as.integer(k/4)

    timeControl <- modl.trControl(initialWindow,horizon,fixedWindow,givenSummary=TRUE)
    tserieDF <- modl.tsToDataFrame(tserie,formula)
    model <- modl.dataMining(Class ~ .,tserieDF,algorithm,timeControl,metric='RMSE',maximize=FALSE)

    winningIndex <- which.min(model$results["RMSE"][,1])
    errors <- c(model$results["RMSE"][winningIndex,1]
                ,model$results["MAE"][winningIndex,1]
                ,model$results["MAPE"][winningIndex,1])
  }
  else stop('Invalid modelling method')

  names(errors) <- c("RMSE","MAE","MAPE")

  #creating the object
  obj <- list()
  class(obj) <- "modl"
  obj$tserie <- tserie
  obj$tserieDF <- tserieDF
  obj$method <- method
  obj$algorithm <- algorithm
  obj$horizon <- horizon
  obj$model <- model
  obj$errors <- errors

  return (obj)
}

#' Automatic ARIMA model
#'
#' Assuming "tserie" is stationary, returns the best arima model
#'
#' @param tserie A ts object.
#' @return ARIMA model.
#' @author Alberto Vico Moreno
#' @export
#' @examples
#' modl.arima(AirPassengers)
modl.arima <- function(tserie){
  return (forecast::auto.arima(tserie,seasonal=FALSE))
}

#' Ts to data frame transformation
#'
#' Transform a ts object into a data frame using the given formula.
#'
#' @param tserie A ts object.
#' @param formula An integer vector. Contains the indexes from the \code{tserie} wich will indicate how to extract the features. The last value will be the class index. Default value: c(1:16). Has to be length 6 minimum.
#' @return the time serie as data frame
#' @author Alberto Vico Moreno
#' @export
#' @examples
#' modl.tsToDataFrame(AirPassengers,formula=c(1,3,4,5,6,7))
#' modl.tsToDataFrame(AirPassengers,formula=c(1:20))
modl.tsToDataFrame <- function(tserie,formula=NULL){
  if(!stats::is.ts(tserie)) stop('Not a ts object')

  k <- 60 #minimum number of observations to make good predictions

  if(is.null(formula)) formula <- c(1:16)
  else if(length(formula)<6) stop('Parameter formula needs to have at least 6 elements')

  numE <- max(formula)-min(formula)+1 #number of features per observation
  n <- length(formula)-1 #portion of the time serie needed to make an observation

  if(length(tserie)<numE) stop('Sample exceeds ts length')

  mat <- data.frame()

  continuar <- TRUE
  i <- 0

  while(continuar){
    aux <- 1 #aux index
    vaux <- NULL #aux vector
    while(aux <= n+1){ #loop to make an observation
      vaux <- c(vaux,tserie[i+formula[aux]])
      aux <- aux+1
    }
    mat <- rbind(mat,vaux) #adding the observation to the data frame
    i <- i+1
    if(i > length(tserie)-numE) continuar <- FALSE #not enough elements left to make another observation
  }

  if(i < k) stop('Not enough observations to make good prediction: you need minimum a 60 length series.\n')

  for(i in 1:(n+1)){
    if (i==n+1) colnames(mat)[i] <- 'Class'
    else colnames(mat)[i] <- paste0('c_',i)
  }

  return (mat)
}

#' Control the splitting to train the data
#'
#' Creates the needed \code{caret::trainControl} object to control the training
#' splitting.
#'
#' We always split using method "timeslice", wich is the better for time series.
#' More information on how this works on \url{http://topepo.github.io/caret/data-splitting.html#data-splitting-for-time-series}.
#'
#' @param initialWindow An integer. The initial number of consecutive values in each training set sample. Default value: 30.
#' @param horizon An integer. The number of consecutive values in test set sample. Default value: 15.
#' @param fixedWindow A logical: if FALSE, the training set always start at the first sample and the training set size will vary over data splits. Default value: TRUE.
#' @param givenSummary A logical. Indicates if it should be used the customized summaryFunction(?trainControl for more info) modl.sumFunction or not. Default is FALSE; this will use default \code{caret} metrics.
#' @return trainControl object
#' @author Alberto Vico Moreno
#' @export
#' @examples
#' modl.trControl(initialWindow=30,horizon=15,fixedWindow=TRUE,givenSummary=TRUE)
modl.trControl <- function(initialWindow,horizon,fixedWindow,givenSummary=FALSE){
  k <- 60 #minimum number of observations to make good predictions

  if(givenSummary)
    return (caret::trainControl(method="timeslice",initialWindow=initialWindow,horizon=horizon,fixedWindow=fixedWindow,summaryFunction = modl.sumFunction))
  else
    return (caret::trainControl(method="timeslice",initialWindow=initialWindow,horizon=horizon,fixedWindow=fixedWindow))
}

#' Train the data
#'
#' Train the time serie(as data frame) to build the model.
#'
#' @param form A formula of the form y ~ x1 + x2 + ...
#' @param tserieDF Data frame.
#' @param algorithm A string. Algorithm to perform the training. Full list at \url{http://topepo.github.io/caret/train-models-by-tag.html}. Only regression types allowed.
#' @param timeControl trainControl object.
#' @param metric A string. Specifies what summary metric will be used to select the optimal model. Possible values in \code{caret} are "RMSE" and "Rsquared". "RMSE" set as default. If you used a custom summaryFunction(see ?trainControl) your metrics will prevail over default.
#' @param maximize A logical. Should the metric be maximized or minimized? Default is FALSE, since that is what makes sense for time series.
#' @return train object
#' @author Alberto Vico Moreno
#' @export
#' @examples
#' \donttest{
#' modl.dataMining(form=Class ~ .,
#'  tserieDF=modl.tsToDataFrame(AirPassengers,formula=c(1:20)),
#'  algorithm='rpart',
#'  timeControl=modl.trControl(initialWindow=30,horizon=15,fixedWindow=TRUE))
#' }
modl.dataMining <- function(form,tserieDF,algorithm,timeControl,metric="RMSE",maximize=FALSE){
  return (caret::train(form,data=tserieDF,method=algorithm,trControl=timeControl,metric=metric,maximize=maximize))
}

#sets custom metrics to the data mining model (RMSE, MAE, MAPE)
modl.sumFunction <- function(data,lev = NULL,model = NULL){
  rmse <- caret::RMSE(data$obs,data$pred)
  mae <- Metrics::mae(data$obs,data$pred)
  mape <- TSPred::MAPE(data$obs,data$pred)

  output <- c(rmse,mae,mape)
  names(output) <- c("RMSE","MAE","MAPE")

  return (output)
}

##generic functions

#' Generic function
#'
#' Summary of object modl
#' @param object \code{prep} object
#' @param ... ignored
#' @export
#' @examples
#' summary(modl(prep(AirPassengers)))
summary.modl <- function(object,...){
  cat("Fitted model object\n")

  if(object$method=='arima') cat("~Modelling method: ARIMA\n")
  else if(object$method=='dataMining') cat ("~Modelling method: Data mining\n~Algorithm: ",object$algorithm,"\n")

  cat("~Original time serie: \n")
  utils::str(object$tserie)

  cat("~Original time serie as data frame")
  utils::str(object$tserieDF)

  cat("~Model (best tuning parameters): \n")
  if(object$method=='dataMining') print(object$model$bestTune)
  else if(object$method=='arima') print(forecast::arimaorder(object$model))

  cat("~Error measures: \n")
  print(object$errors)
}

#' Generic function
#'
#' Prints object modl
#' @param x \code{prep} object
#' @param ... ignored
#' @export
#' @examples
#' print(modl(prep(AirPassengers)))
print.modl <- function(x,...){
  cat("Fitted model object\n")
  cat("Class: modl\n\n")
  cat("Attributes: \n\n")
  cat("$method: ",x$method,"\n\n")
  if(!is.null(x$algorithm)) cat("$algorithm: ",x$algorithm,"\n\n")
  cat("$tserie: \n")
  print(x$tserie)
  cat("\n")
  cat("$tserieDF: \n")
  utils::str(x$tserieDF)
  cat("\n")
  cat("$model: \n")
  print(x$model)
  cat("\n")
  cat("$errors: \n")
  print(x$errors)
}
