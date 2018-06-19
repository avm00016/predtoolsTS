#' Predictions
#'
#' Performs predictions over a trained model.
#'
#' Predicts future values over a "modl" object which can be ARIMA or data mining, and returns
#' the predictions. Data mining predictions start right after the last value
#' contained in the training data, so they overlap with the end of the original.
#'
#' The object contains only two time series: the original one
#' and the predictions. You can just set these series aswell.
#'
#'
#' @param model A \code{modl} object. Contains the trained model we want to predict with.
#' @param n.ahead Number of values to predict ahead of the end of the original time serie. Default value is 20. Must ve lower than 100.
#' @param tserie A \code{ts} object.
#' @param predictions A \code{ts} object.
#' @return A list is returned of class \code{pred} containing:
#' \item{tserie}{Original time serie.}
#' \item{predictions}{Time serie with the predictions.}
#' @author Alberto Vico Moreno
#' @seealso{
#' \code{\link{modl}}
#' \code{\link{pred.arima}},
#' \code{\link{pred.dataMining}},
#' \code{\link{pred.compareModels}}
#' }
#' @export
#' @examples
#' prediction <- pred(model=modl(prep(AirPassengers)),n.ahead=25)
#' pred(tserie=prediction$tserie, predictions=prediction$predictions)
pred <- function(model=NULL,n.ahead=20,tserie=NULL,predictions=NULL){
  if(!is.null(model)){
    if(!methods::is(model,"modl")) stop("Not a modl object")

    predictions <- NULL

    if(model$method=='arima') predictions <- pred.arima(model$model,n.ahead=n.ahead)

    else if(model$method=='dataMining') predictions <- pred.dataMining(model,n.ahead=n.ahead)

    #creating the object
    obj <- list()
    class(obj) <- "pred"
    obj$tserie <- model$tserie
    obj$predictions <- predictions

    return (obj)
  }else{
    if(!stats::is.ts(tserie)) stop('Tserie not a ts object')
    if(!stats::is.ts(predictions)) stop('Predictions not a ts object')

    obj <- list()
    class(obj) <- "pred"
    obj$tserie <- tserie
    obj$predictions <- predictions

    return (obj)
  }
}

#' Predicts for ARIMA
#'
#' Performs predictions over an ARIMA model using the \code{stats::predict} function.
#'
#' @param model An ARIMA model.
#' @param n.ahead Number of values to predict.
#' @return A \code{ts} object containing the predictions.
#' @author Alberto Vico Moreno
#' @export
#' @examples
#' pred.arima(forecast::auto.arima(prep(AirPassengers)$tserie),n.ahead=30)
pred.arima <- function(model,n.ahead){
  if(!methods::is(model,"Arima")) stop("Not an Arima object")
  if(n.ahead > 100) stop("n.ahead must be lower than 100")

  return (stats::predict(model,n.ahead=n.ahead,se.fit=FALSE))
}

#' Predicts for data mining methods
#'
#' Performs predictions over a data mining model using the \code{caret::predict.train} function.
#'
#' @param model A \code{modl} object.
#' @param n.ahead Number of values to predict.
#' @return A \code{ts} object containing the predictions.
#' @author Alberto Vico Moreno
#' @export
#' @examples
#' m <- modl(prep(AirPassengers),method='dataMining',algorithm='rpart')
#' pred.dataMining(m,n.ahead=15)
pred.dataMining <- function(model,n.ahead){
  if(!methods::is(model,"modl")) stop("Not a modl object")
  if(model$method != "dataMining") stop("model method has to be dataMining.")
  if(n.ahead > 100) stop("n.ahead must be lower than 100")

  data <- model$tserieDF #original tserie as data frame
  index <- length(data[,1]) - model$horizon #index of the last observation contained into the training data

  obs <- data[index,] #last possible observation from the real serie

  predictions <- NULL

  for(i in 1:(n.ahead+model$horizon)){
    for(j in 1:(ncol(obs))){
      colnames(obs)[j] <- paste0('c_',j) #naming the features
    }
    predictions[i] <- caret::predict.train(model$model,obs) #store the prediction
    obs <- data.frame(c(obs[nrow(obs),2:ncol(obs)],predictions[i])) #calculate the new observation with the last predicted value
  }

  predictionsTS <- stats::ts(predictions,start=end(model$tserie)-c(0,model$horizon),frequency=frequency(model$tserie)) #time serie for the predictions

  return (predictionsTS)
}

#' Compare different predictions
#'
#' Plots the original time serie along with 2-5 predictive models.
#'
#' This function aims to ease the comparation between different predictive models
#' by plotting them into the same graphic.
#'
#' @param originalTS A \code{ts} object
#' @param p_1 A \code{ts} object
#' @param p_2 A \code{ts} object
#' @param p_3 A \code{ts} object. Default is NULL.
#' @param p_4 A \code{ts} object. Default is NULL.
#' @param p_5 A \code{ts} object. Default is NULL.
#' @param legendNames String vector with the names for the legend. Has to be same length as number of time series we are plotting(including the original one). Default is NULL.
#' @param colors Vector with the colors. Has to be same length as number of time series we are plotting(including the original one). Default is NULL.
#' @param legend A logical. Do we want a legend? Default is TRUE.
#' @param legendPosition A string with the position of the legend (bottomright, topright, ...). Default is NULL.
#' @param yAxis A string. Name for the y axis. "Values" as default.
#' @param title A string. Title for the plot. Default is "Predictions".
#' @author Alberto Vico Moreno
#' @import stats
#' @export
#' @examples
#' \donttest{
#' data(AirPassengers)
#' #pre-processing
#' p <- prep(AirPassengers)
#' #modelling
#' arima.modl <- modl(p)
#' cart.modl <- modl(p,method='dataMining',algorithm='rpart')
#' #predicting
#' arima.pred <- pred(arima.modl,n.ahead=30)
#' cart.pred <- pred(cart.modl,n.ahead=45)
#' #post-processing
#' arima.pred <- postp(arima.pred,p)
#' cart.pred <- postp(cart.pred,p)
#' #visual comparison
#' pred.compareModels(AirPassengers,arima.pred$predictions,cart.pred$predictions
#' ,legendNames=c('AirPassengers','ARIMA','CART'),yAxis='Passengers',legendPosition = 'topleft')
#' }
pred.compareModels <- function(originalTS,p_1,p_2,p_3=NULL,p_4=NULL,p_5=NULL,
                               legendNames=NULL,colors=NULL,legend=TRUE,legendPosition=NULL,yAxis="Values",title="Predictions"){

  if(!is.ts(originalTS) || !is.ts(p_1) || !is.ts(p_2)) stop('Not a ts object')

  n <- 2 #number of prediction time series
  maxX <- max(end(p_1)[1],end(p_2)[1]) #maximum value of x axis (year)
  if(!is.null(p_3)) {
    if(!is.ts(p_3)) stop('Not a ts object')
    if(end(p_3)[1] > maxX) maxX <- end(p_3)[1]
    n <- n+1
  }
  if(!is.null(p_4)) {
    if(!is.ts(p_4)) stop('Not a ts object')
    if(end(p_4)[1] > maxX) maxX <- end(p_4)[1]
    n <- n+1
  }
  if(!is.null(p_5)) {
    if(!is.ts(p_5)) stop('Not a ts object')
    if(end(p_5)[1] > maxX) maxX <- end(p_5)[1]
    n <- n+1
  }
  maxX <- maxX+1 #set the maximum a year after the last detected

  minY <- min(originalTS,p_1,p_2,p_3,p_4,p_5) #minimum value of y axis
  maxY <- max(originalTS,p_1,p_2,p_3,p_4,p_5) #maximum value of y axis

  xlim <- c(start(originalTS)[1],maxX)
  ylim <- c(minY,maxY)

  if(is.null(colors)){ #set colors
    colors <- c('black','green','blue','darkgoldenrod1','aquamarine3','darkorchid' )
    colors <- c(colors[1:(n+1)])
  }else if(length(colors) != (n+1)) stop('Vector "colors" wrong size. Has to contain colors for every time serie (including the original one)')

  if(is.null(legendNames)){ #set legend names
    legendNames <- c('Original TS','Predictions','Predictions','Predictions','Predictions','Predictions')
    legendNames <- c(legendNames[1:(n+1)])
  }else if(length(legendNames) != (n+1)) stop('Vector "legendNames" wrong size. Has to contain legend names for every time serie (including the original one)')

  graphics::plot(originalTS,xlim=xlim,ylim=ylim,col=colors[1],ylab=yAxis) #plot original ts
  title(main=paste0(title))

  #draw the predictions
  graphics::lines(p_1,col=colors[2])
  graphics::lines(p_2,col=colors[3])
  colorCounter <- 4
  if(!is.null(p_3)) {
    graphics::lines(p_3,col=colors[colorCounter])
    colorCounter <- colorCounter+1
  }
  if(!is.null(p_4)) {
    graphics::lines(p_4,col=colors[colorCounter])
    colorCounter <- colorCounter+1
  }
  if(!is.null(p_5)) {
    graphics::lines(p_5,col=colors[colorCounter])
    colorCounter <- colorCounter+1
  }

  #draw the legend
  if(legend){
    if(is.null(legendPosition))
      legendPosition <- "bottomright"
    legend(legendPosition,lty=c(1,1),col=colors,legend=legendNames)
  }
}


#generic functions

#' Generic function
#'
#' Plots object prep
#' @param x \code{pred} object
#' @param ylab ylab
#' @param main main
#' @param ... ignored
#' @export
#' @examples
#' plot(pred(modl(prep(AirPassengers))))
plot.pred <- function(x,ylab="Values",main="Predictions",...){

  maxX <- end(x$predictions) #maximum value of x axis
  minY <- min(x$tserie,x$predictions) #minimum value of y axis
  maxY <- max(x$tserie,x$predictions) #maximum value of y axis

  xlim <- c(start(x$tserie)[1],maxX[1])
  ylim <- c(minY,maxY)

  colors <- c('black','green')

  graphics::plot(x$tserie,xlim=xlim,ylim=ylim,col=colors[1],ylab=ylab,main=main) #plot original ts

  graphics::lines(x$predictions,col=colors[2]) #draw the predictions
}

#' Generic function
#'
#' Summary of object pred
#' @param object \code{prep} object
#' @param ... ignored
#' @export
#' @examples
#' summary(pred(modl(prep(AirPassengers))))
summary.pred <- function(object,...){
  cat("Predicted time serie object\n\n")
  cat("~Original time serie:\n")
  utils::str(object$tserie)
  cat("\n")
  cat("~Predictions (n=",length(object$predictions),"):\n")
  utils::str(object$predictions)
}

#' Generic function
#'
#' Prints object pred
#' @param x \code{prep} object
#' @param ... ignored
#' @export
#' @examples
#' print(pred(modl(prep(AirPassengers))))
print.pred <- function(x,...){
  cat("Predicted time serie object\n\n")
  cat("Class: pred\n\n")
  cat("Attributes: \n")
  cat("$tserie: \n")
  print(x$tserie)
  cat("\n")
  cat("$predictions: \n")
  print(x$predictions)
}
