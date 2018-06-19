#' Automatic pre-preprocessing
#'
#' This function performs pre-processing on a time series object(ts) to treat
#' heterocedasticity, trend and seasonality in order to make the serie stationary.
#'
#' Returns an object \code{prep} which stores all data needed to undo the changes later on.
#'
#' This function provides an automatic way of pre-processing based on unit root tests, but
#' this is not the perfect way to do it. You should always check manually if
#' the given time serie is actually stationary, and modify the parameters according
#' to your thoughts.
#'
#' @param tserie A ts object.
#' @param homogenize.method A string. Current methods available are "log" and "boxcox". Method "log" is set as default. If you don't want to perform this transformation, set method as "none".
#' @param detrend.method A string. Current methods available are "differencing" and "sfsm". Method "differencing" is set as default. If you don't want to perform this transformation, set method as "none".
#' @param nd A number. Number of differences you want to apply to the "differencing" detrending method. As default its value is NULL, which means nd will be calculated internally.
#' @param deseason.method A string. Current methods available are "differencing". Method "differencing" is set as default. If you don't want to perform this transformation, set method as "none".
#' @param nsd A number. Number of seasonal differences you want to apply to the "differencing" deseasoning method. As default its value is NULL, which means nsd will be calculated internally.
#' @param detrend.first A boolean. TRUE if detrending method is applied first, then deseasoning. FALSE if deseasoning method is applied first. Default is TRUE.
#' @return A list is returned of class \code{prep} containing:
#' \item{tserie}{Processed ts object.}
#' \item{homogenize.method}{Method used for homogenizing.}
#' \item{detrend.method}{Method used for detrending.}
#' \item{nd}{Number of differences used on detrending through differencing.}
#' \item{firstvalues}{First \code{nd} values of the original series.}
#' \item{deseason.method}{Method used for deseasoning.}
#' \item{nsd}{Number of seasonal differences used on deseasoning through differencing.}
#' \item{firstseasons}{First \code{nsd} seasons of the original series.}
#' \item{detrend.first}{Processed ts object}
#' \item{means}{Vector of means used in "sfsm" detrending method.}
#' \item{lambda}{Coefficient used in "boxcox" transformation.}
#' \item{start}{Start of the original time serie.}
#' \item{length}{Length of the original time serie.}
#' @author Alberto Vico Moreno
#' @seealso{
#' \code{\link{prep.homogenize.log}},
#' \code{\link{prep.homogenize.boxcox}},
#' \code{\link{prep.detrend.differencing}},
#' \code{\link{prep.detrend.sfsm}},
#' \code{\link{prep.deseason.differencing}},
#' \code{\link{prep.check.acf}},
#' \code{\link{prep.check.adf}}
#' }
#' @export
#' @references \url{https://www.otexts.org/fpp/8/1}
#' @examples
#' prep(AirPassengers)
#' prep(AirPassengers,homogenize.method='boxcox',detrend.method='none')
prep <- function(tserie,homogenize.method='log'
                      ,detrend.method='differencing',nd=NULL
                      ,deseason.method='differencing',nsd=NULL
                      ,detrend.first=TRUE){
  if(!stats::is.ts(tserie)) stop('Not a ts object')
  if(stats::frequency(tserie) == 1) deseason.method='none'

  newts <- tserie
  means <- NULL
  lambda <- NULL
  firstvalues <- NULL
  firstseasons <- NULL
  start <- start(tserie)
  length <- length(tserie)

  #homogenizing
  if(homogenize.method=='log') newts <- prep.homogenize.log(newts)
  else if(homogenize.method=='boxcox'){
    bc <- prep.homogenize.boxcox(newts)
    newts <- bc[[1]]
    lambda <- bc[[2]]
  }
  else if(homogenize.method!='none') stop('Invalid homogenizing method')

  #detrending and deseasoning
  if(detrend.first==FALSE){

    if(deseason.method=='differencing'){
      dsts <- prep.deseason.differencing(newts,nsd)
      if(dsts[[2]]>0){
        newts <- dsts[[1]]
        nsd <- dsts[[2]]
        firstseasons <- dsts[[3]]
      }else deseason.method='none'
    }else if(deseason.method!='none') stop('Invalid deseasoning method')

    if(detrend.method=='differencing'){
      dtts <- prep.detrend.differencing(newts,nd)
      if(dtts[[2]]>0){
        newts <- dtts[[1]]
        nd <- dtts[[2]]
        firstvalues <- dtts[[3]]
      }else detrend.method='none'

    }else if(detrend.method=='sfsm'){
      dtts <- prep.detrend.sfsm(newts)
      newts <- dtts[[1]]
      means <- dtts[[2]]
    }else if(detrend.method!='none') stop('Invalid detrending method')

  }else{

    if(detrend.method=='differencing'){
      dtts <- prep.detrend.differencing(newts,nd)
      if(dtts[[2]]>0){
        newts <- dtts[[1]]
        nd <- dtts[[2]]
        firstvalues <- dtts[[3]]
      }else detrend.method='none'

    }else if(detrend.method=='sfsm'){
      dtts <- prep.detrend.sfsm(newts)
      newts <- dtts[[1]]
      means <- dtts[[2]]
    }else if(detrend.method!='none') stop('Invalid detrending method')

    if(deseason.method=='differencing'){
      dsts <- prep.deseason.differencing(newts,nsd)
      if(dsts[[2]]>0){
        newts <- dsts[[1]]
        nsd <- dsts[[2]]
        firstseasons <- dsts[[3]]
      }else deseason.method='none'
    }else if(deseason.method!='none') stop('Invalid deseasoning method')
  }

  #creating the object
  obj <- list()
  class(obj) <- "prep"
  obj$tserie <- newts
  obj$homogenize.method <- homogenize.method
  obj$detrend.method <- detrend.method
  obj$nd <- nd
  obj$firstvalues <- firstvalues
  obj$deseason.method <- deseason.method
  obj$nsd <- nsd
  obj$firstseasons <- firstseasons
  obj$detrend.first <- detrend.first
  obj$means <- means
  obj$lambda <- lambda
  obj$start <- start
  obj$length <- length

  return (obj)
}

#' Logarithmic transformation
#'
#' Performs a logarithmic transformation to a time serie.
#'
#' @param tserie a \code{ts} object
#' @return \code{ts} object with transformed time serie
#' @export
#' @examples
#' prep.homogenize.log(AirPassengers)
prep.homogenize.log <- function(tserie){
  if(!stats::is.ts(tserie)) stop('Not a ts object')

  return (log(tserie))
}

#' Box-Cox transformation
#'
#' Performs a Box-Cox transformation to a time serie.
#'
#' @param tserie a \code{ts} object
#' @return A list is returned containing:
#' \item{boxcox}{Transformed ts object.}
#' \item{lambda}{Lambda value.}
#' @references Box-Cox transformation: \url{https://en.wikipedia.org/wiki/Power_transform#Box.E2.80.93Cox_transformation}
#' @export
#' @examples
#' prep.homogenize.log(AirPassengers)
prep.homogenize.boxcox <- function(tserie){
  if(!stats::is.ts(tserie)) stop('Not a ts object')
  lambda <- forecast::BoxCox.lambda(tserie)
  return (list(boxcox=forecast::BoxCox(tserie,lambda),lambda=lambda))
}

#' Detrend with differencing method
#'
#' Performs differencing with lag=1.
#'
#' If no number of differences is specified, the function will make an estimation
#' of the number of differences needed based on unit root test provided by \code{forecast::ndiffs}
#'
#' @param tserie a \code{ts} object
#' @param nd number of differences to apply. As default its value is NULL; in this case, the function will perform an automatic estimation of \code{nd}.
#' @return A list is returned containing:
#' \item{tserie}{Transformed ts object.}
#' \item{nd}{Number of differencies applied.}
#' \item{firstvalues}{Lost values after differencing.}
#' @export
#' @examples
#' prep.detrend.differencing(AirPassengers)
#' prep.detrend.differencing(AirPassengers,nd=2)
prep.detrend.differencing <- function(tserie,nd=NULL){
  if(!stats::is.ts(tserie)) stop('Not a ts object')

  firstvalues <- NULL

  if(is.null(nd)){
    nd <- forecast::ndiffs(tserie)
    if(nd > 0){
      firstvalues <- tserie[1:nd]
      tserie <- diff(tserie,differences=nd)
    }
    return (list(tserie=tserie,nd=nd,firstvalues=firstvalues))
  }else{
    if(nd > 0) firstvalues <- tserie[1:nd]
    return (list(tserie=tserie <- diff(tserie,differences=nd),nd=nd,firsvalues=firstvalues))
  }
}

#' Detrend with "substracting full-season means" method
#'
#' Performs "substracting full-season means" method to go for a totally automatic
#' approach.
#'
#' Under this detrending scheme, a series is first split into segments. The length
#' of the segments is equal to the length of seasonality(12 for monthly).
#' The mean of the historical observations within each of these segments is substacted
#' from every historical observation in the segment.
#' To get the detrended serie we do:
#' \code{ds = xi - m}
#' Being \code{xi} the actual values on the time series and \code{m} the mean of the segment of \code{xi}
#'
#' @param tserie a \code{ts} object
#' @return A list is returned containing:
#' \item{tserie}{Transformed ts object.}
#' \item{means}{Vector containing the historical means.}
#' @export
#' @examples
#' prep.detrend.sfsm(AirPassengers)
prep.detrend.sfsm <- function(tserie){
  if(!stats::is.ts(tserie)) stop('Not a ts object')

  first <- start(tserie)[2] #first observation index of the time serie
  index <- 1 #general index
  means <- NULL #means vector
  cont <- TRUE

  while(cont){
    aux <- NULL

    while(index <= length(tserie)){ #values of a season
      aux <- c(aux,tserie[index])
      index <- index + 1
      if((index+first-2) %% frequency(tserie) == 0) break
    }

    means <- c(means,mean(aux)) #create means vector

    if(index > length(tserie)) cont <- FALSE
  }

  mindex <- 1
  index <- 1

  while(index <= length(tserie)){
    tserie[index] <- tserie[index]-means[mindex]
    index <- index+1
    if((index+first-2) %% frequency(tserie) == 0) mindex <- mindex+1
  }

  return (list(tserie=tserie,means=means))
}

#' Deseason with differencing method
#'
#' Performs differencing with lag=frequency.
#'
#' If no number of differences is specified, the function will make an estimation
#' of the number of differences needed based on unit root test provided by \code{forecast::nsdiffs}
#'
#' @param tserie a \code{ts} object
#' @param nsd number of seasonal differences to apply. As default its value is NULL; in this case, the function will perform an automatic estimation of \code{nsd}.
#' @return A list is returned containing:
#' \item{tserie}{Transformed ts object.}
#' \item{nsd}{Number of seasonal differencies applied.}
#' \item{firstseasons}{Lost values after differencing.}
#' @export
#' @examples
#' prep.deseason.differencing(AirPassengers)
#' prep.deseason.differencing(AirPassengers,nsd=2)
prep.deseason.differencing <- function(tserie, nsd=NULL){
  if(!stats::is.ts(tserie)) stop('Not a ts object')

  firstseasons <- NULL

  if(is.null(nsd)){
    nsd <- forecast::nsdiffs(tserie)
    if(nsd > 0){
      firstseasons <- tserie[1:(frequency(tserie)*nsd)]
      tserie <- diff(tserie,lag=frequency(tserie),differences=nsd)
    }
    return (list(tserie=tserie,nsd=nsd,firstseasons=firstseasons))
  }else{
    if(nsd > 0) firstseasons <- tserie[1:(frequency(tserie)*nsd)]
    return (list(tserie=diff(tserie,lag=frequency(tserie),differences=nsd),nsd=nsd,firstseasons=firstseasons))
  }
}


#' Autocorrelation function
#'
#' Plots the autocorrelation function to check stationarity
#'
#' For a stationary time series, the ACF will drop to zero
#' relatively quickly, while the ACF of non-stationary data decreases slowly.
#' Also, for non-stationary data, the value is often large and positive.
#' @param tserie a \code{ts} or a \code{prep} object
#' @export
#' @examples
#' prep.check.acf(AirPassengers)
#' prep.check.acf(prep(AirPassengers))
prep.check.acf <- function(tserie){
  if(stats::is.ts(tserie)) stats::acf(tserie)
  else if(class(tserie)=='prep') stats::acf(tserie$tserie)
  else stop('Not a ts/prep object')
}

#' Augmented Dickey-Fuller test
#'
#' Performs ADF test just as another tool to check stationarity.
#'
#' Shows the results of an ADF test. A p-value<0.05 suggests the data is stationary.
#'
#' @param tserie a \code{ts} or a \code{prep} object
#' @export
#' @examples
#' prep.check.adf(AirPassengers)
#' prep.check.adf(prep(AirPassengers))
prep.check.adf <- function(tserie){
  if(stats::is.ts(tserie)) return (tseries::adf.test(tserie,alternative="stationary"))
  else if(class(tserie)=='prep') return (tseries::adf.test(tserie$tserie,alternative="stationary"))
  else stop('Not a ts/prep object')
}

#generic functions

#' Generic function
#'
#' Plots object prep
#' @param x \code{prep} object
#' @param ylab ylab
#' @param xlab xlab
#' @param ... ignored
#' @export
#' @examples
#' plot(prep(AirPassengers),ylab="Stationary AisPassengers")
plot.prep <- function(x,ylab="Preprocessed time serie",xlab="",...){
  graphics::plot(x$tserie,ylab=ylab,xlab=xlab)
}

#' Generic function
#'
#' Summary of object prep
#' @param object \code{prep} object
#' @param ... ignored
#' @export
#' @examples
#' summary(prep(AirPassengers))
summary.prep <- function(object,...){
  cat("Preprocessed time series object\n")

  if(object$homogenize.method=='log') cat("~Homogenizing method: logarithmic transformation\n")
  else if(object$homogenize.method=='boxcox') cat ("~Homogenizing method: Box-Cox transformation with lambda=",object$lambda,"\n")
  else cat("~Transformation not applied\n")

  if(object$detrend.method=='differencing') {
    cat("~Detrending method: differencing\n Number of differences: ",object$nd,"\n")
    cat("First original values: ",object$firstvalues,"\n")
  }
  else if(object$detrend.method=='sfsm'){
    cat("~Detrending method: substracting full-season means\n Season means: \n")
    utils::str(object$means)
  }
  else cat("~No detrending performed\n")

  if(object$deseason.method=='differencing') {
    cat("~Deseason method: differencing\n Number of seasonal differences: ",object$nsd,"\n")
    cat("First original seasons: ",object$firstseasons,"\n")
  }
  else cat("~No deseason performed\n")

  if(object$deseason.method!='none' && object$detrend.method!='none'){
    if(object$detrend.first==TRUE) cat("~Detrending applied before deseasoning\n")
    else cat("~Detrending applied after deseasoning\n")
  }

  cat("~Original serie start: ",object$start,"\n")
  cat("~Original serie length: ",object$length,"\n")

  cat("~Preprocessed time serie:\n")
  utils::str(object$tserie)
}

#' Generic function
#'
#' Prints object prep
#' @param x \code{prep} object
#' @param ... ignored
#' @export
#' @examples
#' print(prep(AirPassengers))
print.prep <- function(x,...){
  cat("Preprocessed time series object\n\n")
  cat("Class: prep\n\n")
  cat("Attributes: \n")
  cat("$homogenize.method: ",x$homogenize.method,"\n")
  if(!is.null(x$lambda)) cat("$lambda: ",x$lambda,"\n")
  cat("$detrend.method: ",x$detrend.method,"\n")
  if(!is.null(x$nd)) cat("$nd: ",x$nd,"\n")
  if(!is.null(x$firstvalues)) cat("$firstvalues: ",x$firstvalues,"\n")
  if(!is.null(x$means)) cat("$means: ",x$means,"\n")
  cat("$deseason.method: ",x$deseason.method,"\n")
  if(!is.null(x$nsd)) cat("$nsd: ",x$nsd,"\n")
  if(!is.null(x$firstseasons)) cat("$firstseasons: ",x$firstseasons,"\n")
  cat("$detrend.first: ",x$detrend.first,"\n")
  cat("$start: ",x$start,"\n")
  cat("$length: ",x$length,"\n")
  cat("$tserie: \n")
  print(x$tserie)
}
