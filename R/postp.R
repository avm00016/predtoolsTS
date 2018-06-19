#' Post-processing of pre-processed data
#'
#' Using the \code{prep} data we undo the changes on a \code{pred} object.
#'
#' @param prd A \code{pred} object.
#' @param pre A \code{prep} object.
#' @return A \code{pred} object with reverted transformations.
#' @author Alberto Vico Moreno
#' @seealso{
#' \code{\link{pred}}
#' \code{\link{prep}},
#' \code{\link{postp.homogenize.log}},
#' \code{\link{postp.homogenize.boxcox}},
#' \code{\link{postp.detrend.differencing}},
#' \code{\link{postp.detrend.sfsm}},
#' \code{\link{postp.deseason.differencing}}
#' }
#' @import stats
#' @export
#' @examples
#' preprocess <- prep(AirPassengers)
#' prediction <- pred(modl(preprocess),n.ahead=30)
#' postp.prediction <- postp(prediction,preprocess)
postp <- function(prd,pre){
  if(!methods::is(prd,"pred")) stop("Not a pred object")
  if(!methods::is(pre,"prep")) stop("Not a prep object")

  u <- c(window(prd$tserie,end=(start(prd$predictions)-c(0,1))),prd$predictions) #union of both time series as a vector, leaving the test part out from the original time serie
  tserie <- prd$tserie #original time serie

  #detrend and deseason post-processing
  if(pre$detrend.first==TRUE){

    if(pre$deseason.method!='none')
      if(pre$deseason.method=='differencing'){
        u <- postp.deseason.differencing(u,nsd=pre$nsd,firstseasons=pre$firstseasons,frequency=frequency(prd$tserie))
        tserie <- postp.deseason.differencing(tserie,nsd=pre$nsd,firstseasons=pre$firstseasons,frequency=frequency(prd$tserie))
      }

    if(pre$detrend.method!='none')
      if(pre$detrend.method=='differencing'){
        u <- postp.detrend.differencing(u,nd=pre$nd,firstvalues=pre$firstvalues)
        tserie <- postp.detrend.differencing(tserie,nd=pre$nd,firstvalues=pre$firstvalues)
      }
      else if(pre$detrend.method=='sfsm'){
        u <- postp.detrend.sfsm(u,pre$means,start(prd$tserie),frequency(prd$tserie))
        tserie <- postp.detrend.sfsm(tserie,pre$means,start(prd$tserie),frequency(prd$tserie))

      }

  }else{

    if(pre$detrend.method!='none')
      if(pre$detrend.method=='differencing'){
        u <- postp.detrend.differencing(u,nd=pre$nd,firstvalues=pre$firstvalues)
        tserie <- postp.detrend.differencing(tserie,nd=pre$nd,firstvalues=pre$firstvalues)
      }
      else if(pre$detrend.method=='sfsm'){
        u <- postp.detrend.sfsm(u,pre$means,start(prd$tserie),frequency(prd$tserie))
        tserie <- postp.detrend.sfsm(tserie,pre$means,start(prd$tserie),frequency(prd$tserie))
      }
    if(pre$deseason.method!='none')
      if(pre$deseason.method=='differencing'){
        u <- postp.deseason.differencing(u,nsd=pre$nsd,firstseasons=pre$firstseasons,frequency=frequency(prd$tserie))
        tserie <- postp.deseason.differencing(tserie,nsd=pre$nsd,firstseasons=pre$firstseasons,frequency=frequency(prd$tserie))
      }
  }

  #homogenize post-processing
  if(pre$homogenize.method=='log'){
    u <- postp.homogenize.log(u)
    tserie <- postp.homogenize.log(tserie)
  }
  else if(pre$homogenize.method=='boxcox'){
    u <- postp.homogenize.boxcox(u,pre$lambda)
    tserie <- postp.homogenize.boxcox(tserie,pre$lambda)
  }

  predictions <- stats::window(ts(u,start=start(tserie),frequency=frequency(tserie)),start=start(prd$predictions))

  return (pred(tserie=tserie,predictions=predictions))
}

#' Undo logarithmic transformation
#'
#' Uses exponent to reverse the logarithm
#'
#' @param tserie A \code{ts} object.
#' @return A \code{ts} object.
#' @author Alberto Vico Moreno
#' @export
#' @examples
#' postp.homogenize.log(prep.homogenize.log(AirPassengers))
postp.homogenize.log <- function(tserie){
  return (exp(tserie))
}

#' Undo Box-Cox transformation
#'
#' @param tserie A \code{ts} object.
#' @param lambda A numeric.
#' @return A \code{ts} object.
#' @author Alberto Vico Moreno
#' @export
#' @examples
#' p <- prep.homogenize.boxcox(AirPassengers)
#' postp.homogenize.boxcox(p$tserie,p$lambda)
postp.homogenize.boxcox <- function(tserie,lambda){
  if(lambda==0) return (exp(tserie))
  else (lambda*tserie + 1)^(1/lambda)
}

#' Undo detrend(differencing)
#'
#' Uses inverse differences to revert the changes
#'
#' @param tserie A \code{ts} object.
#' @param nd Number of differences.
#' @param firstvalues Values lost on the original differences
#' @return A \code{ts} object.
#' @author Alberto Vico Moreno
#' @export
#' @examples
#' p <- prep.detrend.differencing(AirPassengers)
#' postp.detrend.differencing(p$tserie,p$nd,p$firstvalues)
postp.detrend.differencing <- function(tserie,nd,firstvalues){
  return (stats::diffinv(tserie,lag=1,differences=nd,firstvalues))
}

#' Undo detrend(substracting full-means method)
#'
#' @param tserie A \code{ts} object.
#' @param means A numeric vector.
#' @param start Start of original time serie
#' @param frequency Frequency of the original time serie
#' @return A \code{ts} object.
#' @author Alberto Vico Moreno
#' @export
#' @import stats
#' @examples
#' p <- prep.detrend.sfsm(AirPassengers)
#' postp.detrend.sfsm(p$tserie,p$means,start(AirPassengers),frequency(AirPassengers))
postp.detrend.sfsm <- function(tserie,means,start,frequency){
  first <- start[2] #first observation index of the time serie
  index <- 1 #general index
  mindex <- 1 #means index
  cont <- TRUE
  newts <- tserie

  while(cont){

    #start building new time serie
    for(i in 1:frequency){
      newts[i+index-1] <- tserie[i+index-1] + means[mindex]
      if(i+first-1 %% frequency == 0) {
        index<-index-(frequency-i)
        break
      }
    }

    mindex <- mindex+1
    index <- index+frequency

    if(mindex > length(means)) cont <- FALSE
  }

  cont <- TRUE

  if(index<length(tserie)){
    #here starts the predicted values aproximation
    while(cont){
      means[mindex] <- mean(means[mindex-1],means[mindex-2]) #new means are the average of the last two

      for(i in 1:frequency){
        newts[index] <- tserie[index]+means[mindex]
        index <- index+1
        if(index>length(tserie)){
          cont <- FALSE
          break
        }
      }
      mindex <- mindex+1
    }
  }

  return (newts)
}

#' Undo deseason(differencing)
#'
#' Uses inverse seasonal differences to reverse the changes
#'
#' @param tserie A \code{ts} object.
#' @param nsd Number of seasonal differences.
#' @param firstseasons Values lost on the original differences
#' @param frequency Frequency of the original time serie
#' @return A \code{ts} object.
#' @author Alberto Vico Moreno
#' @export
#' @examples
#' p <- prep.deseason.differencing(AirPassengers)
#' postp.deseason.differencing(p$tserie,p$nsd,p$firstseasons,frequency(AirPassengers))
postp.deseason.differencing <- function(tserie,nsd,firstseasons,frequency){
  return (stats::diffinv(tserie,lag=frequency,differences=nsd,firstseasons))
}

