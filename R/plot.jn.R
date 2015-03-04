#' Plot Function For Johnson-Neyman
#'
#' Plot function for objects of class \code{"jn"}.
#'
#' @param x an object of class \code{"jn"}.
#' @param xlab a title for the x axis (character).
#' @param ylab a title for the y axis (character).
#' @param xlim coordinates range for x axis (numeric vector).
#' @param ylim coordinates range for y axis (numeric vector)
#' @param \dots Additional arguments (not supported yet).
#'
#' @return none
#'
#' @method plot jn
#' @S3method plot jn
#' @examples
#' \dontrun{
#' myModel <- lm('DV ~ IV + MOD', data=someData)
#' jnresults <- jn(myModel, dv='DV', iv='IV', mod='MOD')
#' plot(jnresults)
#' }
#' @rdname plot.jn
#' @export

plot.jn <- function(x, xlab='', ylab='', xlim=0, ylim=0, ...){
  if(is(x,'lm')) {
    x <- jn(x)
  }

  if(is(x,'jn')) {
    if(length(x$error) > 0){
      cat('Error: ',x$error)
    } else {
      #default auto scaling
      if(missing(ylim)){
        ylim=c(floor(min(x$llci)),ceiling(max(x$ulci)))
      }
      if(missing(xlim)){
        xlim=c(floor(min(x$x)),ceiling(max(x$x)))
      }
      plot(x=x$x, y=x$y, ylim=ylim, xlim=xlim, ylab="" ,xlab="", type="l")

      #default axis labels
      if(ylab==''){
        ylab=paste("Conditional Effect of", x$iv, "on", x$dv)
      }
      if(xlab==''){
        xlab=x$mod
      }

      cat(paste('Values of', x$mod, 'indicated by the shaded region\n'))
      srvect <- vector()
      srvect = rbind(srvect, unlist(lapply(x$rsll, '[[', 1)))
      srvect = rbind(srvect, unlist(lapply(x$rsul, '[[', 1)))
      rownames(srvect) <- c('Lower Bound:','Upper Bound:')
      print(srvect)
      for(i in 1:nrow(x$signintervals))
      {
        ypolycoords <- x$x[x$signintervals[i,1]:x$signintervals[i,2]]
        polygon(c(ypolycoords,rev(ypolycoords)),c(x$llci[x$signintervals[i,1]:x$signintervals[i,2]],x$ulci[x$signintervals[i,2]:x$signintervals[i,1]]),col = "grey50", density = c(10, 20), angle = c(-45, 45), border = FALSE)
      }
      lines(x=x$x,y=x$y,lwd=10)
      lines(x=x$x, y=x$ulci, col="grey50",lty=2,lwd=8)
      lines(x=x$x, y=x$llci, col="grey50",lty=2,lwd=8)
      axis(side=1, lwd=10)
      mtext(side=1, xlab, line=2,font=2)
      axis(side=2, lwd=10)
      mtext(side=2, ylab, line=2,font=2)
      abline(a=0, b=0)
    }
  }
}
