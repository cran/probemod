#' Print Function For Johnson-Neyman
#'
#' Prints function for objects of class \code{"jn"}
#'
#' @param x an object of class \code{"jn"}.
#' @param \dots Additional arguments (not supported yet).
#' @return none
#'
#' @examples
#' \dontrun{
#' myModel <- lm('DV ~ IV + MOD', data=someData)
#' jnresults <- jn(myModel, dv='DV', iv='IV', mod='MOD')
#' jnresults
#' }
#' @rdname print.jn
#' @export

print.jn <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("Conditional effects of ", x$iv, " on ", x$dv, " at values of ", x$mod ,"\n")
  showrect <- cbind(round(x$x, digits=4), round(x$y, digits=4), round(x$se, digits=4), round(x$t, digits=4), round(x$p, digits=4), round(x$llci, digits=4), round(x$ulci, digits=4))
  colnames(showrect) <- c(x$iv,'effect','se','t','p','llci','ulci')
  print(showrect)
}
