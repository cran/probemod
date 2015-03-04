#' Johnson-Neyman Technique
#'
#' Probe moderation effect using the Johnson-Neyman technique
#'
#' @param model regression model (lm).
#' @param dv dependent variable (character).
#' @param iv independent variable (character).
#' @param mod moderator variable(s) (character or character vector).
#' @param alpha alpha level to use (numeric).
#'
#' @return A list with the elements
#'
#' @examples
#' \dontrun{
#' myModel <- lm('DV ~ IV + MOD', data=someData)
#' jnresults <- jn(myModel, dv='DV', iv='IV', mod='MOD')
#' jnresults <- jn(myModel, dv='DV', iv='IV', mod='MOD', alpha=.01)
#' plot(jnresults)
#' }
#' @rdname jn
#'
#' @references
#' Spiller, S. A., Fitzsimons, G. J., Lynch, J. G., Jr, & McClelland, G. H. (2013). Spotlights, floodlights, and the magic number zero: Simple effects tests in moderated regression. Journal of Marketing Research, 50(2), 277-288.
#'
#' Bauer, D. J., & Curran, P. J. (2005). Probing interactions in fixed and multilevel regression: Inferential and graphical techniques. Multivariate Behavioral Research, 40(3), 373-400.
#'
#' @export

jn <- function(model,dv,iv,mod,alpha) UseMethod("jn")
jn <- function(model,dv,iv,mod,alpha =.05){
  #run checks on params
  if(!is(model,'lm')) {
    stop('this method currently supports only the lm object')
  }
  if(!is.character(dv) || !is.character(iv) || !is.character(mod)){
    stop('params dv/iv/mod incorrectly specified')
  }
  if(!is.numeric(alpha)){
    stop('param alpha has to be numeric')
  }

  jnret <- list()
  if(is(model,'lm')) {
    data <- model$model
    jnret$iv <- iv
    jnret$dv <- dv
    jnret$mod <- mod
    jnret$n <- nrow(model$model)
    tcrit=qt(p=alpha/2,df=model$df,lower.tail = FALSE)
    if(length(jnret$mod) > 1){
      jnret$error <- 'The current version supports only one moderator'
    } else {
      beta.hat <- coef(model)
      cov <- vcov(model)
      jnret$x <- seq(min(data[[jnret$mod]]), max(data[[jnret$mod]]), length.out = 1000)
      interactionterm <- paste(jnret$iv,jnret$mod,sep=':')

      jna <- tcrit^2*cov[interactionterm, interactionterm] - beta.hat[interactionterm]^2
      jnb <- 2 * (tcrit^2*cov[jnret$iv, interactionterm] - beta.hat[iv]*beta.hat[interactionterm])
      jnc <- tcrit^2*cov[jnret$iv, jnret$iv] - beta.hat[iv]^2
      jnret$rsul <- list()
      jnret$rsll <- list()
      rs <- c( (-jnb - sqrt(jnb^2 - 4*jna*jnc)) / (2 * jna), (-jnb + sqrt(jnb^2 - 4*jna*jnc)) / (2 * jna))
      jnret$rsul$x <- max(rs)
      jnret$rsul$y <- beta.hat[jnret$iv] + beta.hat[interactionterm]*jnret$rsul$x
      jnret$rsul$yse <- sqrt(cov[jnret$iv, jnret$iv] + jnret$rsul$x^2*cov[interactionterm, interactionterm] + 2*jnret$rsul$x*cov[jnret$iv, interactionterm])
      jnret$rsul$t <- jnret$rsul$y/jnret$rsul$yse
      jnret$rsul$p <- round(2*pt(-abs(jnret$rsul$t),df=model$df),digits=4)
      jnret$rsll$x <- min(rs)
      jnret$rsll$y <- beta.hat[jnret$iv] + beta.hat[interactionterm]*jnret$rsll$x
      jnret$rsll$yse <- sqrt(cov[jnret$iv, jnret$iv] + jnret$rsll$x^2*cov[interactionterm, interactionterm] + 2*jnret$rsll$x*cov[jnret$iv, interactionterm])
      jnret$rsll$t <- jnret$rsul$y/jnret$rsll$yse
      jnret$rsll$p <- round(2*pt(-abs(jnret$rsll$t),df=model$df),digits=4)

      jnret$y <- beta.hat[jnret$iv] + beta.hat[interactionterm]*jnret$x
      jnret$yse <- sqrt(cov[jnret$iv, jnret$iv] + jnret$x^2*cov[interactionterm, interactionterm] + 2*jnret$x*cov[jnret$iv, interactionterm])
      jnret$t <- cbind(t=jnret$y/jnret$yse)
      jnret$p=round(2*pt(-abs(jnret$t),df=model$df),digits=4)
      jnret$ulci <- jnret$y + tcrit*jnret$yse
      jnret$llci <- jnret$y - tcrit*jnret$yse

      jnret$signintervals <- vector()
      trackstartpoly <- 0
      trackstartvalpoly <- vector()
      for(i in 1:length(jnret$x)) {
        if(jnret$ulci[i] >= 0 && jnret$llci[i] >= 0) {
          if(trackstartpoly == 0) {
            trackstartpoly <- 1
            trackstartvalpoly <- i
          }
        }
        else if(jnret$ulci[i] <= 0 && jnret$llci[i] <= 0) {
          if(trackstartpoly == 0) {
            trackstartpoly <- 1
            trackstartvalpoly <- i
          }
        }
        else if(trackstartpoly == 1){
          trackstartpoly <- 0
          jnret$signintervals <- rbind(jnret$signintervals,c(trackstartvalpoly,i-1))
          trackstartvalpoly <- vector()
        }
      }
      #wrap up poly calulations
      if(trackstartpoly == 1){
        trackstartpoly <- 0
        jnret$signintervals <- rbind(jnret$signintervals,c(trackstartvalpoly,i))
        trackstartvalspoly <- vector()
      }
    }
  }
  jnret$call <- match.call()
  class(jnret) <- "jn"
  jnret
}
