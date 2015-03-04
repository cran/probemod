#' Pick-A-Point Technique
#'
#' Probe moderation effect using the Pick-A-Point technique
#'
#' @param model regression model (lm).
#' @param dv dependent variable (character).
#' @param iv independent variable (character).
#' @param mod moderator variable(s) (character or character vector).
#' @param points list of points to test for each moderators (list).
#' @param method method to use. Possible values are: \code{"meansd", "percentiles"}, \code{centered="meansd"} by default.
#' @param alpha alpha level to use (numeric).
#'
#' @return A list with the elements
#'
#' @examples
#' \dontrun{
#' myModel <- lm('dv ~ iv + mod', data=someData)
#' pickapoint(myModel, dv='DV', iv='IV', mod='MOD')
#' pickapoint(myModel, dv='DV', iv='IV', mod='MOD', alpha=.01)
#' pickapoint(myModel, dv='DV', iv='IV', mod='MOD', method='percentiles')
#' pickapoint(myModel, dv='DV', iv='IV', mod='MOD', points=c(1,2,3))
#' }
#' @rdname pickapoint
#'
#' @references
#' Spiller, S. A., Fitzsimons, G. J., Lynch, J. G., Jr, & McClelland, G. H. (2013). Spotlights, floodlights, and the magic number zero: Simple effects tests in moderated regression. Journal of Marketing Research, 50(2), 277-288.
#'
#' Aiken, L. S., & West, S. G. (1991). Multiple regression: Testing and interpreting interactions. Thousand Oaks, CA: Sage Publications.
#'
#' @export

pickapoint <- function(model,dv,iv,mod,points,method,alpha) UseMethod("pickapoint")
pickapoint <- function(model,dv,iv,mod,points,method="meansd",alpha=.05){
  #run checks on params
  if(!is(model,'lm')) {
    stop('this method currently supports only the lm object')
  }
  if(!is.character(dv) || !is.character(iv) || !is.character(mod)){
    stop('params dv/iv/mod incorrectly specified')
  }
  if(!missing(points)) {
    if(!is(points,'list')){
      stop('param points has to be a list')
    }
  } else{
    points <- list()
  }
  if(!(method == 'meansd' || method == 'percentiles')) {
    stop('param method has to be either "meansd" or "percentiles"')
  }
  if(!is.numeric(alpha)){
    stop('param alpha has to be numeric')
  }

  qrange <- c(.10,.25,.50,.75,.90)
  papret <- list()

  data <- model$model
  papret$method <- method
  papret$iv <- iv
  papret$dv <- dv
  papret$mod <- mod
  beta.hat <- coef(model)
  tcrit=qt(p=alpha/2,df=model$df,lower.tail = FALSE)
  cov <- vcov(model)
  modpprefvals <- vector('list')
  modppmatrix <- vector('list')

  interactionterms <- vector()
  for (i in 1:length(mod)){
    if(length(points[[mod[i]]]) > 0) {
      modpprefvals[[mod[i]]] <- points[[mod[i]]]
    } else if(length(unique(data[[mod[i]]])) == 2){
      #binary treatment for binary variables
      modpprefvals[[mod[i]]] <- unique(data[[mod[i]]])
    } else if (method == 'percentiles'){
      #calculate values for different percentiles
      modpprefvals[[mod[i]]] <- quantile(data[[mod[i]]], qrange)
    } else{
      #calculate values for mean +/-1sd
      modpprefvals[[mod[i]]] <- c(mean(data[[mod[i]]]) - sd(data[[mod[i]]]), mean(data[[mod[i]]]), mean(data[[mod[i]]]) + sd(data[[mod[i]]]) )
    }
    modppmatrix[[mod[i]]] <- c(1:length(modpprefvals[[mod[i]]]))
    interactionterms <- c(interactionterms, paste(papret$iv,papret$mod[i],sep=':'))
  }
  modpp <- expand.grid(modppmatrix)
  modvals <- vector()
  for(i in 1:ncol(modpp)){
    modvals <- cbind(modvals,modpprefvals[[i]][modpp[,i]])
  }

  papret$outcome <- vector()
  for(i in 1:nrow(modvals)){
    cy <- beta.hat[papret$iv]
    cyse <- cov[papret$iv, papret$iv]
    for(j in 1:ncol(modvals))
    {
      cy <- cy + (beta.hat[interactionterms[j]] * modvals[i,j])
      cyse <- cyse + (modvals[i,j]^2*cov[interactionterms[j], interactionterms[j]] + 2*modvals[i,j]*cov[papret$iv, interactionterms[j]])
    }
    cyse <- sqrt(cyse)
    ct <- cy/cyse
    cp <- 2*pt(-abs(ct),df=model$df)
    culci <- cy + tcrit*cyse
    cllci <- cy - tcrit*cyse
    papret$outcome <- rbind(papret$outcome,c(modvals[i,],cy,cyse,ct,cp,cllci,culci))
  }

  colnames(papret$outcome) <- c(mod,'Effect','SE','t','p','llci','ulci')
  rownames(papret$outcome) <- rep('',nrow(papret$outcome))
  papret$call <- match.call()
  class(papret) <- "pickapoint"
  papret
}

