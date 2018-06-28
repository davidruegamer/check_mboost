#########################################################################
########### helper functions ############################################
#########################################################################

# calculate cumsum per factor 
cumsum_grouped <- function(x, fac) {
  # x: numeric vector to be summed over
  # fac: factor of the same length as x indicating the grouping sructure for 
  #       which cumsum is computed sperately. If fac is not a factor, as.factor(fac) is used.
  stopifnot(length(x) == length(fac))
  cusu <- numeric(length = length(x))
  cusu[order(fac)] <- unlist( tapply(x, fac, cumsum) )
  return(cusu)
}

#########################################################################
########### functions which return vector ###############################
#########################################################################

# edf definition particularly used for additive models
edf2 <- function(trhat, trhatsq) 2 * trhat - trhatsq

# corrected AICs for both edf versions
aiccor_edf1 <- function(object, trhat) 
  log(object$risk()[-1]/ 
        (sum(object$`(weights)`[!is.na(fitted(object))]))) + 
  (1 + trhat/(sum(object$`(weights)`[!is.na(fitted(object))]))) / 
  (1 - (trhat + 2)/(sum(object$`(weights)`[!is.na(fitted(object))])))
aiccor_edf2 <- function(object, trhat, trhatsq) log(object$risk()[-1]/
                                               (sum(object$`(weights)`[!is.na(fitted(object))]))) + 
  (1 + edf2(trhat, trhatsq)/(sum(object$`(weights)`[!is.na(fitted(object))])))/
  (1 - (edf2(trhat, trhatsq) + 2)/(sum(object$`(weights)`[!is.na(fitted(object))])))

# classical AICs for both edf versions
aicclas_edf1 <- function(object, trhat) 2 * object$risk()[-1] + 2 * trhat
aicclas_edf2 <- function(object, trhat, trhatsq) 2 * object$risk()[-1] + 
  2 * edf2(trhat, trhatsq)

# Classical BICs for both edf versions
bic_edf1 <- function(object, trhat) 2 * object$risk()[-1] + length(object$response) * trhat
bic_edf2 <- function(object, trhat, trhatsq) 2 * object$risk()[-1] + 
  length(object$response) * edf2(trhat, trhatsq)

# general minimum description length (gMDL) for both edf versions
gmdl_edf1 <- function(object, trhat){ 

  s_edf1 <- function(object) object$risk()[-1]/ 
    ((sum(object$`(weights)`[!is.na(fitted(object))])) - trhat)  
  
  log(s_edf1(object)) + trhat/(sum(object$`(weights)`[!is.na(fitted(object))])) * 
  log((sum(object$response^2) - object$risk()[-1])/(trhat * s_edf1(object)))
  
}
gmdl_edf2 <- function(object, trhat, trhatsq){ 
  
  s_edf2 <- function(object, trhat, trhatsq) 
    object$risk()[-1]/ ((sum(object$`(weights)`[!is.na(fitted(object))])) - 
                          edf2(trhat, trhatsq))
  
  log(s_edf2(object, trhat, trhatsq)) + 
  edf2(trhat, trhatsq)/(sum(object$`(weights)`[!is.na(fitted(object))])) * 
  log((sum(object$response^2) - object$risk()[-1])/
        (edf2(trhat, trhatsq) * s_edf2(object, trhat, trhatsq)))
  
}



# extract cummulated explained risk / learner
## note that the function relies on risk() and, thus, returns the inbag / 
##    out of bag risk depending on boost_control()
## note that in mboost the risk corresponds to the sum of the loss rather than its mean 
##    -> we use the mean
extract_cum_expl_risk <- function(object) {
  risk_diff <- - diff(risk(object))/length(object$response)
  learner_selected <- selected(object)
  explained_risk <- cumsum_grouped(risk_diff, learner_selected)
  attr( explained_risk, "initial_risk") <- risk(object)[1]/length(object$response)
  return( explained_risk )
}

# extract selection path
extract_sel_path <- function(object) {
  learner_selected <- selected(object)
  sels <- integer(length = length(learner_selected))
  sels[order(learner_selected)] <- unlist( sapply( as.numeric(table(learner_selected)), function(n) 1L:n ) )
  attr( sels, "mstop" ) <- mstop(object)
  return(sels)
}

