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
########### functions which return scalar (iteration specific) ##########
#########################################################################

# trace of the squared hat matrix
trhatsq <- function(object) 
  sum(diag(crossprod(attr(hatvalues(object), "hatmatrix"))))

#########################################################################
########### functions which return vector ###############################
#########################################################################

# original definition of edf
edf1 <- function(object) attr(hatvalues(object), "trace")
# edf definition particularly used for additive models
edf2 <- function(object, trhatsq) 2 * edf1(object) - trhatsq

# corrected AICs for both edf versions
aiccor_edf1 <- function(object) log(object$risk()[-1]/ 
                                      (sum(object$`(weights)`[!is.na(fitted(object))]))) + 
  (1 + edf1(object)/(sum(object$`(weights)`[!is.na(fitted(object))]))) / 
  (1 - (edf1(object) + 2)/(sum(object$`(weights)`[!is.na(fitted(object))])))
aiccor_edf2 <- function(object, trhatsq) log(object$risk()[-1]/
                                               (sum(object$`(weights)`[!is.na(fitted(object))]))) + 
  (1 + edf2(object, trhatsq)/(sum(object$`(weights)`[!is.na(fitted(object))])))/
  (1 - (edf2(object, trhatsq) + 2)/(sum(object$`(weights)`[!is.na(fitted(object))])))

# classical AICs for both edf versions
aicclas_edf1 <- function(object) 2 * object$risk()[-1] + 2 * edf1(object)
aicclas_edf2 <- function(object, trhatsq) 2 * object$risk()[-1] + 
  2 * edf2(object, trhatsq)

# Classical BICs for both edf versions
bic_edf1 <- function(object) 2 * object$risk()[-1] + length(object$response) * edf1(object)
bic_edf2 <- function(object, trhatsq) 2 * object$risk()[-1] + 
  length(object$response) * edf2(object, trhatsq)

# general minimum description length (gMDL) for both edf versions
gmdl_edf1 <- function(object){ 

  s_edf1 <- function(object) object$risk()[-1]/ ((sum(object$`(weights)`[!is.na(fitted(object))])) - edf1(object))  
  
  log(s_edf1(object)) + edf1(object)/(sum(object$`(weights)`[!is.na(fitted(object))])) * 
  log((sum(object$response^2) - object$risk()[-1])/(edf1(object) * s_edf1(object)))
  
}
gmdl_edf2 <- function(object, trhatsq){ 
  
  s_edf2 <- function(object, trhatsq) 
    object$risk()[-1]/ ((sum(object$`(weights)`[!is.na(fitted(object))])) - edf2(object, trhatsq))
  
  log(s_edf2(object, trhatsq)) + 
  edf2(object, trhatsq)/(sum(object$`(weights)`[!is.na(fitted(object))])) * 
  log((sum(object$response^2) - object$risk()[-1])/
        (edf2(object, trhatsq) * s_edf2(object, trhatsq)))
  
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

