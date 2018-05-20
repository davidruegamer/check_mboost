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
aiccor_edf1 <- function(object) log(object$risk()[-1]/ sumw) + 
  (1 + edf1(object)/sumw) / (1 - (edf1(object) + 2)/sumw)
aiccor_edf2 <- function(object, trhatsq) log(object$risk()[-1]/sumw) + 
  (1 + edf2(object, trhatsq)/sumw)/
  (1 - (edf2(object, trhatsq) + 2)/sumw)

# classical AICs for both edf versions
aicclas_edf1 <- function(object) 2 * object$risk()[-1] + 2 * edf1(object)
aicclas_edf2 <- function(object, trhatsq) 2 * object$risk()[-1] + 
  2 * edf2(object, trhatsq)

# Classical BICs for both edf versions
bic_edf1 <- function(object) 2 * object$risk()[-1] + n * edf1(object)
bic_edf2 <- function(object, trhatsq) 2 * object$risk()[-1] + 
  n * edf2(object, trhatsq)

# general minimum description length (gMDL) for both edf versions
gmdl_edf1 <- function(object){ 

  s_edf1 <- function(object) object$risk()[-1]/ (sumw - edf1(object))  
  
  log(s_edf1(object)) + edf1(object)/sumw * 
  log((sum(object$response^2) - object$risk()[-1])/(edf1(object) * s_edf1(object)))
  
}
gmdl_edf2 <- function(object, trhatsq){ 
  
  s_edf2 <- function(object, trhatsq) 
    object$risk()[-1]/ (sumw - edf2(object, trhatsq))
  
  log(s_edf2(object, trhatsq)) + 
  edf2(object, trhatsq)/sumw * 
  log((sum(object$response^2) - object$risk()[-1])/
        (edf2(object, trhatsq) * s_edf2(object, trhatsq)))
  
}

