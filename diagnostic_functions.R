#########################################################################
########### functions which return scalar (iteration specific) ##########
#########################################################################

# trace of the squared hat matrix
traceHatsq <- function(object) sum(diag(crossprod(attr(hatvalues(object), "hatmatrix"))))

#########################################################################
########### functions which return vector ###############################
#########################################################################

edf1 <- function(object) attr(hatvalues(object), "trace")
edf2 <- function(object, traceHatSqVec) 2 * edf1(object) - traceHatSqVec
s_edf1 <- function(object) object$risk()[-1]/ (sumw - edf1(object))
s_edf2 <- function(object, traceHatSqVec) 
  object$risk()[-1]/ (sumw - edf2(object, traceHatSqVec))

aiccor_edf1 <- function(object) log(object$risk()[-1]/ sumw) + 
  (1 + edf1(object)/sumw) / (1 - (edf1(object) + 2)/sumw)
aicclas_edf1 <- function(object) 2 * object$risk()[-1] + 2 * edf1(object)
gmdl_edf1 <- function(object) log(s_edf1(object)) + edf1(object)/sumw * 
  log((sum(object$response^2) - object$risk()[-1])/(edf1(object) * s_edf1(object)))
bic_edf1 <- function(object) 2 * object$risk()[-1] + n * edf1(object)

aiccor_edf2 <- function(object, traceHatSqVec) log(object$risk()[-1]/sumw) + 
  (1 + edf2(object, traceHatSqVec)/sumw)/
  (1 - (edf2(object, traceHatSqVec) + 2)/sumw)
aicclas_edf2 <- function(object, traceHatSqVec) 2 * object$risk()[-1] + 
  2 * edf2(object, traceHatSqVec)
gmdl_edf2 <- function(object, traceHatSqVec) log(s_edf2(object, traceHatSqVec)) + 
  edf2(object, traceHatSqVec)/sumw * 
  log((sum(object$response^2) - object$risk()[-1])/
        (edf2(object, traceHatSqVec) * s_edf2(object, traceHatSqVec)))
bic_edf2 <- function(object, traceHatSqVec) 2 * object$risk()[-1] + 
  n * edf2(object, traceHatSqVec)
