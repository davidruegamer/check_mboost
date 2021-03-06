#' @title Diagnostic function for mboost objects
#'
#' @param object an mboost object
#' @param plot logical; whether to show diagnostic plots
#' @param what a character vector specifying what is computed per 
#' default. Additional measures can be defined by \code{FUN_iter} and 
#' \code{FUN_obj}.
#' @param FUN_iter a named list of functions
#' that are applied over the all iterations of the mboost object
#' and returned in the results
#' @param FUN_obj a named list of functions
#' that are applied on the given mboost object and return
#' in the results
#' @param ups if available a list of residual matrices of each iteration
#' as returned by the \code{getUpsilons} function.
#' @param ... further paramters passed to ??
#' @return A data.frame with diagnostic measures
#'
#' @examples 
#' cars.gb <- gamboost(dist ~ speed, data = cars, dfbase = 4,
#'                     control = boost_control(mstop = 50))
#' res <- check_mboost(cars.gb)
#' str(res,1)
#' @export
#' @import mboost
#'
check_mboost <- function(
  object,
  plot = TRUE,
  what = c("inbagrisk", 
           "edf1", "edf2", 
           "AIC1", "AIC2", 
           "AICc1", "AICc2", 
           "BIC1", "BIC2", 
           "gMDL1", "gMDL2",
           "selectionFreq", 
           "cumulativeExplainedRisk"),
  FUN_iter = list(), # list(example = function(object) trhatsq(object)),
  FUN_obj = list(), # list(example2 = function(object) mstop(object)),
  ups = NULL,
  ...
)
{
  
  ####### initial check
  if(!inherits(object, "mboost"))
    stop("object is not an mboost object.")
  
  ####### global definitions
  weights <- model.weights(object)
  n <- length(object$response)
  mstopinit <- mstop(object)
  selcourse <- selected(object)
  
  ####### further checks
  if(any(weights == 0)) 
    warning("zero weights") 
  if(mstopinit == 0)
    stop("Diagnostic tool not meaningful for objects fitted for zero iterations.")
  
  ####### get model degrees of freedom
  
  if(!("glmboost" %in% class(object)) & !("gamboost" %in% class(object)))
  {
    
   if(is.null(ups)) ups <- getUpsilons(object)
   resToHat <- function(mat){
     mat <- -1*mat
     diag(mat) <- diag(mat) + 1
     return(mat)
   }
   traceHat <- sapply(ups[2:(mstop(object)+1)], 
                      function(x) sum(diag(resToHat(x))))
   traceHatHat <- sapply(ups[2:(mstop(object)+1)], 
                         function(x) sum(diag(crossprod(resToHat(x)))))
    
     
  }else{
    
    traceHatHat <- sapply(1:mstopinit, function(m)
      sum(diag(crossprod(attr(hatval(object[m]), "hatmatrix")))))
    traceHat <- attr(hatvalues(object[mstopinit]), "trace")
      
  }
  
  ####### apply iterative functions
  
  default_funs <- list(
    # some more functions here ...,
    inbagrisk = function(object) object$risk()[-1],
    AIC1 = function(object) 
      aicclas_edf1(object, trhat = traceHat),
    AIC2 = function(object) 
      aicclas_edf2(object, trhat = traceHat, trhatsq = traceHatHat),
    AICc1 = function(object) 
      aiccor_edf1(object, trhat = traceHat),
    AICc2 = function(object) 
      aiccor_edf2(object, trhat = traceHat, trhatsq = traceHatHat),
    BIC1 = function(object) 
      bic_edf1(object, trhat = traceHat), 
    BIC2 = function(object) 
      bic_edf2(object, trhat = traceHat, trhatsq = traceHatHat),
    gMDL1 = function(object) 
      gmdl_edf1(object, trhat = traceHat), 
    gMDL2 = function(object) 
      gmdl_edf2(object, trhat = traceHat, trhatsq = traceHatHat),
    selectionFreq = function(object) 
      extract_sel_path(object),
    cumulativeExplainedRisk = function(object)
      extract_cum_expl_risk(object)
  )
  
  # exclude gMDL in case the repsonse is not numeric
  if(!is.numeric(object$response))
  {
    
    what <- setdiff(what, c("gMDL1", "gMDL2"))
    
    
  }
  
  default_funs <- default_funs[names(default_funs) %in% what]
  
  # define place holders
  resSq <- rep(NA, mstopinit)
  deffunres <- funitres <- NULL
  if(length(default_funs) > 0) 
    deffunres <- matrix(NA, nrow = mstopinit, ncol = length(default_funs))
  if(length(FUN_iter) > 0) 
    funitres <- matrix(NA, nrow = mstopinit, ncol = length(FUN_iter))
  
  for(m in mstopinit:1){
    
    resSq[m] <- as.numeric(var(object[m]$resid()))

    if(length(FUN_iter) > 0) 
      funitres[m,] <- sapply(FUN_iter, function(fun) fun(object[m]))
  
  }

  iterfunres <- funitres
  if(!is.null(iterfunres)) names(iterfunres) <- c(names(FUN_iter))
  
  ####### apply vec functions

  invisible(object[mstopinit])
  
  FUN_obj <- 
    c(default_funs,
      FUN_obj)
  
  vecfunres <- sapply(FUN_obj, function(fun) fun(object))
  # check that all functions did return a vector
  if(any(sapply(vecfunres, NCOL) != 1))
    stop("functions specified in the list of FUN_obj must return a vector.")
  vecfunres <- as.data.frame(vecfunres)  
  
  ####### combine results  
  
  res <- data.frame(selection = selcourse)
  if(!is.null(iterfunres)) res <- cbind(res, iterfunres)
  if(!is.null(vecfunres)) res <- cbind(res, vecfunres)
  res$residualVariance <- resSq
  res$edf1 <- traceHat
  res$edf2 = edf2(trhat = traceHat, trhatsq = traceHatHat)

  ####### attributes
    
  attr(res, "dfinit") <- if(class(object)[1] == "glmboost") NA else extract(object, "df")  
  attr(res, "lambda") <- if(class(object)[1] == "glmboost") NA else extract(object, "lambda")
    
  class(res) <- c("check_mboost", "data.frame")
  
  ####### plot
  
  if(plot)
  {
    
    ### do something
    
  }
  
  invisible(res)
    
}
