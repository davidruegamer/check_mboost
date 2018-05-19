#' @title Diagnostic function for mboost objects
#'
#' @param object an mboost object
#' @param plot logical; whether to show diagnostic plots
#' @param FUN_iter a function defining a list of functions
#' that are applied over the all iterations of the mboost object
#' @param FUN_obj a function defining a list of functions
#' that are applied on the given mboost object
#' @param ... further paramters passed to ??
#'
#'
#'
check_mboost <- function(
  object,
  plot = TRUE,
  FUN_iter = function(object) 
    list(fun1 = fun1(object), 
         fun2 = fun2(object)),
  FUN_obj = function(object) 
    list(fun3 = fun3(object), 
         fun4 = fun4(object)),
  ...
)
{
  
  ####### initial check
  if(!inherits(object, "mboost"))
    stop("object is not an mboost object.")
  
  ####### global definitions
  weights <- model.weights(object)
  sumw <- sum(weights[!is.na(fitted(object))])
  n <- length(object$response)
  mstopinit <- mstop(object)
  selcourse <- selected(object)
  
  # copy of the mboost object - necessary when running iterative first?
  object_copy <- object
  
  ####### further checks
  if(any(weights == 0)) 
    warning("zero weights") 
  if(mstopinit == 0)
    stop("Diagnostic tool not meaningful for objects fitted for zero iterations.")
  
  ####### apply iterative functions
  
  iterfunres <- sapply(mstopinit:1, function(m) FUN_iter(object[m]), 
                       simplify = TRUE)
  iterfunres <- (t(iterfunres))[mstopinit:1,]
  
  ####### apply vec functions

  vecfunres <- FUN_obj(object)
  # check that all functions did return a vector
  if(any(sapply(vecfunres, NCOL) != 1))
    stop("functions specified in the list of FUN_obj must return a vector.")
  vecfunres <- do.call("cbind", vecfunres)
    
  ####### combine results  
  
  res <- cbind(selection = selcourse,
               iterfunres,
               vecfunres)   
    
  ####### attributes
    
  attributes(res, "dfinit") <- extract(object, "df")  
  attributes(res, "lambda") <- extract(object, "lambda")
    
  class(res) <- "check_mboost"  
  
  ####### plot
  
  if(plot)
  {
    
    ### do something
    
  }
  
  invisible(res)
    
}