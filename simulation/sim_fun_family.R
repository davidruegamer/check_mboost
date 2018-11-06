
library(FDboost)

### Function to create functional family out of scalar family

scalar_L2 <- Gaussian()

# function calculating trapezoidal rule weights for one function
trapez_weights <- function(t, range) {
  t_diffs <- diff( c(range[1], sort(t), range[2]) )
  t_diffs[-c(1,length(t_diffs))] <- t_diffs[-c(1,length(t_diffs))]/2
  weights_sorted <- t_diffs[-1] + t_diffs[-length(t_diffs)] 
  return(weights_sorted[order(order(t))])
}


# Function converting scalar family to functional family. 
# CAUTION: Might generalize to other families, but only built and tested with respect to Gaussian()
fam_to_fun_fam <- function(object, id, t, range) {
  # make integration weights
  int_weights <- tapply(t, id, trapez_weights, range = range)
  # adjust family object
  environment( object@risk )$loss <- function(y, f) int_weights * environment( object@risk )$loss(y, f)
  object@ngradient <- function(y, f, w = 1) int_weights * object@ngradient(y, f, w = 1)
  object@name <- paste("functional", object@name)
} 

### 

nsim <- 1

# number of grid points for response function
ngrid <- 10
ncoef <- 5
nsamples <- 10
degree <- 3 # for cubic B-splines

data_list <- list()

set.seed(94038)

# function for random spline generation
bs_fun <- function(t) bs( x = t, knots = seq(0, 1, len = ncoef - degree)[c(-1, -(ncoef + degree +1))] , degree = degree, intercept = TRUE, Boundary.knots = c(0,1))
t_dense <- seq(0, 1, len = 100)

## Conduct simulation

for(i in 1:nsim) {
  
  data_list[[i]] <- list()
  # sample grid
  data_list[[i]]$t <- runif(ngrid * nsamples)
  # data_list[[i]]$t <- rep( seq(0,1, len = ngrid), nsamples) # for equal spaced grid
  # sample random spline coefficients
  data_list[[i]]$mean_coefs <- rnorm(ncoef)
  data_list[[i]]$coefs <- matrix( rnorm(ncoef*nsamples, mean = rep(data_list[[i]]$mean_coefs, nsamples), sd = .2 ), nrow = ncoef )
  ## evaluate function
  B_dense <- bs_fun(t_dense)
  data_list[[i]]$mean_y <- B_dense %*% data_list[[i]]$mean_coefs
  # dense samples
  data_list[[i]]$y_dense <- B_dense %*% data_list[[i]]$coefs
  
    # array of design matrices
    B <- aperm( array( bs_fun(data_list[[i]]$t), dim = c(ngrid, nsamples, ncoef)), c(1,3,2))
    data_list[[i]]$y <- matrix( nrow = ngrid, ncol = nsamples )
    for(j in 1:nsamples) data_list[[i]]$y[, j] <- B[,,j] %*% data_list[[i]]$coefs[, j]
    data_list[[i]]$y <- as.vector( data_list[[i]]$y )
  
  data_list[[i]]$ID <- rep(1:nsamples, each = ngrid)
  
  ## Fit usual FDboost with trapez rule for base-learner fit
  model_classic <- FDboost( y ~ 1, timeformula = ~ bbs(t, boundary.knots = c(0,1)), id = ~ ID, numInt = "Riemann", offset = "scalar", data = data_list[[i]])
  system.time((cv_classic <- applyFolds(model_classic)))
  data_list[[i]]$pred_classic <- predict(model_classic[mstop(cv_classic)])
  data_list[[i]]$dense_classic <- predict(model_classic[mstop(cv_classic)], newdata = data.frame(t = t_dense))
  
  ## Fit functional family created from scalar family
  fun_L2 <- with(data_list[[i]], fam_to_fun_fam(scalar_L2, id = ID, t = t, range = c(0,1)))
  model_alternative <- FDboost( y ~ 1, timeformula = ~ bbs(t, boundary.knots = c(0,1)), id = ~ ID, numInt = "equal" , offset = "scalar", data = data_list[[i]])
  system.time((cv_alternative <- applyFolds(model_alternative)))
  data_list[[i]]$pred_alternative <- predict(model_alternative[mstop(cv_alternative)])
  data_list[[i]]$dense_alternative <- predict(model_alternative[mstop(cv_alternative)], newdata = data.frame(t = t_dense))
  
}

## Visualize results

library(magrittr)
library(ggplot2)
library(dplyr)

## True underlying functions
matplot(t_dense, data_list[[i]]$y_dense, t = "l")
lines(t_dense, data_list[[i]]$mean_y, lwd = 3)
## Mark function evaluations
data_list[[i]] %$% points(t, y, col = ID)

pldat <- data_list[[i]] %$% data_frame( y = y, t = t, ID = ID, pred_classic = pred_classic, pred_alternative = pred_alternative )

ggplot(pldat, aes( x = t, group = ID )) + geom_line(aes(y = y), size = .5 ) + geom_point(aes(y = y)) + 
  geom_line( aes(y = pred_classic, col = "classic", group = NULL ) ) + 
  geom_line( aes(y = pred_alternative, col = "alternative", group = NULL ) )

dense_dat <- data_list[[i]] %$% data_frame( t = t_dense, mean_y = as.vector(mean_y), pred_classic = as.vector(dense_classic), pred_alternative = as.vector(dense_alternative) )

ggplot(dense_dat, aes( x = t )) + geom_line(aes(y = mean_y, col = "original")) + 
  geom_line(aes(y = pred_classic, col = "classic")) + 
  geom_line(aes(y = pred_alternative, col = "alternative"))


