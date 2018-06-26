
#### functions calulating the Kullback-Leibler-Divergence (KLD) for different distributions ####

# functions computing the Kullback-Leibler-Divergence for different distributions (assuming they have been correctly specified)
# mu, sigma, nu: true parameters
# mu.hat, sigma.hat, nu.hat: parameter estimates

# Kullback-Leibler-Divergence of an estimated Normal Distribution
Gaussian_KLD <- function(mu, mu.hat, sigma, sigma.hat) { 
  log(sigma.hat) - log(sigma) + .5*(sigma^2-sigma.hat^2) /sigma.hat^2 + .5*(mu-mu.hat)^2 /sigma.hat^2 }
# Kullback-Leibler-Divergence of an estimated Gamma Distribution
Gamma_KLD <- function(mu, mu.hat, sigma, sigma.hat) { 
  lgamma(sigma) - lgamma(sigma.hat) + sigma*log(mu/sigma) - (sigma - sigma.hat)*digamma(sigma.hat) - sigma*log(mu.hat/sigma.hat) - sigma.hat + sigma*mu.hat/mu }
# Kullback-Leibler-Divergence of an estimated Zero Adapted Gamma Distribution
ZAGA_KLD <- function(mu, mu.hat, sigma, sigma.hat, nu, nu.hat) { 
  nu*(log(nu)-log(nu.hat)) + (1-nu)*gammaKLD(mu = mu, mu.hat = mu.hat, sigma = sigma, sigma.hat = sigma.hat) }
# Kullback-Leibler-Divergence of an estimated Binomial Distribution
Binomial_KLD <- function(mu, mu.hat) {
  mu * (log(mu) - log(mu.hat)) + (1-mu) * (log(1-mu) - log(1-mu.hat))
}


#### Method returning the KLD for each observation for an fitting object ####

calc_KLD <- function(object, mu, sigma, nu, family ) {
  UseMethod("calc_KLD", object)
}


# Function calculating the Kullback-Leibler-Divergence for fit-object (given normal distribution assumption)
# mu and sigma are the true values as matrix for FDboost(LSS) and as vector for the others
calc_KLD.mboost <- function(object, mu, sigma = NULL, nu = NULL, family = c("Gaussian", "Binomial", "GaussianLSS", "GammaLSS", "ZAGA")) # "object" corresponds to the fitting object, "data" is a data.frame/list containing the new.dat for prediction and the corresponding original mu and sigma values for each observation
{
  family <- match.arg(family)
  
  if(inherits(object,"mboostLSS"))  
  {
    mu.hat <- predict(object$mu, type = "response")  # add newdata=data as soon as the problem with the newdata in fdboost - prediction is solved!!!!
    sigma.hat <- predict(object$sigma, type = "response")
    if(family == "ZAGA") nu.hat <- predict(object$nu, type = "response")
  } 
  if(!(inherits(object,"mboostLSS")))  
  {
    mu.hat <- predict(object, type = "response")
    sigma.hat <- NULL
  }
  
  # If sigma is assumed to be constant and not modeled, 
  # the constant sigma.hat ist chosen to be the sigma.hat minimizing the joint KLD for the present estimations of mu
  if(is.null(sigma.hat)) {
    if(family == "Gaussian") sigma.hat <- sd(resid(object))
  }
  
  # return the Kullback-Leibler-Divergence
  kld <- switch(family, 
                Gaussian = Gaussian_KLD(mu = mu, mu.hat = mu.hat, sigma = sigma, sigma.hat = sigma.hat),
                GaussianLSS = Gaussian_KLD(mu = mu, mu.hat = mu.hat, sigma = sigma, sigma.hat = sigma.hat),
                Binomial = Binomial_KLD(mu = mu, mu.hat = mu.hat),
                GammaLSS = Gamma_KLD(mu = mu, mu.hat = mu.hat, sigma = sigma, sigma.hat = sigma.hat),
                ZAGA = ZAGA_KLD(mu = mu, mu.hat = mu.hat, sigma = sigma, sigma.hat = sigma.hat, nu = nu, nu.hat = nu.hat)
  )
  
  return(kld)
}


#### Function returning the mean KLD in each boosting iteration for an mboost object ####

trace_KLD <- function(object, mu, sigma, nu, family, grid = 0:mstop(object)) {
  meankld <- vector(length = length(grid))
  for(i in seq_along(grid)) meankld[i] <- mean(calc_KLD(object[i], mu, sigma, nu, family))
  attr(meankld, "grid") <- grid
  return(meankld)
}
