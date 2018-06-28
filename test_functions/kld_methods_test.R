source("kld_methods.R")
library(mboost)


#### Gaussian example ####

mod <- mboost(dist ~ bols(speed), data = cars)

# sample data
mu <- fitted(mod)
sigma <- sd(resid(mod))

set.seed(903248)
cars$dist_sim <- rnorm(length(mu), mu, sigma)

mod1 <- mboost(dist_sim ~ bols(speed), data = cars) 

# compute KLD
boxplot(calc_KLD(mod1, mu = mu, sigma = sigma, family = "Gaussian"))

mean_kld <- trace_KLD(mod1, mu = mu, sigma = sigma, family = "Gaussian")
plot(attr(mean_kld, "grid"), mean_kld, t = "l")
min(mean_kld)


#### Binomial example ####

cars$bin_dist <- with(cars, factor(dist > mean(dist)))

bod <- mboost(bin_dist ~ bols(speed), data = cars, family = Binomial())

# sample data
Pi <- predict(bod, type = "response")

set.seed(903248)
cars$bin_dist_sim <- factor(rbinom(length(Pi), 1, Pi))

bod1 <- mboost(bin_dist_sim ~ bols(speed), data = cars, family = Binomial()) 

# compute KLD
boxplot(calc_KLD(bod1, mu = Pi, family = "Binomial"))

bean_kld <- trace_KLD(bod1, mu = mu, family = "Binomial")
plot(attr(bean_kld, "grid"), bean_kld, t = "l")
min(mean_kld)

