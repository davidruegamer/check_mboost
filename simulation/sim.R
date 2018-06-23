# save results?
save <- FALSE
# number of simulation repetitions
simReps <- 1000
# maximal mstop
fix_mstop <- 500
# number observations
n <- 300


set.seed(201089)
# covariates
x1 <- scale(rnorm(n), scale = F)
x2 <- scale(rnorm(n), scale = F)
x3 <- scale(rnorm(n), scale = F)
x4 <- scale(rnorm(n), scale = F)
x5 <- scale(rnorm(n), scale = F)
x6 <- scale(rnorm(n), scale = F)
x7 <- scale(rnorm(n), scale = F)
x8 <- scale(rnorm(n), scale = F)
x9 <- scale(rnorm(n), scale = F)
x10 <- scale(rnorm(n), scale = F)

x11 <- scale(rnorm(n), scale = F)
x12 <- scale(rnorm(n), scale = F)
x13 <- scale(rnorm(n), scale = F)
x14 <- scale(rnorm(n), scale = F)
x15 <- scale(rnorm(n), scale = F)
x16 <- scale(rnorm(n), scale = F)
x17 <- scale(rnorm(n), scale = F)
x18 <- scale(rnorm(n), scale = F)
x19 <- scale(rnorm(n), scale = F)
x20 <- scale(rnorm(n), scale = F)

# true linear predictor
mu <- as.numeric(sin(2*x1) + 0.5*x2^2 - 2*x11)
mu <- mu - mean(mu)
sdmu <- sd(mu)

# settings
settings <- expand.grid(list(SNR = c(0.5, 2),
                             knots = c("equal", "unequal"), 
                             df = c("equal", "small difference", "different"),
                             mstop = c("kfold", "bootstrap"), 
                             nu = c(0.001, 0.1, 0.5)))



# list of covariates / baselearner
xList <- list(x1,x2,x3,x4,x5,
              x6,x7,x8,x9,x10,
              x11,x12,x13,x14,x15,
              x16,x17,x18,x19,x20)



#################################################################
################## iteration over settings ######################

for(this_set in 1:nrow(settings)){
  
  set.seed(42)
  
  ###### define setting-specific stuff #####
  # std.dev.
  this_sigma = sdmu / settings$SNR[this_set]
  # knots
  if(settings$knots[this_set]=="equal") this_knots <- rep(10, 10) else
    this_knots <- 4:13
  # df
  this_df <- rep(5, 10)
  if(settings$df[this_set]=="small difference")
    this_df <- this_df + rnorm(10,0,0.001)
  # mstop
  this_type <- settings$mstos[this_set]
  # nu
  this_nu <- settings$nu[this_set]
  
  
  # baselearner definition
  bx1 <- bbs(x1, knots = this_knots[1], df = this_df[1])
  bx2 <- bbs(x2, knots = this_knots[2], df = this_df[2])
  bx3 <- bbs(x3, knots = this_knots[3], df = this_df[3])
  bx4 <- bbs(x4, knots = this_knots[4], df = this_df[4])
  bx5 <- bbs(x5, knots = this_knots[5], df = this_df[5])
  bx6 <- bbs(x6, knots = this_knots[6], df = this_df[6])
  bx7 <- bbs(x7, knots = this_knots[7], df = this_df[7])
  bx8 <- bbs(x8, knots = this_knots[8], df = this_df[8])
  bx9 <- bbs(x9, knots = this_knots[9], df = this_df[9])
  bx10 <- bbs(x10, knots = this_knots[10], df = this_df[10])
  
  bx11 <- bols(x11, df = 1, intercept=FALSE)
  bx12 <- bols(x12, df = 1, intercept=FALSE)
  bx13 <- bols(x13, df = 1, intercept=FALSE)
  bx14 <- bols(x14, df = 1, intercept=FALSE)
  bx15 <- bols(x15, df = 1, intercept=FALSE)
  bx16 <- bols(x16, df = 1, intercept=FALSE)
  bx17 <- bols(x17, df = 1, intercept=FALSE)
  bx18 <- bols(x18, df = 1, intercept=FALSE)
  bx19 <- bols(x19, df = 1, intercept=FALSE)
  bx20 <- bols(x20, df = 1, intercept=FALSE)
  
  blList <- list(bx1, bx2, bx3, bx4, bx5,
                 bx6, bx7, bx8, bx9, bx10,
                 bx11, bx12, bx13, bx14, bx15,
                 bx16, bx17, bx18, bx19, bx20)
  
  if(settings$df[this_set]!="different") blList <- blList[1:10]
  
  
  res <- mclapply(1:simReps, function(nr){
    
    set.seed(nr)
    
    y <- as.numeric(scale(mu + rnorm(n, 0, this_sigma), scale=F))
    modOrg <- mboost_fit(blList, offset = 0, response = y, 
                         control = boost_control(nu = this_nu, mstop = fix_mstop))
    
    cvr <- cvrisk(object = modOrg, 
                  folds = cv(model.weights(modOrg), 
                             type = this_type),
                  papply = lapply)
    
    
    cmr <- check_mboost(modOrg)
    res <- list(checkRes = cmr, 
                resampRes = cvr,
                nr = nr)
    
    return(res)
    
  }, mc.cores=64)
  
  # attach setting
  
  res <- list(res, setting = settings[this_set,])
  
  if(save) saveRDS(res, file=paste0("sim_results/result_",this_set,".RDS"))
  
}