#' Run the FastRegLS algorithm to fit the parameters of the MELS model 
#' 
#' 
#'  Some text
#' 
#' @param y Response variable. Should be a vector.
#' @param X Mean model fixed-effects design matrix. First column should be an intercept. 
#' @param Z Within-subject variance model fixed-effects design matrix. First column should be an intercept.
#' @param id Vector of subject IDs. Should be in order from 1 to the number of subjects.
#' @param nIter Number of times to execute the inner loop. Defaults to 1. 
#' @param getSigmaSE Whether or not to return the standard errors of the random location and random scale variances. Defaults to \code{TRUE}. If \code{FALSE}, returns NA in these fields.
#'  
#' @return A list with the following elements:
#'   \item{beta}{Estimate of beta, the mean fixed-effect coefficients}
#'   \item{tau}{Estimate of tau, the WS variance fixed-effect coefficients}
#'   \item{tau_l}{Estiamte of tau_l, the effect of the random location effect on the WS variance}
#'   \item{sigma_sq_nu}{Estimate of random location variance}
#'   \item{sigma_sq_omega}{Estimate of random scale variance}
#'   \item{se_beta}{Standard error of beta}
#'   \item{se_tau}{Standard error of tau}
#'   \item{se_tau_l}{Standard error of tau_l}
#'   \item{se_sigma_sq_nu}{Standard error of random location variance}
#'   \item{se_sigma_sq_omega}{Standard error of random scale variance}
#'   \item{nu}{Posterior means of the random location effects, ordered by subject}
#'   \item{omega}{Posterior means of the random scale effects, ordered by subject}
#'   \item{marginal_variances}{Diagonal of the fitted covariance matrix}
#'    
#' @export "fastregls"


fastregls = function(y, X, Z, id, nIter = 1, getSigmaSE = TRUE) {
  
  k = dim(X)[2]
  p = dim(Z)[2]
  
  id.factor = as.factor(id)
  
  #Create a vector with the group sizes
  groupSizes = tabulate(id, nbins = max(id))
  
  #create a vector with indices of first observation of each subject
  inds = rep(0, max(id))
  inds[1] = 1
  for (i in 2:length(inds)) {
    inds[i] = inds[i-1] + groupSizes[i-1]
  }
  
  #Initialize beta
  
  mod = lmer(y~X-1+(1|id.factor), REML = FALSE)
  betaHat = unname(summary(mod)$coef[,1])
  blups = repVec(ranef(mod)$id$'(Intercept)', groupSizes)
  sigma_sq_nu_hat <- summary(mod)$varcor$id[[1]]
  
  #Initialize alpha
  
  delta = y-X%*%betaHat-blups
  tempResponse = log(delta^2)
  mod = lmer(tempResponse~Z-1+ blups + (1|id.factor), REML = FALSE)
  alphaHat = unname(summary(mod)$coef[,1])[1:p] + c(1.2704, rep(0, p-1))
  
  #Initialize rho
  
  corrHat = unname(summary(mod)$coef[,1])[p+1]
  

  gammaHat = repVec(ranef(mod)$id$'(Intercept)', groupSizes)
  
  sigmaSE = NA
  
  for (h in 1:nIter) {
    
    #Estimate random effects
    
    a = (y-X%*%betaHat)/exp(Z %*% alphaHat + corrHat*blups + gammaHat)
    a = sumVec(a, groupSizes)
    b = ( 1/sigma_sq_nu_hat +   sumVec(1/exp(Z %*% alphaHat + corrHat*blups + gammaHat), groupSizes))^(-1)
    blups = as.vector(a*b)
    blups = repVec(blups, groupSizes)
   
    #Estimate beta
    tempResponse = y - blups
    mod = lm(tempResponse~X-1, weights = exp(-Z%*%alphaHat- corrHat*blups-gammaHat))
    betaHat = unname(summary(mod)$coef[,1])
    betaSE = diag(vcov(mod))^0.5
 

    tempResponse = log((y - blups - X %*% betaHat)^2)
 
    mod = lmer(tempResponse~Z-1+ blups +(1|id.factor), REML = FALSE)
    
    
    #Estimate alpha
    
    alphaHat = unname(summary(mod)$coef[,1])[1:p] + c(1.2704, rep(0, p-1))
    alphaSE = diag(vcov(mod))[1:p]^0.5
    
    #Estimate rho
    
    corrHat = unname(summary(mod)$coef[,1])[p+1]
    corrSE = diag(vcov(mod))[p+1]^0.5
    
    gammaHat = repVec(ranef(mod)$id$'(Intercept)', groupSizes)
    
    #Estimate random scale variance
    
    sigmaSqHat = summary(mod)$varcor$id[[1]]
    
    
    #get random location effect variance 
    #code taken from Ben Bolker's stack overflow post:
    #https://stackoverflow.com/questions/31694812/standard-error-of-variance-component-from-the-output-of-lmer
    
 
    if (getSigmaSE == 1 && h == nIter) {
    dd.ML <- lme4:::devfun2(mod,useSc=TRUE,signames=FALSE)
    vv <- as.data.frame(VarCorr(mod)) 
    pars <- vv[,"sdcor"]
    hh1 <- numDeriv::hessian(dd.ML,pars)
    vv2 <- 2*solve(hh1)  
    #multiply by 2 to get se of variance instead of standard deviation
    sigmaSE = 2*sqrt(diag(vv2))[1] 
    }
  
    
    #Estimate sigma_sq_nu within loop 
    tempResponse = log(blups[inds]^2)
    sigma_sq_nu_hat <- exp(mean(tempResponse) + 1.2704)
    
  }
  
  mod = lmer(y~X-1 + (1|id.factor), REML = FALSE, weights = exp(-Z%*%alphaHat- corrHat*blups-gammaHat))
  betaSE = diag(vcov(mod))^0.5
  betaHat = unname(summary(mod)$coef[,1])
  
  #get random location effect variance 
  #code taken from Ben Bolker's stack overflow post:
  #https://stackoverflow.com/questions/31694812/standard-error-of-variance-component-from-the-output-of-lmer

  sigma_sq_nu_hat_SE = NA
  if (getSigmaSE == 1) {
    sigma_sq_nu_hat <- summary(mod)$varcor$id[[1]]
    dd.ML <- lme4:::devfun2(mod,useSc=TRUE,signames=FALSE)
    vv <- as.data.frame(VarCorr(mod)) 
    pars <- vv[,"sdcor"]
    hh1 <- numDeriv::hessian(dd.ML,pars)
    vv2 <- 2*solve(hh1)  
    #multiply by 2 to get se of variance instead of standard deviation
    sigma_sq_nu_hat_SE = 2*sqrt(diag(vv2))[1]  
  }
  
  
  marginal_variances = exp(Z%*%alphaHat + corrHat*blups + gammaHat)
  
  return(list(beta = betaHat, 
              tau = alphaHat, 
              tau_l = corrHat, 
              sigma_sq_nu = sigma_sq_nu_hat,
              sigma_sq_omega = sigmaSqHat, 
              se_beta = betaSE, 
              se_tau = alphaSE,
              se_sigma_sq_nu = sigma_sq_nu_hat_SE, 
              se_tau_l = corrSE, 
              se_sigma_sq_omega = sigmaSE, 
              nu = blups, 
              omega = gammaHat, 
              marginal_variances = marginal_variances))
}