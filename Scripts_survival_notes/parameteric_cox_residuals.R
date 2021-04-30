library(survival)

#################### CREDITS##################################

# https://blogs2.datall-analyse.nl/2016/02/19/rcode_martingale_residuals_aft_model/

################## Thanks a lot

##Weibull model
#generate artificial data

#parameters
beta1 <- 2
beta2 <- 3
beta3 <- -1.5
interceptT <- 5

#distribution covariates
n <- 300
x1 <- rnorm(n,0,2)
x2 <- rnorm(n,0,2)

#event time
time <- rweibull(n, shape=.7,
                 scale=exp(interceptT+beta1*x1+beta2*x2+beta3*x2^2))

#apply some kind of censoring (here: type II censoring)
#status: 1=death/failure, 0=alive/censored
status <- as.numeric(time<sort(time)[n-50])
maxtime <- sort(time)[(n-50)]
time <- ifelse(status==0, maxtime, time)

#correctly specified Weibull model
mod <- survreg(Surv(time, status)~ x1 + x2 + I(x2^2), dist="weibull")
summary(mod)
#shape parameter
(beta <- 1/mod$scale)

#null model
sr <- survreg(Surv(time, status)~ 1, dist="weibull")

#standardized residuals
linFit <- predict(sr, type="lp")
sderr <- (log(time)-linFit)/sr$scale

#Cox-Snell residuals
#are defined as: -log(S[t]), where S[t] is the survival function of the specified distribution
CoxSnellResidual <- function (standRes, weight=1, dist){
  standRes <- standRes[rep(seq_len(length(standRes)), weight)]
  if (dist=="lognormal") {csr <- -log(1-pnorm(standRes))}
  else if (dist=="weibull") {csr <- -log(exp(-exp(standRes)))}
}

#note: use the same weights as those employed in the survreg function
cxsn <- CoxSnellResidual(standRes=sderr, dist="weibull")

#Martingale residuals
#Status: 1=death/failure, 0=alive/censored
MartingaleResidual <- function (CoxSnellResid, Status, weight=1) {
  delta <- Status[rep(seq_len(length(Status)), weight)]
  martingale <- delta-CoxSnellResid
  data.frame(Martingale=martingale, Status=delta)
}

#note: use the same weights as those employed in the survreg function
mgale <- MartingaleResidual(CoxSnellResid=cxsn, Status=status)

#Martingale residuals range from +1 to -inf.
max(mgale$Martingale)
#for a Weibull model these residuals sum to zero (Collett, p.235)
sum(mgale$Martingale)

#the Martingale residuals plots for x1 and x2 detect the correct
#functional form, but have the wrong (opposite) sign
#for instance, the slope of the plot for x1 is negative, whereas the
#coefficient for x1 in the Weibull model is positive
plotMartingale <- function (covariate, martingaleResidual,
                            nameCovariate="covariate", weight=1) {
  cov <- covariate[rep(seq_len(length(covariate)), weight)]
  plot(cov, martingaleResidual$Martingale, pch=martingaleResidual$Status,
       xlab=nameCovariate, ylab="Martingale Residual")
  abline(h=0, lty=2, col="gray40")
  lines(lowess(cov, martingaleResidual$Martingale, iter=0), col="blue")
  legend("bottomleft", pch=c(1,0), legend=c("event","censored"))
}

#note: use the same weights as those employed in the survreg function
plotMartingale(covariate=x1, nameCovariate="x1", martingaleResidual=mgale)
plotMartingale(covariate=x2, nameCovariate="x2", martingaleResidual=mgale)



##lognormal model
#generate artificial data

#parameters
beta1 <- 2
beta2 <- 3
beta3 <- -1.5
interceptT <- 5

#distribution covariates
n <- 300
x1 <- rnorm(n,0,2)
x2 <- rnorm(n,0,2)

#event time
time <- rlnorm(n, meanlog=interceptT+beta1*x1+beta2*x2+beta3*x2^2, sdlog=2)

#apply some kind of censoring (here: type II censoring)
#status: 1=death/failure, 0=alive/censored
status <- as.numeric(time<sort(time)[n-50])
maxtime <- sort(time)[(n-50)]
time <- ifelse(status==0, maxtime, time)

#correctly specified lognormal model
mod <- survreg(Surv(time, status)~ x1 + x2 + I(x2^2), dist="lognormal")
summary(mod)

#null model
sr <- survreg(Surv(time, status)~ 1, dist="lognormal")

#standardized residuals
linFit <- predict(sr, type="lp")
sderr <- (log(time)-linFit)/sr$scale

#Cox-Snell residuals
#(remember to use the same weights as those employed in the survreg function)
cxsn <- CoxSnellResidual(standRes=sderr, dist="lognormal")

#Martingale residuals
#Status: 1=death/failure, 0=alive/censored
#(remember to use the same weights as those employed in the survreg function)
mgale <- MartingaleResidual(CoxSnellResid=cxsn, Status=status)

#Martingale residuals range from +1 to -inf
max(mgale$Martingale)
#note that for a lognormal model these residuals do not sum to zero
sum(mgale$Martingale)

#the Martingale residuals plots for x1 and x2 detect the correct
#functional form, but have the wrong (opposite) sign
#(remember to use the same weights as those employed in the survreg function)
plotMartingale(covariate=x1, nameCovariate="x1", martingaleResidual=mgale)
plotMartingale(covariate=x2, nameCovariate="x2", martingaleResidual=mgale)