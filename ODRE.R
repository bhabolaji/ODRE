
set.seed(1234)
library(quantoptr)

expit <- function(x){ exp(x)/(1+exp(x)) }; logit <- function(x){ log(x/(1-x)) }

n <-4*1000; nsim <- 500; rateseq <- seq(0.1,0.5,by=0.05); res2 <- NULL
error <- rnorm(n, sd= 0.5)

for (rate in rateseq){
rate <- 0.1
  res <- data.frame(matrix(nrow=nsim,ncol=4))
  colnames(res) <- c("plugin","xl","drl","oracle.drl")
for (i in 1:nsim){
## simulate data
i <- 1
s <- sort(rep(1:4,n/4)); x <- (runif(n,-1,1)); ps <- 0.1 + 0.8*(x>0)
mu0 <- (x <= -.5)*0.5*(x+2)^2 + (x/2+0.875)*(x>-1/2 & x<0) +
  (x>0 & x<.5)*(-5*(x-0.2)^2 +1.075) + (x>.5)*(x+0.125); mu1 <-  mu0 ; tau <- 0.5
a <- rbinom(n,1,ps); y <- a*mu1 + (1-a)*mu0 + rnorm(n,sd=(.2-.1*cos(2*pi*x)))*error
## estimate nuisance functions
df <- data.frame(x=x[s==4],y=y[s==4],a=a[s==4])
##
pihat <- expit( logit(ps) + rnorm(n,mean=1/(n/4)^rate,sd=1/(n/4)^rate))
mu1hat <- predict(smooth.spline(x[a==1 & s==2],y[a==1 & s==2]),x)$y
mu0hat <- predict(smooth.spline(x[a==0 & s==2],y[a==0 & s==2]),x)$y
## construct estimators
plugin <- mu1hat-mu0hat
x1 <- predict(smooth.spline(x[a==1 & s==3],(y-mu0hat)[a==1 & s==3]),x)$y
x0 <- predict(smooth.spline(x[a==0 & s==3],(mu1hat-y)[a==0 & s==3]),x)$y
xl <- pihat*x0 + (1-pihat)*x1
pseudo <- ((a-pihat)/(pihat*(1-pihat)))*(y-a*mu1hat-(1-a)*mu0hat) + mu1hat-mu0hat
drl <- predict(smooth.spline(x[s==3],pseudo[s==3]),x)$y
pseudo.or <- ((a-ps)/(ps*(1-ps)))*(y-a*mu1-(1-a)*mu0) + mu1-mu0
  

oracle.drl <- predict(smooth.spline(x[s==3],pseudo.or[s==3]),x)$y
## save MSEs
res$plugin[i] <- (n/4)*mean((plugin-tau)[s==4]^2)
res$xl[i] <- (n/4)*mean((xl-tau)[s==4]^2)
res$drl[i] <- (n/4)*mean((drl-tau)[s==4]^2)
res$oracle.drl[i] <- (n/4)*mean((oracle.drl-tau)[s==4]^2)
# ATE
res$te <- mean(drl[s==4])

drl_tr <- drl[s==4] >0
# Quantoptr
#fit.test <- IPWE_Qopt(data = df, regimeClass = "a~x",
                      #tau = 0.5, moPropen="a~x", cl.setup=1, it.num=1,pop.size=1000,
                      #s.tol=0.3)
qtres_temp <-IPWE_Mopt(df,regimeClass = a~x,
                    moPropen = a~x)
qtr_tr <- (cbind(1, x[s==4]) %*% qtres_temp$coef.orgn.scale) >0
}
#res2 <- rbind(res2, c(rate, apply(res2,2,mean)))
}

