####### Modification of code from github repo Gaynor "Mediation analysis for common binary outcomes"

## simulate following section 3.1
SimulatePaper <- function(n,k){
  mu_a <- 0.4
  mu_c <- mu_a*0.3
  sig <- 0.75
  beta <- matrix(c(0.1,0.5,0.4),ncol=1)
  theta <- matrix(c(k,0.4,0.5,0.25),ncol=1)
  a <- matrix(rnorm(n,mean=mu_a,sd=sig),ncol=1)
  con <- matrix(rnorm(n,mean=mu_c,sd=sig),ncol=1)
  m <- cbind(1,a,con)%*%beta + rnorm(n,mean=0,sd=sig)
  p <- 1 / (1 + exp(-cbind(1,a,m,con)%*%theta))
  y <- rbinom(n,size=1,prob=p)
  prev <- mean(y)
  return(list(y=y,con=as.vector(con),a=as.vector(a),m=as.vector(m),prev=prev))
}

## from a,m,con,y to xx,mm,co,y data format
ReformatData <- function(dat){
  xx <- matrix(dat$a,ncol=1)
  colnames(xx) <- "xx1"
  mm <- matrix(dat$m,ncol=1)
  colnames(mm) <- "mm1"
  co <- matrix(dat$con,ncol=1)
  colnames(co) <- "co1"
  path <- matrix(1,nrow=1,ncol=1)
  colnames(path) <- colnames(mm)
  rownames(path) <- colnames(xx)
  return(list(y=dat$y,xx=xx,mm=mm,co=co,path=path,family="binomial"))
}


directEffectProposed   <-function(s,t0,t1,t2,t4,c,sigma,b0,b1,b2){ 
  return( oddCalc(pnorm(    (s*t0 + s*t1 + s*t4*c + s*t2*(b0 + b2*c)) / (sqrt(1 + (s*t2*sigma)^2))  )) /
            oddCalc(pnorm(    (s*t0 +        s*t4*c + s*t2*(b0 + b2*c)) / (sqrt(1 + (s*t2*sigma)^2))  ))) }
indirectEffectProposed <-function(s,t0,t1,t2,t4,c,sigma,b0,b1,b2){ 
  return( oddCalc(pnorm(    (s*t0 + s*t1 + s*t4*c + s*t2*(b0 + b1 + b2*c)) / (sqrt(1 + (s*t2*sigma)^2))  )) / 
            oddCalc(pnorm(    (s*t0 + s*t1 + s*t4*c + s*t2*(b0 +      b2*c)) / (sqrt(1 + (s*t2*sigma)^2))  ))) }
oddCalc<-function(p){ return((p/(1-p))) }
mse <- function(sm){ mean(sm$residuals^2) }
erfinv <-  function(x) { qnorm((1 + x) /2) / sqrt(2) }


ComputeProbitApprox <- function(dat){
  y <- dat$y
  a <- dat$a
  m <- dat$m
  con <- dat$con
  cCon<- median(con)
  logitModel <- glm(y ~ a + m  + con, na.action=na.omit, family=binomial(link="logit"))
  probitModel <- glm(y ~ a + m  + con, na.action=na.omit, family=binomial(link="probit"))
  sEst <- median(coef(probitModel)/coef(logitModel), na.rm = TRUE)
  linearModel <- lm(m ~ a + con, na.action=na.omit)
  logitModelCoef <- coef(logitModel); linearModelCoef <- coef(linearModel)
  de <- directEffectProposed(s= sEst,t0=logitModelCoef[1],
                             t1=logitModelCoef[2],t2=logitModelCoef[3],t4=logitModelCoef[4],
                             c=cCon,sigma=sqrt(mse(linearModel)),b0=linearModelCoef[1],b1=linearModelCoef[2],
                             b2=linearModelCoef[3])
  ind <- indirectEffectProposed(s= sEst,t0=logitModelCoef[1],
                             t1=logitModelCoef[2],t2=logitModelCoef[3],t4=logitModelCoef[4],
                             c=cCon,sigma=sqrt(mse(linearModel)),b0=linearModelCoef[1],b1=linearModelCoef[2],
                             b2=linearModelCoef[3])
  return(c(de,ind))
}


