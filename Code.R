### Probabilistic Modelling
### Summer Semester 2020
### Project: An Econimic Evaluation of Emperical Exchange Rate Models


### Library
library(dplyr)
library(readxl)
library(rjags)
library(psych)

### Upload data from excel file
Data_Ex_Rate <-read_excel(path = file.choose())

### Function for finding a predictive variable Xt for MF1, MF2, MF3, FP
pred_data <-function(Data_Rate){
  
  #Transform variable Ym to a time series 
  Ym.ts<-ts(Data_Rate$Ym, start=1979, frequency = 12)  
  
  #Decompose time series and deseasonalize
  Ym.decom <- decompose(Ym.ts, type = "mult")
  Data_Rate$Ym_ds <- as.vector (Ym.ts/Ym.decom$seasonal)
  
  #take logarithmic transformations of the data
  Data_Rate[,2:8] <- log(Data_Rate[,2:8])
  
  #Calculate the monetary fundamentals series zt according to Equation (2) 
  Data_Rate$Z <- (Data_Rate$M-Data_Rate$Mm)-(Data_Rate$Y-Data_Rate$Ym_ds)
  
  #Find xt for MF1 according to Equation (1)
  Data_Rate$MF1_xt <- Data_Rate$Z - Data_Rate$S
  Data_Rate$MF1_xt <- (Data_Rate$MF1_xt/lag(Data_Rate$MF1_xt)-1)*100 
  
  #Estimate linear regression and use residuals as xt for MF2
  lmodel <- lm(S ~ Z, data= Data_Rate)
  Data_Rate$MF2_xt <- -lmodel$residuals
  Data_Rate$MF2_xt <- (Data_Rate$MF2_xt/lag(Data_Rate$MF2_xt)-1)*100 
  
  #Estimate linear regression with simple time trend and use residuals as xt for MF3
  #as.numeric(rownames(Data_Ex_Rate)) is t variable
  ltmodel <- lm(S ~ Z + as.numeric(rownames(Data_Rate)), data=Data_Rate )
  Data_Rate$MF3_xt <- -ltmodel$residuals
  Data_Rate$MF3_xt <- (Data_Rate$MF3_xt/lag(Data_Rate$MF3_xt)-1)*100 
  
  #Find xt for Fp according to Equation (5)
  Data_Rate$FP <- Data_Rate$F - Data_Rate$S
  Data_Rate$FP <- ((Data_Rate$FP-lag(Data_Rate$FP))/lag(Data_Rate$FP))*100 
  Data_Rate$d_S <- (Data_Rate$S - lag(Data_Rate$S))/lag(Data_Rate$S)*100
  
  Data_Rate$FP[which(!is.finite(Data_Rate$FP))]<-0
  Data_Rate$d_S[which(!is.finite(Data_Rate$d_S))]<-0
  Data_Rate$MF1_xt[which(!is.finite(Data_Rate$MF1_xt))]<-0
  Data_Rate$MF2_xt[which(!is.finite(Data_Rate$MF2_xt))]<-0
  Data_Rate$MF3_xt[which(!is.finite(Data_Rate$MF3_xt))]<-0
  
  return(Data_Rate)
}

### Descriptive statistics
describe(Data_Ex_Rate, type = 2)


### Division of dataset into two subset
Data_in_sample <-pred_data(Data_Ex_Rate[1:324,])
Data_out_sample <-pred_data(Data_Ex_Rate[325:nrow(Data_Ex_Rate),])


#### JAGS Modelling Linear Regression in-sample--------

LR_model <-function(y,x,n){
 
 LR_string = " model {
    for (i in 2:n) {   
       y[i] ~ dnorm(mu[i], prec) 
       mu[i] = b[1] + b[2]*xt[i-1] 
     }

    for (j in 1:2) {
       b[j] ~ dnorm(0.0, 1.0/1.0e6)
    }

    prec ~ dgamma(5/2.0, 5*10.0/2.0)
    sig2 = 1.0 / prec
    sig = sqrt(sig2)
    } "

set.seed(72)
data1_jags = list(y=y, n=n, xt=x)
params1 = c("b", "sig")

inits1 = function() {
  inits = list("b"=rnorm(2,0.0,100.0), "prec"=rgamma(1,1.0,1.0))
}
LR_MF = jags.model(textConnection(LR_string), data=data1_jags, inits=inits1, n.chains=3)
return(LR_MF)
}

LR_coef <-function(LR_mf){
   update(LR_mf, 20000) 
   params1 = c("b", "sig")
   LR_sim = coda.samples(model=LR_mf,
                            variable.names=params1,
                            n.iter=100000)
  # Combine multiple chains
  LR_csim = do.call(rbind, LR_sim) 
  
 return (LR_sim)
}

LR_RW_model <-function(y,n){
  
LR_RW_string = " model {
  for (i in 1:n) {   
     y[i] ~ dnorm(b, prec) 
   }
   b ~ dnorm(0.0, 1.0/1.0e6)
   prec ~ dgamma(5/2.0, 5*10.0/2.0)
   sig2 = 1.0 / prec
   sig = sqrt(sig2)
 }"

#Set random seed and identify data that go into JAGS.  
set.seed(72)
data1_jags = list(y=y, n=n) 

params1 = c("b", "sig")

inits1 = function() {
  inits = list("b"=rnorm(1,0.0,100.0), "prec"=rgamma(1,1.0,1.0))
}

LR_RW = jags.model(textConnection(LR_RW_string), data=data1_jags, inits=inits1, n.chains=3)
return(LR_RW)
}


### Implementation of Random Walk model
LR_RW1 <- LR_RW_model(Data_in_sample$d_S[-1],nrow(Data_in_sample)-1)
LR_RW_coef <-LR_coef(LR_RW1)
par(mar=c(1,1,1,1))
plot(LR_RW_coef)
summary(LR_RW_coef)

### Implementation of MF1 model
LR_MF1 <-LR_model(Data_in_sample$d_S[-1],Data_in_sample$MF1_xt[-1],nrow(Data_in_sample)-1)
LR_MF1_coef <-LR_coef(LR_MF1)
par(mar=c(1,1,1,1))
plot(LR_MF1_coef)
summary(LR_MF1_coef)

### Implementation of MF2 model
LR_MF2 <-LR_model(Data_in_sample$d_S[-1],Data_in_sample$MF2_xt[-1],nrow(Data_in_sample)-1)
LR_MF2_coef<-LR_coef(LR_MF2)
par(mar=c(1,1,1,1))
plot(LR_MF2_coef)
summary(LR_MF2_coef)

### Implementation of MF3 model
LR_MF3 <-LR_model(Data_in_sample$d_S[-1],Data_in_sample$MF3_xt[-1],nrow(Data_in_sample)-1)
LR_MF3_coef<-LR_coef(LR_MF3)
par(mar=c(1,1,1,1))
plot(LR_MF3_coef)
summary(LR_MF3_coef)

### Implementation of FP model
LR_FP <-LR_model(Data_in_sample$d_S[-1],Data_in_sample$FP[-1],nrow(Data_in_sample)-1)
LR_FP_coef<-LR_coef(LR_FP)
par(mar=c(1,1,1,1))
plot(LR_FP_coef)
summary(LR_FP_coef)


### The convergence diagnostic-------
gelman.diag(LR_RW_coef)
gelman.diag(LR_MF1_coef)
gelman.diag(LR_MF2_coef)
gelman.diag(LR_MF3_coef)
gelman.diag(LR_FP_coef)

autocorr.diag(LR_RW_coef)
autocorr.diag(LR_MF1_coef)
autocorr.diag(LR_MF2_coef)
autocorr.diag(LR_MF3_coef)
autocorr.diag(LR_FP_coef)

### DIC -----------
dic.samples(LR_RW1, n.iter=20000)
dic.samples(LR_MF1, n.iter=20000)
dic.samples(LR_MF2, n.iter=20000)
dic.samples(LR_MF3, n.iter=20000)
dic.samples(LR_FP,  n.iter=20000)


### Prediction-----
predicted <- function(y,x,z){
  Coef_m <-summary(y)
  Pred_var_mean <- cbind(1,x) %*% Coef_m$statistics[1:2,1]
  Pred_var_up <-cbind(1,x)%*% Coef_m$statistics[1:2,1]+Coef_m$statistics[3,1]
  Pred_var_down <-cbind(1,x)%*% Coef_m$statistics[1:2,1]-Coef_m$statistics[3,1]
  plot(z, type="l", xlab="Month", ylab="Value of exchange rate")
  lines(Pred_var_mean, col="red")
  lines(Pred_var_down, col="green")
  lines(Pred_var_up, col="blue")
}


### Forecast on models
LR_MF1_pred <-predicted(LR_MF1_coef,Data_out_sample$MF1_xt[-1],Data_out_sample$d_S[-1])
LR_MF2_pred <-predicted(LR_MF2_coef,Data_out_sample$MF2_xt[-1],Data_out_sample$d_S[-1])
LR_MF3_pred <-predicted(LR_MF3_coef,Data_out_sample$MF3_xt[-1],Data_out_sample$d_S[-1])
LR_FP_pred <-predicted(LR_FP_coef,Data_out_sample$FP[-1],Data_out_sample$d_S[-1])


### Forecast on Random Walk model
Coef_m <-summary(LR_RW_coef)
x<-rep(1,nrow(Data_out_sample)-1)
Pred_var_mean <- x * Coef_m$statistics[1,1]
Pred_var_up <-x * Coef_m$statistics[1,1]+Coef_m$statistics[2,1]
Pred_var_down <-x * Coef_m$statistics[1,1]-Coef_m$statistics[2,1]
plot(Data_out_sample$d_S[-1], type="l", xlab="Month", ylab="Value of exchange rate")
lines(Pred_var_mean, col="red")
lines(Pred_var_down, col="green")
lines(Pred_var_up, col="blue")














