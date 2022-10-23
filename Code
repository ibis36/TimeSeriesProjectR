#Packages
library("CADFtest")
library(readxl)
library(forecast)
library("fGarch")
library("vars")
library("urca")
library(tsibbledata)

#### Data importation ####
aus_prod <- aus_production[1:196,c(1,6,7)]
attach(aus_prod)
View(aus_prod)

#### TS definition ####
#Elec in levels
elec_ts <- ts(Electricity, frequency = 4, start=c(1956,1))
ts.plot(elec_ts) 

#Elec in log 
logelec_ts <- log(elec_ts)
ts.plot(logelec_ts)
max.lag <- round(sqrt(length(logelec_ts))) #14 = square root of the sample size which can serve for the maximum lag length 
CADFtest(logelec_ts, type="trend", criterion = "BIC", max.lag.y = max.lag)
#H0 not rejected -> trend stochastic

#Elec in log-diff
dlogelec_ts <- diff(log(elec_ts))
ts.plot(dlogelec_ts)
acf(dlogelec_ts)
#Seasonal plot
monthplot(dlogelec_ts)

#Elec in log-diff and seasonal diff
ddlogelec_ts <- diff(dlogelec_ts, lag=4)
ts.plot(ddlogelec_ts)
par(mfrow=c(1,1))
acf(ddlogelec_ts)
pacf(ddlogelec_ts)
CADFtest(ddlogelec_ts, type="drift", criterion = "BIC", max.lag.y = max.lag)
#H0 is rejected > time series stationary
Box.test(ddlogelec_ts, lag = max.lag, type = "Ljung-Box")
#H0 is rejected > time series is not white noise

#### Model definition ####
#change parameters for each model

#SARIMA(1,1,1)(2,1,0)
fit_sarima2 <- arima(logelec_ts, order = c(1,1,1), seasonal = list(order = c(2,1,0)))
fit_sarima2 
abs(fit_sarima2$coef/sqrt(diag(fit_sarima2$var.coef))) 
#Validity :
ts.plot(fit_sarima2$residuals)
acf(fit_sarima2$residuals)
Box.test(fit_sarima2$residuals, lag=max.lag, type="Ljung-Box") 
#H0 not rejected -> white noise -> valid
AIC(fit_sarima2, k=log(196))

#### Forecast ####
#Plot
myforecastMA<-predict(fit_sarima2,n.ahead=10)
expected<-myforecastMA$pred
lower<-myforecastMA$pred-qnorm(0.975)*myforecastMA$se;
upper<-myforecastMA$pred+qnorm(0.975)*myforecastMA$se;
cbind(lower,expected,upper)
plot.ts(logelec_ts,xlim=c(2000,2008),ylim=c(10.7,11.1))
lines(expected,col="red")
lines(lower,col="blue")
lines(upper,col="blue")

#Comparing forecast peformance 
y<-logelec_ts
S=round(0.75*length(y))
h=1
error1.h<-c()
for (i in S:(length(y)-h))
{
  mymodel.sub<-arima(y[1:i], order = c(2,1,2),seasonal=c(2,1,1))
  predict.h<-predict(mymodel.sub,n.ahead=h)$pred[h]
  error1.h<-c(error1.h,y[i+h]-predict.h)
}
error2.h<-c()
for (i in S:(length(y)-h))
{
  mymodel.sub<-arima(y[1:i], order = c(1,1,1),seasonal=c(2,1,0))
  predict.h<-predict(mymodel.sub,n.ahead=h)$pred[h]
  error2.h<-c(error2.h,y[i+h]-predict.h)
}
mean(abs(error1.h))
mean(abs(error2.h)) 
mean(abs(error1.h^2))
mean(abs(error2.h^2)) 
dm.test(error1.h,error2.h,h=h,power=1)

#### Multivariate analysis ####
#Gas in levels
gas_ts <- ts(Gas, frequency = 4, start=c(1956,1))
ts.plot(gas_ts)

#Gas in log
loggas_ts <- log(gas_ts)
ts.plot(loggas_ts)
ts.plot(logelec_ts, loggas_ts, col=c("blue", "orange"))
CADFtest(loggas_ts,type="trend",criterion="BIC",max.lag.y=max.lag)
#H0 is not rejected the trend are stochastic

#Gas in log-diff
dloggas_ts <- diff(log(gas_ts))
ts.plot(dloggas_ts)
ts.plot(dlogelec_ts, dloggas_ts, col=c("black","red"))
#seasonality 
monthplot(dloggas_ts) 

#Gas in log-dif and seasonal-diff
ddloggas_ts <-diff(dloggas_ts, lag=4)
ts.plot(ddloggas_ts)
CADFtest(ddloggas_ts,type="drift",criterion="BIC",max.lag.y=max.lag)
#H0 rejected -> stationary

# Cointegration (Engle-granger test)
fit_ci <- lm(logelec_ts~loggas_ts)
res_fit_ci <- fit_ci$residuals
CADFtest(res_fit_ci,type="drift",criterion="BIC",max.lag.y=max.lag)
# The test stat -2.48 is bigger that the Engle-Granger ADF test stat for one explanatory variable -3.41
# Thus we do not reject H0 of no cointegration and conclude that log(Elec) and log(Gas) are not cointegrated

#### VAR model ####
ddlogdata <- data.frame(ddlogelec_ts, ddloggas_ts)
names(ddlogdata) <- c("ddlogelec", "ddloggas")
attach(ddlogdata)

#Automatic selection
VARselect(ddlogdata,lag.max=10,type="const") #SC=4

#Model
fit_var4<-VAR(ddlogdata,type="const",p=4)
summary(fit_var4)

#Analyse residuals
var4_residuals<-resid(fit_var4)
ts.plot(var4_residuals)
par(mfrow=c(2,1))
acf(var4_residuals[,1])
acf(var4_residuals[,2])
ccf(var4_residuals[,1],var4_residuals[,2])

