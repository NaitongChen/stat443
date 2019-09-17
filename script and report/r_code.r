#Save data file in working directory before running the following

#Data
monthly_data <- read.csv("data.csv",header = TRUE)
hold.month<-monthly_data[c(193:nrow(monthly_data)),]
train.month<-monthly_data[c(1:(nrow(monthly_data)-nrow(hold.month))),]
tsSO2 <- ts(train.month$SO2, frequency = 12, start = c(2001,1))
tsCO <- ts(train.month$CO, frequency = 12, start = c(2001,1))
tsNO2 <- ts(train.month$NO2, frequency = 12, start = c(2001,1))
cormat <- cor(monthly_data[c(2,3,4)])
par(mfrow = c(1,1))
plot(tsCO, main ="Time Series Plot of Carbon Monoxide (CO)", xlab = "Year", ylab = expression(paste("Carbon Monoxide (mg per m"^"3",")")))
par(mfrow = c(1,2))
acf(tsCO, main = "Autocorrelation Plot of CO")
pacf(tsCO, main = "Partial Correlation Plot of CO")

#Persistence
nt<-nrow(train.month)
nh<-nrow(hold.month)
mse<-0
fc<-train.month[nt,2] # CO
zt<-hold.month[1,2]
fcerror<-zt-fc
mse<-mse+fcerror^2
for (i in 2:nh){
  fc<-hold.month[i-1,2]
  zt<-hold.month[i,2]
  fcerror<-zt-fc
  mse<-mse+fcerror^2
}
rmse<-sqrt(mse/nh)
fc.persistence <- c(train.month[nt,2],hold.month[1:nh-1,2])

#Average 
cumsum<-sum(train.month[,2])
mse<-0
fc<-cumsum/nt
zt<-hold.month[1,2]
fcerror<-zt-fc
fc.average <- c(rep(NA,12))
fc.average[1] <- fc
mse<-mse+fcerror^2
for (i in 2:nh){
  cumsum<-cumsum+hold.month[i-1,2]
  fc<-cumsum/(nt+i-1)
  fc.average[i] <- fc
  zt<-hold.month[i,2]
  fcerror<-zt-fc
  mse<-mse+fcerror^2
}
rmse<-sqrt(mse/nh)

#Holt-Winters
hw.co <- HoltWinters(tsCO,seasonal='additive')
#print(sqrt(hw.co$SSE/(nrow(train.month)-12))) # in sample rmse = 0.08503725
par(mfrow = c(1,1))
acf(ts(hw.co$fitted[,1]-hw.co$fitted[,2]-hw.co$fitted[,4], start = c(2001,1)),main = "Autocorrelation Plot of Detrended and \n Deseasonalized CO trend", cex.main = 0.75) 

nt<-nrow(train.month)
nh<-nrow(hold.month)
alpha<-hw.co$alpha
beta<-hw.co$beta
gamma<-hw.co$gamma
lev<-hw.co$coefficients[1]
tre<-hw.co$coefficients[2]
s1<-hw.co$coefficients[3]
s2<-hw.co$coefficients[4]
s3<-hw.co$coefficients[5]
s4<-hw.co$coefficients[6]
s5<-hw.co$coefficients[7]
s6<-hw.co$coefficients[8]
s7<-hw.co$coefficients[9]
s8<-hw.co$coefficients[10]
s9<-hw.co$coefficients[11]
s10<-hw.co$coefficients[12]
s11<-hw.co$coefficients[13]
s12<-hw.co$coefficients[14]

sse<-0
fc<-lev+tre+s1
zt<-hold.month[1,2]
newfc<-fc
fc.hw<-vector()
fc.hw[1]<-fc
fcerror<-zt-fc
sse<-sse+fcerror^2
vprev<-lev
bprev<-tre
for(i in 2:nh){
  # d = 12
  d=12
  if ((i-1)%%d==0 ){s=s12}
  if ((i-1)%%d==1 ){s=s1}
  if ((i-1)%%d==2 ){s=s2}
  if ((i-1)%%d==3 ){s=s3}
  if ((i-1)%%d==4 ){s=s4}
  if ((i-1)%%d==5 ){s=s5}
  if ((i-1)%%d==6 ){s=s6}
  if ((i-1)%%d==7 ){s=s7}
  if ((i-1)%%d==8 ){s=s8}
  if ((i-1)%%d==9 ){s=s9}
  if ((i-1)%%d==10 ){s=s10}
  if ((i-1)%%d==11 ){s=s11}
  if ((i-1)%%d==12 ){s=s12}
  
  if (i%%d==0 ){sn=s11}
  if (i%%d==1 ){sn=s1}
  if (i%%d==2 ){sn=s2}
  if (i%%d==3 ){sn=s3}
  if (i%%d==4 ){sn=s4}
  if (i%%d==5 ){sn=s5}
  if (i%%d==6 ){sn=s6}
  if (i%%d==7 ){sn=s7}
  if (i%%d==8 ){sn=s8}
  if (i%%d==9 ){sn=s9}
  if (i%%d==10 ){sn=s10}
  if (i%%d==11 ){sn=s11}
  if (i%%d==12 ){sn=s12}
  
  vnew<-alpha*(hold.month[i-1,2]-s)+(1-alpha)*(vprev+bprev)
  bnew<-beta*(vnew-vprev)+(1-beta)*bprev
  snew<-gamma*(hold.month[i-1,2]-vnew)+(1-gamma)*s
  fc<-vnew+bnew+sn
  fc.hw[i]<-fc
  zt<-hold.month[i,2]
  fcerror<-zt-fc
  sse<-sse+fcerror^2
  vprev<-vnew
  bprev<-bnew
  if ((i-1)%%d==0 ){s12<-snew}
  if ((i-1)%%d==1 ){s1<-snew}
  if ((i-1)%%d==2 ){s2<-snew}
  if ((i-1)%%d==3 ){s3<-snew}
  if ((i-1)%%d==4 ){s4<-snew}
  if ((i-1)%%d==5 ){s5<-snew}
  if ((i-1)%%d==6 ){s6<-snew}
  if ((i-1)%%d==7 ){s7<-snew}
  if ((i-1)%%d==8 ){s8<-snew}
  if ((i-1)%%d==9 ){s9<-snew}
  if ((i-1)%%d==10 ){s10<-snew}
  if ((i-1)%%d==11 ){s11<-snew}
  if ((i-1)%%d==12 ){s12<-snew}
}
rmse=sqrt(sse/nh)

#AR2
ar2=arima(ts(train.month[,2]),order=c(2,0,0),method = 'CSS')
sqrt(ar2$sigma2) #in-sample:  0.1024577
nt<-nrow(train.month)
nh<-nrow(hold.month)
phivec<-as.vector(c(ar2$coef[1],ar2$coef[2]))
miu<-ar2$coef[3]
y<-c(train.month[nrow(train.month)-1,2],train.month[nrow(train.month),2],hold.month[,2])
length(y) # 14
y<-y-miu
fc.ar2<-vector()
sse<-0
for (i in 1:nh){
  u<-c(y[i+1],y[i])
  fc<-miu+phivec%*% as.vector(u)
  fc.ar2[i]<-fc
  zt<-hold.month[i,2]
  fcerror<-zt-fc
  sse<-sse+fcerror^2
}
rmse=sqrt(sse/nh)

#AR3
ar3=arima(ts(train.month[,2]),order=c(3,0,0),method = 'CSS')
sqrt(ar3$sigma2) #in-sample:0.100132
nt<-nrow(train.month)
nh<-nrow(hold.month)
phivec<-as.vector(c(ar3$coef[1],ar3$coef[2],ar3$coef[3]))
miu<-ar3$coef[4]
y<-c(train.month[nrow(train.month)-2,2],train.month[nrow(train.month)-1,2],train.month[nrow(train.month),2],hold.month[,2])
length(y) # 15
y<-y-miu
fc.ar3<-vector()
sse<-0
for (i in 1:nh){
  u<-c(y[i+2],y[i+1],y[i])
  fc<-miu+phivec%*% as.vector(u)
  fc.ar3[i]<-fc
  zt<-hold.month[i,2]
  fcerror<-zt-fc
  sse<-sse+fcerror^2
}
rmse=sqrt(sse/nh)

#AR1
ar1=arima(ts(train.month[,2]),order=c(1,0,0),method = 'CSS')
sqrt(ar1$sigma2) #in-sample: 0.1146686
nt<-nrow(train.month)
nh<-nrow(hold.month)
phivec<-ar1$coef[1]
miu<-ar1$coef[2]
y<-c(train.month[nrow(train.month),2],hold.month[,2])
length(y) # 13
y<-y-miu
fc.ar1<-vector()
sse<-0
for (i in 1:nh){
  u<-y[i]
  fc<-miu+phivec%*% as.vector(u)
  fc.ar1[i]<-fc
  zt<-hold.month[i,2]
  fcerror<-zt-fc
  sse<-sse+fcerror^2
}
rmse=sqrt(sse/nh) 

###ARMAX
par(mfrow = c(1,2))
ccf(train.month$NO2, train.month$CO, main = expression("Cross Correlation Plot \n of NO2 and CO"), cex.main = 0.8)
ccf(train.month$SO2, train.month$CO, main = expression("Cross Correlation Plot \n of SO2 and CO"), cex.main = 0.8)
train.month4 = train.month[4:nrow(train.month),]
train.month4$NO2lag3 = train.month[1:(nrow(train.month)-3),3]
train.month4$NO2lag2 = train.month[2:(nrow(train.month)-2),3]
train.month4$NO2lag1 = train.month[3:(nrow(train.month)-1),3]
train.month4$SO2lag2 = train.month[2:(nrow(train.month)-2),4]
train.month4$SO2lag1 = train.month[3:(nrow(train.month)-1),4]
armax=arima(ts(train.month4[,2]),order=c(2,0,0), xreg=train.month4[1:nrow(train.month4),c(5,6,7,8,9)], method='CSS')
sqrt(armax$sigma2) # in-sample: 0.09839809
# out-of-sample RMSE
hold =  hold.month
hold$NO2lag3 = c(train.month[((nrow(train.month)-2):nrow(train.month)),3],hold[1:(nrow(hold)-3),3])
hold$NO2lag2 = c(train.month[((nrow(train.month)-1):nrow(train.month)),3],hold[1:(nrow(hold)-2),3])
hold$NO2lag1 = c(train.month[(nrow(train.month)),3],hold[1:(nrow(hold)-1),3])
hold$SO2lag2 = c(train.month[((nrow(train.month)-1):nrow(train.month)),4],hold[1:(nrow(hold)-2),4])
hold$SO2lag1 = c(train.month[(nrow(train.month)),4],hold[1:(nrow(hold)-1),4])
train<-train.month4[,c(2,5,6,7,8,9)]
hold <-hold[,c(2,5,6,7,8,9)]
phivec<-armax$coef[c(1,2)]
betas<-armax$coef[3:8]
y<-rbind(train[nrow(train)-1,],train[nrow(train),],hold)
r<-vector()
for (j in 1:nrow(y)){
  r[j]<-y$CO[j]-c(1,y$NO2lag3[j],y$NO2lag2[j],y$NO2lag1[j],y$SO2lag2[j],y$SO2lag1[j]) %*%betas
}
fcvec<-vector()
sse<-0
for (i in 1:nh){
  fc<- t(c(1,y$NO2lag3[j],y$NO2lag2[j],y$NO2lag1[j],y$SO2lag2[j],y$SO2lag1[j])) %*%betas + t(phivec)%*%c(r[i+1],r[i])
  fcvec[i]<-fc
  zt<-hold[i,1]
  fcerror<-zt-fc
  sse<-sse+fcerror^2
}
rmse=sqrt(sse/nh)

#SARIMA
arima.season=arima(ts(train.month[,2]),order=c(2,0,0),seasonal=list(order=c(0,1,0),period=12),method = "CSS")
sqrt(arima.season$sigma2) # 0.09162692
# out-of-sample RMSE calculation
phivec<-arima.season$coef[c(1,2)]
s<-12
y<-c(train.month[c((nrow(train.month)-13):nrow(train.month)),2],hold.month[,2])
length(y) # p+s+nh =2+12+12=26
fc.sarima<-vector()
sse<-0
for (i in 1:nh){
  u<-c(y[i+13],y[1+12])-c(y[i+1],y[i])
  fc<- y[i+2]+ phivec %*%u
  fc.sarima[i]<-fc
  zt<-hold.month[i,2]
  fcerror<-zt-fc
  sse<-sse+fcerror^2
}
rmse=sqrt(sse/nh)

##Holt Winter
ts.all.co = ts(monthly_data$CO[1:(nrow(monthly_data))],frequency = 12, start = c(2001,1))
hw.co.fc <- HoltWinters(ts.all.co,seasonal='additive')
fc.ho = predict(hw.co.fc,n.ahead = 12)
##seasonal arima
arima.season.fc=arima(ts.all.co,order=c(2,0,0),seasonal=list(order=c(0,1,0),period=12))
fc.season = predict(arima.season.fc,n.ahead = 12)
z = ts(c(ts.all.co,fc.ho),frequency = 12,start=c(2001,1))
par(mfrow = c(1,1))
plot(z,ylab='CO',main='Forecast of CO Level in 2018')
lines(fc.ho,col='red')
lines(fc.season$pred,col='blue')
legend(x = 'topright',0.5,legend = c('CO level 2001-2017','Holt-Winter Forecast','SARIMA Forecast'),col = c("black","red","blue"),pch=15, cex = 0.6)
##zoom in version
z2 =  ts(c(ts.all.co[169:204],fc.ho),frequency = 12,start=c(2015,1))
plot(z2,ylab='CO',main='Forecast of CO Level in 2018')
lines(fc.ho,col='red')
lines(fc.season$pred,col='blue')
legend(x = 'topright',0.5,legend = c('CO level 2015-2017','Holt-Winter Forecast','SARIMA Forecast'),col = c("black","red","blue"),pch=15, cex = 0.6)

#Summary of Results
month.names <- c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December')
results_table <- data.frame(month.names, round(monthly_data[193:204,2],3),round(fc.persistence,3),round(fc.average,3),round(fc.ar2,3),round(fc.hw,3),round(fc.sarima,3))
colnames(results_table) <- c("2017 Months","Holdout", "Persistence", "Average", "AR(2)", "Holt-Winters","SARIMA")
forecast_table <- data.frame(month.names,round(fc.ho,3),round(fc.season$pred,3))
colnames(forecast_table) <- c("2018 Months","Holt-Winters","SARIMA")
