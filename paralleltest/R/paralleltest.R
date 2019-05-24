
#' Title parallel test
#'
#' @param data y: Response value;
#'            x2: Dilution Dose;
#'            x1: the group number - 0 for the standard group, 1,2,3... for the test group
#' @param bd the base number of logarithmic transformation of dilution dose
#' @param bp the base number of logarithelphmic transformation of relative potency
#' @param j the number of test group
#' @param snc  Standard Nominal Concentration
#'
#' @return a list
#' @export
#'
#' @examples
#' result <- parallel_test(data=mydata, bd=2, br=2, j=1, snc=48000000)
#'
parallel_test<-function(data,bd,br,j,snc){

data$logx2<-log(data$x2,bd)

data_stan = subset(data,x1==0)
data_test = subset(data,x1==j)
mydata = rbind(data_stan,data_test)

fitF = lm(y~x1+logx2+x1*logx2+1,data=mydata)
fitR = lm(y~x1+logx2+1,data=mydata)
summary_F<-summary(fitF)$coefficients
summary_R<-summary(fitR)$coefficients


anov_test<-anova(fitR,fitF)


Ns = nrow(data_stan)
SXX_stan = sum(data_stan$logx2^2)-(sum(data_stan$logx2)^2)/Ns
SXY_stan = sum(data_stan$y*data_stan$logx2)-sum(data_stan$logx2)*sum(data_stan$y)/Ns

Nt = nrow(data_test)
SXX_test = sum(data_test$logx2^2)-(sum(data_test$logx2)^2)/Nt
SXY_test = sum(data_test$y*data_test$logx2)-sum(data_test$logx2)*sum(data_test$y)/Nt

b_common = summary(lm(y~x1+logx2+1,data=mydata))$coefficients[3,1]
M = mean(data_test$logx2)-mean(data_stan$logx2)-(mean(data_test$y)-mean(data_stan$y))/b_common

RP = br^M
TP = RP*snc


df = nrow(data_test)+nrow(data_stan)-3
hat_sigmasquare = sum(fitR$residuals^2)/df
t = qt(0.975,df)
g = (hat_sigmasquare*(t^2))/(b_common^2*(SXX_stan+SXX_test))
i = mean(data_test$logx2)-mean(data_stan$logx2)
h = (t/b_common)*sqrt(hat_sigmasquare*((1-g)*(1/Ns+1/Nt)+(M-i)/(SXX_stan+SXX_test)))

M_L = i+(M-i-h)/(1-g)
M_U = i+(M-i+h)/(1-g)

RP_L = br^M_L
RP_U = br^M_U

TP_L = snc*RP_L
TP_U = snc*RP_U


b_stan = summary(lm(y~logx2+1,data=data_stan))$coefficients[2,1]
b_stan_U = confint(lm(y~logx2+1,data=data_stan))[2,2]
b_stan_L = confint(lm(y~logx2+1,data=data_stan))[2,1]

b_test = summary(lm(y~logx2+1,data=data_test))$coefficients[2,1]
b_test_U = confint(lm(y~logx2+1,data=data_test))[2,2]
b_test_L = confint(lm(y~logx2+1,data=data_test))[2,1]

b_common_U = confint(lm(y~x1+logx2+1,data=mydata))[3,2]
b_common_L = confint(lm(y~x1+logx2+1,data=mydata))[3,1]

res_slopes = matrix(c(round(b_stan,3),round(b_stan_L,3),round(b_stan_U,3),round(b_test,3),round(b_test_L,3),round(b_test_U,3),round(b_common,3),round(b_common_L,3),round(b_common_U,3)),byrow=T,nrow=3)
rownames(res_slopes)<-c("Standard Slope","Test Slope","Common Slope")
colnames(res_slopes)<-c("Value","Lower","Upper")

res_Fieller = matrix(c(round(M,3),round(M_L,3),round(M_U,3),round(RP,3),round(RP_L,3),round(RP_U,3),TP,TP_L,TP_U),byrow=T,nrow=3)
rownames(res_Fieller)<-c("Potency","Relative Potency","Test Potency")
colnames(res_Fieller)<-c("Value","Lower","Upper")

res<-list(j,summary_F,summary_R,res_slopes,anov_test,res_Fieller)
names(res)<-c("the group of test","Full Model","Reduce Model","Confidence intervals of Slopes","ANOVA Test","Confidence intervals of Potency,Relative Potency,Test Potency")


par(mfrow=c(1,2))
plot(data_test$logx2,data_test$y,col="Red",pch=17,cex=1,main="FullModel",
    xlab="Log(Dose)",ylab="Cycle threshold",
    xlim=c(0.9*min(mydata$logx2),1.1*max(mydata$logx2)),
    ylim=c(0.9*min(mydata$y),1.1*max(mydata$y)))
lines(data_test$logx2,predict(fitF,data_test),col="Red",lwd=2)

points(data_stan$logx2,data_stan$y,col="Black",pch=15,cex=0.9)
lines(data_stan$logx2,predict(fitF,data_stan),col="Black",lwd=2,lty=3)

plot(data_test$logx2,data_test$y,col="Red",pch=17,cex=1,
    xlab="Log(Dose)",ylab=" ",
    xlim=c(0.9*min(mydata$logx2),1.1*max(mydata$logx2)),
    ylim=c(0.9*min(mydata$y),1.1*max(mydata$y)),main="ReduceModel")
lines(data_test$logx2,predict(fitR,data_test),col="Red",lwd=2)

points(data_stan$logx2,data_stan$y,col="Black",pch=15,cex=0.9)
lines(data_stan$logx2,predict(fitR,data_stan),col="Black",lwd=2,lty=3)

l = legend("topright",c("test","standard"),lty=c(1,3),pch=c(17,15),col=c(2,9),cex=1)

res
}








