d<-read.table("CompareResultsLasso_CompareSparkvsR.csv",header=TRUE,sep=",")
str(d)
d$diff_mse<-d$MSE_R - d$MSE_Spark
d$diff_lambda<-d$lambda_R - d$lambda_Spark

jpeg("Diff_MSE_stnr.jpg")
with(d,plot(signal.to.noise.ratio, diff_mse, pch="+",xlab="Signal-to-noise ratio",ylab="Mean Squared Error in R - Spark"))
dev.off()

with(d,plot(num.vars, diff_mse, pch="+",xlab="Number of features",ylab="Mean Squared Error in R - Spark"))

jpeg("Diff_MSE_observations_per_feature.jpg")
with(d,plot(num.observations/num.vars, diff_mse, pch="+",xlab="Number of observations/features",ylab="Mean Squared Error in R - Spark"))
dev.off()

d$diff_lambda_min<-abs(d$diff_lambda)
d$diff_lambda_min<-ifelse(d$diff_lambda_min<0.1,0.1,d$diff_lambda_min)

jpeg("Diff_MSE_observations_per_feature_bubble_plot.jpg")
with(d, symbols(x=(num.observations/num.vars), y=diff_mse, 
                circles=d$diff_lambda_min, inches=1/3,
                  ann=F, bg="steelblue2", fg=NULL,
                xlab="Number of observations/features",
                ylab="Mean Squared Error in R - Spark"))
dev.off()


with(d,tapply(diff_mse,polynomial.degree,mean))
#1          2 
#-16.14489 -105.55090 

with(d,tapply(diff_mse,polynomial.degree,sd))
#1         2 
#33.61262 130.07721 
