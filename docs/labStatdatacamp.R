############
# Stat Data Camp 2021 - 
# Lab. :  Lasso con modello logistico, modello multinomiale, lasso per SVM, group lasso, fused lasso,  modelli additivi sparsi
############
# b. scarpa
############

#------------------
#------------------
#------------------
#Penalized logistic regression
#------------------
#------------------
#------------------
#
###Internet Ad Data
#This dataset represents a set of possible advertisements on Internet pages. 
#The features encode the geometry of the image (if available) as well as 
#phrases occuring in the URL, the image's URL and alt text, the anchor text, 
#and words occuring near the anchor text. The task is to predict whether 
#an image is an advertisement ("ad") or not ("nonad").
#
#3 continous; others binary; this is the "STANDARD encoding" mentioned in  Kushmerick (1999).
#
#N. Kushmerick (1999) Learning to remove internet advertisements. 
#In Proceedings of the Third Annual Conference on Autonomous Agents. 
#
rm(list=ls())
load("InternetAd.RData")
ls()
str(InternetAd)
attach(InternetAd)
#
dim(x)
class(x)
nnzero(x)
100*nnzero(x)/prod(dim(x))
table(y)
library(glmnet)
fit <- glmnet(x, y, standardize=FALSE, family="binomial")
plot(fit)
#metto un limite inferiore e superiore ai coefficienti
fitpos <- update(fit, lower=-1, upper=1) 
plot(fitpos,ylim=c(-0.3,3))
#
plot(fit,xvar="dev")
cfit <- cv.glmnet(x, y, standardize=FALSE, family="binomial")
plot(cfit)
cfit$lambda.min
cfit$lambda.1se

#
fit <- glmnet(x,y,family="binomial")
plot(fit)
plot(fit, xvar="dev")
cvfit.a=cv.glmnet(x,y, family="binomial")
plot(cvfit.a)

cvfit.a$lambda.min
cvfit.a$lambda.1se


#
cvfit.b=cv.glmnet(x,y, family="binomial", type.measure="auc")
plot(cvfit.b)
cvfit.b$lambda.min
cvfit.b$lambda.1se

cvfit.c=cv.glmnet(x,y, family="binomial", type.measure="class")
plot(cvfit.c)
cvfit.c$lambda.min
cvfit.c$lambda.1se

coef(fit, s=cvfit.a$lambda.min)[coef(fit, cvfit.a$lambda.min)!=0]
length(coef(fit, s=cvfit.a$lambda.min)[coef(fit, cvfit.a$lambda.min)!=0])
length(coef(fit, s=cvfit.b$lambda.min)[coef(fit, cvfit.b$lambda.min)!=0])
length(coef(fit, s=cvfit.c$lambda.min)[coef(fit, cvfit.c$lambda.min)!=0])

detach(InternetAd)
rm(list=ls())
#############
#------------------
#------------------
#------------------
# Lasso for multinomial model
#------------------
#------------------
#------------------

library(glmnet)


#### Multiclass classification
#Consider the study by Ramaswamy et al. (2001) 
#The authors collected 198 tumor samples, spanning 14 common 
#tumor types, which they split into training (n = 144) 
#and testing (n = 54) sets 
#Their goal was to develop a method for classifying the 
#samples based on microarray data only, with the clinical 
#goal of aiding in diagnosis, particularly of metastatic 
#or atypical tumors
#
#
load("ramaswamy2001.RData")
attach(Ramaswamy2001)
table(y)
fit=glmnet(X, y, family="multinomial")
par(mfrow=c(4, 4))
plot(fit,xvar="dev")
cvfit=cv.glmnet(X,y,family="multinomial", nfolds=8)
par(mfrow=c(1, 1))
plot(cvfit) 

pred=predict(cvfit,newx=Ramaswamy2001$X.test,type="class")
table(pred, Ramaswamy2001$y.test)
mean(pred!=Ramaswamy2001$y.test)
#random accuracy: 
#14/(14*14)=1/14
t1<-table(Ramaswamy2001$y.test, as.factor(pred))
sum(rowSums(t1)*colSums(t1)/(NROW(Ramaswamy2001$X.test)^2))
sum(diag(outer(rowSums(t1),colSums(t1))/(NROW(Ramaswamy2001$X.test)^2)))
#
detach(Ramaswamy2001)
rm(list=ls())
#
#------------------
#------------------
#------------------
# Lasso for SVM
#------------------
#------------------
#------------------
#
#- Internet Ad Data
#------------------
#install.packages("sparseSVM")
library(sparseSVM)

###Internet Ad Data
#This dataset represents a set of possible advertisements on 
#Internet pages. 
#The features encode the geometry of the image (if available) 
#as well as 
#phrases occuring in the URL, the image's URL and alt text, 
#the anchor text, 
#and words occuring near the anchor text. The task is to 
#predict whether 
#an image is an advertisement ("ad") or not ("nonad").
#
#3 continous; others binary; this is the "STANDARD encoding" 
#mentioned in the [Kushmerick, 99].)
#One or more of the three continous features are missing in 28% 
#of the instances; missing values should be interpreted as "unknown".
#
#N. Kushmerick. Learning to remove internet advertisements. 
#In Proceedings of the Third Annual Conference on Autonomous 
#Agents. 1999
#
#------
load("InternetAd.RData")
#
x<- InternetAd$x
y <- InternetAd$y

lasso.svm <- sparseSVM(x, y)
plot(lasso.svm, xvar="norm")
plot(lasso.svm, xvar="lambda")
#
cv.svm <- cv.sparseSVM(x,y)
cv.svm$lambda.min
plot(cv.svm )

plot(lasso.svm)
abline(v = log(cv.svm$lambda.min), col = "red", lty = "dashed")

#increase the part around the choice of lambda
plot(lasso.svm, ylim=c(-25,25))
abline(v = log(cv.svm$lambda.min), col = "red", lty = "dashed")

#coefficient with minimum lambda
coef(lasso.svm)[, cv.svm$lambda == cv.svm$lambda.min][coef(lasso.svm)[, cv.svm$lambda == cv.svm$lambda.min]!=0]
#
#
#exercise: compare with sparse logistic
#
rm(list=ls())
############
#------------------
#------------------
#------------------
# Lasso for ligistic - Group lasso
#------------------
#------------------
#------------------

#logistic lasso
#introduction to tidyverse
#install.packages("tidyverse")
library(tidyverse)
#install.packages("caret")
library(caret)
library(glmnet)

load("pima.RData")
#pregnant: Number of times pregnant
#glucose: Plasma glucose concentration (glucose tolerance test)
#pressure: Diastolic blood pressure (mm Hg)
#triceps: Triceps skin fold thickness (mm) 
#####
#####
#insulin: 2-Hour serum insulin (mu U/ml)
#mass: Body mass index (weight in kg/(height in m)\^2)
#pedigree: Diabetes pedigree function
#age:Age (years)
#diabetes: Class variable (test for diabetes)
#
#
PimaIndiansDiabetes2 <- na.omit(PimaIndiansDiabetes2)
#
#
# Inspect the data
sample_n(PimaIndiansDiabetes2, 3) 
#
# Split the data into training and test set
set.seed(123)
training.samples <- PimaIndiansDiabetes2$diabetes %>% 
  createDataPartition(p = 0.8, list = FALSE)  
train.data  <- PimaIndiansDiabetes2[training.samples, ]
test.data <- PimaIndiansDiabetes2[-training.samples, ]
#
# Indicator code for categorical predictor variables
x <- model.matrix(diabetes~., train.data)[,-1]
x.test <- model.matrix(diabetes ~., test.data)[,-1]
#
# Convert the outcome (class) to a numerical variable
y <- ifelse(train.data$diabetes == "pos", 1, 0)
#
#
# Find the best lambda using cross-validation
set.seed(123) 
cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial")
plot(cv.lasso)
plot(cv.lasso, xlim=c(-5,-4), ylim=c(0.90,0.91))
#
cv.lasso$lambda.min
cv.lasso$lambda.1se
coef(cv.lasso, cv.lasso$lambda.min)
coef(cv.lasso, cv.lasso$lambda.1se)
#
#
#----Compute the final model using lambda.min:
# Final model with lambda.min
lasso.model <- glmnet(x, y, alpha = 1, family = "binomial",
                      lambda = cv.lasso$lambda.min)
#
# Make prediction on test data
probabilities <- lasso.model %>% predict(newx = x.test, type="response")
predicted.classes <- as.factor(ifelse(probabilities > 0.5, "pos", "neg"))
#
# Model accuracy
observed.classes <- test.data$diabetes
mean(predicted.classes == observed.classes)
#
table(predicted.classes,observed.classes)
#
confusionMatrix(predicted.classes, observed.classes)    
confusionMatrix(table(predicted.classes, observed.classes))  #another call fo the function

#-----Compute the final model using lambda.1se:
# Final model with lambda.1se
lasso.model <- glmnet(x, y, alpha = 1, family = "binomial",
                      lambda = cv.lasso$lambda.1se)
#
# Make prediction on test data
probabilities <- lasso.model %>% predict(newx = x.test, type="response")
# probabilities1 <- predict(lasso.model, newx = x.test, type="response")
predicted.classes <- ifelse(probabilities > 0.5, "pos", "neg")
#
# Model accuracy rate
observed.classes <- test.data$diabetes
mean(predicted.classes == observed.classes)
table(predicted.classes,observed.classes)
confusionMatrix(table(predicted.classes, observed.classes))  #library caret
#confusionMatrix(as.factor(predicted.classes), observed.classes)  #library caret


#----full logistic model
# Fit the model
full.model <- glm(diabetes ~., data = train.data, family = binomial)
#
# Make predictions
probabilities <- full.model %>% predict(test.data, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, "pos", "neg")
#
# Model accuracy
observed.classes <- test.data$diabetes
mean(predicted.classes == observed.classes)
confusionMatrix(table(predicted.classes, observed.classes))  #library caret
#------------------
#------------------
#------------------
#lasso for SVM
library(sparseSVM)
#
lasso.svm <- sparseSVM(x, y)
plot(lasso.svm, xvar="norm")
plot(lasso.svm, xvar="lambda")
#
cv.svm <- cv.sparseSVM(x,y)
cv.svm$lambda.min
plot(cv.svm )

plot(lasso.svm)
abline(v = log(cv.svm$lambda.min), col = "red", lty = "dashed")

#zoom around the choice of lambda
plot(lasso.svm, ylim=c(-0.01,0.001))
abline(v = log(cv.svm$lambda.min), col = "red", lty = "dashed")

#fit with minimum lambda
fit.svm <- sparseSVM(x, y, alpha = 1)
coef(fit.svm)[, cv.svm$lambda == cv.svm$lambda.min]
coef(fit.svm)[, cv.svm$lambda == cv.svm$lambda.min][coef(fit.svm)[, cv.svm$lambda == cv.svm$lambda.min]!=0]
#
####comparison with other models
#the coefficients for the logistic where
coef(cv.lasso, cv.lasso$lambda.min)
#all together
cbind(svm=coef(fit.svm)[, cv.svm$lambda == cv.svm$lambda.min], coef(cv.lasso, cv.lasso$lambda.min), logistic=coef(full.model))



#------------------
#------------------
#------------------
# elastic net with svm
#------------------
#------------------
#------------------

elnet.svm <- sparseSVM(x, y, alpha=0.5)
plot(elnet.svm, xvar="norm")
plot(elnet.svm, xvar="lambda")

cv.svm.en <- cv.sparseSVM(x,y, alpha=0.5)
cv.svm.en$lambda.min
plot(cv.svm.en )
abline(v = log(cv.svm.en$lambda.min), col = "red", lty = "dashed")
#
plot(elnet.svm)
abline(v = log(cv.svm.en$lambda.min), col = "red", lty = "dashed")
#
#
#stima con lambda minimo
coef(elnet.svm)[, cv.svm.en$lambda == cv.svm.en$lambda.min][coef(elnet.svm)[, cv.svm.en$lambda == cv.svm.en$lambda.min]!=0]
#
#all together
cbind(svm=coef(fit.svm)[, cv.svm$lambda == cv.svm$lambda.min], coef(cv.lasso, cv.lasso$lambda.min), logistic=coef(full.model), coef(elnet.svm)[, cv.svm.en$lambda == cv.svm.en$lambda.min][coef(elnet.svm)[, cv.svm.en$lambda == cv.svm.en$lambda.min]!=0])
#
rm(list=ls())
#------------------
#------------------
#------------------
#--------------Group lasso
#------------------
#------------------
#------------------
#
#Predict the US Unemployment Rate using some information 
#that’s provided by the Fed. These are the CCAR/DFAST 
#regulatory stress testing scenarios that banks need to use 
#to project their financial statement line items. The Fed provides 
#forward looking projections for a variety of different 
#macro-economic indicators. Let's use them to create forward 
#looking unemployment projections.

#
#install.packages("gglasso")
library(gglasso)

hist=read.csv("Historical Data.csv") #50 observations x 15 variables
proj=read.csv("Projections.csv") #13 observations

#install.packages("zoo")
library(zoo) # for Date management --> as.yearqtr
hist=data.frame(Date=as.Date(as.yearqtr(hist[,1])), hist[,-1])

plot(y=hist$Unemployment.Rate, x=hist$Date, main="Unemployment",
   lwd=2,col="slateblue",type="l",  xlab="Time",ylab="Unemployment %")
grid()

names(hist)

#response variable Unemployment
Y=hist[,4]
#Remove Dates and Unemployment from the model matrix 
X=hist[,c(-1,-4)]
is.matrix(X)
X=as.matrix(X)
#
#
#GDP, income and inflation rate as Group1
#All treasury yields and BBB Corporate yield as Group2
#Mortage rate, prime rate, HPI% change and commercial RE% change as Group3
#Dow Jones and VIX as Group4
#
grp=c(1,1,1,2,2,2,2,3,3,3,3,4,4)
?gglasso
fit=gglasso(x=X, y=Y, group=grp, loss='ls')
coef.mat=fit$beta
plot(fit)


#Group1 enters the equation - first model where group 1 is present
g1=max(which(coef.mat[1,]==0))

#Group2 enters the equation - first model where group 2 is present
g2=max(which(coef.mat[4,]==0))

#Group3 enters the equation - first model where group 3 is present
g3=max(which(coef.mat[8,]==0))

#Group4 enters the equation - first model where group 4 is present
g4=max(which(coef.mat[12,]==0))

#Coefficient Plot

plot(fit$b0,main="Coefficient vs Step",
     ylab="Intercept",xlab="Step (decreasing Lambda =>)",
     xlim=c(-1,100),
     ylim=c(5,max(fit$b0)),
     type="l",lwd=4)
grid()

x=c(g1,g2,g3,g4)
y=c(fit$b0[g1],fit$b0[g2],fit$b0[g3],fit$b0[g4])

points(x=x,y=y,pch=13,lwd=2,cex=2,
     xlim=c(-1,100),ylim=c(5,max(fit$b0)),
     xaxt='n',yaxt='n',xlab="",ylab="")

lmda=round(fit$lambda[c(g1,g2,g3,g4)],2)
text(x=x-0.5,y=y+0.1,labels=c("Group1","Group2",
                      "Group3","Group4"),pos=3,cex=0.9)
text(x=x-0.5,y=y-0.1,labels=paste("Lambda\n=",lmda),pos=1,cex=0.8)

#The intercept is not penalized and hence is always present 
#in the regression equation. But as the plot above shows, 
#each group enters the regression equation at a particular 
#value of lambda.
coef.mat[,c(12,13,14)]

#Note how at step = 13 group 3 enters the regression equation 
#and all the variables in group3 have coefficients > 0. 
#Group4 was already in the regression equation starting at 
#step = 1 (see plot)
#Hence, group4 enters the regression first, then group3, 
#then group1 and group2 enters last.
#Let’s get to cross validation and some forward looking projections.

#
#
##Cross Validation
#
set.seed(55)
fit.cv=cv.gglasso(x=X,y=Y,group=grp,nfolds=10)
plot(fit.cv)
#
#lambda.min is at the border
#
set.seed(55)
fit.cv=cv.gglasso(x=X,y=Y,group=grp,nfolds=10, lambda.factor=0.0004)
plot(fit.cv)

#Pick the best Lambda
lmbda=fit.cv$lambda.1se
lmbda1 <- fit.cv$lambda.min
plot(fit)
abline(v=log(lmbda), lty=2, col=2)
abline(v=log(lmbda1), lty=2, col=2)

coefs=coef(object=fit,s=lmbda)
coefs1=coef(object=fit,s=lmbda1)
cbind(coefs,coefs1)



#At best lambda get coefficients and fitted values
plt=cbind(Y,predict(object=fit, newx=X, s=lmbda),predict(object=fit, newx=X, s=lmbda1)) #,type='link'))
matplot(plt,main="Predicted vs Actual", type='l',lwd=2,
        ylab="Unemplyoment %",
        xlab="Time", xlim=c(1,65))
grid()
legend(x=40,y=5,legend=c("Actual","Fitted-1se","Fitted-min", "Predicted"), fill=1:4,bty="n",cex=0.7)
#
#
#Get forward looking projections
X.p=as.matrix(proj[,-c(1,4)]) #Remove Dates and Unemployment from the model matrix 
plt=cbind(proj$Unemployment.Rate, predict(object=fit,
                                newx=X.p, s=lmbda), predict(object=fit,
                                newx=X.p, s=lmbda1)) #,type='link'))
matplot(51:63,plt,main="Predicted vs Fed Projections",type='l',lwd=2,
        ylab="Unemplyoment %",
        xlab="Time", add=T, col=4)
#
#
#-----there is a discontinuity between fitted and projected
#
#-----connection with the previous prediction
#At best lambda get coefficients and fitted values
plt=cbind(Y,predict(object=fit, newx=X, s=lmbda),predict(object=fit, newx=X, s=lmbda1)) #,type='link'))
matplot(plt,main="Predicted vs Actual", type='l',lwd=2,
        ylab="Unemplyoment %",
        xlab="Time", xlim=c(1,65))
grid()
legend(x=40,y=5,legend=c("Actual","Fitted-1se","Fitted-min", "Predicted"), fill=1:4,bty="n",cex=0.7)
#
#
#
#
#Get forward looking projections
X.p1<- rbind(X[nrow(X),],X.p)        
plt=cbind(c(Y[length(Y)],proj$Unemployment.Rate), predict(object=fit,
        		                  newx=X.p1, s=lmbda), predict(object=fit,
                		          newx=X.p1, s=lmbda1)) #,type='link'))
matplot(50:63,plt,main="Predicted vs Fed Projections",type='l',lwd=2,
         ylab="Unemplyoment %",
		 xlab="Time", add=T, col=4)
#        
#   
#the model has all the groups and hence all the variables 
#are in the regression equation. 
############
rm(list=ls())
#------------------
#------------------
#------------------
# Fused Lasso - Lasso for  SVM
#------------------
#------------------
#------------------
#
#-Glioma data
#--------------------
#install.packages("flsa")
library(flsa)
# CGH example 
l2r <- read.table("glioma.txt")$x
#
plot(l2r, pch=16, cex=0.7, xlab='Genome order', ylab=expression(Log[2]*' ratio'))
#

fit <- flsa(l2r)
  plot(l2r, pch=16, cex=0.7, xlab='Genome order', ylab=expression(Log[2]*' ratio'), col='gray', bty='n', las=1)
  lines(1:length(l2r), flsaGetSolution(fit, lambda1=0, lambda2=0), col='slateblue', lwd=1)
 # lasso solution - without "fusion"
  plot(l2r, pch=16, cex=0.7, xlab='Genome order', ylab=expression(Log[2]*' ratio'), col='gray', bty='n', las=1)
   lines(1:length(l2r), flsaGetSolution(fit, lambda1=1, lambda2=0), col='slateblue', lwd=1)
 # lasso and  "fusion"
  plot(l2r, pch=16, cex=0.7, xlab='Genome order', ylab=expression(Log[2]*' ratio'), col='gray', bty='n', las=1)
   lines(1:length(l2r), flsaGetSolution(fit, lambda1=1, lambda2=1), col='slateblue', lwd=1)
#only "fusion"
  plot(l2r, pch=16, cex=0.7, xlab='Genome order', ylab=expression(Log[2]*' ratio'), col='gray', bty='n', las=1)
   lines(1:length(l2r), flsaGetSolution(fit, lambda1=0, lambda2=1), col='slateblue', lwd=1)
#
#
#different choices for l2
par(mfrow=c(2,2))
l2 <- c(0.1, 0.5, 1, 3)
for (i in 1:length(l2)) {
  plot(l2r, pch=16, cex=0.7, xlab='Genome order', ylab=expression(Log[2]*' ratio'), col='gray', bty='n', las=1)
  lines(1:length(l2r), flsaGetSolution(fit, lambda1=0.1, lambda2=l2[i]), col='slateblue', lwd=1)
  abline(h=0, lty=2)
  mtext(bquote(paste(lambda[1] == 0.1, "     ",   lambda[2] == .(l2[i]))))
#  readline()
}
#par(mfrow=c(1,1))
#
#
#different choices for l1
#par(mfrow=c(2,2))
l1 <- c(0.001,0.01, 0.1, 0.5)
for (i in 1:length(l1)) {
  plot(l2r, pch=16, cex=0.7, xlab='Genome order', ylab=expression(Log[2]*' ratio'), col='gray', bty='n', las=1)
  lines(1:length(l2r), flsaGetSolution(fit, lambda1=l1[i], lambda2=3), col='slateblue', lwd=1)
  abline(h=0, lty=2)
  mtext(bquote(paste(lambda[1] ==  .(l1[i]), "     ",   lambda[2] ==3)))
# readline()
}
par(mfrow=c(1,1))


rm(list=ls())
#------------------

#------------------
#------------------
#------------------
#-------GAM
#------------------
#------------------
#------------------

library(splines)
#install.packages("grpreg")
library(grpreg)


#attachData(Scheetz2006)
#save("fData","X","y", file="Scheetz2006.RData")
load("Scheetz2006.RData")
# Restrict to Chr5
ind <- which(fData$Chr==5)
df <- 3
#
HH <- sapply(as.data.frame(X[,ind]),
             ns, simplify="array", df=df)
#
H <- t(apply(HH, 1, rbind))
#
group <- rep(1:length(ind), each=df)
#
#
cvfit <- cv.grpreg(H, y, group)
#
plot(cvfit)
#
plot(cvfit$fit)
   abline(v=cvfit$lambda.min)

Lines <- function(id, fit, l, ...) {
  x <- X[,ind[id]]
  xx <- seq(min(x), max(x), len=101)
  BB <- predict(ns(x, df=df), newx=xx)
  ind <- 1+((id-1)*df+1):((id-1)*df+df)
  yhat <- as.numeric(BB %*% coef(fit, lambda=l)[ind])
  i <- which.min(abs(xx-mean(x)))
  yhat <- yhat  - yhat[i]
  lines(xx, yhat, ...)
  yhat
}


sel.id<- NULL
for (id in 1:length(ind)) {
#for (id in c(4)) {
	plot(X[,ind[id]],y-mean(y), pch=19, cex=.7, xlab="PNISR", ylab="TRIM32", las=1, bty='n')
	foo<-Lines(id, cvfit$fit, cvfit$lambda.min, lwd=3, col=4)
	if(sum(abs(foo))>1) {sel.id<-c(sel.id,id) ;cat(id, "")}
	}


par(mfrow=c(3,3))
for (id in sel.id[1:9]) {
	plot(X[,ind[id]],y-mean(y), pch=19, cex=.7, xlab=id, ylab="TRIM32", las=1, bty='n', ylim=c(-0.4,0.4))
	foo<-Lines(id, cvfit$fit, cvfit$lambda.min, lwd=3, col=4)
}



par(mfrow=c(3,2))
for (id in sel.id[10:5]) {
	plot(X[,ind[id]],y-mean(y), pch=19, cex=.7, xlab=id, ylab="TRIM32", las=1, bty='n', ylim=c(-0.4,0.4))
	foo<-Lines(id, cvfit$fit, cvfit$lambda.min, lwd=3, col=4)
}



