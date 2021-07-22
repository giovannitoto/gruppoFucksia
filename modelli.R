# Modelli con risposta "study_condition"
# - considerando aggregazione a livello di specie

rm(list=ls())
load("dataset_fix.RData")
load("data_fix_family.RData")

# Elimino country
data_fix$country <- NULL

# Unisco due dataset: li andro' a dividere nuovamente in seguito
data_fix <- cbind(data_fix, data_fix_family[,10:74])

s <- data_fix[data_fix$study_name != "FengQ_2015", ]
v <- data_fix[data_fix$study_name == "FengQ_2015", ]

# Elimino study_name
s <- s[, -1]
v <- v[, -1]

# Bilancio il dataset in modo da avere una proporzione 40%/60%
sb1 <- s[s$study_condition=="adenoma", ]
sb0 <- s[s$study_condition=="control", ]
# In pratica, scarto oss. dal gruppo piu' numeroso in modo da avere gruppi bilanciati
# In questo caso sb0 ha troppe oss. => scarto quelle in eccesso
set.seed(496)
idx <- sample(1:nrow(sb0), trunc(nrow(sb1)*0.6/0.4))
sb0 <- sb0[idx, ]
# Sostituisco data_fix con la versione bilanciata
s <- rbind(sb1, sb0)  
rm(list=c("sb1","sb0","idx"))

prop.table(table(s$study_condition))
prop.table(table(v$study_condition))

# ---------------------------------------------------------------------------- #

# Dataset di stima e convalida
set.seed(28)
cb1 <- sample(1:nrow(s), trunc(nrow(s)*2/3))
cb2 <- setdiff(1:nrow(s), cb1)

# Definisco tabella per il confronto dei modelli
tab_confronto <- matrix(NA, nrow = 1, ncol=3)
colnames(tab_confronto) <- c("dataset", "model", "err")
tab_confronto <- as.data.frame(tab_confronto)

# ---------------------------------------------------------------------------- #

# PROBLEMIIII
colnames(s) # 512, 519, 523, 524, 528, 535, 536, 557, 569
colnames(s)[512] <- "Bacillales_unclassified" 
colnames(v)[512] <- "Bacillales_unclassified" 

colnames(s)[519] <- "Candidatus_Gastranaerophilales_unclassified" 
colnames(v)[519] <- "Candidatus_Gastranaerophilales_unclassified" 

colnames(s)[523] <- "Clostridiales_Family_XIII_Incertae_Sedis" 
colnames(v)[523] <- "Clostridiales_Family_XIII_Incertae_Sedis" 

colnames(s)[524] <- "Clostridiales_unclassified" 
colnames(v)[524] <- "Clostridiales_unclassified" 

colnames(s)[528] <- "Corynebacteriales_unclassified" 
colnames(v)[528] <- "Corynebacteriales_unclassified" 

colnames(s)[535] <- "Eukaryota_unclassified" 
colnames(v)[535] <- "Eukaryota_unclassified" 

colnames(s)[536] <- "Firmicutes_unclassified" 
colnames(v)[536] <- "Firmicutes_unclassified" 

colnames(s)[557] <- "Proteobacteria_unclassified" 
colnames(v)[557] <- "Proteobacteria_unclassified" 

colnames(s)[569] <- "Tissierellia_unclassified" 
colnames(v)[569] <- "Tissierellia_unclassified" 

# In alternativa, si puo' fare in un colpo solo
# library(stringr)
# colnames(s) <- str_replace_all(colnames(s), " ", "_")
# colnames(v) <- str_replace_all(colnames(v), " ", "_")


datasets_s <- list(s[,c(1,2,3,7)], s[,c(1,4,5,6, 8:506)], s[,1:506],
                   s[,c(1,4,5,6, 507:571)], s[, c(1:7, 507:571)])
datasets_v <- list(v[,c(1,2,3,7)], v[,c(1,4,5,6, 8:506)], v[,1:506],
                   v[,c(1,4,5,6, 507:571)], v[,c(1:7, 507:571)])
datasets_names <- c("cliniche", "specie - batteri", " specie - cliniche + batteri",
                    "family - batteri", " family - cliniche + batteri")

# ---------------------------------------------------------------------------- #

### MODELLO LOGISTICO

for (j in 1:length(datasets_names)) {
  # Stimo il modello
  fit.model <- glm(study_condition ~ ., family=binomial, data=datasets_s[[j]])
  # Effettuo previsioni
  y.model <- predict(fit.model, newdata=datasets_v[[j]], type="response")
  et.model <- table(y.model>0.5, datasets_v[[j]]$study_condition)
  e.model <- 1 - sum(diag(et.model))/sum(et.model)
  tab_confronto <- rbind(tab_confronto,
                         c(datasets_names[[j]],"Modello logistico", e.model))
  cat("Data", j, "\n")
}
tab_confronto <- tab_confronto[-1,]

# ---------------------------------------------------------------------------- #

# LASSO (alpha=1)
library(glmnet)

for (j in 1:length(datasets_names)) {
  x   <- model.matrix(~ ., data=datasets_s[[j]][,-1])
  x[,colnames(x)!="gendermale"] <- apply(x[,colnames(x)!="gendermale"], 2,
                                         function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y)))
  x.v <- model.matrix(~ ., data=datasets_v[[j]][,-1])
  x.v[,colnames(x)!="gendermale"] <- apply(x.v[,colnames(x)!="gendermale"], 2,
                                           function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y)))
  fit.model <- cv.glmnet(x, datasets_s[[j]]$study_condition,
                         family="binomial", type.measure="class",
                         alpha=1, nfolds=10)
  # Effettuo previsioni
  y.model <- predict(fit.model, newx=x.v, type="response")
  et.model <- table(y.model>0.5, datasets_v[[j]]$study_condition)
  e.model <- 1 - sum(diag(et.model))/sum(et.model)
  tab_confronto <- rbind(tab_confronto,
                         c(datasets_names[[j]],"Lasso", e.model))
  cat("Data", j, "\n")
}


# ---------------------------------------------------------------------------- #

# ELASTIC NET
library(glmnet)

for (HYP in c(0,0.2,0.4,0.6,0.8)) {
  for (j in 1:length(datasets_names)) {
    x   <- model.matrix(~ ., data=datasets_s[[j]][,-1])
    x[,colnames(x)!="gendermale"] <- apply(x[,colnames(x)!="gendermale"], 2,
                                           function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y)))
    x.v <- model.matrix(~ ., data=datasets_v[[j]][,-1])
    x.v[,colnames(x)!="gendermale"] <- apply(x.v[,colnames(x)!="gendermale"], 2,
                                             function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y)))
    fit.model <- cv.glmnet(x, datasets_s[[j]]$study_condition,
                           family="binomial", type.measure="class",
                           alpha=HYP, nfolds=10)
    # Effettuo previsioni
    y.model <- predict(fit.model, newx=x.v, type="response")
    et.model <- table(y.model>0.5, datasets_v[[j]]$study_condition)
    e.model <- 1 - sum(diag(et.model))/sum(et.model)
    tab_confronto <- rbind(tab_confronto,
                           c(datasets_names[[j]],paste("ElasticNet",HYP), e.model))
    cat("Data", j, "alpha", HYP, "\n")
  }
}

# ---------------------------------------------------------------------------- #
# Group lasso
library(gglasso)

# Si può fare solo con il dataset quello con tutto a livello di specie
# Per raggrupppare le covariate posso mettere tutte quelle cliniche assieme (in caso le droppa),
# poi mettere magari assiemem quelle delle reads
# e poi per le specie devo scegliere un livello della tassonomia

j = 3

y <- ifelse(datasets_s[[j]]$study_condition == "adenoma", 1, 0)

x   <- model.matrix(~ ., data=datasets_s[[j]][,-1])
x[,colnames(x)!="gendermale"] <- apply(x[,colnames(x)!="gendermale"], 2,
                                       function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y)))
x.v <- model.matrix(~ ., data=datasets_v[[j]][,-1])
x.v[,colnames(x)!="gendermale"] <- apply(x.v[,colnames(x)!="gendermale"], 2,
                                         function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y)))

colnames(x)[1:20]

dim(rowData(se_fix))
colnames(rowData(se_fix))
rowData(se_fix)$Family

dim(x)
colnames(x)

family <- c(1,1,1,2,2,2,1,rowData(se_fix)$Family) # e poi gli indici delle famigile

##Cross Validation
#
# PENSO CHE BISOGNA STANDARDIZZARE
set.seed(55)
fit.cv = cv.gglasso(x=x, y=y,
                    loss="ls", pred.loss = "misclass",
                    group=family, nfolds=5)
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

# ---------------------------------------------------------------------------- #

# RANDOM FOREST (CV)
library(randomForest)
library(e1071)

set.seed(496)
for (j in 1:length(datasets_names)) {
  MTRY_RANGE <- 3*(1:15)[3*(1:15)<NCOL(datasets_s[[j]])]
  forest.tune <- tune(randomForest, study_condition ~ ., data=datasets_s[[j]],
                      ranges=list(mtry=MTRY_RANGE),
                      tunecontrol=tune.control(cross=5),
                      ntree=50)
  best.mtry <- min(forest.tune$best.parameters, NCOL(datasets_s[[j]]))
  fit.model <- randomForest(study_condition ~ ., data=datasets_s[[j]], mtry=best.mtry)
  # Effettuo previsioni
  y.model <- predict(fit.model, newdata=datasets_v[[j]])
  et.model <- table(y.model, datasets_v[[j]]$study_condition)
  e.model <- 1 - sum(diag(et.model))/sum(et.model)
  tab_confronto <- rbind(tab_confronto,
                         c(datasets_names[[j]],paste("Random Forest - mtry", best.mtry), e.model))
  cat("Data", j, "\n")
  rm(list=c("forest.tune","fit.model","y.model"))
}

# ---------------------------------------------------------------------------- #

### SUPPORT VECTOR MACHINES (SVM)
library(e1071)

set.seed(496)

for (KERNEL in c("sigmoid", "radial")) {
  for (j in 1:length(datasets_names)) {
    fit.svm.cv <- tune(svm, study_condition ~ ., data=datasets_s[[j]],
                       ranges=list(cost=2^(-6:0)),
                       tunecontrol=tune.control(cross=5),
                       kernel=KERNEL)
    best.cost <- fit.svm.cv$best.parameters
    fit.model <- svm(study_condition ~ ., data=datasets_s[[j]],
                     cost=best.cost, probability=T, kernel=KERNEL)
    y.model <- predict(fit.model, newdata=datasets_v[[j]], decision.values=T)
    et.model <- table(y.model, datasets_v[[j]]$study_condition)
    e.model <- 1 - sum(diag(et.model))/sum(et.model)
    tab_confronto <- rbind(tab_confronto,
                           c(datasets_names[[j]], paste("SVM",KERNEL,"- cost",best.cost), e.model))
    cat("Data", j, KERNEL, "\n")
  }
}

tab_confronto

# ---------------------------------------------------------------------------- #

### BOOSTING
library(ada)

set.seed(28)
for (HYP in c(50,100,200)) {
  for (j in 1:length(datasets_names)) {
    fit.model <- ada(study_condition ~ ., data=datasets_s[[j]], iter=HYP)
    # Previsione
    y.model <- predict(fit.model, newdata=datasets_v[[j]])
    et.model <- table(y.model, datasets_v[[j]]$study_condition)
    e.model <- 1 - sum(diag(et.model))/sum(et.model)
    tab_confronto <- rbind(tab_confronto,
                           c(datasets_names[[j]], paste("Boosting",HYP), e.model))
    cat("Data", j, "\n")
  }
  
}

# ---------------------------------------------------------------------------- #

# ALBERO DI CLASSIFICAZIONE
library(tree)

set.seed(28)
for (j in 1:length(datasets_names)) {
  fit.tree <- tree(study_condition ~ ., data=datasets_s[[j]], split="deviance",
                   control=tree.control(nobs=nrow(datasets_s[[j]]),
                                        minsize=1,
                                        mindev=0))
  cv.tree <- cv.tree(fit.tree, FUN=prune.tree, K=5, method="misclass")
  J <- min(cv.tree$size[cv.tree$dev==min(cv.tree$dev)])
  # Stimo modello finale
  final.tree <- prune.tree(fit.tree, best=J)
  y.treecv <- predict(final.tree, newdata=datasets_v[[j]], type="class")
  et.model <- table(y.treecv, datasets_v[[j]]$study_condition)
  e.model <- 1 - sum(diag(et.model))/sum(et.model)
  tab_confronto <- rbind(tab_confronto,
                         c(datasets_names[[j]], paste("Albero di Classificazione",J), e.model))
  cat("Data", j, "best", J, "\n")
}

# ---------------------------------------------------------------------------- #

# MARS
library(polspline)

set.seed(28)
for (j in 1:length(datasets_names)) {
  fit.mars <- polyclass(datasets_s[[j]]$study_condition, datasets_s[[j]][,-1])  
  
  y.mars <- ppolyclass(datasets_v[[j]][,-1], fit.mars)[,2]
  et.mars <- table(y.mars>0.5, datasets_v[[j]]$study_condition)
  e.mars <- 1 - sum(diag(et.mars))/sum(et.mars)
  tab_confronto <- rbind(tab_confronto,
                         c(datasets_names[[j]], paste("MARS"), e.mars))
  cat("Data", j, "\n")
}

# ---------------------------------------------------------------------------- #

# BAGGING con OOB

library(ipred)

set.seed(28)
for (j in 1:length(datasets_names)) {
  # Costruisco ciclo per vedere se il numero di alberi costruito (nbag) e' sufficiente
  err <- matrix(NA, 15, 2)
  for (i in 1:NROW(err)) {
    err[i,] <- c(i*10, bagging(study_condition ~ ., data=datasets_s[[j]],
                               nbagg=i*10, coob=T)$err)
    cat(i*10, " ")
  }
  NBAGG <- err[which.min(err[,2]),1]
  # Stimo modello migliore
  fit.bagging <- bagging(study_condition ~ .,data=datasets_s[[j]],
                         coob=T, nbagg=NBAGG)
  # Previsione
  y.bagging <- predict(fit.bagging, newdata=datasets_v[[j]])
  et.bagging <- table(y.bagging, datasets_v[[j]]$study_condition)
  e.bagging <- 1 - sum(diag(et.bagging))/sum(et.bagging)
  tab_confronto <- rbind(tab_confronto,
                         c(datasets_names[[j]], paste("Bagging",NBAGG), e.bagging))
  cat("Data", j, "nbagg", NBAGG, "\n")
}

# ---------------------------------------------------------------------------- #



