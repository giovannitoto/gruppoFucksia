rm(list = ls())
load("data_fix105_family.RData")
load("dataset_fix2_105.RData")

# Unisco due dataset: li andro' a dividere nuovamente in seguito
data_fix2 <- cbind(data_fix2, data_fix2_family[,11:56])

# ---------------------------------------------------------------------------- #

# Definisco tabella per il confronto dei modelli
tab_confronto <- c("DATA", "MODELLO", "TASSO ERRATA CLASSIFICAZIONE")

tab_confronto <- matrix(NA, nrow = 1, ncol = 3)
colnames(tab_confronto) <- c("dataset", "model", "err")
tab_confronto <- as.data.frame(tab_confronto)

# ---------------------------------------------------------------------------- #

#Guardo i colnames di data_fix2
colnames(data_fix2)
colnames(data_fix2)[320] <- "Bacillales_unclassified"
colnames(data_fix2)[326] <- "Clostridiales_unclassified"
colnames(data_fix2)[335] <- "Eukaryota_unclassified"
colnames(data_fix2)[336] <- "Firmicutes_unclassified"
colnames(data_fix2)[351] <- "Proteobacteria_unclassified"

# ---------------------------------------------------------------------------- #

cv_data <- list(data_fix2[,c(1,2,3,7,8,9,10)], data_fix2[,c(1,4,5,6, 11:315)],
                data_fix2[,1:315], data_fix2[,c(1,4,5,6, 316:361)],
                data_fix2[,c(1:10, 316:361)])

cv_data_names <- c("cliniche", "specie - batteri", " specie - cliniche + batteri",
                    "family - batteri", " family - cliniche + batteri")

# ---------------------------------------------------------------------------- #

# K-fold Cross-Validation
FOLDS <- 5
set.seed(2021-07-22)
folds <- sample(cut(seq(1,nrow(data_fix2)), breaks=FOLDS, labels=FALSE))

'cv_errors <- matrix(NA, nrow = 1, ncol = FOLDS+3)
colnames(cv_errors) <- c("dataset", "model", "mean_err", paste("err", 1:FOLDS, sep = "_"))
cv_errors <- as.data.frame(cv_errors)'

cv_err_glm <- matrix(NA, nrow = 1, ncol = FOLDS+3)
colnames(cv_err_glm) <- c("dataset", "model", "mean_err", paste("err", 1:FOLDS, sep = "_"))
cv_err_glm <- as.data.frame(cv_err_glm)

cv_err_ppr <- cv_err_cart <- cv_err_bagging <- cv_err_rf <- cv_err_boosting <- cv_err_svm <- cv_err_glm

cv_err_elasticnet <- matrix(NA, nrow = 1, ncol = FOLDS+4)
colnames(cv_err_elasticnet) <- c("dataset", "alpha", "model", "mean_err", paste("err", 1:FOLDS, sep = "_"))
cv_err_elasticnet <- as.data.frame(cv_err_elasticnet)

for (k in 1:FOLDS) { print(table(cv_data[[1]][which(folds==k, arr.ind=TRUE), 1])) }

if (file.exists("Errori_prev_FengQ.RData")) load("Errori_prev_FengQ.RData")

# ---------------------------------------------------------------------------- #

  # MODELLO LOGISTICO
  for(j in 1:length(cv_data)) {
    err <- NULL
    for(i in 1:FOLDS){
      testIndexes <- which(folds==i, arr.ind=TRUE)
      # Stimo il modello
      fit.model <- glm(study_condition ~ ., family=binomial, data=cv_data[[j]][-testIndexes, ])
      # Calcolo tasso di errata classificazione
      pred <- predict(fit.model, newdata=cv_data[[j]][testIndexes, ])
      tmp_tab <- table(pred, cv_data[[j]][testIndexes, ]$study_condition)
      err <- c(err, 1 - sum(diag(tmp_tab))/sum(tmp_tab))
      cat("Data", j, "Fold", i, "\n")
    }
    cv_err_glm <- rbind(cv_err_glm,
                       c(cv_data_names[j], "Modello logistico", mean(err), err))
  }
  
  cv_err_glm <- cv_err_glm[-1,]
  cv_err_glm
  # ---------------------------------------------------------------------------- #
  
  # PROJECTION PURSUIT REGRESSION (PPR)
  
  set.seed(28)
  for(j in 1:length(cv_data)) {
    for (HYP in c(2,4,6,8,10)) {
      err <- NULL
      for(i in 1:FOLDS){
        testIndexes <- which(folds==i, arr.ind=TRUE)
        # Stimo il modello
        fit.model <- ppr(as.numeric(cv_data[[j]][-testIndexes, 1] == "adenoma") ~ .,
                         data=cv_data[[j]][-testIndexes,-1], nterms=HYP)
        # Calcolo tasso di errata classificazione
        pred <- predict(fit.model, newdata=cv_data[[j]][testIndexes, ]) > 0.5
        tmp_tab <- table(pred, cv_data[[j]][testIndexes, ]$study_condition)
        err <- c(err, 1 - sum(diag(tmp_tab))/sum(tmp_tab))
        cat("Data", j, "Fold", i, "nterms", HYP, "\n")
      }
      cv_err_ppr <- rbind(cv_err_ppr,
                         c(cv_data_names[j], paste("PPR",HYP), mean(err), err))
    }
  } #non va oltre il primo dataset (quello clinico)
  cv_err_ppr <- cv_err_ppr[-1,]
  cv_err_ppr
  # ---------------------------------------------------------------------------- #
  
  # ALBERO DI CLASSIFICAZIONE
  library(tree)
  
  for(j in 1:length(cv_data)) {
    for (HYP in c(2,3,5,6,7,8,9,10,20,30,40,50,70,90,100)) {
      err <- NULL
      # Salto se il numero di variabili e' troppo piccolo
      #if (HYP>=NCOL(cv_data[[j]])) { next }
      for(i in 1:FOLDS){
        testIndexes <- which(folds==i, arr.ind=TRUE)
        # Stimo il modello
        fit.model <- tree(study_condition ~ ., data=cv_data[[j]][-testIndexes, ], split="deviance",
                          control=tree.control(nobs=nrow(cv_data[[j]][-testIndexes, ]),
                                               minsize=1,
                                               mindev=0))
        fit.model <- prune.tree(fit.model, best=HYP)
        # Calcolo tasso di errata classificazione
        pred <- predict(fit.model, newdata=cv_data[[j]][testIndexes, ], type="class")
        tmp_tab <- table(pred, cv_data[[j]][testIndexes, ]$study_condition)
        err <- c(err, 1 - sum(diag(tmp_tab))/sum(tmp_tab))
        cat("Data", j, "Fold", i, "J", HYP, "\n")
      }
      cv_err_cart <- rbind(cv_err_cart,
                         c(cv_data_names[j], paste("Albero di Classificazione", HYP), mean(err), err))
    }
  }
  cv_err_cart <- cv_err_cart[-1,]
  cv_err_cart
  
  'i = 5
  testIndexes <- which(folds==i, arr.ind=TRUE)
  dim(cv_data[[4]][-testIndexes, ])'
  
  # ---------------------------------------------------------------------------- #
  
  # BAGGING
  library(ipred)
  
  set.seed(496)
  for(j in 1:length(cv_data)) {
    for (HYP in 10*(5:15)) {
      err <- NULL
      for(i in 1:FOLDS){
        testIndexes <- which(folds==i, arr.ind=TRUE)
        # Stimo il modello
        fit.model <- bagging(study_condition ~ .,
                             data=cv_data[[j]][-testIndexes, ], nbagg=HYP)
        # Calcolo tasso di errata classificazione
        pred <- predict(fit.model, newdata=cv_data[[j]][testIndexes, ])
        tmp_tab <- table(pred, cv_data[[j]][testIndexes, ]$study_condition)
        err <- c(err, 1 - sum(diag(tmp_tab))/sum(tmp_tab))
        cat("Data", j, "Fold", i, "nbagg", HYP, "\n")
      }
      cv_err_bagging <- rbind(cv_err_bagging,
                         c(cv_data_names[j], paste("Bagging",HYP), mean(err), err))
    }
  }
  cv_err_bagging <- cv_err_bagging[-1,]
  cv_err_bagging
  # ---------------------------------------------------------------------------- #
  
  # RANDOM FOREST + BAGGING
  library(randomForest)
  
  set.seed(28)
  for(j in 1:length(cv_data)) {
    for (HYP in c(1,3,5,10,20,30,40,50,75,100,150,200)) {
      err <- NULL
      # Se MTRY e' troppo grande, salto il passaggio
      if (HYP>=NCOL(cv_data[[j]])) { next }
      for(i in 1:FOLDS){
        testIndexes <- which(folds==i, arr.ind=TRUE)
        # Stimo il modello
        fit.model <- randomForest(study_condition ~., data=cv_data[[j]][-testIndexes, ],
                                  nodesize=1, mtry=HYP, ntree=500, importance = T)
        # Calcolo tasso di errata classificazione
        pred <- predict(fit.model, newdata=cv_data[[j]][testIndexes, ])
        tmp_tab <- table(pred, cv_data[[j]][testIndexes, ]$study_condition)
        err <- c(err, 1 - sum(diag(tmp_tab))/sum(tmp_tab))
        cat("Data", j, "Fold", i, "mtry", HYP, "\n")
      }
      cv_err_rf <- rbind(cv_err_rf,
                         c(cv_data_names[j], paste("RandomForest",HYP), mean(err), err))
    }
  }
  cv_err_rf <- cv_err_rf[-1,]
  cv_err_rf
  
  # ---------------------------------------------------------------------------- #
  
  # BOOSTING
  library(ada)
  
  set.seed(496)
  for(j in 1:length(cv_data)) {
    err <- NULL
    for(i in 1:FOLDS){
      testIndexes <- which(folds==i, arr.ind=TRUE)
      # Stimo il modello
      fit.model <- ada(study_condition ~ ., data=cv_data[[j]][-testIndexes, ], iter=100)
      # Calcolo tasso di errata classificazione
      pred <- predict(fit.model, newdata=cv_data[[j]][testIndexes, ])
      tmp_tab <- table(pred, cv_data[[j]][testIndexes, ]$study_condition)
      err <- c(err, 1 - sum(diag(tmp_tab))/sum(tmp_tab))
      cat("Data", j, "Fold", i, "\n")
    }
    cv_err_boosting <- rbind(cv_err_boosting,
                       c(cv_data_names[j], "Boosting", mean(err), err))
  }
  cv_err_boosting <- cv_err_boosting[-1,]
  cv_err_boosting
  
  # ---------------------------------------------------------------------------- #
  
  # SUPPORT VECTOR MACHINES - RADIALE
  library(e1071)
  KERNEL <- "radial"
  
  set.seed(28)
  for(j in 1:length(cv_data)) {
    for (HYP in 0.5*(1:20)) {
      err <- NULL
      # Se MTRY e' troppo grande, salto il passaggio
      for(i in 1:FOLDS){
        testIndexes <- which(folds==i, arr.ind=TRUE)
        # Stimo il modello
        fit.model <- svm(study_condition ~ ., data=cv_data[[j]][-testIndexes, ],
                         cost=HYP, kernel=KERNEL)
        # Calcolo tasso di errata classificazione
        pred <- predict(fit.model, newdata=cv_data[[j]][testIndexes, ], decision.values=TRUE)
        tmp_tab <- table(pred, cv_data[[j]][testIndexes, ]$study_condition)
        err <- c(err, 1 - sum(diag(tmp_tab))/sum(tmp_tab))
        cat("Data", j, "Fold", i, "mtry", HYP, "\n")
      }
      cv_err_svm <- rbind(cv_err_svm,
                         c(cv_data_names[j], paste("SVM",KERNEL,HYP), mean(err), err))
    }
  }
  cv_err_svm <- cv_err_svm[-1,]
  cv_err_svm
  # ---------------------------------------------------------------------------- #
  
  # SUPPORT VECTOR MACHINES - SIGMOIDE
  library(e1071)
  KERNEL <- "sigmoid"
  
  set.seed(28)
  for(j in 1:length(cv_data)) {
    for (HYP in 0.5*(1:20)) {
      err <- NULL
      # Se MTRY e' troppo grande, salto il passaggio
      for(i in 1:FOLDS){
        testIndexes <- which(folds==i, arr.ind=TRUE)
        # Stimo il modello
        fit.model <- svm(study_condition ~ ., data=cv_data[[j]][-testIndexes, ],
                         cost=HYP, kernel=KERNEL)
        # Calcolo tasso di errata classificazione
        pred <- predict(fit.model, newdata=cv_data[[j]][testIndexes, ], decision.values=TRUE)
        tmp_tab <- table(pred, cv_data[[j]][testIndexes, ]$study_condition)
        err <- c(err, 1 - sum(diag(tmp_tab))/sum(tmp_tab))
        cat("Data", j, "Fold", i, "mtry", HYP, "\n")
      }
      cv_err_svm <- rbind(cv_err_svm,
                         c(cv_data_names[j], paste("SVM",KERNEL,HYP), mean(err), err))
    }
  }
  cv_err_svm[101:200,]
  
  # ---------------------------------------------------------------------------- #
  
  # ELASTIC NET
  library(glmnet)
  
  set.seed(28)
  for (ALPHA in c(0,0.2,0.4,0.6,0.8,1)) {
    for(j in 1:length(cv_data)) {
      for (HYP in 10^seq(3,-3,length=30)) {
        err <- NULL
        for(i in 1:FOLDS){
          testIndexes <- which(folds==i, arr.ind=TRUE)
          x   <- model.matrix(~ ., data=cv_data[[j]][-testIndexes,-1])
          x[,colnames(x)!="gendermale"] <- apply(x[,colnames(x)!="gendermale"], 2,
                                                 function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y)))
          x.v <- model.matrix(~ ., data=cv_data[[j]][testIndexes,-1])
          x.v[,colnames(x)!="gendermale"] <- apply(x.v[,colnames(x)!="gendermale"], 2,
                                                   function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y)))
          # Stimo il modello
          fit.lasso <- glmnet(x, cv_data[[j]][-testIndexes,1], family="binomial",
                              lambda=HYP, alpha=ALPHA)
          # Calcolo tasso di errata classificazione
          pred <- predict(fit.lasso, newx=x.v, type="response")
          tmp_tab <- table(pred, cv_data[[j]][testIndexes, ]$study_condition)
          err <- c(err, 1 - sum(diag(tmp_tab))/sum(tmp_tab))
          cat("Data", j, "Fold", i, "alpha", ALPHA, "hyp", HYP, "\n")
        }
        cv_err_elasticnet <- rbind(cv_err_elasticnet,
                                   c(cv_data_names[j], ALPHA, paste("ElasticNet",HYP), mean(err), err))
      }
    }
  }
  cv_err_elasticnet <- cv_err_elasticnet[-1,]
  cv_err_elasticnet

# ---------------------------------------------------------------------------- #

#Controllo errori

#GLM
cv_err_glm[,1:3] #fa tutto schifo

#PPR
cv_err_ppr[,1:3] #con HYP 2 mean_err = 0.523809523809524 -> solo cliniche (da errore con il resto)

#CART

#cliniche
cv_err_cart[cv_err_cart$dataset == "cliniche",][which.min(cv_err_cart
                      [cv_err_cart$dataset == "cliniche","mean_err"]),1:3] #0.2667
#specie - batteri
cv_err_cart[cv_err_cart$dataset == "specie - batteri",][which.min(cv_err_cart
                      [cv_err_cart$dataset == "specie - batteri","mean_err"]),1:3] #0.5333
#specie - cliniche + batteri
cv_err_cart[cv_err_cart$dataset == " specie - cliniche + batteri",][which.min(cv_err_cart
                      [cv_err_cart$dataset == " specie - cliniche + batteri","mean_err"]),1:3] #0.5047
#family - batteri
cv_err_cart[cv_err_cart$dataset == "family - batteri",][which.min(cv_err_cart
                      [cv_err_cart$dataset == "family - batteri","mean_err"]),1:3] #0.5048
# family - cliniche + batteri
cv_err_cart[cv_err_cart$dataset == " family - cliniche + batteri",][which.min(cv_err_cart
                      [cv_err_cart$dataset == " family - cliniche + batteri","mean_err"]),1:3] #0.4952

#BAGGING

#cliniche
cv_err_bagging[cv_err_bagging$dataset == "cliniche",][which.min(cv_err_bagging
                      [cv_err_bagging$dataset == "cliniche","mean_err"]),1:3] #0.2095
#specie - batteri
cv_err_bagging[cv_err_bagging$dataset == "specie - batteri",][which.min(cv_err_bagging
                      [cv_err_bagging$dataset == "specie - batteri","mean_err"]),1:3] #0.4381
#specie - cliniche + batteri
cv_err_bagging[cv_err_bagging$dataset == " specie - cliniche + batteri",][which.min(cv_err_bagging
                      [cv_err_bagging$dataset == " specie - cliniche + batteri","mean_err"]),1:3] #0.4381
#family - batteri
cv_err_bagging[cv_err_bagging$dataset == "family - batteri",][which.min(cv_err_bagging
                      [cv_err_bagging$dataset == "family - batteri","mean_err"]),1:3] #0.4381
# family - cliniche + batteri
cv_err_bagging[cv_err_bagging$dataset == " family - cliniche + batteri",][which.min(cv_err_bagging
                      [cv_err_bagging$dataset == " family - cliniche + batteri","mean_err"]),1:3] #0.3714

#RANDOM FOREST

#cliniche
cv_err_rf[cv_err_rf$dataset == "cliniche",][which.min(cv_err_rf
                      [cv_err_rf$dataset == "cliniche","mean_err"]),1:3] #0.2667
#specie - batteri
cv_err_rf[cv_err_rf$dataset == "specie - batteri",][which.min(cv_err_rf
                      [cv_err_rf$dataset == "specie - batteri","mean_err"]),1:3] #0.4381
#specie - cliniche + batteri
cv_err_rf[cv_err_rf$dataset == " specie - cliniche + batteri",][which.min(cv_err_rf
                      [cv_err_rf$dataset == " specie - cliniche + batteri","mean_err"]),1:3] #0.4286
#family - batteri
cv_err_rf[cv_err_rf$dataset == "family - batteri",][which.min(cv_err_rf
                      [cv_err_rf$dataset == "family - batteri","mean_err"]),1:3] #0.4476
# family - cliniche + batteri
cv_err_rf[cv_err_rf$dataset == " family - cliniche + batteri",][which.min(cv_err_rf
                      [cv_err_rf$dataset == " family - cliniche + batteri","mean_err"]),1:3] #0.4190

#BOOSTING

#cliniche
cv_err_boosting[cv_err_boosting$dataset == "cliniche",][which.min(cv_err_boosting
                      [cv_err_boosting$dataset == "cliniche","mean_err"]),1:3] #0.3048
#specie - batteri
cv_err_boosting[cv_err_boosting$dataset == "specie - batteri",][which.min(cv_err_boosting
                      [cv_err_boosting$dataset == "specie - batteri","mean_err"]),1:3] #0.5333
#specie - cliniche + batteri
cv_err_boosting[cv_err_boosting$dataset == " specie - cliniche + batteri",][which.min(cv_err_boosting
                      [cv_err_boosting$dataset == " specie - cliniche + batteri","mean_err"]),1:3] #0.4190
#family - batteri
cv_err_boosting[cv_err_boosting$dataset == "family - batteri",][which.min(cv_err_boosting
                      [cv_err_boosting$dataset == "family - batteri","mean_err"]),1:3] #0.4286
# family - cliniche + batteri
cv_err_boosting[cv_err_boosting$dataset == " family - cliniche + batteri",][which.min(cv_err_boosting
                      [cv_err_boosting$dataset == " family - cliniche + batteri","mean_err"]),1:3] #0.4381
#SVM

#cliniche
cv_err_svm[cv_err_svm$dataset == "cliniche",][which.min(cv_err_svm
                      [cv_err_svm$dataset == "cliniche","mean_err"]),1:3] #0.3714
#specie - batteri
cv_err_svm[cv_err_svm$dataset == "specie - batteri",][which.min(cv_err_svm
                      [cv_err_svm$dataset == "specie - batteri","mean_err"]),1:3] #0.4476
#specie - cliniche + batteri
cv_err_svm[cv_err_svm$dataset == " specie - cliniche + batteri",][which.min(cv_err_svm
                      [cv_err_svm$dataset == " specie - cliniche + batteri","mean_err"]),1:3] #0.4476
#family - batteri
cv_err_svm[cv_err_svm$dataset == "family - batteri",][which.min(cv_err_svm
                      [cv_err_svm$dataset == "family - batteri","mean_err"]),1:3] #0.4476
# family - cliniche + batteri
cv_err_svm[cv_err_svm$dataset == " family - cliniche + batteri",][which.min(cv_err_svm
                      [cv_err_svm$dataset == " family - cliniche + batteri","mean_err"]),1:3] #0.4476
#RIDGE: alpha = 0
#cliniche
cv_err_elasticnet[cv_err_elasticnet$dataset == "cliniche" & cv_err_elasticnet$alpha == "0",][
  which.min(cv_err_elasticnet[cv_err_elasticnet$dataset == "cliniche" &
                                cv_err_elasticnet$alpha == "0","mean_err"]),1:4] #0.9429
#specie - batteri
cv_err_elasticnet[cv_err_elasticnet$dataset == "specie - batteri" & cv_err_elasticnet$alpha == "0",][
  which.min(cv_err_elasticnet[cv_err_elasticnet$dataset == "specie - batteri" &
                                cv_err_elasticnet$alpha == "0","mean_err"]),1:4] #0.9619
#specie - cliniche + batteri
cv_err_elasticnet[cv_err_elasticnet$dataset == " specie - cliniche + batteri" & cv_err_elasticnet$alpha == "0",][
  which.min(cv_err_elasticnet[cv_err_elasticnet$dataset == " specie - cliniche + batteri" &
                                cv_err_elasticnet$alpha == "0","mean_err"]),1:4] #0.9619
#family - batteri
cv_err_elasticnet[cv_err_elasticnet$dataset == "family - batteri" & cv_err_elasticnet$alpha == "0",][
  which.min(cv_err_elasticnet[cv_err_elasticnet$dataset == "family - batteri" &
                                cv_err_elasticnet$alpha == "0","mean_err"]),1:4] #0.9333
# family - cliniche + batteri
cv_err_elasticnet[cv_err_elasticnet$dataset == " family - cliniche + batteri" & cv_err_elasticnet$alpha == "0",][
  which.min(cv_err_elasticnet[cv_err_elasticnet$dataset == " family - cliniche + batteri" &
                                cv_err_elasticnet$alpha == "0","mean_err"]),1:4] #0.9429
#LASSO: alpha = 1
#cliniche
cv_err_elasticnet[cv_err_elasticnet$dataset == "cliniche" & cv_err_elasticnet$alpha == "1",][
  which.min(cv_err_elasticnet[cv_err_elasticnet$dataset == "cliniche" &
                                cv_err_elasticnet$alpha == "1","mean_err"]),1:4] #0.5524
#specie - batteri
cv_err_elasticnet[cv_err_elasticnet$dataset == "specie - batteri" & cv_err_elasticnet$alpha == "1",][
  which.min(cv_err_elasticnet[cv_err_elasticnet$dataset == "specie - batteri" &
                                cv_err_elasticnet$alpha == "1","mean_err"]),1:4] #0.5524
#specie - cliniche + batteri
cv_err_elasticnet[cv_err_elasticnet$dataset == " specie - cliniche + batteri" & cv_err_elasticnet$alpha == "1",][
  which.min(cv_err_elasticnet[cv_err_elasticnet$dataset == " specie - cliniche + batteri" &
                                cv_err_elasticnet$alpha == "1","mean_err"]),1:4] #0.5524
#family - batteri
cv_err_elasticnet[cv_err_elasticnet$dataset == "family - batteri" & cv_err_elasticnet$alpha == "1",][
  which.min(cv_err_elasticnet[cv_err_elasticnet$dataset == "family - batteri" &
                                cv_err_elasticnet$alpha == "1","mean_err"]),1:4] #0.5524
# family - cliniche + batteri
cv_err_elasticnet[cv_err_elasticnet$dataset == " family - cliniche + batteri" & cv_err_elasticnet$alpha == "1",][
  which.min(cv_err_elasticnet[cv_err_elasticnet$dataset == " family - cliniche + batteri" &
                                cv_err_elasticnet$alpha == "1","mean_err"]),1:4] #0.5523
#ELASTIC NET

#cliniche
cv_err_elasticnet[cv_err_elasticnet$dataset == "cliniche",][which.min(cv_err_elasticnet
                      [cv_err_elasticnet$dataset == "cliniche","mean_err"]),1:4] #0.5524
#specie - batteri
cv_err_elasticnet[cv_err_elasticnet$dataset == "specie - batteri",][which.min(cv_err_elasticnet
                      [cv_err_elasticnet$dataset == "specie - batteri","mean_err"]),1:4] #0.5524
#specie - cliniche + batteri
cv_err_elasticnet[cv_err_elasticnet$dataset == " specie - cliniche + batteri",][which.min(cv_err_elasticnet
                      [cv_err_elasticnet$dataset == " specie - cliniche + batteri","mean_err"]),1:4] #0.5524
#family - batteri
cv_err_elasticnet[cv_err_elasticnet$dataset == "family - batteri",][which.min(cv_err_elasticnet
                      [cv_err_elasticnet$dataset == "family - batteri","mean_err"]),1:4] #0.5524
# family - cliniche + batteri
cv_err_elasticnet[cv_err_elasticnet$dataset == " family - cliniche + batteri",][which.min(cv_err_elasticnet
                      [cv_err_elasticnet$dataset == " family - cliniche + batteri","mean_err"]),1:4] #0.5524


# ---------------------------------------------------------------------------- #

#Salvo risultati intermedi

save(cv_err_glm, cv_err_ppr, cv_err_cart, cv_err_bagging, cv_err_rf, cv_err_boosting, cv_err_svm,
     cv_err_elasticnet, file = "Errori_prev_FengQ.RData")
# ---------------------------------------------------------------------------- #