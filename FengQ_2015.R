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

cv_data <- list(data_fix2[,c(1,2,3,7,8,9,10)], data_fix2[,c(1,4,5,6, 11:315)],
                data_fix2[,1:315], data_fix2[,c(1,4,5,6, 316:361)],
                data_fix2[,c(1:10, 316:361)])

cv_data_names <- c("specie - cliniche", "specie - batteri", " specie - cliniche + batteri",
                    "family - batteri", " family - cliniche + batteri")



# ---------------------------------------------------------------------------- #

# K-fold Cross-Validation
FOLDS <- 5

folds <- sample(cut(seq(1,nrow(data_fix2)), breaks=FOLDS, labels=FALSE))

cv_errors <- matrix(NULL, nrow = 1, ncol = FOLDS+3)
colnames(cv_errors) <- c("dataset", "model", "mean_err", paste("err", 1:FOLDS, sep = "_"))
cv_errors <- as.data.frame(cv_errors)


for (k in 1:FOLDS) { print(table(cv_data[[1]][which(folds==k, arr.ind=TRUE), 1])) }


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
  cv_errors <- rbind(cv_errors,
                     c(cv_data_names[j], "Modello logistico", mean(err), err))
}


cv_errors
# ---------------------------------------------------------------------------- #

# PROJECTION PURSUIT REGRESSION (PPR)

set.seed(28)
for(j in 1:length(cv_data)) {
  for (HYP in c(2,4,6,8,10)) {
    err <- NULL
    for(i in 1:FOLDS){
      testIndexes <- which(folds==i, arr.ind=TRUE)
      # Stimo il modello
      fit.model <- ppr(as.numeric(study_condition=="adenoma") ~ .,
                       data=cv_data[[j]][-testIndexes,-1], nterms=HYP)
      # Calcolo tasso di errata classificazione
      pred <- predict(fit.model, newdata=cv_data[[j]][testIndexes, ]) > 0.5
      tmp_tab <- table(pred, cv_data[[j]][testIndexes, ]$study_condition)
      err <- c(err, 1 - sum(diag(tmp_tab))/sum(tmp_tab))
      cat("Data", j, "Fold", i, "nterms", HYP, "\n")
    }
    cv_errors <- rbind(cv_errors,
                       c(cv_data_names[j], paste("PPR",HYP), mean(err), err))
  }
}

# ---------------------------------------------------------------------------- #

# ALBERO DI CLASSIFICAZIONE
library(tree)

for(j in 1:length(cv_data)) {
  for (HYP in c(2,3,5,6)) {
    err <- NULL
    # Salto se il numero di variabili e' troppo piccolo
    if (HYP>=NCOL(cv_data[[j]])) { next }
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
    cv_errors <- rbind(cv_errors,
                       c(cv_data_names[j], paste("Albero di Classificazione", HYP), mean(err), err))
  }
}

cv_errors 

i = 5
testIndexes <- which(folds==i, arr.ind=TRUE)
dim(cv_data[[4]][-testIndexes, ])

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
    cv_errors <- rbind(cv_errors,
                       c(cv_data_names[j], paste("Bagging",HYP), mean(err), err))
  }
}

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
    cv_errors <- rbind(cv_errors,
                       c(cv_data_names[j], paste("RandomForest",HYP), mean(err), err))
  }
}

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
  cv_errors <- rbind(cv_errors,
                     c(cv_data_names[j], "Boosting", mean(err), err))
}
cv_errors
# ---------------------------------------------------------------------------- #

# SUPPORT VECTOR MACHINES - RADIALE
library(e1071)
KERNEL <- "radial"

set.seed(28)
for(j in 1:length(cv_data)) {
  for (HYP in c(1,2,3,4,5,6)) {
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
    cv_errors <- rbind(cv_errors,
                       c(cv_data_names[j], paste("SVM",KERNEL,HYP), mean(err), err))
  }
}

# ---------------------------------------------------------------------------- #

# SUPPORT VECTOR MACHINES - SIGMOIDE
library(e1071)
KERNEL <- "sigmoid"

set.seed(28)
for(j in 1:length(cv_data)) {
  for (HYP in c(1,2,3,4,5,6)) {
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
    cv_errors <- rbind(cv_errors,
                       c(cv_data_names[j], paste("SVM",KERNEL,HYP), mean(err), err))
  }
}

# ---------------------------------------------------------------------------- #











# ---------------------------------------------------------------------------- #