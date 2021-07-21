# Modelli con risposta "study_condition"
# - considerando aggregazione a livello di specie

rm(list=ls())
load("dataset_fix.RData")

stima <- data_fix[data_fix$study_name != "FengQ_2015", ]
verifica <- data_fix[data_fix$study_name == "FengQ_2015", ]

# Bilancio il dataset in modo da avere una proporzione 40%/60%
stima_fix1 <- stima[stima$study_condition=="adenoma", ]
stima_fix0 <- stima[stima$study_condition=="control", ]
# In pratica, scarto oss. dal gruppo piu' numeroso in modo da avere gruppi bilanciati
# In questo caso sb0 ha troppe oss. => scarto quelle in eccesso
set.seed(28)
idx <- sample(1:nrow(stima_fix0), trunc(nrow(stima_fix1)*0.6/0.4))
stima_fix0 <- stima_fix0[idx, ]
# Sostituisco data_fix con la versione bilanciata
stima <- rbind(stima_fix1, stima_fix0)  
prop.table(table(stima$study_condition))
rm(list=c("stima_fix1","stima_fix0","idx"))

prop.table(table(stima$study_condition))
# adenoma control 
# 0.4     0.6 
prop.table(table(verifica$study_condition))

# elimino study_name
stima <- stima[, -1]
verifica <- verifica[, -1]
# ---------------------------------------------------------------------------- #

# K-fold Cross-Validation
FOLDS <- 5

folds <- sample(cut(seq(1,nrow(data_fix)), breaks=FOLDS, labels=FALSE))

cv_errors <- matrix(NA, nrow = 1, ncol = FOLDS+3)
colnames(cv_errors) <- c("dataset", "model", "mean_err", paste("err", 1:FOLDS, sep = "_"))
cv_errors <- as.data.frame(cv_errors)

cv_data <- list(data_fix[,c(1,2,3,4,8)], data_fix[,c(1,5,6,7,9:507)], data_fix)
cv_data_names <- c("cliniche", "batteri", "cliniche + batteri")


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
  for (HYP in c(2,5,10,15,20,30,40,50,60,70,80,90,100)) {
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