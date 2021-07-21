# Modelli con risposta "study_condition"
# - considerando aggregazione a livello di specie

rm(list=ls())
load("dataset_fix_fixcl.RData")

# Tolgo variabili collineari
data_fix <- data_fix[,-c(419+8, 420+8, 493+8)]

# ---------------------------------------------------------------------------- #

# K-fold Cross-Validation
FOLDS <- 10

folds <- sample(cut(seq(1,nrow(data_fix)), breaks=FOLDS, labels=FALSE))

cv_errors <- matrix(NA, nrow = 1, ncol = FOLDS+3)
colnames(cv_errors) <- c("dataset", "model", "mean_err", paste("err", 1:FOLDS, sep = "_"))
cv_errors <- as.data.frame(cv_errors)

cv_data <- list(data_fix[,c(1,2,3,4,8)], data_fix[,c(1,5,6,7,9:507)], data_fix)
cv_data_names <- c("cliniche", "batteri", "cliniche + batteri")

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
  cv_errors <- rbind(cv_errors, c(cv_data_names[j], "Modello logistico", mean(err), err))
}

# ---------------------------------------------------------------------------- #

# ALBERO DI CLASSIFICAZIONE
library(tree)

for(j in 1:length(cv_data)) {
  err <- NULL
  for (J in c(2,5,10,15,20,30,40,50,60,70,80,90,100)) {
    # Salto se il numero di variabili e' troppo piccolo
    if (J>=NCOL(cv_data[[j]])) { next }
    for(i in 1:FOLDS){
      testIndexes <- which(folds==i, arr.ind=TRUE)
      # Stimo il modello
      fit.model <- tree(study_condition ~ ., data=cv_data[[j]][-testIndexes, ], split="deviance",
                        control=tree.control(nobs=nrow(cv_data[[j]][-testIndexes, ]),
                                             minsize=1,
                                             mindev=0))
      fit.model <- prune.tree(fit.model, best=J)
      # Calcolo tasso di errata classificazione
      pred <- predict(fit.model, newdata=cv_data[[j]][testIndexes, ], type="class")
      tmp_tab <- table(pred, cv_data[[j]][testIndexes, ]$study_condition)
      err <- c(err, 1 - sum(diag(tmp_tab))/sum(tmp_tab))
      cat("Data", j, "Fold", i, "J", J, "\n")
    }
    cv_errors <- rbind(cv_errors, c(cv_data_names[j], paste("Albero di Classificazione", J), mean(err), err))
  }
}

# ---------------------------------------------------------------------------- #

# RANDOM FOREST + BAGGING CON OTTIMIZZAZIONE
library(randomForest)

for(j in 1:length(cv_data)) {
  err <- NULL
  for (MTRY in c(5,10,20,30,40,50,60,70,80,90,100,200,300,400)) {
    # Se MTRY e' troppo grande, salto il passaggio
    if (MTRY>=NCOL(cv_data[[j]])) { next }
    for(i in 1:FOLDS){
      testIndexes <- which(folds==i, arr.ind=TRUE)
      set.seed(28)
      rf1 <- randomForest(study_condition ~., data=cv_data[[j]][-testIndexes, ], nodesize=1, 
                          mtry=3, ntree=500, importance = T)
      # Calcolo tasso di errata classificazione
      pred <- predict(rf1, newdata=cv_data[[j]][testIndexes, ])
      tmp_tab <- table(pred, cv_data[[j]][testIndexes, ]$study_condition)
      err <- c(err, 1 - sum(diag(tmp_tab))/sum(tmp_tab))
      cat("Data", j, "Fold", i, "mtry", MTRY, "\n")
    }
    cv_errors <- rbind(cv_errors, c(cv_data_names[j], paste("RandomForest -",MTRY), mean(err), err))
  }
}

# ---------------------------------------------------------------------------- #

