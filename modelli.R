# Modelli con risposta "study_condition"
# - considerando aggregazione a livello di specie

load("dataset_fix_fixcl.RData")

# Tolgo varaibili collineari
data_fix <- data_fix[,-c(419+8, 420+8, 493+8)]
data_fix_cl <- data_fix_cl[,-c(419+5, 420+5, 493+5)]

prop.table(table(data_fix$study_condition))

str(data_fix)

# ---------------------------------------------------------------------------- #
# K-fold Cross-Validation
FOLDS <- 10

folds <- sample(cut(seq(1,nrow(data_fix)), breaks=FOLDS, labels=FALSE))

cv_errors <- matrix(NA, nrow = 1, ncol = FOLDS+3)
colnames(cv_errors) <- c("dataset", "model", "mean_err", paste("err", 1:FOLDS, sep = "_"))
cv_errors <- as.data.frame(cv_errors)

cv_data <- list(data_fix, data_fix[,-c(9:ncol(data_fix))], data_fix[,c(1,9:ncol(data_fix))],
                data_fix_cl, data_fix_cl[,-c(6:ncol(data_fix_cl))], data_fix_cl[,c(1,6:ncol(data_fix_cl))])
cv_data_names <- c("data_fix", "data_fix_nobact", "data_fix_bact", "data_cl", "data_cl_nobact", "data_cl_bact")

#########################################################################################
#########################################################################################
# RANDOM FOREST + BAGGING
library(randomForest)

# SENZA OTTIMIZZAZIONE
# ---------------------------------------------------------------------------- #
for(j in 1:length(cv_data)){
  err <- NULL
  for(i in 1:FOLDS){
    testIndexes <- which(folds==i, arr.ind=TRUE)
    set.seed(123)
    rf1 <- randomForest(study_condition ~., data=cv_data[[j]][-testIndexes, ], nodesize=1, 
                        mtry=3, ntree=500, importance = T)
    # Calcolo tasso di errata classificazione
    pred <- predict(rf1, newdata=cv_data[[j]][testIndexes, ])
    tmp_tab <- table(pred, cv_data[[j]][testIndexes, ]$study_condition)
    err <- c(err, 1 - sum(diag(tmp_tab))/sum(tmp_tab))
    cat("Data", j, "Fold", i, "\n")
  }
  cv_errors <- rbind(cv_errors, c(cv_data_names[j], "RandomForest", mean(err), err))
}

cv_errors

# ---------------------------------------------------------------------------- #

