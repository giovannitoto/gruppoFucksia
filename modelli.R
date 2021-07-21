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
set.seed(28)
idx <- sample(1:nrow(sb0), trunc(nrow(sb1)*0.6/0.4))
sb0 <- sb0[idx, ]
# Sostituisco data_fix con la versione bilanciata
s <- rbind(sb1, sb0)  
prop.table(table(s$study_condition))
rm(list=c("sb1","sb0","idx"))

prop.table(table(s$study_condition))
prop.table(table(v$study_condition))

# ---------------------------------------------------------------------------- #

# Dataset di stima e convalida
set.seed(28)
cb1 <- sample(1:nrow(s), trunc(nrow(s)*2/3))
cb2 <- setdiff(1:nrow(s), cb1)

# Definisco tabella per il confronto dei modelli
tab_confronto <- c("DATA", "MODELLO", "TASSO ERRATA CLASSIFICAZIONE")

tab_confronto <- matrix(NA, nrow = 1, ncol=3)
colnames(tab_confronto) <- c("dataset", "model", "err")
tab_confronto <- as.data.frame(tab_confronto)

# ---------------------------------------------------------------------------- #

datasets_s <- list(s[,c(1,2,3,7)], s[,c(1,4,5,6, 8:506)], s[,1:506],
                   s[,c(1,4,5,6, 507:571)], s[, c(1:7, 506:571)])
datasets_v <- list(v[,c(1,2,3,7)], v[,c(1,4,5,6, 8:506)], v[,1:506],
                   v[,c(1,4,5,6, 507:571)], v[,c(1:7, 506:571)])
datasets_names <- c("specie - cliniche", "specie - batteri", " specie - cliniche + batteri",
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
  #x[,colnames(x)!="gendermale"] <- scale(x[,colnames(x)!="gendermale"])
  x.v <- model.matrix(~ ., data=datasets_v[[j]][,-1])
  #x.v[,colnames(x)!="gendermale"] <- scale(x.v[,colnames(x)!="gendermale"])
  
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

for (HYP in c(0,0.25,0.5,0.75)) {
  for (j in 1:length(datasets_names)) {
    x   <- model.matrix(~ ., data=datasets_s[[j]][,-1])
    #x[,colnames(x)!="gendermale"] <- scale(x[,colnames(x)!="gendermale"])
    x.v <- model.matrix(~ ., data=datasets_v[[j]][,-1])
    #x.v[,colnames(x)!="gendermale"] <- scale(x.v[,colnames(x)!="gendermale"])
    
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

# RANDOM FOREST (CV)
library(randomForest)
library(e1071)

set.seed(496)
for (j in 1:length(datasets_names)) {
  forest.tune <- tune(randomForest, study_condition ~ ., data=datasets_s[[j]],
                      ranges=list(mtry=2*(2:6)),
                      tunecontrol=tune.control(cross=5),
                      ntree=100)
  best.mtry <- min(forest.tune$best.parameters, NCOL(s))
  fit.model <- randomForest(study_condition ~ ., data=datasets_s[[j]], mtry=best.mtry)
  # Effettuo previsioni
  y.model <- predict(fit.model, newdata=datasets_v[[j]])
  et.model <- table(y.model, datasets_v[[j]]$study_condition)
  e.model <- 1 - sum(diag(et.model))/sum(et.model)
  tab_confronto <- rbind(tab_confronto,
                         c(datasets_names[[j]],paste("Random Forest - mtry", best.mtry), e.model))
  cat("Data", j, "\n")
}




# ---------------------------------------------------------------------------- #















