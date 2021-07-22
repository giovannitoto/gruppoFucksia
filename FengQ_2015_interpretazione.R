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

source("lift-roc-tab.R")

# ---------------------------------------------------------------------------- #

### ALBERO DI CLASSIFICAZIONE
library(tree)

# Voglio dividere finche' non avro' due punti per foglia: mi servono opzioni di controllo
# In particolare:
# -    nobs : numero di oss. del dataset
# - minsize : numero massimo di oss. per ogni foglia
# -  mindev : la devianza minima affinche' ci sia un ulteriore suddivisione
fit.tree <- tree(study_condition ~ ., data=cv_data[[1]], split="deviance",
                 control=tree.control(nobs=nrow(cv_data[[1]]),
                                      minsize=1,
                                      mindev=0))
# Stimo modello finale
final.tree <- prune.tree(fit.tree, best=8)

# Rappresentazione grafica
plot(final.tree, type="uniform")
text(final.tree, pretty=1, cex=0.75, label="yprob")

# ---------------------------------------------------------------------------- #

# RANDOM FOREST
library(randomForest)

# Parametri:
# - nodesize specifica quanto fare crescere ciascun albero
# - mtry specifica il numero di variabili per nodo
# - xtest specifica le covariate di un insieme di verifica
# - ytest specifica la risposta di un insieme di verifica

set.seed(28)
fit.forest <- randomForest(study_condition ~ ., data=cv_data[[1]], nodesize=1, mtry=3)
# Grafico degli errori ottenuti con la tecnica out-of-bag
# -  nero : errore complessivo
# - rosso : errore su quelli che non se ne vanno
# - verde : errore su quelli che se ne vanno
plot(fit.forest, main="OOB")

# Grafico per importanza delle var. nella random forest
varImpPlot(fit.forest)

# ---------------------------------------------------------------------------- #