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








