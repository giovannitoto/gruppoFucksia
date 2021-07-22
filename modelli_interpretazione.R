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
datasets_names <- c("cliniche", "specie: batteri", "specie: cliniche+batteri",
                    "family: batteri", "family: cliniche+batteri")

# ---------------------------------------------------------------------------- #

source("lift-roc-tab.R")

# ---------------------------------------------------------------------------- #

# RANDOM FOREST
library(randomForest)

# Parametri:
# - nodesize specifica quanto fare crescere ciascun albero
# - mtry specifica il numero di variabili per nodo
# - xtest specifica le covariate di un insieme di verifica
# - ytest specifica la risposta di un insieme di verifica

# Seleziono mtry e ristimo foresta sull'intero dataset di stima
set.seed(28)
fit.forest <- randomForest(study_condition ~ ., data=datasets_s[[4]],
                           nodesize=1, mtry=39,
                           xtest=datasets_v[[4]][,-1], ytest=datasets_v[[4]][,1])

# Grafico per importanza delle var. nella random forest
varImpPlot(fit.forest)

# fit.forest$importance


# Previsione
et.forest <- tabella.sommario(fit.forest$test$predicted, datasets_v[[4]]$study_condition)
e.forest <- 1 - sum(diag(et.forest))/sum(et.forest)
a.forest <- lift.roc(fit.forest$test$votes[,2],
                     as.numeric(datasets_v[[4]]$study_condition=="control"), type="crude")

# ---------------------------------------------------------------------------- #

### BOOSTING
library(ada)

fit.boosting <- ada(study_condition ~ ., data=datasets_s[[4]], iter=200,
                    test.x=datasets_v[[4]][,-1],
                    test.y=datasets_v[[4]][,1])

# Grafico dell'importanza delle variabili
varplot(fit.boosting)

# Previsione
y.boosting <- predict(fit.boosting, newdata=datasets_v[[4]])
et.boosting <- tabella.sommario(y.boosting, datasets_v[[4]]$study_condition)
y.boosting.prob <- predict(fit.boosting, newdata=v, type="prob")
a.boosting <- lift.roc(y.boosting.prob[,2],
                       as.numeric(datasets_v[[4]]$study_condition=="control"), type="crude")


# ---------------------------------------------------------------------------- #
