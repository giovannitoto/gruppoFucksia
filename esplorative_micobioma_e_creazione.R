# Analisi esplorrativa sul microbioma + modelli che sfruttano la tassonomia
load("dataset_fix.RData")

# Guardo la tassonomia
colnames(rowData(se_fix))
dim(assay(se_fix))

microb <- as.data.frame(cbind(assay(se_fix), rowData(se_fix)))
str(microb) # prime 627 colonne sono le percentuali, poi le variabili delle specie
microb$Kingdom <- as.factor(microb$Kingdom)
microb$Phylum <- as.factor(microb$Phylum)
microb$Class <- as.factor(microb$Class)
microb$Order <- as.factor(microb$Order)
microb$Family <- as.factor(microb$Family)
microb$Genus <- as.factor(microb$Genus)
microb$Species <- as.factor(microb$Species)

# Kingdom
table(microb$Kingdom)
livelli <- levels(microb$Kingdom)
sum_sogg <- matrix(NA, nrow = length(livelli), ncol = 627)
for(i in 1:length(livelli)){
  subs <- microb[microb$Kingdom == livelli[i], 1:627]
  sum_sogg[i,] <- colSums(subs)
}
sum_sogg <- as.data.frame(t(sum_sogg))
colnames(sum_sogg) <- livelli
boxplot(sum_sogg, outline = F)

# Non ha troppo senso perché ho solo batteri

# Phylum
table(microb$Phylum)
livelli <- levels(microb$Phylum)
sum_sogg <- matrix(NA, nrow = length(livelli), ncol = 627)
for(i in 1:length(livelli)){
  subs <- microb[microb$Phylum == livelli[i], 1:627]
  sum_sogg[i,] <- colSums(subs)
}
sum_sogg <- as.data.frame(t(sum_sogg))
colnames(sum_sogg) <- livelli
boxplot(sum_sogg, outline = F)

# Class
table(microb$Class)
dim(table(microb$Class))
livelli <- levels(microb$Class)
sum_sogg <- matrix(NA, nrow = length(livelli), ncol = 627)
for(i in 1:length(livelli)){
  subs <- microb[microb$Class == livelli[i], 1:627]
  sum_sogg[i,] <- colSums(subs)
}
sum_sogg <- as.data.frame(t(sum_sogg))
colnames(sum_sogg) <- livelli
boxplot(sum_sogg, outline = F)

# Order
table(microb$Order)
dim(table(microb$Order))
livelli <- levels(microb$Order)
sum_sogg <- matrix(NA, nrow = length(livelli), ncol = 627)
for(i in 1:length(livelli)){
  subs <- microb[microb$Order == livelli[i], 1:627]
  sum_sogg[i,] <- colSums(subs)
}
sum_sogg <- as.data.frame(t(sum_sogg))
colnames(sum_sogg) <- livelli
boxplot(sum_sogg, outline = F)

# Family
table(microb$Family)
dim(table(microb$Family))
livelli <- levels(microb$Family)
sum_sogg <- matrix(NA, nrow = length(livelli), ncol = 627)
for(i in 1:length(livelli)){
  subs <- microb[microb$Family == livelli[i], 1:627]
  sum_sogg[i,] <- colSums(subs)
}
sum_sogg <- as.data.frame(t(sum_sogg))
colnames(sum_sogg) <- livelli
boxplot(sum_sogg, outline = F)

# Genus
table(microb$Genus)
dim(table(microb$Genus))
livelli <- levels(microb$Genus)
sum_sogg <- matrix(NA, nrow = length(livelli), ncol = 627)
for(i in 1:length(livelli)){
  subs <- microb[microb$Genus == livelli[i], 1:627]
  sum_sogg[i,] <- colSums(subs)
}
sum_sogg <- as.data.frame(t(sum_sogg))
colnames(sum_sogg) <- livelli
boxplot(sum_sogg, outline = F)


##### ACCORPO PER FAMILY su DATI SOTTOINSIEME
################### Dati sottinsieme -> accorpo per FAMILY
load("dataset_fix2_317.RData")

microb3 <- as.data.frame(cbind(assay(se_fix3), rowData(se_fix3)))
str(microb3) # prime 627 colonne sono le percentuali, poi le variabili delle specie
microb3$Kingdom <- as.factor(microb3$Kingdom)
microb3$Phylum <- as.factor(microb3$Phylum)
microb3$Class <- as.factor(microb3$Class)
microb3$Order <- as.factor(microb3$Order)
microb3$Family <- as.factor(microb3$Family)
microb3$Genus <- as.factor(microb3$Genus)
microb3$Species <- as.factor(microb3$Species)

table(microb3$Family)
dim(table(microb3$Family))
livelli <- levels(microb3$Family)
sum_sogg <- matrix(NA, nrow = length(livelli), ncol = 317)
for(i in 1:length(livelli)){
  subs <- microb3[microb3$Family == livelli[i], 1:317]
  sum_sogg[i,] <- colSums(subs)
}
sum_sogg <- as.data.frame(t(sum_sogg))
colnames(sum_sogg) <- livelli
str(sum_sogg) # ho 59 specie

str(data_fix3)
data_fix3_family <- data_fix3[,-c(10:407)]
dim(data_fix3_family)
data_fix3_family <- cbind(data_fix3_family, sum_sogg)
str(data_fix3_family)

save(data_fix3_family, file = "data_fix317_family.RData")

##### ACCORPO PER FAMILY su DATI TOTALI
################### Dati totali -> accorpo per FAMILY
load("dataset_fix.RData")

microb <- as.data.frame(cbind(assay(se_fix), rowData(se_fix)))
str(microb) # prime 627 colonne sono le percentuali, poi le variabili delle specie
microb$Kingdom <- as.factor(microb$Kingdom)
microb$Phylum <- as.factor(microb$Phylum)
microb$Class <- as.factor(microb$Class)
microb$Order <- as.factor(microb$Order)
microb$Family <- as.factor(microb$Family)
microb$Genus <- as.factor(microb$Genus)
microb$Species <- as.factor(microb$Species)

table(microb$Family)
dim(table(microb$Family))
livelli <- levels(microb$Family)
sum_sogg <- matrix(NA, nrow = length(livelli), ncol = 627)
for(i in 1:length(livelli)){
  subs <- microb[microb$Family == livelli[i], 1:627]
  sum_sogg[i,] <- colSums(subs)
}
sum_sogg <- as.data.frame(t(sum_sogg))
colnames(sum_sogg) <- livelli
str(sum_sogg) # ho 65 specie

str(data_fix)
data_fix_family <- data_fix[,-c(9:507)] # tolgo le specie
dim(data_fix_family)
data_fix_family <- cbind(data_fix_family, sum_sogg)
str(data_fix_family)

save(data_fix_family, file = "data_fix_family.RData")
