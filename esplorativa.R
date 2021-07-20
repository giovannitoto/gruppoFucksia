load("adenoma.rda")

head(rowData(se))

# Oggetto generale
se
dim(se)
# le righe sono i soggetti, le colonne le specie di batteri (vero?)

# Guardo colData, ovvero le caratteristiche degli individui
dim(colData(se))
colData(se)

#Creo una copia del dataset
se_fix <- se
'Guardo le covariate una per una'

# study_name
table(is.na(se$study_name))
table(se$study_name) # abbiamo 5 studi

#-> non utile per prevedere la risposta
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "study_name")]

# subject_id
table(is.na(se$subject_id))
length(unique(se$subject_id)) # tutti i soggetti osservati una sola volta

#-> non utile per prevedere la risposta
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "subject_id")]

# body_site
table(is.na(se$body_site))
table(se$body_site) # sono tutti campioni di feci

#-> non utile per prevedere la risposta
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "body_site")]

# antibiotics_current_use
table(is.na(se$antibiotics_current_use)) # ha tanti NA, da droppare
table(se$antibiotics_current_use)

#-> non utile per prevedere la risposta
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "antibiotics_current_use")]

# study_condition -> variabile risposta
table(is.na(se$study_condition))
table(se$study_condition)
prop.table(table(se$study_condition))
str(se$study_condition)
se_fix$study_condition <- as.factor(se_fix$study_condition)

# disease -> altra possibile variabile risposta, assolutamente non una esplicativa
table(is.na(se$disease))
table(se$disease)
t(table(se$disease, se$study_condition))

#-> non è una variabile esplicativa clinica
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "disease")]

# age
table(is.na(se$age))
summary(se$age)
hist(se$age)

boxplot(se$age ~ se$study_condition)

# age_category
table(is.na(se$age_category))
table(se$age_category)
str(se$age_category)

#-> in relazione deterministica con age
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "age_category")]

# gender
table(is.na(se$gender))
prop.table(table(se$gender))
str(se$gender)
se_fix$gender <- as.factor(se_fix$gender)

# country
table(is.na(se$country))
prop.table(table(se$country))
str(se$country)
se_fix$country <- as.factor(se_fix$country)

# non_westernized
table(is.na(se$non_westernized))
prop.table(table(se$non_westernized)) # inultile

#-> non utile per prevedere la risposta
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "non_westernized")]

# sequencing_platform
table(is.na(se$sequencing_platform))
prop.table(table(se$sequencing_platform)) # inultile

#-> non utile per prevedere la risposta
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "sequencing_platform")]

# DNA_extraction_kit (intanto la toglierei)
table(is.na(se$DNA_extraction_kit)) # metà NA
prop.table(table(se$DNA_extraction_kit)) # inultile
table(se$DNA_extraction_kit, se$study_name, useNA = "always") #solo uno studio non rileva il tipo di extraction
#si può pensare di tenerla solo se si stratifica per studio e si considerano gli microbiomi

#-> per ora si toglie
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "DNA_extraction_kit")]

# PMID (inutile, articolo PUBMED)
table(is.na(se$PMID))
prop.table(table(se$PMID)) # da togliere

#-> non utile per spiegare la risposta
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "PMID")]

# number_reads
table(is.na(se$number_reads))
summary(se$number_reads)
hist(se$number_reads, nclass = 50)


# number_bases (numero di basi di tutte le reads)
table(is.na(se$number_bases)) 
summary(se$number_bases)
hist(se$number_bases, nclass = 50)

cor(se$number_reads, se$number_bases)

#Si toglie anche dall'altro dataset perché è fortemente correlato con number_bases
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "number_bases")]


# minimum_read_length
table(is.na(se$minimum_read_length)) 
summary(se$minimum_read_length)
table(se$minimum_read_length)

# median_read_length
table(is.na(se$median_read_length)) 
summary(se$median_read_length)
table(se$median_read_length)

cor(se$median_read_length, se$minimum_read_length)
cor(se$median_read_length, se$number_reads)


# NCBI_accession -> codice univoco per una sequenza
table(is.na(se$NCBI_accession)) # troppi NA

#Si toglie dai dataset
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "NCBI_accession")]

# curator (da togliere)
table(is.na(se$curator)) 
table(se$curator)
levels(as.factor(se$curator))

#Si toglie dai dataset
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "curator")]

# BMI (occhio, 7 NA!)
table(is.na(se$BMI)) 
summary(se$BMI)
table(se$BMI)

#dato che sono solo 7 gli NA, togliamo le osservazioni anziché l'intera variabile
se_fix <- se_fix[,!is.na(se_fix$BMI)]

# location (da togliere, troppi NA)
table(is.na(se$location)) 
table(se$location)

#Si toglie dai dataset
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "location")]


# disease_subtype (ha tanti NA, capire perché parlavano di usarla per classificare)
table(is.na(se$disease_subtype))
table(se$disease_subtype, se$study_condition, useNA = "always")
table(se$disease_subtype, se$disease, useNA = "always")

table(se$disease_subtype, se$country, useNA = "always")
table(se$disease_subtype, se$study_name, useNA = "always")

#Sicuramente non è una variabile esplicativa, per ora la togliamo da entrambi i dataset:
#Valuteremo se usarla come possibile variabile risposta
#Si toglie dai dataset
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "disease_subtype")]

# alcohol
table(is.na(se$alcohol))
table(se$alcohol) # ha solo una modalità, non è una variabile
table(se$alcohol, se$study_condition)

#Si toglie dai dataset
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "alcohol")]

# triglycerides, da togliere
table(is.na(se$triglycerides))
summary(se$triglycerides)
table(is.na(se$triglycerides), is.na(se$hdl))
table(is.na(se$triglycerides), is.na(se$ldl))

#Si toglie dai dataset
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "triglycerides")]

# hdl, colesterolo buono, troppi NA
table(is.na(se$hdl))
table(is.na(se$hdl), se$study_name)

#Si toglie dai dataset
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "hdl")]

# ldl, colesterolo cattivo, troppi NA
table(is.na(se$ldl))
table(is.na(se$ldl), is.na(se$hdl))

#Si toglie dai dataset
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "ldl")]

# hba1c, emoglobina glicata
table(is.na(se$hba1c))
summary(se$hba1c)
table(is.na(se$hba1c), is.na(se$hdl))

#Si toglie dai dataset
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "hba1c")]

# smoker
table(is.na(se$smoker))
table(se$smoker)

#Si toglie dai dataset perché ha solo una modalità e tanti NA
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "smoker")]

# ever_smoker
table(is.na(se$ever_smoker))
table(se$ever_smoker)

#Si toglie dai dataset perché ha solo una modalità e tanti NA
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "ever_smoker")]

# fobt -> presenza di sangue occulto nelle feci
table(is.na(se$fobt))
table(se$fobt)

#Si toglie dai dataset
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "fobt")]


# brinkman_index, tant
table(is.na(se$brinkman_index))
table(se$brinkman_index)

table(is.na(se$brinkman_index), se$study_name, useNA = "always") #solo in YachidaS_2019

#Per ora si toglie da dataset
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "brinkman_index")]

# alcohol_numeric
table(is.na(se$alcohol_numeric))
summary(se$alcohol_numeric)

table(is.na(se$alcohol_numeric), se$study_name, useNA = "always") #solo in YachidaS_2019

#Per ora si toglie da dataset
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "alcohol_numeric")]

# Pacchetto table1 per le analisi esplorative (dove si becca? -> dal sui gitHub)

# PMID -> codice PUBMED, con il codice si ottiene il paper pubblicato sui dati

'----------------------------Esplorativa Tassonometria--------------------------------------------'

colnames(rowData(se))
tmp <- rowData(se)

#Kingdom
table(is.na(tmp$Kingdom))
table(tmp$Kingdom)

#Phylum
table(is.na(tmp$Phylum))
sort(table(tmp$Phylum), decreasing = T) #ci sono le prime quattro modalità alte, le altre si
#potrebbero accorpare
length(table(tmp$Phylum))
sum(sort(table(tmp$Phylum), decreasing = T)[1:4])/ sum(table(tmp$Phylum)) 

#Class
table(is.na(tmp$Class))
sort(table(tmp$Class), decreasing = T)
length(table(tmp$Class))
sum(sort(table(tmp$Class), decreasing = T)[1:6])/ sum(table(tmp$Class)) 
#Anche qui forse ha senso accorpare

#Order
table(is.na(tmp$Order))
length(table(tmp$Order))
sort(table(tmp$Order), decreasing = T)
sum(sort(table(tmp$Order), decreasing = T)[1:7])/ sum(table(tmp$Order)) 

#Family
table(is.na(tmp$Family))
length(table(tmp$Family))
sort(table(tmp$Family), decreasing = T)
sum(sort(table(tmp$Family), decreasing = T)[1:14])/ sum(table(tmp$Family)) 
#qui è più difficile accorpare

#Genus
table(is.na(tmp$Genus))
length(table(tmp$Genus))
sort(table(tmp$Genus), decreasing = T)
sum(sort(table(tmp$Genus), decreasing = T)[1:20])/ sum(table(tmp$Genus)) 

#Species
table(is.na(tmp$Species))
length(table(tmp$Species))
sort(table(tmp$Species), decreasing = T)
#tutte modalità diverse

'----------------------------------Esplorativa Assay-------------------------------------------------'
#assay: relative abundance
dim(assay(se))
#rownames(assay(se))

#togliamo gli zeri
#se <- se[rowMaxs(assay(se)) != 0,]
dim(se[rowMaxs(assay(se)) != 0,])
dim(se[rowMaxs(assay(se)) > 0.1,]) #Questo è il filtro che teniamo
se_fix <- se_fix[rowMaxs(assay(se)) > 0.1,]

cor_fix <- cor(t(assay(se_fix))) - diag(1,502)
row.names(cor_fix) <- NULL
which(cor_fix == 1, arr.ind = T)
which(cor_fix >= 0.9, arr.ind = T)

col <- colorRampPalette(c("blue", "white", "red"))(200)
heatmap(cor_fix, col = col, symm = T)
heatmap(cor_fix, col = col, symm = T, Colv = NA, Rowv = NA)

#Creazione data.frame
dim(colData(se_fix))
dim(assay(se_fix))
data_fix <- as.data.frame(cbind(colData(se_fix), t(assay(se_fix))))

save(data_fix, file = "dataset_fix_fixcl.RData")

# 1: $ study_condition -> risposta
# 2-4 8 : cliniche
#  5-7: collegate a microbioma
# 9-510: proporzioni microbioma