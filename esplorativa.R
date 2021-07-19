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
se_fix_clinical <- se
se_fix <- se
'Guardo le covariate una per una'

# study_name
table(is.na(se$study_name))
table(se$study_name) # abbiamo 5 studi

#-> non utile per prevedere la risposta
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "study_name")]
colData(se_fix_clinical) <- colData(se_fix_clinical)[,which(colnames(colData(se_fix_clinical)) != "study_name")]

# subject_id
table(is.na(se$subject_id))
length(unique(se$subject_id)) # tutti i soggetti osservati una sola volta

#-> non utile per prevedere la risposta
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "subject_id")]
colData(se_fix_clinical) <- colData(se_fix_clinical)[,which(colnames(colData(se_fix_clinical)) != "subject_id")]

# body_site
table(is.na(se$body_site))
table(se$body_site) # sono tutti campioni di feci

#-> non utile per prevedere la risposta
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "body_site")]
colData(se_fix_clinical) <- colData(se_fix_clinical)[,which(colnames(colData(se_fix_clinical)) != "body_site")]

# antibiotics_current_use
table(is.na(se$antibiotics_current_use)) # ha tanti NA, da droppare
table(se$antibiotics_current_use)

#-> non utile per prevedere la risposta
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "antibiotics_current_use")]
colData(se_fix_clinical) <- colData(se_fix_clinical)[,which(colnames(colData(se_fix_clinical)) != "antibiotics_current_use")]

# study_condition -> variabile risposta
table(is.na(se$study_condition))
table(se$study_condition)
prop.table(table(se$study_condition))
str(se$study_condition)
se_fix$study_condition <- as.factor(se_fix$study_condition)
se_fix_clinical$study_condition <- as.factor(se_fix_clinical$study_condition)

# disease -> altra possibile variabile risposta, assolutamente non una esplicativa
table(is.na(se$disease))
table(se$disease)
t(table(se$disease, se$study_condition))

#-> non è una variabile esplicativa clinica
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "disease")]
colData(se_fix_clinical) <- colData(se_fix_clinical)[,which(colnames(colData(se_fix_clinical)) != "disease")]

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
colData(se_fix_clinical) <- colData(se_fix_clinical)[,which(colnames(colData(se_fix_clinical)) != "age_category")]

# gender
table(is.na(se$gender))
prop.table(table(se$gender))
str(se$gender)
se_fix$gender <- as.factor(se_fix$gender)
se_fix_clinical$gender <- as.factor(se_fix_clinical$gender)

# country
table(is.na(se$country))
prop.table(table(se$country))
str(se$country)
se_fix$country <- as.factor(se_fix$country)
se_fix_clinical$country <- as.factor(se_fix_clinical$country)

# non_westernized
table(is.na(se$non_westernized))
prop.table(table(se$non_westernized)) # inultile

#-> non utile per prevedere la risposta
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "non_westernized")]
colData(se_fix_clinical) <- colData(se_fix_clinical)[,which(colnames(colData(se_fix_clinical)) != "non_westernized")]

# sequencing_platform
table(is.na(se$sequencing_platform))
prop.table(table(se$sequencing_platform)) # inultile

#-> non utile per prevedere la risposta
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "sequencing_platform")]
colData(se_fix_clinical) <- colData(se_fix_clinical)[,which(colnames(colData(se_fix_clinical)) != "sequencing_platform")]

# DNA_extraction_kit (intanto la toglierei)
table(is.na(se$DNA_extraction_kit)) # metà NA
prop.table(table(se$DNA_extraction_kit)) # inultile
table(se$DNA_extraction_kit, se$study_name, useNA = "always") #solo uno studio non rileva il tipo di extraction
#si può pensare di tenerla solo se si stratifica per studio e si considerano gli microbiomi

#-> per ora si toglie
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "DNA_extraction_kit")]
colData(se_fix_clinical) <- colData(se_fix_clinical)[,which(colnames(colData(se_fix_clinical)) != "DNA_extraction_kit")]

# PMID (inutile, articolo PUBMED)
table(is.na(se$PMID))
prop.table(table(se$PMID)) # da togliere

#-> non utile per spiegare la risposta
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "PMID")]
colData(se_fix_clinical) <- colData(se_fix_clinical)[,which(colnames(colData(se_fix_clinical)) != "PMID")]

# number_reads
table(is.na(se$number_reads))
summary(se$number_reads)
hist(se$number_reads, nclass = 50)

#Si toglie dal dateset clinico
colData(se_fix_clinical) <- colData(se_fix_clinical)[,which(colnames(colData(se_fix_clinical)) != "number_reads")]

# number_bases (numero di basi di tutte le reads)
table(is.na(se$number_bases)) 
summary(se$number_bases)
hist(se$number_bases, nclass = 50)

cor(se$number_reads, se$number_bases)

#Si toglie dal dateset clinico
colData(se_fix_clinical) <- colData(se_fix_clinical)[,which(colnames(colData(se_fix_clinical)) != "number_bases")]
#Si toglie anche dall'altro dataset perché è fortemente correlato con number_bases
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "number_reads")]


# minimum_read_length
table(is.na(se$minimum_read_length)) 
summary(se$minimum_read_length)
table(se$minimum_read_length)

#Si toglie dal dateset clinico
colData(se_fix_clinical) <- colData(se_fix_clinical)[,which(colnames(colData(se_fix_clinical)) != "minimum_read_length")]

# median_read_length
table(is.na(se$median_read_length)) 
summary(se$median_read_length)
table(se$median_read_length)

cor(se$median_read_length, se$minimum_read_length)
cor(se$median_read_length, se$number_reads)

#Si toglie dal dateset clinico
colData(se_fix_clinical) <- colData(se_fix_clinical)[,which(colnames(colData(se_fix_clinical)) != "median_read_length")]

# NCBI_accession -> codice univoco per una sequenza
table(is.na(se$NCBI_accession)) # troppi NA

#Si toglie dai dataset
colData(se_fix_clinical) <- colData(se_fix_clinical)[,which(colnames(colData(se_fix_clinical)) != "NCBI_accession")]
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "NCBI_accession")]

# curator (da togliere)
table(is.na(se$curator)) 
table(se$curator)
levels(as.factor(se$curator))

#Si toglie dai dataset
colData(se_fix_clinical) <- colData(se_fix_clinical)[,which(colnames(colData(se_fix_clinical)) != "curator")]
colData(se_fix) <- colData(se_fix)[,which(colnames(colData(se_fix)) != "curator")]

# BMI (occhio, 7 NA!)
table(is.na(se$BMI)) 
summary(se$BMI)
table(se$BMI)

#dato che sono solo 7 gli NA, togliamo le osservazioni anziché l'intera variabile
se_fix <- se_fix[,!is.na(se_fix$BMI)]
se_fix_clinical <- se_fix_clinical[,!is.na(se_fix_clinical$BMI)]

# location (da togliere, troppi NA)
table(is.na(se$location)) 
table(se$location)

# disease_subtype (ha tanti NA, capire perché parlavano di usarla per classificare)
table(is.na(se$disease_subtype))
table(se$disease_subtype, se$study_condition, useNA = "always")
table(se$disease_subtype, se$disease, useNA = "always")

table(se$disease_subtype, se$country, useNA = "always")

# alcohol
table(is.na(se$alcohol))
table(se$alcohol)
table(se$alcohol, se$study_condition)

# triglycerides, da togliere
table(is.na(se$triglycerides))
summary(se$triglycerides)

# hdl, colesterolo buono, troppi NA
table(is.na(se$hdl))

# ldl, colesterolo cattivo, troppi NA
table(is.na(se$ldl))

# hba1c, emoglobina glicata
table(is.na(se$hba1c))
summary(se$hba1c)

# smoker
table(is.na(se$smoker))
table(se$smoker)

# ever_smoker
table(is.na(se$ever_smoker))
table(se$ever_smoker)

# fobt
table(is.na(se$fobt))
table(se$fobt)

# brinkman_index, tant
table(is.na(se$brinkman_index))
table(se$brinkman_index)

table(se$brinkman_index, se$study_name, useNA = "always")

# alcohol_numeric
table(is.na(se$alcohol_numeric))
summary(se$alcohol_numeric)


# Pacchetto table1 per le analisi esplorative (dove si becca? -> dal sui gitHub)

# PMID -> codice PUBMED, con il codice si ottiene il paper pubblicato sui dati

'----------------------------------------------------------------------------------'
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

'---------------------------------------------------------------------------------'
#assay: relative abundance
dim(assay(se))
rownames(assay(se))

#togliamo gli zeri
se <- se[rowMaxs(assay(se)) != 0,]
dim(se[rowMaxs(assay(se)) > 0.1,])
rel_abu <- assay(se)
row_max <- rowMaxs(rel_abu)
head(row_max)
boxplot(row_max)
summary(row_max) #si può stabilire una soglia

