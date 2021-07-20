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
se_fix2 <- se

'Guardo le covariate una per una'

# study_name
table(is.na(se$study_name))
table(se$study_name) # abbiamo 5 studi

# subject_id
table(is.na(se$subject_id))
length(unique(se$subject_id)) # tutti i soggetti osservati una sola volta

# body_site
table(is.na(se$body_site))
table(se$body_site) # sono tutti campioni di feci (non ha modalità)

# antibiotics_current_use
table(is.na(se$antibiotics_current_use)) # ha tanti NA, da droppare -> ha solo una modalità, non serve mai
table(se$antibiotics_current_use)

# study_condition -> variabile risposta
table(is.na(se$study_condition))
table(se$study_condition)
prop.table(table(se$study_condition))
str(se$study_condition) #deve diventare un factor

# disease -> altra possibile variabile risposta, assolutamente non una esplicativa
table(is.na(se$disease))
table(se$disease)
t(table(se$disease, se$study_condition))

# age
table(is.na(se$age))
summary(se$age)
hist(se$age)

boxplot(se$age ~ se$study_condition)

# age_category
table(is.na(se$age_category))
table(se$age_category)
str(se$age_category) #deve diventare un fattore (se tenuta), ma è in relazione deterministica con age (si preferisce age)

# gender
table(is.na(se$gender))
prop.table(table(se$gender))
str(se$gender) #deve diventare un factor

# country
table(is.na(se$country))
prop.table(table(se$country))
str(se$country) #deve diventare un factor

# non_westernized
table(is.na(se$non_westernized))
prop.table(table(se$non_westernized)) # inultile, ha solo una modalità

# sequencing_platform
table(is.na(se$sequencing_platform))
prop.table(table(se$sequencing_platform)) # inultile, ha solo una modalità

# DNA_extraction_kit (intanto la toglierei)
table(is.na(se$DNA_extraction_kit)) # 318 NA e 316 osservazioni
prop.table(table(se$DNA_extraction_kit)) # inultile
table(se$DNA_extraction_kit, se$study_name, useNA = "always") #solo uno studio non rileva il tipo di extraction
#si può pensare di tenerla solo se si stratifica per studio e si considerano gli microbiomi

# PMID (inutile, articolo PUBMED)
table(is.na(se$PMID))
prop.table(table(se$PMID)) # da togliere

# number_reads
table(is.na(se$number_reads))
summary(se$number_reads)
hist(se$number_reads, nclass = 50)
hist(log(se$number_reads), nclass = 50)

# number_bases (numero di basi di tutte le reads)
table(is.na(se$number_bases))
summary(se$number_bases)
hist(se$number_bases, nclass = 50)
hist(log(se$number_bases), nclass = 50)

cor(se$number_reads, se$number_bases) #cor 0.91

# minimum_read_length
table(is.na(se$minimum_read_length)) 
summary(se$minimum_read_length)
hist(se$minimum_read_length) #solo poche modalità
table(se$minimum_read_length)

# median_read_length
table(is.na(se$median_read_length)) 
summary(se$median_read_length)
hist(se$median_read_length)
table(se$median_read_length)

cor(se$median_read_length, se$minimum_read_length)
cor(se$median_read_length, se$number_reads)

# NCBI_accession -> codice univoco per una sequenza
table(is.na(se$NCBI_accession)) # 369 NA e 265 osservazioni

# curator (da togliere)
table(is.na(se$curator)) 
table(se$curator) #completamente inutile

# BMI (occhio, 7 NA!)
table(is.na(se$BMI)) #7 NA e 627 osservazioni
summary(se$BMI)
table(se$BMI)

# location (da togliere, troppi NA)
table(is.na(se$location)) #580 NA e 54 osservazioni
table(se$location)

# disease_subtype (ha tanti NA, capire perché parlavano di usarla per classificare)
table(is.na(se$disease_subtype)) #518 NA e 116 osservazioni
table(se$disease_subtype, se$study_condition, useNA = "always")
table(se$disease_subtype, se$disease, useNA = "always")
table(se$disease_subtype, se$country, useNA = "always")
table(se$disease_subtype, se$study_name, useNA = "always")

#Sicuramente non è una variabile esplicativa, per ora la togliamo da entrambi i dataset:
#Valuteremo se usarla come possibile variabile risposta

# alcohol
table(is.na(se$alcohol))
table(se$alcohol) # ha solo una modalità, non è una variabile
table(se$alcohol, se$study_condition)

# triglycerides, da togliere
table(is.na(se$triglycerides)) #528 NA e 106 osservazioni
summary(se$triglycerides)
table(is.na(se$triglycerides), is.na(se$hdl))
table(is.na(se$triglycerides), is.na(se$ldl))

# hdl, colesterolo buono, troppi NA
table(is.na(se$hdl)) #529 NA e 105 osservazioni
table(is.na(se$hdl), se$study_name)

# ldl, colesterolo cattivo, troppi NA
table(is.na(se$ldl)) #529 NA e 105 osservazioni
table(is.na(se$ldl), is.na(se$hdl))

# hba1c, emoglobina glicata
table(is.na(se$hba1c)) #562 NA e 72 osservazioni
summary(se$hba1c)
table(is.na(se$hba1c), is.na(se$hdl))

# smoker
table(is.na(se$smoker))
table(se$smoker) #ha solo una modalità, non è una variabile

# ever_smoker
table(is.na(se$ever_smoker))
table(se$ever_smoker) #ha solo una modalità non è una variabile

# fobt -> presenza di sangue occulto nelle feci
table(is.na(se$fobt)) #478 NA e 156 osservazioni
table(se$fobt)

# brinkman_index, tant
table(is.na(se$brinkman_index)) #317 NA e 317 osservazioni
table(se$brinkman_index)
table(is.na(se$brinkman_index), se$study_name, useNA = "always") #solo in YachidaS_2019

# alcohol_numeric
table(is.na(se$alcohol_numeric)) #317 NA e 317 osservazioni
summary(se$alcohol_numeric)
table(is.na(se$alcohol_numeric), se$study_name, useNA = "always") #solo in YachidaS_2019

#Provo a tenere solo i 116 valori osservati di disease_subtype e vediamo quanti altri NA mi restano:
se_fix2 <- se[,!is.na(se$disease_subtype)]

#Guardo le variabili che avevano NA
table(is.na(se_fix2$DNA_extraction_kit)) #Nessun NA
table(is.na(se_fix2$NCBI_accession)) #27 NA
table(is.na(se_fix2$BMI)) #2 NA
table(is.na(se_fix2$location)) #116 NA
table(is.na(se_fix2$disease_subtype)) #Nessun NA
table(is.na(se_fix2$triglycerides)) #69 NA
table(is.na(se_fix2$hdl)) #69 NA
table(is.na(se_fix2$ldl)) #69 NA
table(is.na(se_fix2$hba1c)) #81 NA
table(is.na(se_fix2$fobt)) #74 NA
table(is.na(se_fix2$brinkman_index)) #116 NA
table(is.na(se_fix2$alcohol_numeric)) #116 NA

#Non va bene, prendo le osservazioni non NA a partire da  un'altra variabile
#Provo a tenere solo i 105 valori osservati di hdl e vediamo quanti altri NA mi restano:
se_fix2 <- se[,!is.na(se$hdl)]

#Guardo le variabili che avevano NA
table(is.na(se_fix2$DNA_extraction_kit)) #Nessun NA
table(is.na(se_fix2$NCBI_accession)) #Nessun NA
table(is.na(se_fix2$BMI)) #Nessun NA
table(is.na(se_fix2$location)) #105 NA
table(is.na(se_fix2$triglycerides)) #Nessun NA
table(is.na(se_fix2$hdl)) #Nessun NA
table(is.na(se_fix2$ldl)) #Nessun NA
table(is.na(se_fix2$hba1c)) #33 NA
table(is.na(se_fix2$fobt)) #105 NA
table(is.na(se_fix2$brinkman_index)) #105 NA
table(is.na(se_fix2$alcohol_numeric)) #105 NA

#Meglio, ma provo ancora:
#tengo solo le 156 osservazioni di fobt e vediamo:
se_fix2 <- se[,!is.na(se$fobt)]

#Guardo le variabili che avevano NA
table(is.na(se_fix2$DNA_extraction_kit)) #Nessun NA
table(is.na(se_fix2$NCBI_accession)) #Nessun NA
table(is.na(se_fix2$BMI)) #2 NA
table(is.na(se_fix2$location)) #102 NA
table(is.na(se_fix2$triglycerides)) #156 NA
table(is.na(se_fix2$hdl)) #156 NA
table(is.na(se_fix2$ldl)) #156 NA
table(is.na(se_fix2$hba1c)) #156 NA
table(is.na(se_fix2$fobt)) #Nessun NA
table(is.na(se_fix2$brinkman_index)) #156 NA
table(is.na(se_fix2$alcohol_numeric)) #156 NA

#Non bene, provo con il brinkman_index:
#Tengo le 317 osservazioni e vediamo:
se_fix2 <- se[,!is.na(se$brinkman_index)]

#Guardo le variabili che avevano NA
table(is.na(se_fix2$DNA_extraction_kit)) #317 NA
table(is.na(se_fix2$NCBI_accession)) #317 NA
table(is.na(se_fix2$BMI)) #Nessun NA
table(is.na(se_fix2$location)) #317 NA
table(is.na(se_fix2$triglycerides)) #317 NA
table(is.na(se_fix2$hdl)) #317 NA
table(is.na(se_fix2$ldl)) #317 NA
table(is.na(se_fix2$hba1c)) #317 NA
table(is.na(se_fix2$fobt)) #317 NA
table(is.na(se_fix2$brinkman_index)) #Nessun NA
table(is.na(se_fix2$alcohol_numeric)) #Nessun NA

#Conclusione facciamo due subset: uno con le 105 osservazioni a partire da hdl e uno con 317 con le righe di brinkman index
#Tengo il subset costruito con le 105 osservazioni di hdl e droppo le variabili inutili:
se_fix2 <- se[,!is.na(se$hdl)]

'---------------------------Tolgo le variabili una alla volta----------------------------------------'
# study_name -> non utile per prevedere la risposta
colData(se_fix2) <- colData(se_fix2)[,which(colnames(colData(se_fix2)) != "study_name")]

# subject_id -> non utile per prevedere la risposta
colData(se_fix2) <- colData(se_fix2)[,which(colnames(colData(se_fix2)) != "subject_id")]

# body_site -> non utile per prevedere la risposta (ha solo una modalità)
colData(se_fix2) <- colData(se_fix2)[,which(colnames(colData(se_fix2)) != "body_site")]

# antibiotics_current_use -> non utile per prevedere la risposta (ha solo una modalità)
colData(se_fix2) <- colData(se_fix2)[,which(colnames(colData(se_fix2)) != "antibiotics_current_use")]

# study_condition -> variabile risposta
table(se_fix2$study_condition)
str(se_fix2$study_condition)
se_fix2$study_condition <- as.factor(se_fix2$study_condition)

# disease -> possibile variabile risposta, assolutamente non una esplicativa
table(is.na(se_fix2$disease))
table(se_fix2$disease)
str(se_fix2$disease)
se_fix2$disease <- as.factor(se_fix2$disease)

# age -> è ok così

# age_category -> in relazione deterministica con age va tolta
colData(se_fix2) <- colData(se_fix2)[,which(colnames(colData(se_fix2)) != "age_category")]

# gender
str(se_fix2$gender)
se_fix2$gender <- as.factor(se_fix2$gender)
levels(se_fix2$gender)

# country
prop.table(table(se_fix2$country)) #ora ha solo una modalità, va tolta
colData(se_fix2) <- colData(se_fix2)[,which(colnames(colData(se_fix2)) != "country")]

# non_westernized -> non utile per prevedere la risposta (ha solo una modalità)
colData(se_fix2) <- colData(se_fix2)[,which(colnames(colData(se_fix2)) != "non_westernized")]

# sequencing_platform -> non utile per prevedere la risposta (ha solo una modalità)
colData(se_fix2) <- colData(se_fix2)[,which(colnames(colData(se_fix2)) != "sequencing_platform")]


# DNA_extraction_kit (intanto la toglierei)
table(is.na(se_fix2$DNA_extraction_kit))
prop.table(table(se$DNA_extraction_kit)) # inultile, va tolta
colData(se_fix2) <- colData(se_fix2)[,which(colnames(colData(se_fix2)) != "DNA_extraction_kit")]


# PMID (inutile, articolo PUBMED)-> non utile per spiegare la risposta
colData(se_fix2) <- colData(se_fix2)[,which(colnames(colData(se_fix2)) != "PMID")]

# number_reads -> ok così

# number_bases (numero di basi di tutte le reads)
cor(se_fix2$number_reads, se_fix2$number_bases)
#Troppo correlato, si toglie bases
colData(se_fix2) <- colData(se_fix2)[,which(colnames(colData(se_fix2)) != "number_bases")]


# minimum_read_length -> ok così

# median_read_length -> ok così

# NCBI_accession -> codice univoco per una sequenza (non informativo per la previsione)
colData(se_fix2) <- colData(se_fix2)[,which(colnames(colData(se_fix2)) != "NCBI_accession")]

# curator (da togliere) -> non utile per la previsione
colData(se_fix2) <- colData(se_fix2)[,which(colnames(colData(se_fix2)) != "curator")]

# BMI -> ok così
table(is.na(se_fix2$BMI)) #nessun NA

# location (da togliere, troppi NA)
table(is.na(se_fix2$location)) #Tutta NA
colData(se_fix2) <- colData(se_fix2)[,which(colnames(colData(se_fix2)) != "location")]


# disease_subtype (ha tanti NA, capire perché parlavano di usarla per classificare)
table(is.na(se_fix2$disease_subtype)) #Ha ancora NA
table(se_fix2$disease_subtype) #Ha solo una modalità, non è più una variabile
colData(se_fix2) <- colData(se_fix2)[,which(colnames(colData(se_fix2)) != "disease_subtype")]

# alcohol
table(is.na(se_fix2$alcohol)) #tutta NA -> si toglie
colData(se_fix2) <- colData(se_fix2)[,which(colnames(colData(se_fix2)) != "alcohol")]

# triglycerides, da togliere
table(is.na(se_fix2$triglycerides))
summary(se_fix2$triglycerides) #ora è ok
hist(se_fix2$triglycerides, nclass = 50)

# hdl, colesterolo buono
table(is.na(se_fix2$hdl)) #ora è ok
hist(se_fix2$hdl, nclass = 50)

# ldl, colesterolo cattivo,
table(is.na(se_fix2$ldl)) #ora è ok
hist(se_fix2$ldl, nclass = 50)

# hba1c, emoglobina glicata
table(is.na(se_fix2$hba1c)) # 33 NA, tocca perdere troppe osservazioni: tolgo la variabile
colData(se_fix2) <- colData(se_fix2)[,which(colnames(colData(se_fix2)) != "hba1c")]
'--------------------SONO ARRIVATO QUI 20/07/21--------------------'
# smoker
table(is.na(se_fix2$smoker))
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

'-----------------Fine--------------------'

#Tengo il subset costruito con le 317 osservazioni di brinkman_index e droppo le variabili inutili:
se_fix3 <- se[,!is.na(se$brinkman_index)]

'----------------------------------Esplorativa Assay-------------------------------------------------'
#assay: relative abundance
dim(assay(se))
#rownames(assay(se))

#togliamo gli zeri
#se <- se[rowMaxs(assay(se)) != 0,]
dim(se[rowMaxs(assay(se)) != 0,])
dim(se[rowMaxs(assay(se)) > 0.1,]) #Questo è il filtro che teniamo
se_fix <- se_fix[rowMaxs(assay(se)) > 0.1,]
se_fix_clinical <- se_fix_clinical[rowMaxs(assay(se)) > 0.1,]

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
data_fix_cl <- as.data.frame(cbind(colData(se_fix_clinical), t(assay(se_fix_clinical))))

save(data_fix, data_fix_cl, file = "dataset_fix_fixcl.RData")
