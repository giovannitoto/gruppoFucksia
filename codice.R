load("adenoma.rda")

head(rowData(se))

# Oggetto generale
se
dim(se)
# le righe sono i soggetti, le colonne le specie di batteri (vero?)

# Guardo colData, ovvero le caratteristiche degli individui
dim(colData(se))
colData(se)

# study_name
table(is.na(se$study_name))
table(se$study_name) # abbiamo 5 studi

# subject_id
table(is.na(se$subject_id))
length(unique(se$subject_id)) # tutti i soggetti osservati una sola volta

# body_site
table(is.na(se$body_site))
table(se$body_site) # sono tutti campioni di feci

# antibiotics_current_use
table(is.na(se$antibiotics_current_use)) # ha tanti NA, da droppare
table(se$antibiotics_current_use)

# study_condition
table(is.na(se$study_condition))
table(se$study_condition)
prop.table(table(se$study_condition))

# disease
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

# gender
table(is.na(se$gender))
prop.table(table(se$gender))

# country
table(is.na(se$country))
prop.table(table(se$country))

# non_westernized
table(is.na(se$non_westernized))
prop.table(table(se$non_westernized)) # inultile

# sequencing_platform
table(is.na(se$sequencing_platform))
prop.table(table(se$sequencing_platform)) # inultile

# DNA_extraction_kit (da buttare)
table(is.na(se$DNA_extraction_kit)) # metà NA
prop.table(table(se$DNA_extraction_kit)) # inultile

# PMID (inutile, articolo PUBMED)
table(is.na(se$PMID))
prop.table(table(se$PMID)) # da togliere

# number_reads
table(is.na(se$number_reads))
summary(se$number_reads)
hist(se$number_reads, nclass = 50)

# number_bases (numero di basi di tutte le reads)
table(is.na(se$number_bases)) 
summary(se$number_bases)
hist(se$number_bases, nclass = 50)

cor(se$number_reads, se$number_bases)

# minimum_read_length
table(is.na(se$minimum_read_length)) 
summary(se$minimum_read_length)
table(se$minimum_read_length)

# median_read_length
table(is.na(se$median_read_length)) 
summary(se$median_read_length)
table(se$median_read_length)

cor(se$median_read_length, se$minimum_read_length)

# NCBI_accession
table(is.na(se$NCBI_accession)) # troppi NA

# curator (da togliere)
table(is.na(se$curator)) 
table(se$curator)
levels(as.factor(se$curator))

# BMI (occhio, 7 NA!)
table(is.na(se$BMI)) 
summary(se$BMI)
table(se$BMI)

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