# gruppoFucksia - Diario di bordo

## Obiettivi
- Stimare un classificazione sulla variabile ```study_condition``` utilizzando solo la proporzione delle varie specie di batteri nei soggetti.
- Stimare un classificazione sulla variabile ```study_condition``` utilizzando sia la proporzione delle varie specie di batteri nei soggetti sia le variabili relative a ogni soggetto e/o ogni specie di batterio.
- Applicare *group lasso* considerando gruppi di specie di batteri, in accordo con la tassonomia.
- Potrebbe essere interessante provare a ricostruire la tassonomia dei batteri a partire dalle altre variabili.

## Giorno 1 - 19/07/2021
Nel file adenoma.rda abbiamo un oggetto ```se``` contenente una serie di dataset che fanno riferimento a 875 specie di batteri e 634 soggetti.

### Analisi preliminari di ```colData(se)```
```colData(se)``` contiene informazioni relative ai 634 soggetti; per ognuno di essi si hanno 33 variabili:
1. ```study_name```: 5 studi
2. ```subject_id```: 634 individui distinti, quindi si ha una osservazione per ogni soggetto.
3. ```body_site``` ha un'unica modalità che è stool.
4. ```antibiotics_current_use``` ha molti valori nulli (529 su 634), tutti i rimanenti sono 'no', quindi buttiamo la variabile.
5. ```study_condition```: 209 malati (33%) e 429 controlli (67%).
6. ```disease``` contiene un eventuali malattie aggiuntive del soggetto: da notare che anche i controlli possono avere delle malattie; la modalità adenoma fa riferimento ai malati che non hanno ulteriori malattie; sarebbe molto utile fare la tabella di frequenze congiunte delle variabili ```study_condition``` e ```disease```.
7. ```age``` va da 21 a 89 anni, concentrata tra i 60 e 70, la mediana è 64. Facendo dei boxplot condizionatamente a ```study_condition```, si osserva la stessa mediana ma maggior variabilità nei controlli: ha senso poiché sono il doppio.
8. ```age_category``` è una versione aggregata di ```age```.
9. ```gender``` presenta due generi: 278 femmine (44%) e 356 maschi (56%).
10. ```country``` ha 6 modalità: CAN (25), USA (29), ITA (51), FRA (103), AUT (108), JPN (318).
11. ```non_westernized``` ha un'unica modalità, quindi la droppiamo.
12. ```sequencing_platform``` ha un'unica modalità, quindi la droppiamo.
13. ```DNA_extraction_kit``` indica come è stato estratto il DNA: per circa metà dei soggetti non è riportato il metodo (tutte e solo i soggetti dello studio Yachida_2019), ergo buttiamo la variabile.
14. ```PMID``` è un identificatore del dataset PUBMED, quindi lo droppiamo; è utile solo per ricavare la pubblicazione di riferimento.
15. ```number_reads``` è il numeri di frammenti di DNA disponibili per il soggetto.
16. ```number_bases``` è il numero di basi per il soggetto; si noti che un read è formato da una sequenza di basi.
17. ```minimum_read_length``` è un conteggio che indica la lunghezza minima osservata dei segmenti del soggetto; in realtà tutte le reads hanno la stessa lunghezza, tuttavia il macchinario che le sequenzia le spezza prima se l'errore di codifica che si forma "leggendo" le basi risulta troppo alto. 
18. ```median_read_length``` è la mediana della lunghezza osservata dei segmenti; si nota una correlazione bassa (0.32) con la lunghezza minima ```minimum_read_length```.
19. ```NCBI_accession``` è piena di valori nulli, ergo la buttiamo.
20. ```curator``` ha due modalità; indica i curatori della pubblicazione? Inutile, quindi da buttare.
21. ```BMI``` (Body Mass Index) ha 7 valori nulli: da capire come gestire la cosa.
22. ```location``` ha quasi esclusivamente valori nulli, ergo la droppiamo.
23. ```disease_subtype``` è una variabile risposta alternativa che permette di effettuare analisi più approfondite passando da classificazione a multiclassificazione; non può essere usata come variabile esplicativa poiché contiene risultati clinici non immediatti da ottenere.
24. ```alcohol``` indica se abbiamo a che fare con un alcolizzato; pieno di valori nulli.
25. ```triglycerides``` è una variabile numerica piena di valori nulli (528), ergo la buttiamo.
26. ```hdl``` è piena di valori nulli (529), ergo la buttiamo.
27. ```ldl``` è piena di valori nulli (529), ergo la buttiamo; verosimilmente coincidono con quelle della variabile ```hdl```.
28. ```hba1c``` è una variabile numerica che fa riferimento all'emoglobina glicata, ovvero l'emoglobina col glucosio; è piena di valori nulli (562), ergo le buttiamo.
29. ```smoker``` indica se abbiamo a che fare con un fumatore; pieno di valori nulli.
30. ```ever_smoker``` indica se il soggetto ha mai fumato; piena di valori nulli.
31. ```fobt``` indica la presenza di sangue occulto nelle feci; piena di valori nulli (478); unico aspetto positivo è che abbiamo sia dei SI sia dei NO.
32. ```brinkman_index``` è una misura di esposizione al fumo di sigaretta; si hanno 317 valori nulli, i valori rimanenti sono interi tra 0 e 3200; è stata rilevata solo nello studio YachidaS_2019.
33. ```alcohol_numeric``` è indicatore di quanto si imbriaga di solito; si hanno 317 valori nulli; si  noti che questa variabile e ```brinkman_index``` sono state raccolte nello studio Yachida_2019.

### Gestione delle variabili in ```colData(se)```
Costruiamo due dataset: il primo, detto *se_fix_clinical*, contiene solo variabili cliniche e il secondo, *se_fix*, contiene tutte le variabili disponibili e utili.
1. ```NO NO``` ```study_name```: inutile.
2. ```NO NO``` ```subject_id```: tutti distinti.
3. ```NO NO``` ```body_site```: unica modalità.
4. ```NO NO``` ```antibiotics_current_use```: troppi NA.
5. ```SI SI``` ```study_condition```: variabile risposta.
6. ```NO NO``` ```disease```: non utile per predere ```study_condition```; potrebbe essere utile per costruire una var. risposta muldimodale.
7. ```SI SI``` ```age```: variabile numerica.
8. ```NO NO``` ```age_category```: informazione già contenuta in ```age```.
9. ```SI SI``` ```gender```: fattore con 2 modalità.
10. ```SI SI``` ```country```: fattore con 6 modalità.
11. ```NO NO``` ```non_westernized```: unica modalità.
12. ```NO NO``` ```sequencing_platform```: unica modalità.
13. ```NO NO``` ```DNA_extraction_kit```: troppi NA **(inutile?)**.
14. ```NO NO``` ```PMID```: inutile.
15. ```NO SI``` ```number_reads```: variabile quantitativa.
16. ```NO NO``` ```number_bases```: correlazione superiore a 0.9 con ```number_basis```.
17. ```NO SI``` ```minimum_read_length```: variabile quantitativa. 
18. ```NO SI``` ```median_read_length```: variabile quantitativa.
19. ```NO NO``` ```NCBI_accession```: troppi NA.
20. ```NO NO``` ```curator```: inutile.
21. ```SI SI``` ```BMI```: variabile quantitativa. **7 valori mancanti => rimuoviamo soggetti corrispondenti**
22. ```NO NO``` ```location```: troppi NA.
23. ```NO NO``` ```disease_subtype```: non rilevante.
24. ```NO NO``` ```alcohol```: troppi NA e unica modalità.
25. ```NO NO``` ```triglycerides```: troppi NA. **da valutare in un dataset ridotto**
26. ```NO NO``` ```hdl```: troppi NA. **da valutare in un dataset ridotto**
27. ```NO NO``` ```ldl```: troppi NA. **da valutare in un dataset ridotto**
28. ```NO NO``` ```hba1c```: troppi NA. **da valutare in un dataset ridotto**
29. ```NO NO``` ```smoker```: troppi NA e unica modalità.
30. ```NO NO``` ```ever_smoker```: troppi NA e unica modalità.
31. ```NO NO``` ```fobt```: troppi NA. **da valutare in un dataset ridotto**
32. ```NO NO``` ```brinkman_index```: troppi NA. **da valutare in un dataset ridotto**
33. ```NO NO``` ```alcohol_numeric```: troppi NA. **da valutare in un dataset ridotto**


### Domande
- ```disease_subtype``` ha un sacco di valori nulli: corrispondono a mancate osservazioni oppure ad altro?
- **```alcohol```, ```smooker``` e ```ever_smoker``` hanno 48 NO e rimanenti NA, gli NA corrispondono ad astenuti oppure SI?** Non ha senso farsi tanti problemi: buttiamo tutto via.
- La variabile ```DNA_extraction_kit```, che indica il metodo di estrazione del DNA, è utile? Ha senso utilizzarla solo se si considerano anche le proporzioni delle specie di batteri nei soggetti, giusto?
- Vogliamo utilizzare solo informazioni relative ad esami semplici da effettuare oppure vogliamo considerare tutta l'informazione disponibile? Più nello specifico, voglio utilizzare variabili ottenute in seguito alla sequenziazione del DNA.
- Meglio eliminare soggetti a causa di pochi valori mancanti oppure fare imputazione in qualche modo?
 

## Giorno 2 - 20/07/2021
Concludiamo le analisi esplorative e gestione delle variabili in ```colData(se)```: i due dataset includono rispettivamente 6 (*se_fix_clinical*) e 9 (*se_fix*) variabili. Passiamo quindi ad un'analisi delle proporzioni delle specie di batteri nei soggetti, ovvero ci concentriamo ora sulla trasposta di ```rowData(se)```.

### Analisi preliminari di ```rowData(se)```
```rowData(se)``` contiene la tassonomia delle 875 specie di batteri osservate; per ognuna di esse si hanno 7 variabili:
1. ```Kingdom``` ha tre modalità: Archaea (7), Bacteria (863) e Eukaryota (5); fattore estremamente sbilanciato.
2. ```Phylum``` ha 17 modalità di cui le prime 4 rappresetano il 94% del dataset (prime nel senso con maggior frequenza); **ha senso accorpare le rimanenti in "altro"?**
3. ```Class``` ha 32 modalità di cui le prime 3 rappresentano il 55% del dataset, le prime 5 il 76%, le prime 6 l'81%; **ha senso accorpare le rimanenti in "altro"? Come?**
4. ```Order``` ha 50 modalità di cui le prime 4 rappresetano il 60% del dataset, le prima 7 il 69%; **ha senso accorpare le rimanenti in "altro"?**
5. ```Family``` ha 97 modalità con frequenza distribuita abbastanza equamente; non è così intuitivo fare accorpamenti vari.
6. ```Genus``` ha 256 modalità.
7. ```Species``` ha 875 modalità, una per ogni specie di batterio.

### Filtraggio e analisi preliminari in ```t(assay(se))```
```assay(se)``` è una matrice *specie x soggetto* il cui elemento *(i,j)* indica la proporzione della specie batterio *i* sul totale delle specie di batteri nel soggetto *j*; si ha quindi che ogni riga somma a 1. Come prima cosa, rimuoviamo le specie di batteri che non compaiono in nessun soggetto, passando così da 875 a 789 specie di batteri. A questo punto si rimuovono le specie che presentano una proporzione molto bassa in tutti i soggetti: per adesso si seleziona 0.1% (0.001) come soglia, passando da 789 a 502 specie.

Osserviamo che le specie 418, 419 e 420 hanno correlazione pari a 1 e quindi sono collineari; analogamente le specie 250 e 493 sono collineari.

# Creazione dei dataset
Generiamo due dataset *soggetto x variabile* come segue:
- ```data_fix_cl <- as.data.frame(cbind(colData(se_fix_clinical), t(assay(se_fix_clinical))))```  \[627x(5+502)\]
- ```data_fix <- as.data.frame(cbind(colData(se_fix), t(assay(se_fix))))```  \[627x(8+502)\]

e li salviamo nel file `dataset_fix_fixcl.Rdata`.

Si noti che a questo punto non sono ancora state rimosse le variabili collineari citate nella sezione precedente.

### Classificazione di study_condition a livello di specie
Prima di applicari i modelli, rimuoviamo le variabili collineari, ovvero le variabili in posizione 419+(5|8), 420+(5|8), 493+(5|8).
Visto che fornisce risultati migliori, decidiamo di fare *Convalida Incrociata* (*CV*).

Consideriamo 6 dataset:
1. data_fix: var. cliniche + var.sui batteri + proporzioni
2. data_fix[,-c(9:ncol(data_fix))]: var. cliniche + var. sui batteri
3. data_fix[,c(1,9:ncol(data_fix))]: proporzioni
4. data_fix_cl: var. cliniche + proporzioni
5. data_fix_cl[,-c(6:ncol(data_fix_cl))]: var. cliniche
6. data_fix_cl[,c(1,6:ncol(data_fix_cl))]: proporzioni


"data_fix", "data_fix_nobact", "data_fix_bact", "data_cl", "data_cl_nobact", "data_cl_bact"

cv_data <- list(data_fix, )


### Note
- Osserviamo che le variabili ```triglycerides```, ```hdl```, ```ldl``` hanno valori nulli sui stessi soggetti; inoltre sono disponibili solo nello studio FengQ_2015; potrebbe aver senso considerare solo i soggetti dello studio FengQ_2015?
- ```brinkman_index``` e ```alcohol_numeric``` è disponibile solo nello studio YachidaS_2019.



## Giorno 3 - 21/07/2021


## Giorno 4 - 22/07/2021


## Giorno 5 - 23/07/2021

## Da fare
- Leggere almeno gli abstract degli articoli contenuti in ```colData(se)$PMID```.
- Finire la gestione delle variabili in ```colData(se)```.
- Valutare in quali soggetti si ha un valore nullo della variabile antibiotics_current_use.




