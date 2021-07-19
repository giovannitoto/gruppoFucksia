# gruppoFucksia - Diario di bordo

## Giorno 1 - 19/07/2021

### Analisi preliminari del dataset
Nel file adenoma.rda abbiamo un oggetto ```se``` contenente una serie di dataset che fanno riferimento a 875 specie di batteri e 634 soggetti.

```colData(se)``` contiene informazioni relative ai 634 soggetti; per ognuno di essi si hanno 33 variabili:
1. ```study_condition```: 5 studi
2. ```subject_id```: 634 individui distinti, quindi si ha una osservazione per ogni soggetto.
3. ```body_site``` ha un'unica modalità che è stool.
4. ```antibiotics_current_use```  ha molti valori nulli (529 su 634), tutti i rimanenti sono 'no', quindi buttiamo la variabile.
5. ```study_condition```: 209 malati (33%) e 429 controlli (67%)
6. ```disease``` contiene un eventuali malattie aggiuntive del soggetto: da notare che anche i controlli possono avere delle malattie; la modalità adenoma fa riferimento ai malati che non hanno ulteriori malattie; sarebbe molto utile fare la tabella di frequenze congiunte delle variabili ```study_condition``` e ```disease```.
7. ```age``` va da 21 a 89 anni, concentrata tra i 60 e 70, la mediana è 64. Facendo dei boxplot condizionatamente a study_condition, si osserva la stessa mediana ma maggior variabilità nei controlli: ha senso poiché sono il doppio.
8. ```age_category``` è una versione aggregata di ```age```.
9. ```gender``` presenta due generi: 278 femmine (44%) e 356 maschi (56%).
10. ```country``` ha 6 modalità: CAN (25), USA (29), ITA (51), FRA (103), AUT (108), JPN (318).
11. ```non_westernized``` ha un'unica modalità, quindi la droppiamo.
12. ```sequencing_platform``` ha un'unica modalità, quindi la droppiamo.
13. ```DNA_extraction_kit``` indica come è stato estratto il DNA: per circa metà dei soggetti non è riportato il metodo (valore mancante), ergo buttiamo la variabile.
14. ```PMID``` è un identificatore del dataset PUBMED, quindi lo droppiamo; è utile solo per ricavare la pubblicazione di riferimento.
15. ```number_reads``` è il numeri di frammenti di DNA disponibili per il soggetto.
16. ```number_bases``` è il numero di basi per il soggetto; si noti che un read è formato da una sequenza di basi.
17. ```minimum_read_length``` è un conteggio che indica la lunghezza minima osservata dei segmenti del soggetto; in realtà tutte le reads hanno la stessa lunghezza, tuttavia il macchinario che le sequenzia le spezza prima se l'errore di codifica che si forma "leggendo" le basi risulta troppo alto. 
18. ```median_read_length``` è la mediana della lunghezza osservata dei segmenti; si nota una correlazione bassa (0.32) con la lunghezza minima ```minimum_read_length```.
19. ```NCBI_accession``` è piena di valori nulli, ergo la buttiamo.
20. ```curator``` ha due modalità; indica i curatori della pubblicazione? Inutile, quindi da buttare.
21. ```BMI``` (Body Mass Index) ha 7 valori nulli: da capire come gestire la cosa.
22. ```location``` ha quasi esclusivamente valori nulli, ergo la droppiamo.
23. ```disease_subtype``` è una variabile risposta alternativa che permette di effettuare analisi più approfondite passando da classificazione a multiclassificazione.
24. ```alcohol``` indica se abbiamo a che fare con un alcolizzato; pieno di valori nulli.
25. ```triglycerides``` è una variabile numerica piena di valori nulli (528), ergo la buttiamo.
26. ```hdl``` è piena di valori nulli (529), ergo la buttiamo.
27. ```ldl``` è piena di valori nulli (529), ergo la buttiamo; verosimilmente coincidono con quelle della variabile ```hdl```.
28. ```hba1c``` è una variabile numerica che fa riferimento all'emoglobina glicata, ovvero l'emoglobina col glucosio; è piena di valori nulli (562), ergo le buttiamo.
29. ```smoker``` indica se abbiamo a che fare con un fumatore; pieno di valori nulli.
30. ```ever_smoker``` indica se il soggetto ha mai fumato; piena di valori nulli.
31. ```fobt``` indica la presenza di sangue occulto nelle feci; piena di valori nulli (478); unico aspetto positivo è che abbiamo sia dei SI sia dei NO.
32. ```brinkman_index``` è una misura di esposizione al fumo di sigaretta; si hanno 317 valori nulli, i valori rimanenti sono interi tra 0 e 3200.
33. ```alcohol_numeric``` è indicatore di quanto si imbriaga di solito; si hanno 317 valori nulli; si  noti che questa variabile e ```brinkman_index``` sono state raccolte nello studio Yachida_2019.


```rowData(se)``` contiene la tassonomia delle 875 specie di batteri osservate; per ognuna di esse si hanno 7 variabili:
1. ```Kingdom``` ha tre modalità: Archaea (7), Bacteria (863) e Eukaryota (5); fattore estremamente sbilanciato.
2. ```Phylum``` ha 17 modalità di cui le prime 4 rappresetano il 94% del dataset (prime nel senso con maggior frequenza); **ha senso accorpare le rimanenti in "altro"?**
3. ```Class``` ha 32 modalità di cui le prime 3 rappresentano il 55% del dataset, le prime 5 il 76%, le prime 6 l'81%; **ha senso accorpare le rimanenti in "altro"? Come?**
4. ```Order``` ha 50 modalità di cui le prime 4 rappresetano il 60% del dataset, le prima 7 il 69%; **ha senso accorpare le rimanenti in "altro"?**
5. ```Family``` ha 97 modalità con frequenza distribuita abbastanza equamente; non è così intuitivo fare accorpamenti vari.
6. ```Genus``` ha 256 modalità.
7. ```Species``` ha 875 modalità, una per ogni specie di batterio.


### Domande
- ```disease_subtype``` ha un sacco di valori nulli: corrispondono a mancate osservazioni oppure ad altro?
- ```alcohol``` ha 48 NO e rimanenti NA, gli NA corrispondono ad astenuti oppure SI?
- ```smooker``` ha 48 NO e rimanenti NA, gli NA corrispondono ad astenuti oppure SI?
- ```ever_smoker``` ha 48 NO e rimanenti NA, gli NA corrispondono ad astenuti oppure SI?


### Chicche
- Potrebbe essere interessante provare a ricostruire la tassonomia dei batteri a partire dalle altre variabili.


## Giorno 2 - 20/07/2021


## Giorno 3 - 21/07/2021


## Giorno 4 - 22/07/2021


## Giorno 5 - 23/07/2021






