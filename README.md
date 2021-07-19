# gruppoFucksia - Diario di bordo

## Giorno 1 - 19/07/2021 mattina

### Struttura del dataset
Nel file adenoma.rda abbiamo un oggetto ```se``` contenente una serie di dataset che fanno riferimento a 875 specie di batteri e 634 soggetti.

\[```colData(se)```\] Per ogni soggetto si hanno 32 variabili:
- 5 studi
- X individui distinti, quindi una osservazione per ogni soggetto
- ```antibiotics_current_use```  ha molti valori nulli (529 su 634), tutti i rimanenti sono nulli, quindi la variabile si butta.
- 209 malati (33%) e 429 controlli (67%)
- ```disease``` contiene un eventuali malattie aggiuntive del soggetto: da notare che anche i controlli possono avere delle malattie; la modalità adenoma fa riferimento ai malati che non hanno ulteriori malattie. E' molto utile fare la tabella di frequenze congiunte delle variabili ```study_condition``` e ```disease```.
- ```eta``` va da 21 a 89 anni, concentrata tra i 60 e 70, la mediana è 64. Facendo dei boxplot condizionatamente a study_condition, si osserva la stessa mediana ma maggior variabilità nei controlli: ha senso poiché sono il doppio.
- Due generi: 278 femmine (44%) e 356 maschi (56%)
- 6 nazioni: CAN (25), USA (29), ITA (51), FRA (103), AUT (108), JPN (318)
- ```DNA_extraction_kit``` indica come è stato estratto il DNA: per circa metà dei soggetti non è riportato il metodo (valore mancante), ergo buttiamo la variabile.
- PMID è un identificatore, quindi lo droppiamo.
- number_reads è il numeri di frammenti di DNA disponibili per il soggetto.

\[```rowData(se)```\] Per ogni specie di batterio si hanno:
-




