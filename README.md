# TREE-BASED-CLASSIFICATION-METHODS
• Thesis on tree-based machine learning classification algorithms such as random forest, boosting and bagging.

• Application of these algorithms to classify tumors as benign or malignant with 98% accuracy through adaboost.

• Application of these algorithms for comparison with parametric methods such as logistic regression, naive Bayes and linear discriminant analysis in which the random forest had a 
  higher accuracy of 6%.

## Analisi Breast Cancer rivista (`breast_cancer_analysis.R`)

Versione migliorata e ampliata della prima parte pratica (il codice originale della tesi resta in `code`).

**Miglioramenti rispetto al codice originale**
- percorso del file configurabile e colonne inutili (`id`, colonna vuota `X`) rimosse per nome invece che per posizione;
- divisione training/test stratificata 70/30 (`createDataPartition`) e classe positiva = **maligno**, così sensibilità e specificità hanno il significato clinico corretto;
- tutti gli iperparametri scelti in 10-fold cross-validation ripetuta 3 volte sugli stessi fold: `cp` dell'albero (regola 1-SE), `mtry` della random forest, numero di alberi e profondità del boosting (prima fissati a 5000 e 4);
- niente più sovrascrittura di `data_cancer$diagnosis` in 0/1 per il boosting;
- codice duplicato (grafici dell'errore OOB, valutazione dei modelli) sostituito da funzioni;
- confronto finale con accuratezza + intervallo di confidenza 95%, sensibilità, specificità, AUC e curve ROC, al posto del boxplot di tre numeri.

**Analisi nuove**
- analisi esplorativa: matrice di correlazione e PCA;
- confronto robusto dei modelli sui 30 fold della CV con test appaiati (un solo test set di ~170 casi non basta a distinguere i metodi);
- confronto dell'importanza delle variabili tra i quattro modelli;
- interpretazione del boosting: grafici di dipendenza parziale e dipendenza congiunta delle due variabili principali;
- scelta della soglia di decisione in base al costo di un tumore maligno non riconosciuto (falso negativo);
- recursive feature elimination: quante variabili servono davvero per classificare bene.

Pacchetti R: `caret`, `rpart`, `randomForest`, `gbm`, `pROC`, `ggplot2`. I grafici vengono salvati in `figure_breast_cancer/`. Tempo di esecuzione: circa 5 minuti.
