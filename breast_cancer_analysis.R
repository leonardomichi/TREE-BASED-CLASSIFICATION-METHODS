# =============================================================================
#  Metodi di classificazione ad albero
#  Parte pratica 1 - Breast Cancer Wisconsin (Diagnostic)
#
#  Versione rivista e ampliata dell'analisi della tesi:
#    1. analisi esplorativa (correlazioni, PCA)
#    2. albero di classificazione con potatura cost-complexity (regola 1-SE)
#    3. bagging, random forest e boosting con tuning in cross-validation
#    4. confronto sul test set (accuratezza, sensibilita', specificita', AUC, ROC)
#    5. confronto robusto con 10-fold CV ripetuta e test appaiati
#    6. NUOVO: interpretazione del boosting (dipendenza parziale e interazioni)
#    7. NUOVO: scelta della soglia di decisione in base al costo dei falsi negativi
#    8. NUOVO: quante variabili servono davvero? (recursive feature elimination)
#
#  Dataset: "Breast Cancer Wisconsin (Diagnostic)" (UCI / Kaggle), 569 tumori,
#  30 variabili numeriche, diagnosi M = maligno, B = benigno.
# =============================================================================


# ---- 0. Impostazioni ---------------------------------------------------------

percorso_dati   <- "data_cancer.csv"        # modificare con il percorso del file
salva_figure    <- TRUE                     # salva i grafici in formato .png
cartella_figure <- "figure_breast_cancer"

# costi usati nella scelta della soglia (sezione 7): un tumore maligno non
# riconosciuto (falso negativo) e' considerato 5 volte piu' grave di un falso
# allarme (falso positivo)
costo_FN <- 5
costo_FP <- 1

pacchetti <- c("caret", "rpart", "randomForest", "gbm", "pROC", "ggplot2")
mancanti  <- pacchetti[!vapply(pacchetti, requireNamespace, logical(1), quietly = TRUE)]
if (length(mancanti) > 0) install.packages(mancanti)
invisible(lapply(pacchetti, library, character.only = TRUE))

col_M       <- "#E69F00"   # maligno
col_B       <- "#56B4E9"   # benigno
col_blu     <- "#1E77E8"
col_arancio <- "#E56E45"
col_modelli <- c("Albero" = "#999999", "Bagging" = "#56B4E9",
                 "Random forest" = "#1E77E8", "Boosting" = "#E56E45")
theme_set(theme_minimal(base_size = 12))

if (salva_figure) dir.create(cartella_figure, showWarnings = FALSE)

# mostra un grafico ggplot e, se richiesto, lo salva su file
mostra <- function(grafico, nome, larghezza = 8, altezza = 5) {
  print(grafico)
  if (salva_figure) {
    ggsave(file.path(cartella_figure, paste0(nome, ".png")), grafico,
           width = larghezza, height = altezza, dpi = 150, bg = "white")
  }
  invisible(grafico)
}

# stessa cosa per i grafici di base R (es. il disegno dell'albero)
mostra_base <- function(disegna, nome, larghezza = 9, altezza = 6) {
  disegna()
  if (salva_figure) {
    png(file.path(cartella_figure, paste0(nome, ".png")),
        width = larghezza, height = altezza, units = "in", res = 150)
    disegna()
    dev.off()
  }
}


# ---- 1. Dati -----------------------------------------------------------------

carica_dati <- function(percorso) {
  dati <- read.csv(percorso, header = TRUE, sep = ",", dec = ".")
  # si eliminano l'identificativo e le colonne completamente vuote (nel file
  # Kaggle la virgola finale crea una colonna "X" di soli NA). Selezionare per
  # nome e' piu' sicuro che togliere le colonne 1 e 33 per posizione.
  dati <- dati[, names(dati) != "id" & colSums(!is.na(dati)) > 0]
  # "M" come primo livello: e' la classe positiva per caret e per le metriche
  dati$diagnosis <- factor(dati$diagnosis, levels = c("M", "B"))
  stopifnot(!anyNA(dati))
  dati
}

data_cancer <- carica_dati(percorso_dati)
variabili   <- setdiff(names(data_cancer), "diagnosis")
p           <- length(variabili)

dim(data_cancer)
table(data_cancer$diagnosis)
round(prop.table(table(data_cancer$diagnosis)), 3)


# ---- 2. Analisi esplorativa --------------------------------------------------

# 2a. Distribuzione delle 10 caratteristiche "medie" per diagnosi (standardizzate)
var_mean <- grep("_mean$", variabili, value = TRUE)
lungo <- data.frame(
  diagnosis = rep(data_cancer$diagnosis, times = length(var_mean)),
  variabile = factor(rep(sub("_mean$", "", var_mean), each = nrow(data_cancer)),
                     levels = sub("_mean$", "", var_mean)),
  valore    = as.vector(scale(data_cancer[, var_mean]))
)
g <- ggplot(lungo, aes(variabile, valore, fill = diagnosis)) +
  geom_boxplot(outlier.size = 0.6, alpha = 0.85) +
  scale_fill_manual(values = c(M = col_M, B = col_B),
                    labels = c(M = "Maligno", B = "Benigno")) +
  labs(x = NULL, y = "Valore standardizzato", fill = "Diagnosi",
       title = "Caratteristiche medie dei nuclei cellulari per diagnosi") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
mostra(g, "01_boxplot_variabili")

# 2b. Matrice di correlazione: molte variabili sono quasi collineari
#     (raggio, perimetro e area hanno correlazione ~0.99). Gli alberi non ne
#     soffrono, ma la collinearita' influenza le misure di importanza.
corr <- cor(data_cancer[, variabili])
corr_lungo <- data.frame(
  v1 = factor(rep(variabili, times = p), levels = variabili),
  v2 = factor(rep(variabili, each = p), levels = rev(variabili)),
  r  = as.vector(corr)
)
g <- ggplot(corr_lungo, aes(v1, v2, fill = r)) +
  geom_tile() +
  scale_fill_gradient2(low = col_blu, mid = "white", high = col_arancio,
                       limits = c(-1, 1)) +
  labs(x = NULL, y = NULL, fill = "r", title = "Matrice di correlazione") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        axis.text.y = element_text(size = 7))
mostra(g, "02_correlazioni", 8, 7.5)

coppie <- which(abs(corr) > 0.95 & upper.tri(corr), arr.ind = TRUE)
data.frame(var1 = rownames(corr)[coppie[, 1]], var2 = colnames(corr)[coppie[, 2]],
           r = round(corr[coppie], 3))

# 2c. Analisi delle componenti principali: le due classi sono gia' ben separate
#     lungo la prima componente, il che spiega le accuratezze elevate.
pca <- prcomp(data_cancer[, variabili], scale. = TRUE)
var_spiegata <- summary(pca)$importance["Proportion of Variance", ]
round(cumsum(var_spiegata)[1:10], 3)
g <- ggplot(data.frame(pca$x[, 1:2], diagnosis = data_cancer$diagnosis),
            aes(PC1, PC2, colour = diagnosis)) +
  geom_point(alpha = 0.7) +
  scale_colour_manual(values = c(M = col_M, B = col_B),
                      labels = c(M = "Maligno", B = "Benigno")) +
  labs(x = sprintf("PC1 (%.1f%% della varianza)", 100 * var_spiegata[1]),
       y = sprintf("PC2 (%.1f%% della varianza)", 100 * var_spiegata[2]),
       colour = "Diagnosi", title = "Proiezione sulle prime due componenti principali")
mostra(g, "03_pca", 7, 5)


# ---- 3. Training set e test set ----------------------------------------------

# Divisione stratificata 70/30: le proporzioni M/B sono le stesse nei due insiemi
set.seed(3)
idx_train <- createDataPartition(data_cancer$diagnosis, p = 0.7, list = FALSE)
train <- data_cancer[idx_train, ]
test  <- data_cancer[-idx_train, ]
rbind(train = table(train$diagnosis), test = table(test$diagnosis))

# Tutti i modelli vengono scelti e confrontati sugli STESSI fold di una
# 10-fold cross-validation ripetuta 3 volte, usando l'AUC come criterio.
set.seed(10)
fold <- createMultiFolds(train$diagnosis, k = 10, times = 3)
controllo <- trainControl(method = "repeatedcv", number = 10, repeats = 3,
                          index = fold, classProbs = TRUE,
                          summaryFunction = twoClassSummary,
                          savePredictions = "final")


# ---- 4. Albero di classificazione --------------------------------------------

# 4a. Albero completo (cp = 0) e cross-validation interna di rpart
set.seed(3)
albero_completo <- rpart(diagnosis ~ ., data = train, method = "class",
                         control = rpart.control(cp = 0, xval = 10))
tabella_cp <- as.data.frame(albero_completo$cptable)
tabella_cp

# Regola 1-SE: il piu' piccolo albero con errore CV entro una deviazione
# standard dal minimo (piu' semplice e quasi altrettanto accurato)
i_min    <- which.min(tabella_cp$xerror)
soglia1se <- tabella_cp$xerror[i_min] + tabella_cp$xstd[i_min]
cp_1se   <- tabella_cp$CP[which(tabella_cp$xerror <= soglia1se)[1]]

g <- ggplot(tabella_cp, aes(nsplit + 1, xerror)) +
  geom_hline(yintercept = soglia1se, linetype = "dashed", colour = col_arancio) +
  geom_errorbar(aes(ymin = xerror - xstd, ymax = xerror + xstd),
                width = 0.2, colour = "grey60") +
  geom_line(colour = col_blu) +
  geom_point(colour = col_blu, size = 2) +
  labs(x = "Numero di nodi terminali", y = "Errore relativo in cross-validation",
       title = "Potatura cost-complexity",
       subtitle = "Linea tratteggiata: soglia della regola 1-SE")
mostra(g, "04_cv_potatura", 7, 4.5)

# 4b. Albero potato, valutato sugli stessi fold degli altri modelli
set.seed(3)
fit_albero <- train(diagnosis ~ ., data = train, method = "rpart",
                    metric = "ROC", trControl = controllo,
                    tuneGrid = data.frame(cp = cp_1se))
albero_potato <- fit_albero$finalModel
albero_potato

mostra_base(function() {
  par(mfrow = c(1, 1), xpd = NA)
  plot(albero_potato, col = "grey", lwd = 3, uniform = TRUE, margin = 0.08)
  text(albero_potato, pretty = 1, col = col_blu, cex = 0.8, use.n = TRUE)
}, "05_albero_potato")


# ---- 5. Bagging e random forest ----------------------------------------------

# Bagging = random forest con mtry = p (tutte le variabili a ogni split)
set.seed(4)
fit_bag <- train(diagnosis ~ ., data = train, method = "rf", metric = "ROC",
                 trControl = controllo, tuneGrid = data.frame(mtry = p),
                 ntree = 500, importance = TRUE)

# Random forest: mtry scelto in cross-validation invece di essere fissato
set.seed(1)
fit_rf <- train(diagnosis ~ ., data = train, method = "rf", metric = "ROC",
                trControl = controllo,
                tuneGrid = data.frame(mtry = c(1, 2, 3, 4, 5, 6, 8, 10, 15, 20)),
                ntree = 500, importance = TRUE)
fit_rf$bestTune

g <- ggplot(fit_rf$results, aes(mtry, ROC)) +
  geom_errorbar(aes(ymin = ROC - ROCSD / sqrt(30), ymax = ROC + ROCSD / sqrt(30)),
                width = 0.3, colour = "grey60") +
  geom_line(colour = col_blu) + geom_point(colour = col_blu, size = 2) +
  geom_vline(xintercept = fit_rf$bestTune$mtry, linetype = "dashed",
             colour = col_arancio) +
  labs(x = "mtry (variabili candidate a ogni split)", y = "AUC in cross-validation",
       title = "Random forest: scelta di mtry", subtitle = "Barre: +/- 1 errore standard")
mostra(g, "06_rf_mtry", 7, 4.5)

# Errore out-of-bag in funzione del numero di alberi
errore_oob <- function(modello, nome) {
  err <- modello$err.rate
  data.frame(Modello = nome,
             Alberi  = rep(seq_len(nrow(err)), times = ncol(err)),
             Tipo    = factor(rep(colnames(err), each = nrow(err)),
                              levels = c("OOB", "M", "B")),
             Errore  = as.vector(err))
}
oob <- rbind(errore_oob(fit_bag$finalModel, "Bagging"),
             errore_oob(fit_rf$finalModel, "Random forest"))
g <- ggplot(oob, aes(Alberi, Errore, colour = Tipo)) +
  geom_line() +
  facet_wrap(~ Modello) +
  scale_colour_manual(values = c(OOB = "black", M = col_M, B = col_B),
                      labels = c(OOB = "OOB complessivo", M = "Maligni", B = "Benigni")) +
  labs(x = "Numero di alberi", y = "Errore out-of-bag", colour = NULL,
       title = "Errore out-of-bag: bagging vs random forest")
mostra(g, "07_errore_oob", 9, 4.5)


# ---- 6. Boosting -------------------------------------------------------------

# Il numero di alberi e la profondita' vengono scelti in cross-validation
# (nella versione originale erano fissati a 5000 e 4). Non serve piu' convertire
# la diagnosi in 0/1 sovrascrivendo il data frame: se ne occupa caret.
griglia_gbm <- expand.grid(n.trees = seq(100, 3000, by = 100),
                           interaction.depth = c(1, 2, 4),
                           shrinkage = 0.01,
                           n.minobsinnode = 10)
set.seed(1)
fit_gbm <- train(diagnosis ~ ., data = train, method = "gbm", metric = "ROC",
                 trControl = controllo, tuneGrid = griglia_gbm, verbose = FALSE)
fit_gbm$bestTune

g <- ggplot(fit_gbm$results,
            aes(n.trees, ROC, colour = factor(interaction.depth))) +
  geom_line() +
  scale_colour_manual(values = c("1" = col_B, "2" = col_blu, "4" = col_arancio)) +
  labs(x = "Numero di alberi", y = "AUC in cross-validation",
       colour = "Profondita'", title = "Boosting: scelta degli iperparametri")
mostra(g, "08_gbm_tuning", 7, 4.5)


# ---- 7. Confronto sul test set -----------------------------------------------

modelli <- list("Albero" = fit_albero, "Bagging" = fit_bag,
                "Random forest" = fit_rf, "Boosting" = fit_gbm)

# probabilita' stimata di tumore maligno per ogni osservazione del test set
prob_test <- lapply(modelli, function(m) predict(m, test, type = "prob")[, "M"])

classifica <- function(prob, soglia = 0.5) {
  factor(ifelse(prob >= soglia, "M", "B"), levels = c("M", "B"))
}

matrici <- lapply(prob_test, function(pr)
  confusionMatrix(classifica(pr), test$diagnosis, positive = "M"))
curve_roc <- lapply(prob_test, function(pr)
  roc(test$diagnosis, pr, levels = c("B", "M"), direction = "<", quiet = TRUE))

matrici

risultati_test <- do.call(rbind, lapply(names(modelli), function(nome) {
  cm <- matrici[[nome]]
  data.frame(Modello        = nome,
             Accuratezza    = cm$overall[["Accuracy"]],
             IC_inf         = cm$overall[["AccuracyLower"]],
             IC_sup         = cm$overall[["AccuracyUpper"]],
             Sensibilita    = cm$byClass[["Sensitivity"]],
             Specificita    = cm$byClass[["Specificity"]],
             AUC            = as.numeric(auc(curve_roc[[nome]])),
             Falsi_negativi = cm$table["B", "M"],
             Falsi_positivi = cm$table["M", "B"])
}))
risultati_test$Modello <- factor(risultati_test$Modello, levels = names(modelli))
print(risultati_test, digits = 3)

# Accuratezza con intervallo di confidenza al 95% (al posto del boxplot di tre
# numeri della versione originale)
g <- ggplot(risultati_test, aes(Accuratezza, Modello, colour = Modello)) +
  geom_errorbarh(aes(xmin = IC_inf, xmax = IC_sup), height = 0.2) +
  geom_point(size = 3) +
  scale_colour_manual(values = col_modelli, guide = "none") +
  labs(x = "Accuratezza sul test set (IC 95%)", y = NULL,
       title = "Accuratezza dei modelli sul test set")
mostra(g, "09_accuratezza_test", 7, 3.5)

g <- ggroc(curve_roc, legacy.axes = TRUE, linewidth = 0.9) +
  geom_abline(linetype = "dotted", colour = "grey50") +
  scale_colour_manual(values = col_modelli,
                      labels = sprintf("%s (AUC = %.3f)", names(modelli),
                                       risultati_test$AUC)) +
  coord_equal() +
  labs(x = "1 - Specificita'", y = "Sensibilita'", colour = NULL,
       title = "Curve ROC sul test set") +
  theme(legend.position = c(0.65, 0.2))
mostra(g, "10_curve_roc", 6, 6)


# ---- 8. Confronto robusto in cross-validation --------------------------------

# Con un test set di ~170 osservazioni gli intervalli di confidenza si
# sovrappongono: un singolo split non basta per stabilire quale metodo sia
# migliore. Sfruttando i 30 fold comuni si ottiene una distribuzione delle
# prestazioni per ogni modello e si possono fare confronti appaiati.
ricampionamenti <- resamples(modelli)
summary(ricampionamenti)

valori <- ricampionamenti$values
metriche <- c(ROC = "AUC", Sens = "Sensibilita'", Spec = "Specificita'")
cv_lungo <- do.call(rbind, lapply(names(modelli), function(nome)
  do.call(rbind, lapply(names(metriche), function(m)
    data.frame(Modello = nome, Metrica = metriche[[m]],
               Valore = valori[[paste(nome, m, sep = "~")]])))))
cv_lungo$Modello <- factor(cv_lungo$Modello, levels = names(modelli))
cv_lungo$Metrica <- factor(cv_lungo$Metrica, levels = metriche)

g <- ggplot(cv_lungo, aes(Modello, Valore, fill = Modello)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.8) +
  facet_wrap(~ Metrica, scales = "free_y") +
  scale_fill_manual(values = col_modelli, guide = "none") +
  labs(x = NULL, y = NULL,
       title = "Prestazioni su 10-fold cross-validation ripetuta 3 volte") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
mostra(g, "11_confronto_cv", 10, 4.5)

# differenze appaiate tra modelli (t-test con correzione di Bonferroni)
differenze <- diff(ricampionamenti)
summary(differenze)


# ---- 9. Importanza delle variabili -------------------------------------------

importanza <- do.call(rbind, lapply(names(modelli), function(nome) {
  imp <- varImp(modelli[[nome]], scale = TRUE)$importance
  data.frame(Modello = nome, Variabile = rownames(imp), Importanza = imp[, 1])
}))
importanza$Modello <- factor(importanza$Modello, levels = names(modelli))

# le 12 variabili piu' importanti in media sui quattro modelli
media_imp <- sort(tapply(importanza$Importanza, importanza$Variabile, mean),
                  decreasing = TRUE)
top_var <- names(media_imp)[1:12]

g <- ggplot(subset(importanza, Variabile %in% top_var),
            aes(Importanza, factor(Variabile, levels = rev(top_var)),
                colour = Modello)) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_colour_manual(values = col_modelli) +
  labs(x = "Importanza (scala 0-100)", y = NULL, colour = NULL,
       title = "Importanza delle variabili nei diversi modelli")
mostra(g, "12_importanza", 8, 5)


# ---- 10. NUOVO: come decide il boosting? -------------------------------------

# Grafici di dipendenza parziale: effetto marginale di una variabile sulla
# probabilita' stimata di tumore maligno, mediando sulle altre variabili.
modello_gbm <- fit_gbm$finalModel
alberi_gbm  <- fit_gbm$bestTune$n.trees
influenza   <- summary(modello_gbm, n.trees = alberi_gbm, plotit = FALSE)
influenza[1:10, ]
var_pdp <- as.character(influenza$var[1:4])

# caret codifica la prima classe ("M") come 1: la risposta e' P(maligno)
pdp <- do.call(rbind, lapply(var_pdp, function(v) {
  griglia <- plot(modello_gbm, i.var = v, n.trees = alberi_gbm,
                  type = "response", return.grid = TRUE)
  data.frame(Variabile = v, x = griglia[[v]], prob = griglia$y)
}))
pdp$Variabile <- factor(pdp$Variabile, levels = var_pdp)
tacche <- do.call(rbind, lapply(var_pdp, function(v)
  data.frame(Variabile = factor(v, levels = var_pdp), x = train[[v]])))

g <- ggplot(pdp, aes(x, prob)) +
  geom_rug(data = tacche, aes(x = x), inherit.aes = FALSE,
           alpha = 0.3, colour = "grey40") +
  geom_line(colour = col_arancio, linewidth = 1) +
  facet_wrap(~ Variabile, scales = "free_x") +
  labs(x = NULL, y = "P(maligno) - dipendenza parziale",
       title = "Boosting: effetto delle variabili piu' influenti")
mostra(g, "13_dipendenza_parziale", 8, 6)

# Dipendenza parziale congiunta delle due variabili piu' influenti: mostra la
# "zona di rischio" individuata dal modello e l'eventuale interazione
griglia2 <- plot(modello_gbm, i.var = var_pdp[1:2], n.trees = alberi_gbm,
                 type = "response", return.grid = TRUE)
names(griglia2)[1:2] <- c("v1", "v2")
g <- ggplot(griglia2, aes(v1, v2, fill = y)) +
  geom_raster(interpolate = TRUE) +
  geom_point(data = train, aes(.data[[var_pdp[1]]], .data[[var_pdp[2]]],
                               shape = diagnosis),
             inherit.aes = FALSE, alpha = 0.6, size = 1.3) +
  scale_fill_gradient2(low = col_B, mid = "white", high = col_M, midpoint = 0.5,
                       limits = c(0, 1)) +
  scale_shape_manual(values = c(M = 17, B = 1),
                     labels = c(M = "Maligno", B = "Benigno")) +
  labs(x = var_pdp[1], y = var_pdp[2], fill = "P(maligno)", shape = "Diagnosi",
       title = "Dipendenza parziale congiunta")
mostra(g, "14_dipendenza_congiunta", 7, 5.5)

# Forza delle interazioni (H di Friedman) tra le coppie delle variabili principali.
# Se la CV sceglie alberi con un solo split (stump) il modello e' additivo:
# le interazioni sono nulle per costruzione e la superficie congiunta qui sopra
# e' la semplice somma degli effetti delle due variabili.
if (fit_gbm$bestTune$interaction.depth > 1) {
  coppie_int <- combn(var_pdp, 2)
  interazioni <- data.frame(
    var1 = coppie_int[1, ], var2 = coppie_int[2, ],
    H = apply(coppie_int, 2, function(cc)
      interact.gbm(modello_gbm, data = train, i.var = cc, n.trees = alberi_gbm))
  )
  print(interazioni[order(-interazioni$H), ])
} else {
  cat("Il boosting scelto usa stump (interaction.depth = 1): modello additivo,",
      "nessuna interazione tra variabili.\n")
}


# ---- 11. NUOVO: scelta della soglia di decisione -----------------------------

# Con la soglia 0.5 un falso negativo (tumore maligno classificato benigno) e un
# falso positivo pesano allo stesso modo. In ambito diagnostico il primo errore
# e' molto piu' grave: si sceglie la soglia che minimizza il costo atteso
#     costo = costo_FN * FN + costo_FP * FP
# usando le previsioni out-of-fold della cross-validation (NON il test set, che
# resta riservato alla verifica finale).
auc_cv   <- sapply(modelli, function(m) max(m$results$ROC))
migliore <- names(which.max(auc_cv))
migliore

oof <- modelli[[migliore]]$pred
oof <- aggregate(M ~ rowIndex + obs, data = oof, FUN = mean)  # media sulle 3 ripetizioni

soglie <- seq(0.01, 0.99, by = 0.01)
curva_soglie <- do.call(rbind, lapply(soglie, function(s) {
  pred <- classifica(oof$M, s)
  FN <- sum(pred == "B" & oof$obs == "M")
  FP <- sum(pred == "M" & oof$obs == "B")
  data.frame(Soglia       = s,
             Sensibilita  = mean(pred[oof$obs == "M"] == "M"),
             Specificita  = mean(pred[oof$obs == "B"] == "B"),
             Costo        = (costo_FN * FN + costo_FP * FP) / nrow(oof))
}))
soglia_ottima <- curva_soglie$Soglia[which.min(curva_soglie$Costo)]
soglia_ottima

curva_lunga <- data.frame(
  Soglia = rep(curva_soglie$Soglia, 3),
  Misura = factor(rep(c("Sensibilita'", "Specificita'", "Costo medio"),
                      each = nrow(curva_soglie)),
                  levels = c("Sensibilita'", "Specificita'", "Costo medio")),
  Valore = c(curva_soglie$Sensibilita, curva_soglie$Specificita, curva_soglie$Costo)
)
g <- ggplot(curva_lunga, aes(Soglia, Valore, colour = Misura)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 0.5, linetype = "dotted", colour = "grey40") +
  geom_vline(xintercept = soglia_ottima, linetype = "dashed", colour = col_arancio) +
  scale_colour_manual(values = c(col_M, col_B, "black")) +
  labs(x = "Soglia su P(maligno)", y = NULL, colour = NULL,
       title = sprintf("%s: scelta della soglia (costo FN = %g x costo FP)",
                       migliore, costo_FN / costo_FP),
       subtitle = sprintf("Tratteggio: soglia ottima = %.2f; puntini: soglia 0.5",
                          soglia_ottima))
mostra(g, "15_scelta_soglia", 8, 4.5)

# Verifica sul test set: soglia 0.5 contro soglia ottima
cm_05     <- confusionMatrix(classifica(prob_test[[migliore]], 0.5),
                             test$diagnosis, positive = "M")
cm_ottima <- confusionMatrix(classifica(prob_test[[migliore]], soglia_ottima),
                             test$diagnosis, positive = "M")
cm_05$table
cm_ottima$table
confronto_soglie <- data.frame(
  Soglia         = c(0.5, soglia_ottima),
  Sensibilita    = c(cm_05$byClass[["Sensitivity"]], cm_ottima$byClass[["Sensitivity"]]),
  Specificita    = c(cm_05$byClass[["Specificity"]], cm_ottima$byClass[["Specificity"]]),
  Falsi_negativi = c(cm_05$table["B", "M"], cm_ottima$table["B", "M"]),
  Falsi_positivi = c(cm_05$table["M", "B"], cm_ottima$table["M", "B"])
)
print(confronto_soglie, digits = 3)


# ---- 12. NUOVO: quante variabili servono davvero? ----------------------------

# Recursive feature elimination con random forest: per ogni fold si calcola
# l'importanza delle variabili e si rimuovono progressivamente le meno utili.
# Ripetere la selezione dentro ogni fold evita una stima troppo ottimistica.
funzioni_rfe <- rfFuncs
funzioni_rfe$summary <- twoClassSummary
set.seed(20)
fit_rfe <- rfe(x = train[, variabili], y = train$diagnosis,
               sizes = c(1:8, 10, 12, 15, 20, 25),
               metric = "ROC",
               rfeControl = rfeControl(functions = funzioni_rfe,
                                       method = "repeatedcv", number = 10,
                                       repeats = 3, index = fold))
fit_rfe

# il modello piu' piccolo con AUC entro l'1% da quella massima
n_parsimonioso <- pickSizeTolerance(fit_rfe$results, metric = "ROC",
                                    tol = 1, maximize = TRUE)
var_parsimoniose <- fit_rfe$optVariables[seq_len(n_parsimonioso)]
n_parsimonioso
var_parsimoniose

g <- ggplot(fit_rfe$results, aes(Variables, ROC)) +
  geom_errorbar(aes(ymin = ROC - ROCSD / sqrt(30), ymax = ROC + ROCSD / sqrt(30)),
                width = 0.4, colour = "grey60") +
  geom_line(colour = col_blu) + geom_point(colour = col_blu, size = 2) +
  geom_vline(xintercept = n_parsimonioso, linetype = "dashed", colour = col_arancio) +
  labs(x = "Numero di variabili", y = "AUC in cross-validation",
       title = "Recursive feature elimination (random forest)",
       subtitle = sprintf("Tratteggio: modello parsimonioso con %d variabili",
                          n_parsimonioso))
mostra(g, "16_selezione_variabili", 7, 4.5)

# random forest con le sole variabili selezionate, valutata sul test set
set.seed(1)
rf_ridotta <- randomForest(x = train[, var_parsimoniose, drop = FALSE],
                           y = train$diagnosis, ntree = 500)
prob_ridotta <- predict(rf_ridotta, test[, var_parsimoniose, drop = FALSE],
                        type = "prob")[, "M"]
cm_ridotta <- confusionMatrix(classifica(prob_ridotta), test$diagnosis,
                              positive = "M")
cm_ridotta
auc(roc(test$diagnosis, prob_ridotta, levels = c("B", "M"), direction = "<",
        quiet = TRUE))


# ---- 13. Riepilogo -----------------------------------------------------------

cat("\n================ RIEPILOGO ================\n")
cat("Risultati sul test set (soglia 0.5):\n")
print(risultati_test[, c("Modello", "Accuratezza", "Sensibilita", "Specificita",
                         "AUC", "Falsi_negativi")], digits = 3, row.names = FALSE)
cat(sprintf("\nAUC media in CV: %s\n",
            paste(sprintf("%s %.4f", names(auc_cv), auc_cv), collapse = ", ")))
cat(sprintf("Soglia ottima per %s (costo FN = %g): %.2f -> falsi negativi sul test set: %d -> %d\n",
            migliore, costo_FN, soglia_ottima,
            confronto_soglie$Falsi_negativi[1], confronto_soglie$Falsi_negativi[2]))
cat(sprintf("Random forest con %d variabili su %d: accuratezza %.3f sul test set\n",
            n_parsimonioso, p, cm_ridotta$overall[["Accuracy"]]))
