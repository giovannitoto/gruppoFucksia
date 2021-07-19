load("adenoma.rda")

head(rowData(se))

# Oggetto generale
se
dim(se)
# le righe sono i soggetti, le colonne le specie di batteri (vero?)

# Guardo colData, ovvero le caratteristiche degli individui
dim(colData(se))
colData(se)

#assay: relative abundance
dim(assay(se))
rownames(assay(se))

#togliamo gli zeri
dim(se[rowMaxs(assay(se)) > 0.1,])
se <- se[rowMaxs(assay(se)) > 0.1,]
rel_abu <- assay(se)
row_max <- rowMaxs(rel_abu)
head(row_max)
boxplot(row_max)
summary(row_max) #si può stabilire una soglia

#Si prova una Random Forest per vedere se ha senso mettere una soglia per filtrare le specie di batteri
#creo il dataset
dim(colData(se))
dim(assay(se))
idx <- rowMaxs(assay(se)) > 0.1
s <- as.data.frame(cbind(colData(se), t(assay(se))))
dim(s)

### RANDOM FOREST (CV)
library(randomForest)
library(e1071)

set.seed(28)
forest.tune <- tune(randomForest, study_condition ~ ., data=s,
                    ranges=list(mtry=2*(2:4)),
                    tunecontrol=tune.control(cross=5),
                    ntree=100, na.action = na.roughfix)
best.mtry <- min(forest.tune$best.parameters, NCOL(s))
# Grafico
plot(forest.tune)
abline(v=best.mtry, lty=2, col=2)

final.forest <- randomForest(y ~ ., data=s, mtry=best.mtry)

# Grafico per importanza delle var. nella random forest
varImpPlot(fit.forest)

# Previsione
y.forest.cv <- predict(final.forest, newdata=v)
et.forest.cv <- tabella.sommario(y.forest.cv, v$y)
e.forest.cv <- 1 - sum(diag(et.forest.cv))/sum(et.forest.cv)
tab_confronto <- rbind(tab_confronto,c("Random Forest (CV)", e.forest.cv))
y.forest.cv.prob <- predict(final.forest, newdata=v, type="prob")
a.forest.cv <- lift.roc(y.forest.cv.prob[,2], ynum.v, type="crude")