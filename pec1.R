library("SummarizedExperiment")
library(dplyr)
library(tibble)
library(bookdown)
library(purrr)
library(ggplot2)
library(nortest)
library(bestNormalize)
library(PCAtools)
library("UpSetR")
library(naniar)
library(pheatmap)


datos<-read.csv("DataValues_S013.csv")
#Eliminamos esa columna añadida como el número de líneas
datos <- datos[, -1] 
mdatos<-read.csv("DataInfo_S013.csv")
mdatos <- mdatos[, -1] 


#Guardamos el rowData
info_samples <- datos[, 1:5]
#Modificamos nuestros datos para crear nuestro objeto 
#SummarizedExperiment
datos <- datos[, -c(1:5)]
mdatos <- mdatos[-c(1:5), ]

#transformacion necesaria a matriz
matDatos <- as.matrix(datos)


#Creamos nuestros ids
id = paste0("SUBJ_",info_samples$SUBJECTS,"_",
            gsub(" ", "_", info_samples$SURGERY)
            , "_","G:",info_samples$Group)
#Lo añadimos a nuestro rowData
info_samples <- info_samples %>%
  mutate(id = id)

#Lo establecemos como nombre de nuestras filas
info_samples<- column_to_rownames(info_samples,'id')

#Establecemos las variables categóricas como factores
info_samples <- info_samples %>%
  mutate(SURGERY = as.factor(SURGERY),
         GENDER = as.factor(GENDER),
         Group = as.factor(Group))

#Finalmente creamos nuestro objeto
sum_ex <- SummarizedExperiment(assays = list(raw=matDatos),
                             rowData=info_samples, 
                             colData=mdatos)

info_samples %>%
  #Esto crea una lista donde cada elemento es la tabla
  map(~ table(.))

#frecuencia relativa porncentual
info_samples %>%
  map(~ prop.table(table(.))*100)

#breve resumen estadístico de variables numéricas
#Como tenemos tantas variables, para mostrar el output en el documento, 
#se seleccionan únicamentelas 10 primeras
summary(datos[1:10])

#visualizacion de datos faltantes
datos2 <- datos %>%
  mutate(id = id)

datos2<- column_to_rownames(datos2,'id')

gg_miss_upset(as.data.frame(t(datos2)),text.scale = 0.3,
              nsets=n_var_miss(as.data.frame(t(datos2))))
n_var_miss(t(datos2))
gg_miss_case(as.data.frame(t(datos2)))

#correcion de datos faltantes
rawValues <- assay(sum_ex, "raw")
NA_clean <- rawValues
NA_clean[is.na(NA_clean)] <- 1
assays(sum_ex)$NA_clean <- NA_clean

#normalizacion 
log <- log(assay(sum_ex, "NA_clean"))
assays(sum_ex)$log <- log
save(sum_ex, file = "PEC1_SE.Rda")

#Análisis de componentes principales (PCA)
pca_res <- pca(t(assay(sum_ex, "NA_clean")), metadata = rowData(sum_ex))
pairsplot(pca_res)
biplot(pca_res,
       colby = 'SURGERY', colkey = c('by pass' = 'forestgreen', 'tubular' = 'purple'),
       colLegendTitle = 'SurgeryType',
       # encircle config
       encircle = TRUE,
       encircleFill = TRUE,
       hline = 0, vline = c(-25, 0, 25),
       legendPosition = 'top', legendLabSize = 16, legendIconSize = 8.0)
##correlaciones
eigencorplot(pca_res,
             components = getComponents(pca_res),
             metavars = colnames(info_samples),
             col = c('white', 'lightblue', 'turquoise', 'blue', 'darkblue'),
             cexCorval = 0.25,
             cexLabX= 0.5,
             cexLabY= 0.6,
             fontCorval = 2,
             posLab = 'all',
             rotLabX = 45,
             scale = TRUE,
             main = bquote(PCA ~ Pearson ~ r^2 ~ correlates),
             plotRsquared = TRUE,
             corFUN = 'pearson',
             corUSE = 'pairwise.complete.obs',
             corMultipleTestCorrection = 'BH',
             signifSymbols = c('****', '***', '**', '*', ''),
             signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1))
#diagrama de cajas
#para nque no se corten los nombres en el eje y
par(mar = c(5, 10, 4, 2))
boxplot(t(assay(sum_ex, "log")),
        xlab = "Concentración",
        horizontal = TRUE,  las = 2, main= "Sujeto/Concentración",
        col=palette.colors())


