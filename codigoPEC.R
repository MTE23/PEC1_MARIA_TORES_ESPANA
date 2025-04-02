#--------------------------------PREPROCESAMIENTO----------------------------------------------------------------

#Seleccionamos y cargamos el dataset selecionado 
#que se obtuvo al descargar del repositorio de github mencionado en el enunciado de la PEC1.
datos<-read.csv("DataValues_S013.csv")
##Eliminamos esa columna añadida como el número de líneas
datos <- datos[, -1] 
mdatos<-read.csv("DataInfo_S013.csv")
mdatos <- mdatos[, -1] 

#Si exploramos nuestros datos nos datos cuenta de que las cinco primeras columnas corresponden al 
#`rowData` (osea la información de las muestras) de forma que modificamos nuestros datos acorde.
##Guardamos el rowData
info_samples <- datos[, 1:5]
##Modificamos nuestros datos para crear nuestro objeto 
##SummarizedExperiment
datos <- datos[, -c(1:5)]
mdatos <- mdatos[-c(1:5), ]

#Además vemos que los mdatos están en formato `data.frame`.
#Esto al algo que deberemos cambiar ya que los objetos de la clase `SummarizedExperiment` 
#requieren matrices numéricas y dadas como `list`.

datos_mat <- as.matrix(datos)

#Además, los objetos de la clase `SummarizedExperiment` requieren que las filas de los datos (`assay()`) y las filas de la información de las muestras (`colData()`) coincidan en los nombres y que los datos categóricos sean factores.
##Creamos nuestros ids
id = paste0("SUBJ_",info_samples$SUBJECTS,"_",
            gsub(" ", "_", info_samples$SURGERY)
            , "_","G:",info_samples$Group)
##Lo añadimos a nuestro rowData
info_samples <- info_samples %>%
  mutate(id = id)

##Lo establecemos como nombre de nuestras filas
info_samples<- column_to_rownames(info_samples,'id')

##Establecemos las variables categóricas como factores
info_samples <- info_samples %>%
  mutate(SURGERY = as.factor(SURGERY),
         GENDER = as.factor(GENDER),
         Group = as.factor(Group))
##Finalmente creamos nuestro objeto
sum_ex <- SummarizedExperiment(assays = list(raw=datos_mat),
                               rowData=info_samples, 
                               colData=mdatos)

#--------------------------------ANALISIS EXPLORATORIO----------------------------------------------------------------
#Tabla de frecuencias de variables
info_samples %>%
  ##Esto crea una lista donde cada elemento es la tabla
  map(~ table(.))

#frecuencia relativa porcentual.
info_samples[c(2:5)] %>%
  map(~ prop.table(table(.))*100)

##Como tenemos tantas variables, para mostrar el output en el documento, 
##se seleccionan únicamentelas 10 primeras
summary(datos[1:10])


datos2 <- datos %>%
  mutate(id = id)
datos2<- column_to_rownames(datos2,'id')
gg_miss_upset(as.data.frame(t(datos2)),text.scale = 0.3,
              nsets=n_var_miss(as.data.frame(t(datos2))))

#--------------------------------TRATAMIENTO DE VALORES FALTANTES----------------------------------------------------------------
#Como todos tienen en algún momento valores faltantes, la opción de quedarnos solo con aquellos que
#estén completos se descarta. Por tanto entonces lo que hacemos es designar 1 a los valores faltantes.

NA_clean <- assay(sum_ex, "raw")
NA_clean[is.na(NA_clean)] <- 1
assays(sum_ex)$NA_clean <- NA_clean
#--------------------------------LOGARITMOS Y GUARDADO EN FORMATO BINARIO----------------------------------------------------------------
#aplicamos logaritmos a nuestros datos sin datos faltantes
log <- log(assay(sum_ex, "NA_clean"))
#la guardamos en nuestro contenedor
assays(sum_ex)$log <- log
#guardamos
save(sum_ex, file = "PEC1_SE.Rda")

#--------------------------------PCA----------------------------------------------------------------
#Análisis de componentes principales (PCA)
pca_res <- pca(t(assay(sum_ex, "NA_clean")), metadata = rowData(sum_ex))
#correlacion entre los diferentes componentes principales
pairsplot(pca_res)
#gráfico de correlación de los componentes principales con las metavariables de las muestras.  El gráfico 
#mostrará entonces la relación entre cada componente y las variables 
#de nuestras metavariables. Se calcula la correlación de Pearson, con corrección por pruebas múltiples
eigencorplot(pca_res,
             components = getComponents(pca_res),
             metavars = colnames(info_samples),
             #paleta de colores
             col = c('white', 'lightblue', 'turquoise', 'blue', 'darkblue'),
             cexCorval = 0.25,
             #tamaño de las fuentes
             cexLabX= 0.5,
             cexLabY= 0.6,
             fontCorval = 2,
             posLab = 'all',
             rotLabX = 45,
             scale = TRUE,
             #leyenda
             main = bquote(PCA ~ Pearson ~ r^2 ~ correlates),
             plotRsquared = TRUE,
             corFUN = 'pearson',
             corUSE = 'pairwise.complete.obs',
             corMultipleTestCorrection = 'BH',
             #para mostrar alta correlación
             signifSymbols = c('****', '***', '**', '*', ''),
             signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1))
#para visualizar las muestras y las variables en el espacio de los componentes principales
biplot(pca_res,
       #variable por la que agrupamos y los colores por los que queremos determinar los grupos
       colby = 'GENDER', colkey = c('M' = 'forestgreen', 'F' = 'purple'),
       colLegendTitle = 'Gender',
       #para que se cree el coloreado de los grupos
       encircle = TRUE,
       encircleFill = TRUE,
       hline = 0, vline = c(-25, 0, 25),
       #posicion de la leyenda y tamaños
       legendPosition = 'top', legendLabSize = 16, legendIconSize = 8.0)

#para evitar que los nombres de los sujetos se salgan del margen
par(mar = c(5, 10, 4, 2))
boxplot(t(assay(sum_ex, "log")),
        xlab = "Concentración",
        horizontal = TRUE,  las = 2, main= "Sujeto/Concentración",
        col=palette.colors())
