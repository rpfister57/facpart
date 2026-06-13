## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 5,      # Width in inches
  fig.height = 5,     # Height in inches
  dpi = 300,          # Resolution (DPI)
  dev = "svg"         # Use PNG for sharp rendering
)


## ----gutt91-------------------------------------------------------------------
library(facpart)

Kor <- gutt91$gutt91_cor
Facets <- gutt91$gutt91_var
round(Kor, 2)
Facets



## ----MDS----------------------------------------------------------------------
Kor_D <- smacof::sim2diss(Kor, method = "corr", to.dist = TRUE)

gutt91_mds <- smacof::mds(Kor_D, type = "ordinal")
gutt91_mds


## ----MDSplot------------------------------------------------------------------
plot(gutt91_mds)

plot(gutt91_mds$conf, 
     asp = 1, las = 1)
text(gutt91_mds$conf, labels = Facets$Material, 
     cex = 0.8, pos = 3)
abline(v=0); abline(h=0)


## ----angularPartition1--------------------------------------------------------
plot(gutt91_mds$conf, 
     asp = 1, las = 1)
text(gutt91_mds$conf, labels = Facets$Material, 
     cex = 0.8, pos = 3)

angularPartition(crd = gutt91_mds$conf,
                 group = Facets$Material, add = TRUE)

