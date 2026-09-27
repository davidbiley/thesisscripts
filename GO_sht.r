library(tidyverse)
library(clusterProfiler)
library(pathview)
library(org.At.tair.db)

agilist <- read_tsv("project/fvesca_sht_agis.txt", col_names=FALSE) %>%
  head(1000)
agilist <- array(unlist(agilist))

GOenriched <- enrichGO(gene = agilist,
                       OrgDb = org.At.tair.db,
                       keyType = "TAIR",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)

png(filename="GO_enrichment_fvesca_sht.png")
dotplot(GOenriched)
dev.off()

sessionInfo()
