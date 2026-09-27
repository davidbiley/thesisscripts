library("BiocManager")
BiocManager::install("tidyverse")
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
library(tidyverse)
library(clusterProfiler)
library(pathview)
library(org.At.tair.db)

aramapping <- read_tsv("project/Athaliana_447_Araport11.protein__v__GCF_002806865.2_ASM280686v2_translated_cds.tsv") %>%
  separate_rows(GCF_002806865.2_ASM280686v2_translated_cds, sep=", ") %>%
  separate_rows(Athaliana_447_Araport11.protein,sep=", ") %>%
  mutate(arabidopsis=gsub("\\.\\d+$","",Athaliana_447_Araport11.protein))%>%
  mutate(protID=str_extract(GCF_002806865.2_ASM280686v2_translated_cds,"XP_[0-9]{9}.[0-9]")) %>%
  dplyr::select(protID,arabidopsis)

degs <- read_tsv("project/noag_cucu_uv_normal_results_table.tsv") %>%
  dplyr::select(protID=realID,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj)

fullmapping <- aramapping %>%
  left_join(degs, by="protID")%>%
  drop_na()

significantGO <- fullmapping %>%
  filter(abs(log2FoldChange) >= 1 & padj <= 0.05) %>%
  distinct(arabidopsis) %>%
  pull(arabidopsis)

foldchange <- fullmapping %>%
  group_by(arabidopsis) %>%
  slice_max(order_by = abs(log2FoldChange), n = 1) %>%
  ungroup() %>%
  filter(!is.na(log2FoldChange) & abs(log2FoldChange) >= 1 & padj < 0.05) %>%
  {setNames(.$log2FoldChange, .$arabidopsis)}

pathview(gene.data  = foldchange,
         pathway.id = "ath00195",
         species    = "ath",
         gene.idtype = "tair",
         same.layer = FALSE,
         limit       = list(gene = 4),
         low         = "gold",
         high        = "firebrick1",)


GOenriched <- enrichGO(gene = significantGO,
                       OrgDb = org.At.tair.db,
                       keyType = "TAIR",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)

png(filename="GO_enrichment_cpepo_uvn.png")
dotplot(GOenriched)
dev.off()
df <- as.data.frame(GOenriched)
write.table(df,"cpepo_uvn_GOs.tsv",sep="\t")

