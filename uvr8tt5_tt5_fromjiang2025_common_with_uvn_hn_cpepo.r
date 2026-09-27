library("BiocManager")
BiocManager::install("tidyverse")
library(tidyverse)

aramapping <- read_tsv("project/Athaliana_447_Araport11.protein__v__GCF_002806865.2_ASM280686v2_translated_cds.tsv") %>%
  separate_rows(GCF_002806865.2_ASM280686v2_translated_cds, sep=", ") %>%
  separate_rows(Athaliana_447_Araport11.protein,sep=", ") %>%
  mutate(arabidopsis=gsub("\\.\\d+$","",Athaliana_447_Araport11.protein))%>%
  mutate(realID=str_extract(GCF_002806865.2_ASM280686v2_translated_cds,"XP_[0-9]{9}.[0-9]")) %>%
  select(realID,arabidopsis)

udegs <- read_tsv("project/noag_cucu_uv_normal_results_table.tsv") %>%
  select(realID,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,CU2,CU3,CU4,CH1, CH3, CH4,CN1,CN3,CN4)

hdegs <- read_tsv("project/noag_cucu_high_normal_results_table.tsv") %>%
  select(realID,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,CU2,CU3,CU4,CH1, CH3, CH4,CN1,CN3,CN4)

htdegs <- read_tsv("project/noag_cucu_high_normal_results_table.tsv") %>%
  select(realID,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj)

table <- htdegs %>% left_join(udegs, by="realID") %>%
  select(realID,log2FoldChange.x,padj.x,padj.y,CU2,CU3,CU4,CH1, CH3, CH4,CN1,CN3,CN4)

annos <- read_tsv("project/Arabidopsis_thaliana_araport11.anno.txt") %>%
  mutate(arabidopsis= str_extract(arabidopsis,"AT[0-9]G[0-9]{5}")) %>%
  mutate(Anno= gsub(".https:.*","",Anno))

uvr8tt5 <- read_csv("project/DEGs_uvr8_tt5.csv")

tt5 <- read_csv("project/DEGs_tt5.csv")

fullmapping <- aramapping %>%
  left_join(udegs, by="realID") %>%
  left_join(annos, by="arabidopsis") %>%
  group_by(realID) %>%
  filter(arabidopsis %in% uvr8tt5$arabidopsis,padj<=0.05) %>%
  group_by(realID) %>% slice_head() %>% ungroup() %>%
  select(realID,padj,Anno,CU2,CU3,CU4,CH1, CH3, CH4,CN1,CN3,CN4) %>%
  drop_na()

print(fullmapping, n=100)

write.table(fullmapping, "uvr8tt5_fromjiang2025_common_with_uvn_cpepo.tsv", sep="\t" )
