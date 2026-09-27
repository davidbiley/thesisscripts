BiocManager::install("DESeq2")
BiocManager::install("tximport")
BiocManager::install("tidyverse")
BiocManager::install("ggplot2")
BiocManager::install("ggrepel")
library("DESeq2")
library("tximport")
library("tidyverse")
library("ggplot2")
library("ggrepel")

csvfile <- "project/metadatacucu2.csv"
sampleTable <- read.csv(csvfile, row.names=1, sep="\t")%>%
  rev()
sampleTable <- sampleTable[rev(rownames(sampleTable)),]
sampleTable$group <- as.factor( sampleTable$group )

summary(sampleTable)

samples <- rownames(sampleTable)
files <- file.path("project", samples, "abundance.h5")
names(files) <- rownames(sampleTable)
txi <- tximport(files, type = "kallisto", txOut = TRUE)

ddsMat <- DESeqDataSetFromTximport(txi, colData=sampleTable, design= ~ group )
nrow(ddsMat)

dds <- ddsMat[ rowSums(counts(ddsMat)) > 1, ]
nrow(dds)

rld <- rlog(dds)
ramp <- 1:2/2
cols <- c( rgb(ramp, 0, 0), rgb(0, ramp, 0), rgb(ramp, 0, ramp), rgb(ramp, 0, ramp) )

pcaData <- plotPCA(rld, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = group, label = name)) +
  geom_point(size = 3) +
  geom_text_repel() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw()

ggsave("cpepopca.png",width = 5, height = 6, units = "in",plot = get_last_plot())

dds <- DESeq(dds)
resuvh <- results(dds, contrast = c("group", "uv", "high"), alpha=.05)
resuvn <- results(dds, contrast = c("group", "uv", "normal"), alpha=.05)
reshn <- results(dds, contrast = c("group", "high", "normal"), alpha=.05)

tpm_matrix <- as.data.frame(txi$abundance)
tpm_matrix <- rownames_to_column(tpm_matrix, "geneID") %>%
  mutate(realID=str_extract(geneID,"XP_[0-9]{9}.[0-9]"))
cfulltable <- as.data.frame(resuvn)
cfulltable <- rownames_to_column(cfulltable, "geneID") %>%
  mutate(realID=str_extract(geneID,"XP_[0-9]{9}.[0-9]")) %>%
  left_join(tpm_matrix, by="realID") %>%
  select(realID, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj,CU2,CU3,CU4,CH1, CH3, CH4,CN1,CN3,CN4)

write.table( cfulltable, "noag_cucu_high_normal_results_table.tsv", sep="\t" )

cpepocaropath <- "project/cpepoca_kipes.txt"
cpepoflavpath <- "project/cpepofl_kipes.txt"
cpepomybpath <- "project/mybs_cpepo.tsv"
cpepomybs <- read_tsv(cpepomybpath) %>% mutate(realID=NewMYB)
cpepocarotenoids <- read_tsv(cpepocaropath) %>% mutate(realID=str_extract(ID,"XP_[0-9]{9}.[0-9]")) 
cpepoflavonoids <- read_tsv(cpepoflavpath) %>% mutate(realID=str_extract(ID,"XP_[0-9]{9}.[0-9]"))
cpepokipes <- rbind(cpepocarotenoids,cpepoflavonoids)

cfulltable <- cfulltable %>% 
  left_join(cpepokipes, by="realID") %>%
  left_join(cpepomybs, by="realID") %>%
  select(realID,Gene,ConservedResidues,RefMYB,ConservedRegions,log2FoldChange,padj,CU2,CU3,CU4,CH1, CH3, CH4,CN1,CN3,CN4)

gois <- c("XP_023517857.1", "XP_023545723.1", "XP_023537989.1",	"XP_023524405.1",
          "XP_023554180.1",	"XP_023518138.1",	"XP_023526786.1",	"XP_023530127.1",
          "XP_023541696.1", "XP_023553842.1", "XP_023541702.1",	"XP_023547794.1",
          "XP_023513811.1",	"XP_023526659.1", "XP_023548447.1", "XP_023523754.1",
          "XP_023526588.1", "XP_023523755.1", "XP_023528106.1",	"XP_023532357.1",
          "XP_023553243.1", "XP_023544822.1", "XP_023543404.1", "XP_023517487.1", "XP_023517444.1")

cfulltable <- cfulltable %>%filter(ConservedResidues == 100 | realID %in% gois | RefMYB!="<NA>" ) %>% filter(padj < .0501 )

write.table( cfulltable, "noag_cpepo_DGE_kipes_selected_uvn.tsv", sep="\t" )

sessionInfo()

