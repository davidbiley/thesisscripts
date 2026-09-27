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

csvfile <- "project/metadatafraga2.tsv"
sampleTable <- read.csv(csvfile, row.names=1, sep="\t")
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

ggsave("fvescapca.png",width = 5, height = 6, units = "in",plot = get_last_plot())

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
  left_join(tpm_matrix, by="realID") #%>%
  select(realID, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj,FU1, FU2, FU4,FH1,FH3,FH4,FN1,FN3,FN4) %>% filter(padj < .05 )

write.table( cfulltable, "noag_fraga_high_normal_results_table.tsv", sep="\t" )

cpepocaropath <- "project/fvescaca_kipes.txt"
cpepoflavpath <- "project/fvescafl_kipes.txt"
cpepomybpath <- "project/mybs_fvesca.tsv"
cpepomybs <- read_tsv(cpepomybpath) %>% mutate(realID=NewMYB)
cpepocarotenoids <- read_tsv(cpepocaropath) %>% mutate(realID=str_extract(ID,"XP_[0-9]{9}.[0-9]")) 
cpepoflavonoids <- read_tsv(cpepoflavpath) %>% mutate(realID=str_extract(ID,"XP_[0-9]{9}.[0-9]"))
cpepokipes <- rbind(cpepocarotenoids,cpepoflavonoids)

cfulltable <- cfulltable %>% 
  left_join(cpepokipes, by="realID") %>%
  left_join(cpepomybs, by="realID") %>%
  select(realID,Gene,ConservedResidues,RefMYB,ConservedRegions,log2FoldChange,padj,FU1, FU2, FU4,FH1,FH3,FH4,FN1,FN3,FN4)

gois <- c("XP_004288153.1", "XP_004291469.1", "XP_004297046.1",
          "XP_004309213.1",	"XP_004287841.1",	"XP_004301156.1", "XP_004294698.1",
          "XP_004310213.1", "XP_004310071.1", "XP_004299737.1", "XP_004290595.1",
          "XP_004295153.1", "XP_004294955.1", "XP_004293649.1", "XP_004291518.1",
          "XP_004293648.1", "XP_004303231.1", "XP_004303244.1", "XP_004288323.1",
          "XP_004287992.1", "XP_004304049.1", "XP_004307159.1", "XP_004289513.1",
          "XP_004288495.1", "XP_004288495.1", "XP_004296150.1")

maybes <- c("XP_004291858.1", "XP_004306543.1", "XP_004306542.1", "XP_004307451.1",
           "XP_004287814.1", "XP_004290837.1", "XP_004301809.1")

tpms <- tpm_matrix %>%
  left_join(cpepokipes, by="realID") %>%
  left_join(cpepomybs, by="realID") %>%
  filter(realID %in% maybes) %>%
  select(realID,Gene,RefMYB,FU1, FU2, FU4,FH1,FH3,FH4,FN1,FN3,FN4)
write.table(tpms,"fvesca_maybes.tsv",sep="\t")

cfulltable <- cfulltable %>%filter(ConservedResidues == 100 | realID %in% gois | RefMYB!="<NA>" ) %>% filter(padj < .05 )

write.table( cfulltable, "noag_fvesca_DGE_kipes_selected_hn.tsv", sep="\t" )

sessionInfo()

