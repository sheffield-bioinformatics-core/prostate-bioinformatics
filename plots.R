sampleinfo <- read.delim("meta_data/sampleInfo_corrected.txt")
raw <- read_tsv("raw_data/selected_prostate_tcga_raw.tsv")
genes <- pull(raw,X)
cts <- as.matrix(raw[,-1])
rownames(cts) <- genes
dds <- DESeqDataSetFromMatrix(countData = cts, 
                              colData = sampleinfo,
                              design = ~Status)

p <- sampleinfo %>% 
  mutate(LibrarySize = colSums(assay(dds))/1e6) %>% 
  ggplot(aes(x = ID,y=LibrarySize)) + geom_col(fill = "steelblue") + ylab("Millions of Reads") + theme(axis.text.x = element_text(angle=90))

ggsave(p, filename = "images/challenge1.png")

p <- plotPCA(vst(dds), intgroup="gleason_score",returnData = TRUE) %>% mutate(LibrarySize= colSums(assay(dds))/1e6) %>% 
  ggplot(aes(x = PC1, y=PC2,size=LibrarySize,col=gleason_score)) + geom_point()
ggsave(p, filename="images/challenge4.png")


de_gleason <- readRDS("Robjects/de_gleason.rds")

results_gleason <- as.data.frame(results(de_gleason, contrast=c("gleason_score","G9","G6"))) %>% 
  rownames_to_column("GeneID")

results_ordered <- arrange(results_gleason, padj)
top_genes <- results_ordered$GeneID[1:75]


pheatmap(assay(vsd)[top_genes,],
         annotation_col = sampleInfo,scale = "row",
         labels_col = dds$ID)