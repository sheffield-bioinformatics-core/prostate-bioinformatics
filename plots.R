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
