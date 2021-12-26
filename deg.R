library(tximport)
library(DESeq2)
library(tidyverse)
library("GSEABase")
library(clusterProfiler)
library(RColorBrewer)
library(pheatmap)
library(wesanderson)
library(org.Sc.sgd.db)
library(GOstats)
library(edgeR)
library(treemap)

# read in the file from url
samples <- read_csv("Sample_data.csv")
# look at the first 6 lines
samples

tx2gene_map <- read_tsv("GCA_000146045.2_R64_tx2gene.tsv")

txi <- tximport(files = samples$quant_file, type = "salmon", tx2gene = tx2gene_map)

colnames(txi$counts) <- samples$sample

dds <- DESeqDataSetFromTximport(txi = txi, 
                                colData = samples, 
                                design = ~condition)

dds$condition <- relevel(dds$condition, ref = "wt")

dds <- DESeq(dds)

res <- results(dds)
head(res)

summary(res, alpha = 0.05)

rld <- rlog(dds, blind=TRUE)
head(assay(rld), 3)

vsd <- vst(dds, blind = TRUE)
head(assay(vsd), 3)

data1 <- plotPCA(rld, returnData=TRUE)
data1$group<-gsub(" : ","_",as.character(data1$group))
percentVar1 <- round(100 * attr(data1, "percentVar"))

PCA<-ggplot(data1, aes(PC1, PC2, color = condition))+ theme_bw()+
  geom_point(size=9, alpha = 0.8) + scale_colour_manual(values = c("#44aabb","#bbbbbb"))+
  xlab(paste0("PC1: ",percentVar1[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar1[2],"% variance")) +
  theme(text = element_text(size=20)) + ggtitle("rlog PCA")
PCA

data1 <- plotPCA(vsd, returnData=TRUE)
data1$group<-gsub(" : ","_",as.character(data1$group))
percentVar1 <- round(100 * attr(data1, "percentVar"))

PCA<-ggplot(data1, aes(PC1, PC2, color = condition))+ theme_bw()+
  geom_point(size=9, alpha = 0.8) + scale_colour_manual(values = c("#44aabb","#bbbbbb"))+
  xlab(paste0("PC1: ",percentVar1[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar1[2],"% variance")) +
  theme(text = element_text(size=20)) + ggtitle("vst PCA")
PCA

df <- as.data.frame(colData(rld)[,c("condition", "sample")])

mat_colors1<-list(sample = brewer.pal(12, "Paired")[0:6])
names(mat_colors1$sample)<- df$sample

mat_colors <- list(condition = brewer.pal(12, "Paired")[7:8])
names(mat_colors$condition) <- c("wt", "snf2")

genes <- order(res$padj)[1:1000]

pheatmap(assay(rld)[genes, ], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, annotation_colors = c(mat_colors1, mat_colors), fontsize = 12)

df <- as.data.frame(colData(vsd)[,c("condition", "sample")])

pheatmap(assay(vsd)[genes, ], cluster_rows=TRUE, show_rownames=FALSE, show_colnames = FALSE,
         cluster_cols=FALSE, annotation_col=df, annotation_colors = c(mat_colors1, mat_colors), fontsize = 12)

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

summary(res, alpha = 0.05)

GO_df = toTable(org.Sc.sgdGO)
head(GO_df)

goframeData = data.frame(GO_df$go_id, GO_df$Evidence, GO_df$systematic_name)
names(goframeData) = c("GO", "Evidence", "gene_id")
goframeData$GO <- as.character(goframeData$GO)
goframeData$Evidence <- as.character(goframeData$Evidence)
goframeData$gene_id <- as.character((goframeData$gene_id))
head(goframeData)

DEres <- as.data.frame(res)

upFrame <- DEres[DEres$padj <= 0.05 & DEres$log2FoldChange >= 1,]
upFrame <- upFrame[complete.cases(upFrame),]
dim(upFrame)

universe = Lkeys(org.Sc.sgdGO)
genes = row.names(upFrame)
params <- GSEAGOHyperGParams(name="Yeast, wt v. snf2", geneSetCollection= gsc, geneIds = genes, universeGeneIds = universe, ontology = "BP", pvalueCutoff = 0.05, conditional = TRUE, testDirection = "over")
Over <- hyperGTest(params)
Over

SC <- toTable(org.Sc.sgdGENENAME)

DEres <- as.data.frame(res)
DEres$systematic_name<- row.names(DEres)
DEres<- merge(DEres, SC, by="systematic_name")

gene <- DEres$gene_name

gene.df <- bitr(gene, fromType = "GENENAME",
                toType = c("ENTREZID", "ENSEMBL"),
                OrgDb = org.Sc.sgd.db)

names(gene.df)<- c("gene_name", "EntrezID", "ENSEMBL" )

go.df <- bitr(gene, fromType = "GENENAME",
              toType = c("GO"),
              OrgDb = org.Sc.sgd.db)

names(go.df)<- c("gene_name", "GO", "evidence", "ontology")

DEres<- merge(DEres, gene.df, by = "gene_name" )

go_dict <- merge(go.df, gene.df, by= "gene_name")

BPdict <- go_dict[go_dict$ontology == "BP", ]

BPdict<- BPdict[,c(2,6)]

upFrame <- DEres[DEres$padj <= 0.05 & DEres$log2FoldChange > 0.0,]
upFrame <- upFrame[complete.cases(upFrame),]

universe<- DEres$systematic_name
gene <- upFrame$systematic_name

enrchdUP_GO<-enricher(gene, pAdjustMethod = "none", universe, pvalueCutoff = 0.05,TERM2GENE= BPdict)

upBP_GOs <- data.frame(enrchdUP_GO)

upBP_GOs

upBP_GOs$ID


revigo.names <- c("term_ID","description","freqInDbPercent","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0016192","vesicle-mediated transport",1.085,0.567,0.000,"vesicle-mediated transport"),
                     c("GO:0019344","cysteine biosynthetic process",0.131,0.502,0.000,"cysteine biosynthesis"),
                     c("GO:0006333","chromatin assembly or disassembly",0.120,0.236,0.311,"cysteine biosynthesis"),
                     c("GO:0042407","cristae formation",0.013,0.294,0.120,"cysteine biosynthesis"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

treemap(
  stuff,
  index = c("representative","description"),
  vSize = "uniqueness",
  type = "categorical",
  vColor = "representative",
  palette = wes_palette("BottleRocket2"),
  title = "REVIGO Gene Ontology treemap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",     # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

UpEnsemble <- upFrame$ENSEMBL[!is.na(upFrame$ENSEMBL)]
write.csv(UpEnsemble, "upEnsembleIDs.csv", row.names = FALSE, col.names = FALSE, quote = FALSE )