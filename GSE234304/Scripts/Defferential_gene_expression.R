#loading required packages
library(tidyverse)
library(DESeq2)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(rio)

#Importing count data and metadata
GSE223119 <- import("Data/GSE223119_raw_counts_GRCh38.p13_NCBI.tsv")
metadata <- import("Data/SraRunTable.csv")

#formatting count data and metadata for differential expression analysis
count_data <- GSE223119 |> column_to_rownames("GeneID") |> as.matrix()

coldata <- metadata[,c("Sample Name","tissue_pair")] |> 
  column_to_rownames("Sample Name") |> 
  dplyr::rename("condition" = "tissue_pair")

coldata$condition <- coldata$condition |>
  as.character() |>
  (\(x) gsub("tissue pair: para-carcinoma tissue.*", "paracancerous", x))() |>
  (\(x) gsub("tissue pair: colorectal carcinoma tissue.*", "cancerous", x))() |>
  as.factor()

# Match coldata with count matrix
coldata <- coldata |> 
  filter(row.names(coldata) %in% colnames(count_data)) |> 
  arrange(match(row.names(coldata),colnames(count_data)))

# Create DESEq2 data set object
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = coldata,
  design = ~condition
)

dds <- dds[rowSums(counts(dds) >= 20) >= 10, ]
dds$condition <- relevel(dds$condition, ref="paracancerous")
dds <- DESeq(dds)

result <- results(dds)
result$GeneID <- rownames(result)

write.csv(result, "outputs/DESeq2_result.csv",
          row.names = FALSE)

#Annotating Genes
res <- import("outputs/DESeq2_result.csv")
res$Gene_Symbol <- mapIds(org.Hs.eg.db,
                          keys = as.character(res$GeneID),
                          column = "SYMBOL",
                          keytype = "ENTREZID",
                          multiVals = "first")

res$Gene_Description <- mapIds(org.Hs.eg.db,
                               keys = as.character(res$GeneID),
                               column = "GENENAME",
                               keytype = "ENTREZID",
                               multiVals = "first")

export(res, "outputs/Annotated_result.csv")


#Finding key genes
EnhancedVolcano(key_genes,
                lab = key_genes$Gene_Symbol,
                x = 'log2FoldChange',
                y = 'padj')

key_genes <- res |> 
  filter(padj < 0.05 & abs(log2FoldChange) >= 2) |> 
  drop_na(Gene_Symbol)

up <- key_genes |> 
  arrange(desc(log2FoldChange)) |> 
  slice_head(n = 30)

down <- key_genes |> 
  arrange(log2FoldChange) |> 
  slice_head(n = 30)

top <- rbind(up, down)

export(top, "outputs/top_genes.csv")
