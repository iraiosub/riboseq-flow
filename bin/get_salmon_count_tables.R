
# also intall readr
library(tximport) # tximeta?
library(DESeq2)



files <- file.path(dir, "salmon", samples$run, "quant.sf")
names(files) <- paste0("sample", 1:6)


txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
head(txi.salmon$counts)

sampleTable <- data.frame(condition = factor(rep(c("A", "B"), each = 3)))
rownames(sampleTable) <- colnames(txi$counts)

dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)


# files <- file.path(dir, "salmon", samples$run, "quant.sf.gz")
#  Typically, abundance is provided by the quantification tools as TPM (transcripts-per-million), 
# while the counts are estimated counts (possibly fractional), and the "length" matrix contains the effective gene lengths. 
# The "length" matrix can be used to generate an offset matrix for downstream gene-level differential analysis of count matrices, as shown below.


# DESeq2, is to use the gene-level estimated counts from the quantification tools, and additionally to use the transcript-level abundance estimates to calculate a gene-level
# offset that corrects for changes to the average transcript length across samples.