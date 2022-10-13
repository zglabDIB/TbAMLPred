print('Load Libraries ...')

library(DESeq2)


print('Reading Associate Files ...')

countdata1 <- read.table("ASSOCIATED_FILES/cd34_PB-AML_matrix.txt", header=TRUE, sep =  "\t", row.names= 1)
countdata2 <- read.table("ASSOCIATED_FILES/cd34_Hl60_matrix.txt", header=TRUE, sep =  "\t", row.names= 1)

print('Reading Input Files ...')
countdata3 <- read.table("INPUT_FILES/Input_Count.txt", header=TRUE, sep =  "\t", row.names= 1)
countdata  <- cbind(countdata1,countdata2,countdata3)
countdata <- as.matrix(countdata)

condition <- factor(c(rep("SAMP", dim(countdata)[2])))


print('Normalizing Data Files ...')

coldata <- data.frame(row.names=colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~1)
dds <- estimateSizeFactors(dds)

resdata1 <- as.data.frame(counts(dds, normalized=TRUE), by="row.names", sort=FALSE)
resdata1=resdata1[c(36:dim(countdata)[2])]

resdata <- cbind(Gene = rownames(resdata1), resdata1)
rownames(resdata) <- NULL
head(resdata)
colnames(resdata)

write.table(resdata, file="INPUT_FILES/Input_Normalized_Count.txt", sep =  "\t", row.names = FALSE, quote =F)






