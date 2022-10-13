print('Load Libraries ...')

library(affy)
library(hgu133plus2.db)
library(frma)


pd=read.AnnotatedDataFrame("INPUT_FILES/phenoData.txt",header=TRUE,row.names=1)

print('Read Affymetrix Array Files ...')

raw_data = ReadAffy(filenames=rownames(pData(pd)))


print('Normalize Data using fRMA ...')


eset_data=frma(raw_data)
eset_data

write.exprs(eset_data,file='INPUT_FILES/Input_INTENSITIES.txt',sep='\t')


