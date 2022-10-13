# TbAMLPred
TbAMLPred classifies the AML samples from the normal ones. This model has been trained using normal CD34+ and AML samples extracted from bone marrow tissues. Our system is modeled in a specialized manner such that TbAMLPred is able to perform prediction task using both peripheral blood and bone marrow also. It supports both microarray and RNA-Seq data. <br />
## Pre-requisite ##
The input provided by user should be Microarray or RNA-Seq data corresponding to CD34 cells isolated from peripheral blood/bone marrow sample of the patient. 
### Input Data format ###
* __Microarray:__ Our prediction algorithm can process microarray data belonging to two Affymetrix platforms;<br/>
  *	GPL570 [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array (54675 probes)
  *	GPL96	[HG-U133A] Affymetrix Human Genome U133A Array (22283 probes)
* __RNA-Seq Data:__ The RNA-seq data needs to be pre-processed to generate raw count values which shall be provided as input to the prediction software. Following are the steps needed for the pre-processing: <br/>
Perform the below steps to generate __ALL_Input_Count.txt__ raw count file  
1.	Quality checking and adapter trimming with ___Cutadapt___  <br/>
___cutadapt -q 35 -m 15 -a Read_1_adapter_sequence -A Read_2_adapter_sequence -o File_1_trimmed.fastq -p File_2_trimmed.fastqFile_1.fastqFile2_2.fastq___
2.	Generate alignment file with spliced aligner like ___Hisat2___ <br/>
(Please use __hg38 GRCh38.p13.genome.fa__ and __gencode.annotation.gtf__ files from http://gencodegenes.org. Follow hisat protocol to build index files with genome fasta) <br/>
___hisat2 -p 8 -q -x /directory_to_hg38_genome_fasta_index_files/GRCh38.p13.genome -1 File_1_trimmed.fastq -2 File_2_trimmed.fastq -S File.sam___<br/>
___samtools view -bS -@ 8File.sam>File_unsorted.bam___<br/>
___samtools sort -@ 8File_unsorted.bam-o File.bam___<br/>
_(Please note that command may vary from version to version, so please check accordingly. Other spliced aligner like STAR can also be used)_<br/>
3.	Generate raw count file with ___FeatureCounts___<br/>
___featureCounts -p -a gencode.annotation.gtf-o File_Count.txt -T 32 -t exon -M -g gene_id File.bam___<br/>
___grep -v "#" File_Count.txt | cut -f1,7| sed "s/.bam//g" >ALL_Input_Count.txt___<br/><br/>
__Software:__ For executing the codes for Microarray or RNA-Seq data, both R and python are required to be installed in local machines.<br/><br/>
__R programming environment:__
    *	___Bioconductor___ <br/>
    *	___affy___ <br/>
    *	___hgu133plus2.db___ <br/>
    *	___frma___ <br/>
    *	___DESeq2___ <br/><br/><br/>
__Python programming Environment:__
    *	___Anaconda software packages for python v3.8___ <br/>
    *	___pandas___ <br/>
    *	___numpy___ <br/>
    *	___Scikit-Learn version 0.23.2___ <br/>
    *	___matplotlib___ <br/>
    *	___pickle___ <br/>
    *	___umap-learn version 0.5.1___ <br/>

