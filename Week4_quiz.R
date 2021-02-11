#Q1:What fraction of reads in yeastRNASeq FASTQ file has an A nucleotide in the 5th base of the read?

library(yeastRNASeq)
library(ShortRead)
fastqFilePath <- system.file("reads", "wt_1_f.fastq.gz", package = "yeastRNASeq")  
reads1 <- readFastq(fastqFilePath)

sreads1 <- sread(reads1)

base5 <- as.character(subseq(sreads1, start = 5, width = 1))    #alternately use alphabetByCycle(sreads1)[ , 5]

fraction <- sum(base5 == "A")/length(base5)
  
fraction       #output: 0.363841


#Q2:What is the average numeric quality value of the 5th base of these reads?

qualsints <- as(quality(reads1), "matrix")

dim(qualsints)          

mean(qualsints[, 5])         #Output: 28.93346 ;5th col is 5th cycle/base qual score


#Q3:We will focus on the interval from 800,000 to 801,000 on yeast chromosome 13.
#In this interval, how many reads are duplicated by position? In  leeBamViews experiment data package

library(leeBamViews)

bamFilePath <- system.file("bam", "isowt5_13e.bam", package="leeBamViews")

bamfile <- BamFile(bamFilePath)

gr=GRanges("Scchr13:800000-801000")

params <- ScanBamParam(which = gr, what = scanBamWhat())

aln <- scanBam(bamfile, param = params)

tl <- table(aln$`Scchr13:800000-801000`$pos)

sum(tl[tl>1])                  #$pos>1, output:129


#Q4: What is the average number of reads across the 8 samples falling in this interval?

bpaths <- list.files(system.file("bam", package="leeBamViews"), pattern = "bam$", full=TRUE)

gr4=GRanges("Scchr13:807762-808068")    #novel transcribed regions in yeast Scchr13:807762-808068.

bviews <- BamViews(bamPaths = bpaths, bamRanges = gr4)   #views object for bam files in gr4 ranges

aln4 <- scanBam(bviews)

q4 <- sapply(aln4, function(x) {length(x$`Scchr13:807762-808068`$pos) })   #number of reads in the region for each file (each list in aln4)

mean(q4)         #output: 90.25


#Q5: 

library(oligo)
library(GEOquery)

getGEOSuppFiles("GSE38792")

list.files("GSE38792")

untar("GSE38792/GSE38792_RAW.tar", exdir = "GSE38792/CEL")

list.files("GSE38792/CEL")

celfiles <- list.files("GSE38792/CEL", full = TRUE)

rawData <- read.celfiles(celfiles)

#clean up phenotype info for rawData:

filename <- sampleNames(rawData)
pData(rawData)$filename <- filename
sampleNames <- sub(".*_", "", filename)
sampleNames <- sub(".CEL.gz$", "", sampleNames)
sampleNames(rawData) <- sampleNames
pData(rawData)$group <- ifelse(grepl("^OSA", sampleNames(rawData)), "OSA", "Control")
pData(rawData)

#data normalization with rma()
normData <- rma(rawData)
normData                        #an ExpressionSet
boxplot(normData)

mean(exprs(normData)["8149273", normData$group == "Control"])    #output: 7.02183, col 1-8 are control group


#Q6:Use the limma package to fit a two group comparison between the control group and the OSA group, 
#and borrow strength across the genes using eBayes(). Include all 18 samples in the model fit.
#Question: What is the absolute value of the log foldchange (logFC) of the gene with the lowest P.value

library(limma)

normData$group <- factor(normData$group)
design <- model.matrix(~ normData$group)
fit <- lmFit(normData, design)
fit <- eBayes(fit)
topTable(fit)              #output: 0.7126484; abs value of logFC for gene with lowest p.value is the 1st row of top table


#Q7: differentially expressed genes with adj P value < 0.05
top <- topTable(fit)
top[top$adj.P.Val < 0.05, ]       #output : 0 genes 

#Q8: 

library(minfiData)
library(GEOquery)
minfiData::RGsetEx

head(pData(RGsetEx))

preproc8 <- preprocessFunnorm(RGsetEx)

opensea <- preproc8[getIslandStatus(preproc8) == "OpenSea"]

beta_opensea <- getBeta(opensea)

meanbetasample <- colMeans(beta_opensea)

vector <- as.vector(meanbetasample)
normal_mean <- sum(vector[c(1,2,5)])/3
cancer_mean <- sum(vector[c(3,4,6)])/3
q8_ans <- normal_mean - cancer_mean          #output: .08863657; incorrect ans

##lauras_means <- c(0.69666, 0.70895, 0.60899, 0.64245, 0.70595, 0.60626)
#normal_mean <- sum(lauras_means[c(1,2,5)])/3
#cancer_mean <- sum(lauras_means[c(3,4,6)])/3
#q8_ans <- normal_mean - cancer_mean          #output: .08462; correct ans using Laura's means, something has changed in packages, logic is correct


#Q9:How many of these DNase hypersensitive sites contain one or more CpGs on the 450k array?
  
library(AnnotationHub)

ah <- AnnotationHub()
ah_hs <- subset(ah, species == "Homo sapiens")
query(ah_hs, c("Caco2", "DNase"))
DNase_hypercaco2 <- ah_hs[["AH22442"]]

rangesCpg <- granges(preproc8)
q9_ans <- subsetByOverlaps(DNaseHyperCaco2, rangesCpg)     #ans : 40151 ranges in q9_ans

  
#Q10: 

library(DESeq2)  
library(zebrafishRNASeq)
data("zfGenes")
head(zfGenes)

countData = as.matrix(zfGenes)
countData = subset(countData,substr(rownames(countData),1,4) != "ERCC")
condition = factor(substr(colnames(countData), 1, 3))
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)
dds <- DESeq(dds)                                   #tried betaPrior = TRUE but still got incorrect ans
res <- results(dds, independentFiltering = FALSE)
   
res <- res[order(res$padj),]
length(which(res$padj <= 0.05))     #output: 91, correct ans is 87. 
