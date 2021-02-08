##This is week 2 assignment code for Bioconductor course 

##Q1
#load BSGenome package as genomes have there own class (because they are big)
library(BSgenome)   

allgenomes <- available.genomes()   #lists which genomes are currently available

grep("Hsapiens", allgenomes, value = TRUE) #grep for human

BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")   #install the needed genome if not already installed, takes time

library("BSgenome.Hsapiens.UCSC.hg19")

#GC content of chr22 in the hg19 build of human genome
gc <- letterFrequency(Hsapiens$chr22, "GC")   #G|C freq
totalwithoutN <- letterFrequency(Hsapiens$chr22, "ATGC")   #excluding Ns
q1_ans <- gc/totalwithoutN
q1_ans    #output is 0.4798807

##Q2: what's the mean GC content of H3K27me3 "narrowPeak" regions from Epigenomics
##Roadmap from the H1 stem cell line on chr22
library(AnnotationHub)
ah <- AnnotationHub()
ah_hs <- subset(h, species == "Homo sapiens")
query(ah_hs, c("E003", "H3K27me3"))
histonemod <- ah_hs[["AH29892"]]  #picked the narrowpeak record from the above query results

#keep peaks for chr22 only as asked in question2
histone_22 <- keepSeqlevels(histonemod, "chr22", pruning.mode = "coarse")
histoneViews <- Views(Hsapiens, histone_22)  #gives Views in Hsapiens for the histonemods in chr22(histone_22) 
gcHistone <- letterFrequency(histoneViews, "GC", as.prob = TRUE)  #gc frequency in narrowpeaks in chr22
q2_ans <- sum(gcHistone)/length(gcHistone)    #mean GC content, output: 0.528866


##Q3:What is the correlation between GC content and “signalValue” of these regions (on chr22)
q3 <- cor(gcHistone, histone_22$signalValue)    #correlation
q3    #output: .004467924

##Q4:what is the correlation between the “signalValue” of the “narrowPeak” 
##regions and the average “fc.signal” across the same regions? (doing for chr22)

ahub.bw <- subset(ah, rdataclass == "BigWigFile" & species == "Homo sapiens")  #look for bigwig files
q <- query(ahub.bw, c("E003", "H3K27me3") )   #query for H1 cell line, E003, H3K27me3 histone mod

bw <- q[[1]]       #took a looong time to download, so check before running this command if this is needed again, 
                    #this will download the fc.signal bigwig file, which was 1st query result, took forever 

chr22 <- GRanges(seqnames = "chr22", ranges = IRanges(start = 1, end = seqlengths(histone_22)))  #GRanges object of length chr22

out.chr22 <- import(bw, which = chr22, as = "Rle")    #import from bw file, only the chr22 GRanges which we created above, do it as Rle

fcsignal_mean <- aggregate(out.chr22$chr22, histone_22, FUN = mean)

signalvalue_histone22 <- histone_22$signalValue

q4_ans <- cor(fcsignal_mean, signalvalue_histone22)   #output: 0.9149614

##Q5: How many bases on chr22 have an fc.signal greater than or equal to 1? Use objects from Q4
fc_above1 <- slice(out.chr22$chr22, 1)    #gives output as views, for values of fc.signal rle chr22 greater than or equal to 1
q5_ans <- sum(width(fc_above1))           #sum of number of bases on chr22 from above views object; output: 10914671

##Q6:Identify the regions of the genome where the fc.signal in E003 is 0.5 or lower and the signal in E055 is 2 or higher.

qu6 <- query(ahub.bw, c("E055", "H3K27me3"))    #query for E055, histone modf, in bigwig files
bw_e055 <- qu6[[1]]                   #download the required fc.signal bigwig for E055, took a long time to download

fcsig_chr22 <- import(bw_e055, which = chr22, as = "Rle")   #import from bw_e005 file, only the chr22 GRanges which we created above, do it as Rle

fc_upper2 <- slice(fcsig_chr22$chr22, lower = 2)     #slice above obtained rle with lower bound at 2 (includes 2 default)

fc_lower0.5 <- slice(out.chr22$chr22, upper = 0.5)  #slice the fcsignal from chr22 of e003 rle with upper bound at 0.5 (includes 2 default)

ir_e003 <- as(fc_lower0.5, "IRanges")       #need to convert Views to IRanges (or Granges) before for using intersect()
ir_e055 <- as(fc_upper2, "IRanges")

q6 <- intersect(ir_e003, ir_e055)           #Iranges object, number of ranges is the number of regions satisfying condition in ques
sum(width(q6))                  #number of bases; output:1869937

##Q7: What is the average observed-to-expected ratio of CpG dinucleotides for CpG Islands on chromosome 22?
cpg <- query(ah_hs, c("CpG", "hg19"))
cpg <- cpg[["AH5086"]]        #CpG islands GRanges for hg19 human

library("BSgenome.Hsapiens.UCSC.hg19")

cpg_22 <- keepSeqlevels(cpg, "chr22", pruning.mode = "coarse")

cpgViews <- Views(Hsapiens, cpg_22)

dnfreq <- dinucleotideFrequency(cpgViews)
cgfreq_obs <- dnfreq[, "CG"]

alphafreq <- alphabetFrequency(cpgViews)

expectedfreq <- (alphafreq[, "C"] * alphafreq[, "G"]) / width(cpgViews)

obstoexpec <- cgfreq_obs/expectedfreq

mean(obstoexpec)      #output: 0.8340929

##Q8: How many TATA boxes are there on chr 22 of build hg19 of the human genome? Search both reverse and forward strand i.e reverse complement of the pattern
#Use DNAString(pattern), in order to use reverseComplement(); output for below expression: 27263
countPattern("TATAAA", Hsapiens$chr22) + countPattern(reverseComplement(DNAString("TATAAA")), Hsapiens$chr22)


##Q9:How many promoters of transcripts on chromosome 22 containing 
##a coding sequence, contains a TATA box on the same strand as the transcript?  
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

#get codons by tx_id from txdb
cds <- cdsBy(txdb, by="tx")

cds_22 <- keepSeqlevels(cds, "chr22", pruning.mode = "coarse")    #keep chr22 only, GRangesList object

namescds <- names(cds_22)     #get names of the GRanges lists, needed for next operation

trans_cds <- subset(transcripts(txdb), tx_id %in% namescds)     #get transcripts with tx_id in namescds from above
proms <- promoters(trans_cds, upstream = 900, downstream = 100)      #promoters for transcripts on chr22, containing a cds

tata_hs <- vmatchPattern("TATAAA", Hsapiens)    #gives a GRanges object

tata_hs22 <- keepSeqlevels(tata_hs, "chr22", pruning.mode = "coarse")

x <- subsetByOverlaps(proms, tata_hs22)     #gives GRanges object #output: 193 (countOverlaps() gives different 
                                            #output requiring sum(output>0) to get right answer )

##Q10:How many bases on chr22 are part of more than one promoter of a coding sequence?

library("BSgenome.Hsapiens.UCSC.hg19")

cover <- coverage(proms)      #coverage() gives rle object for all chr, each chr rle has counts of the 
                              #number of ranges that cover each base   

sum(width(slice(cover$chr22, 2)))     #output: 306920. width is the bases in a range, slice at 2 will give bases in  
                                      # a range with more than or equal 2 promoter




