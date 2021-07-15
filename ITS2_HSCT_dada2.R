getwd()
list.files()
setwd("/change/to/directory/for/ITS2")

library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

#path 1
path <- "/change/to/directory/for/ITS2"  ## CHANGE to the directory containing the fastq files.

fnFs <- sort(list.files(path, pattern = "_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq.gz", full.names = TRUE))

#ITS2 primers
FWD <- "GTGARTCATCGAATCTTT"  
REV <- "GATATGCTTAAGTTCAGCGGGT"  

#different orients/compliments
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)


#filter to remove N's
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) 
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

# remove hsct_a7 and hsct_e4 as no reads post processing for HSCT
fnFs <- fnFs[!(grepl("hsct_a7_", fnFs) | grepl("hsct_e4_", fnFs))]
fnRs <- fnRs[!(grepl("hsct_a7_", fnRs) | grepl("hsct_e4_", fnRs))]
fnFs.filtN <- fnFs.filtN[!(grepl("hsct_a7_", fnFs.filtN) | grepl("hsct_e4_", fnFs.filtN))]
fnRs.filtN <- fnRs.filtN[!(grepl("hsct_a7_", fnRs.filtN) | grepl("hsct_e4_", fnRs.filtN))]



primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
for(i in 1:length(fnFs.filtN)){
  if(gsub("/change/to/directory/for/ITS2/filtN/","", fnFs.filtN[[i]]) %in% list.files("/change/to/directory/for/ITS2/filtN/")){
    print(rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[i]]), 
                FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[i]]), 
                REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[i]]), 
                REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[i]]))) 
  }
}

#Manually go through quality plots
plotQualityProfile(fnFs[1:12])
plotQualityProfile(fnRs[1:12])

#import cutadapt
cutadapt <- "/change/todirectory/for/cutadapt"
system2(cutadapt, args = "--version")

#set path for output
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

#set up FWD/REV + complements
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

#check to ensure removed
for(i in 1:length(fnFs.cut)){
  if(gsub("/change/to/directory/for/ITS2/cutadapt/","", fnFs.cut[[i]]) %in% list.files("/change/to/directory/for/ITS2/cutadapt/")){
    print(rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[i]]), 
                FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[i]]), 
                REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[i]]), 
                REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[i]]))) 
  }
}

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2.fastq.gz", full.names = TRUE))


# Extract sample names, ENSURE filenames have correct format with delimiters:
get.sample.name <- function(fname) strsplit(basename(fname), "_R")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)




#set path for filtered
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))

#perform trimming
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)


head(out)

plotQualityProfile(filtFs[1:12])
plotQualityProfile(filtRs[1:12])


#learn errors
set.seed(101)
errF <- learnErrors(filtFs, multithread = TRUE)
set.seed(101)
errR <- learnErrors(filtRs, multithread = TRUE)

#check errors
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

#derep
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#run dada2
set.seed(101)
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
set.seed(101)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

#merge
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)

#remove chimeras
set.seed(101)
seqtab_final <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE)

#check read attrition
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab_final))
colnames(track) <- c("in", "filt", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track
write.table(track, "track_ITS2_HSCT_reads.txt")

#### Assign taxonomy with unite

#note: even with 'set.seed' there can be some very minor variation between runs of assignTaxonomy
unite.ref <- "/change/to/directory/for/sh_general_release_s_10.05.2021/sh_general_release_dynamic_s_10.05.2021.fasta"  # CHANGE ME to location on your machine
set.seed(101)
taxa <- assignTaxonomy(seqtab_final, unite.ref, multithread = TRUE, tryRC = TRUE)



#import mapping file and select data for this project

#mapping file - change path and metadata file name to match - please contact florent.malard@inserm.fr for this
mapping_file <- read.table("/change/to/directory/for/mapping_file.txt", header=T)

#remove "ITS2_" from sample names
rownames(seqtab_final) <- gsub("ITS2_", "", rownames(seqtab_final))

#remove post HSCT samples
mapping_file <- mapping_file[mapping_file$time_point == "pre",]

#subset sequence table and mapping file to include only samples present in both
seqtab_final <- seqtab_final[rownames(seqtab_final) %in% rownames(mapping_file),]
mapping_file <- mapping_file[rownames(mapping_file) %in% rownames(seqtab_final),]

#ensure all present in both
all(rownames(seqtab_final) %in% rownames(mapping_file))
all(rownames(mapping_file) %in% rownames(seqtab_final))

#match order and ensure correct
mapping_file <- mapping_file[match(rownames(seqtab_final), rownames(mapping_file)),]
all(rownames(mapping_file) == rownames(seqtab_final))



###Make phyloseq object
library(phyloseq)
tax.mat <- tax_table(taxa)

ps <- phyloseq(tax_table(tax.mat),
               sample_data(mapping_file),
               otu_table(seqtab_final, taxa_are_rows = FALSE))

sample_data(ps)$Anon <- paste0("S_", sample_data(ps)$Anon)
ps.fung <- ps

#uncomment to save
#save(ps.fung, file = "ps.fung.dada2.RData")
