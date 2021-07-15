

########################

library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

#path 1
path <- "/change/to/directory/for/16S"  ## CHANGE ME to the directory containing the fastq files.

fnFs <- sort(list.files(path, pattern = "_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq.gz", full.names = TRUE))

fns <- sort(list.files(path, pattern=".fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fns), ".fastq",1), `[`, 1)

#inspect quality files throughout
plotQualityProfile(fnFs[1:12])
plotQualityProfile(fnRs[1:12])


#create directory for filtered reads
dir.create(paste(path, "filtered_use_qual_scores", sep = '/'))
filt_path <- file.path(path, "filtered_use_qual_scores") 
filtFs <- file.path(filt_path, paste0(sample.names[grepl("R1", sample.names)], "_filtFs.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names[grepl("R2", sample.names)], "_filtRs.fastq.gz"))


#filter files according to parameters

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(10, 10), truncLen=c(230, 230),
              maxN=0, maxEE=c(2,2), truncQ=2, multithread = T,
              compress=TRUE)

#learn errors
set.seed(101)
errF <- learnErrors(filtFs, multithread=FALSE)
set.seed(101)
errR <- learnErrors(filtRs, multithread=FALSE)

#plot errors
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

#dereplicate
derepFs <- derepFastq(filtFs)
sam.names <- sapply(strsplit(basename(filtFs), "_R"), `[`, 1)
names(derepFs) <- sam.names


derepRs <- derepFastq(filtRs)
names(derepRs) <- sam.names

#run dada2
set.seed(101)
ddF <- dada(derepFs, err=errF, multithread = TRUE)

set.seed(101)
ddR <- dada(derepRs, err=errR, multithread = TRUE)

#merger reads
mergers <- mergePairs(ddF, derepFs, ddR, derepRs)
seqtab <- makeSequenceTable(mergers)

#remove chimeras
set.seed(101)
seqtab_final <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE)

#check read counts
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(ddF, getN), sapply(ddR, getN), sapply(mergers, getN), rowSums(seqtab_final))
colnames(track) <- c("in", "filtered","denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sam.names
track

write.table(track, "track_16S_HSCT_reads.txt")



#assign taxonomy
set.seed(101)
taxa <- assignTaxonomy(seqtab_final, "/change/to/directory/for/silva_138/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
set.seed(101)
taxa <- addSpecies(taxa, "/change/to/directory/for/silva_138/silva_species_assignment_v138.1.fa.gz")








#import mapping file - change path and metadata file name to match - please contact florent.malard@inserm.fr for this
mapping_file <- read.table("/change/to/directory/for/mapping_file.txt", header=T)

#remove "16S_" from sample names
rownames(seqtab_final) <- gsub("16S_", "", rownames(seqtab_final))
                               
#remove samples post HSCT
mapping_file <- mapping_file[mapping_file$time_point == "pre",]



seqtab_final <- seqtab_final[rownames(seqtab_final) %in% rownames(mapping_file),]
mapping_file <- mapping_file[rownames(mapping_file) %in% rownames(seqtab_final),]

all(rownames(seqtab_final) %in% rownames(mapping_file))
all(rownames(mapping_file) %in% rownames(seqtab_final))

mapping_file <- mapping_file[match(rownames(seqtab_final), rownames(mapping_file)),]
all(rownames(mapping_file) == rownames(seqtab_final))



###Make phyloseq object
library(phyloseq)
tax.mat <- tax_table(taxa)

ps <- phyloseq(tax_table(tax.mat),
               sample_data(mapping_file),
               otu_table(seqtab_final, taxa_are_rows = FALSE))

sample_data(ps)$Anon <- paste0("S_", sample_data(ps)$Anon)
ps.bact <- ps

#uncomment to save
#save(ps.bact, file = "ps.bact.dada2.RData")






