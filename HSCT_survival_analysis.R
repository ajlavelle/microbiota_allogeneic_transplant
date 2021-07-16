#functions for Survival analysis

#assign taxa names to otu table from taxa slot in phyloseq object
assign_taxa_names <- function(x, ps.x=NULL, tax_level="Phylum"){
  #check that ps in correct orientation
  if(taxa_are_rows(ps.x)){
    stop("Taxa are rows - please re-orientate phyloseq object.")
  }
  if(!(all(colnames(x) == rownames(tax_table(ps.x))))){
    stop("Error! Taxa names do not match.")
  }
  #assign taxa name at lowest available level
  tax.tab <- tax_table(ps.x)
  s.vec <- tax.tab[,7]
  temp.vec <- is.na(s.vec) | grepl("unidentified|unclassified|ambiguous|metagenome|uncultured", ignore.case = TRUE, s.vec)
  if(TRUE %in% temp.vec){
    s.vec[temp.vec] <- tax.tab[,6][temp.vec]
  }
  temp.vec <- is.na(s.vec) | grepl("unidentified|unclassified|ambiguous|metagenome|uncultured", ignore.case = TRUE, s.vec)
  if(TRUE %in% temp.vec){
    s.vec[temp.vec] <- tax.tab[,5][temp.vec]
  }
  temp.vec <- is.na(s.vec) | grepl("unidentified|unclassified|ambiguous|metagenome|uncultured", ignore.case = TRUE, s.vec)
  if(TRUE %in% temp.vec){
    s.vec[temp.vec] <- tax.tab[,4][temp.vec]
  }
  temp.vec <- is.na(s.vec) | grepl("unidentified|unclassified|ambiguous|metagenome|uncultured", ignore.case = TRUE, s.vec)
  if(TRUE %in% temp.vec){
    s.vec[temp.vec] <- tax.tab[,3][temp.vec]
  }
  temp.vec <- is.na(s.vec) | grepl("unidentified|unclassified|ambiguous|metagenome|uncultured", ignore.case = TRUE, s.vec)
  if(TRUE %in% temp.vec){
    s.vec[temp.vec] <- tax.tab[,2][temp.vec]
  }
  temp.vec <- is.na(s.vec) | grepl("unidentified|unclassified|ambiguous|metagenome|uncultured", ignore.case = TRUE, s.vec)
  if(TRUE %in% temp.vec){
    s.vec[temp.vec] <- tax.tab[,1][temp.vec]
  }
  #add ASV and unique number to prevent duplicate column names
  s.vec <- paste(s.vec, paste("ASV", seq(1,length(s.vec),1), sep = "."), sep = ".")
  s.vec <- gsub("/", "_", s.vec)
  s.vec <- gsub("-", "_", s.vec)
  #assign names to otu table
  colnames(x) <- s.vec
  otus <- as.data.frame(as.matrix(x))
  #extract matching phylum-level assignments for plotting
  phylum.vec <- tax.tab[,tax_level]
  rownames(phylum.vec) <- NULL
  phylum.vec <- as.vector(phylum.vec)
  #return as list
  return.list <- vector("list",2)
  names(return.list) <- c("otus", "phyla")
  return.list$otus <- otus
  return.list$phyla <- phylum.vec
  return(return.list)
}


#convert a normalised ('closed') count table to binary based on median
make_binary <- function(x){
  if(!(is.data.frame(x))){
    x <- is.data.frame(x)
  }
  x.bin <- apply(x, 2, function(y) y > median(y))
  x.bin[x.bin == TRUE] <- "High"
  x.bin[x.bin == FALSE] <- "Low"
  return(x.bin)
}

#wrap other functions to output binary count table with relevant clinical variables and binary diversity index (Shannon)
binary_taxa_surv <- function(ps.x, agglom=NULL, read_depth=NULL, prev=NULL, surv_cols=NULL, tax_level="Phylum"){
  #keep count data
  otus.orig <- otu_table(ps.x)
  surv_data.orig <- sample_data(ps.x)
  #agglomerate
  if(!(is.null(agglom))){
    ps.x.agglom <- tax_glom(ps.x, taxrank = agglom, NArm = FALSE)
  } else {
    ps.x.agglom <- ps.x
  }
  #extract otu table
  otus <- otu_table(ps.x.agglom)
  
  #assign taxa names
  named.df <- assign_taxa_names(otus, ps.x = ps.x.agglom, tax_level = tax_level)
  otus <- named.df$otus
  phyla.list <- named.df$phyla
  #make otu table df
  otus <- as.data.frame(as.matrix(otus))
  
  #keep samples with reads > read_depth
  otus.filt <- otus[rowSums(otus) > read_depth,]
  
  #normalise data
  otus.filt <- sweep(otus.filt, 1, rowSums(otus.filt), "/")
  
  #filter otus  with > prev of zeros
  filt.vec <- (colSums(otus.filt == 0)/nrow(otus.filt)) <= prev
  
  #filter otu table and phylum vector
  otus.filt <- otus.filt[,filt.vec]
  phyla.filt <- phyla.list[filt.vec]
  
  #convert otu table to binary based on median count
  otus.bin <- make_binary(otus.filt)
  #remove samples from original count table that did not pass read threshold
  otus.orig <- otus.orig[rownames(otus.orig) %in% rownames(otus.bin),]
  #calculate Shannon diversity
  shannon <- vegan::diversity(otus.orig, index = "shannon")
  #convert Shannon diversity to binary
  shannon.bin <- unlist(lapply(shannon, function (x) if(x > median(shannon))("High") else {"Low"}))
  #remove samples from metadata that did not pass read count threshold
  surv_data <- surv_data.orig[rownames(surv_data.orig) %in% rownames(otus.bin),]
  #ensure sample names in all tables match, then combine select metadata columns, binary otu table and binary Shannon diversity
  if(all(rownames(otus.bin) == rownames(surv_data)) & all(rownames(otus.bin) == rownames(otus.orig))){
    otus.bin.surv <- data.frame(surv_data[,surv_cols], otus.bin, shannon.bin)
  } else {
    print("Error! Mismatch in sample names.")
  }
  #return data and phylum-level assignments as a list
  return.list <- vector("list",2)
  names(return.list) <- c("otus", "phyla")
  return.list$otus <- otus.bin.surv
  return.list$phyla <- phyla.filt
  return(return.list)
}







### Import libraries ###
library(phyloseq)
library(vegan)
library(survival)
library(survminer)



#import data
load("ps.bact.dada2.RData")
load("ps.fung.dada2.RData")

#filter ASVs not assigned to a Phylum
ps.bact <- subset_taxa(ps.bact, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
ps.fung <- subset_taxa(ps.fung, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))


#remove duplicates
ps.f.start <- subset_samples(ps.fung, !(rownames(sample_data(ps.fung)) %in% c("hsct_g2", "hsct_k6", "hsct_p6", "hsct_t5")))
ps.b.start <- subset_samples(ps.bact, !(rownames(sample_data(ps.bact)) %in% c("hsct_e4", "hsct_f8", "hsct_k6", "hsct_p6", "hsct_t5")))

#######################################################################

                  ######## Univariate ########

#######################################################################


### Fungal ###

#Fungal data, genus level
surv.f.list <- binary_taxa_surv(ps.f.start, agglom = "Genus", read_depth = 200, prev = 0.9, surv_cols = c("OS_months", "OS_1_death"))
surv.f.bin <- surv.f.list$otus

#perform univariate logrank
os.pvals <- c()
for(i in 3:ncol(surv.f.bin)){
  if(length(unique(surv.f.bin[,i])) > 1){
    x <- survdiff(Surv(surv.f.bin[,1], surv.f.bin[,2], type = "right") ~ surv.f.bin[,i], data = surv.f.bin)
    p_value <- 1 - pchisq(x$chisq, length(x$n) - 1)
    os.pvals <- c(os.pvals, p_value)
  }
}

#make table
OS.f.table <- data.frame(Variable = colnames(surv.f.bin)[3:ncol(surv.f.bin)], Phylum=c(surv.f.list$phyla,NA) , 
                         OS_p.value = as.numeric(os.pvals), 
                         OS_BH.adj = p.adjust(as.numeric(os.pvals), method = "fdr"))

#view table
OS.f.table[OS.f.table$OS_BH.adj < 0.2 & OS.f.table$OS_p.value < 0.05,]

### Figure 1B
#plot FDR-corrected p-values for logrank
OS <- OS.f.table[order(OS.f.table$OS_BH.adj),]
OS$Variable <- as.character(OS$Variable)
OS$Variable <- gsub(".ASV.\\d+","", OS$Variable)
OS$Variable <- gsub("shannon.bin","Shannon diversity", OS$Variable)
OS$Variable <- factor(OS$Variable, levels = rev(as.character(OS$Variable)))

ggplot(OS, aes(x=-1*log10(OS_BH.adj),y=Variable, fill=Phylum)) + geom_bar(stat="identity") + theme_classic() +
  geom_vline(xintercept = -1*log10(0.05), colour="red", lty=2) + 
  geom_vline(xintercept = -1*log10(0.1), colour="grey", lty=2) +
  geom_vline(xintercept = -1*log10(0.2), colour="grey") +
  xlab("-log10(FDR)") + ylab("Genus") + scale_fill_manual(values = c("darkorange","lightblue"), na.value="grey47") +
  ggtitle("Log-rank test overall survival - fungal genera") + xlim(c(0,1.9))



### Figure 1C - Candida genus
#plot KM
fit<- survfit(Surv(OS_months, OS_1_death, type = "right") ~ g__Candida.ASV.4, data = surv.f.bin)
ggsurvplot(fit, pval = TRUE, conf.int = TRUE, risk.table = TRUE, palette = c("red", "blue"))


#agglomerate at species level - set tax_level to genus
surv.f.list <- binary_taxa_surv(ps.f.start, agglom = "Species", read_depth = 200, prev = 1, surv_cols = c("OS_months", "OS_1_death"), tax_level = "Genus")
surv.f.bin <- surv.f.list$otus

#evaluate Candida species
surv.candida.bin <- surv.f.bin[,c(TRUE, TRUE, grepl("candida", ignore.case = T, surv.f.list$phyla), FALSE)]
#filter Candida species >90% low - as median split => 0
surv.candida.bin <- surv.candida.bin[,colSums(surv.candida.bin == "Low")/nrow(surv.candida.bin) <= 0.9]

#univariate test
os.pvals <- c()
for(i in 3:ncol(surv.candida.bin)){
  if(length(unique(surv.candida.bin[,i])) > 1){
    x <- survdiff(Surv(surv.candida.bin[,1], surv.candida.bin[,2], type = "right") ~ surv.candida.bin[,i], data = surv.candida.bin)
    p_value <- 1 - pchisq(x$chisq, length(x$n) - 1)
    os.pvals <- c(os.pvals, p_value)
  }
}

#make table
OS.f.table <- data.frame(Variable = colnames(surv.candida.bin)[3:ncol(surv.candida.bin)], 
                         OS_p.value = as.numeric(os.pvals), 
                         OS_BH.adj = p.adjust(as.numeric(os.pvals), method = "fdr"))
#view table
OS.f.table[OS.f.table$OS_BH.adj < 0.2 & OS.f.table$OS_p.value < 0.05,]

### Figure 1D - Candida albicans
#plot KM
fit<- survfit(Surv(OS_months, OS_1_death, type = "right") ~ s__albicans.ASV.4, data = surv.candida.bin)
ggsurvplot(fit, pval = TRUE, conf.int = TRUE, risk.table = TRUE, palette = c("red", "blue"))
















#Disease-free survival
surv.f.list <- binary_taxa_surv(ps.f.start, agglom = "Genus", read_depth = 200, prev = 0.9, surv_cols = c("DFS_months", "DFS_1_relapse_or_death"))
surv.f.bin <- surv.f.list$otus

#univariate
DFS.pvals <- c()
for(i in 3:ncol(surv.f.bin)){
  if(length(unique(surv.f.bin[,i])) > 1){
    x <- survdiff(Surv(surv.f.bin[,1], surv.f.bin[,2], type = "right") ~ surv.f.bin[,i], data = surv.f.bin)
    p_value <- 1 - pchisq(x$chisq, length(x$n) - 1)
    DFS.pvals <- c(DFS.pvals, p_value)
  }
}
#table
DFS.f.table <- data.frame(Variable = colnames(surv.f.bin)[3:ncol(surv.f.bin)], Phylum=c(surv.f.list$phyla,NA) , 
                         DFS_p.value = as.numeric(DFS.pvals), 
                         DFS_BH.adj = p.adjust(as.numeric(DFS.pvals), method = "fdr"))
#view
DFS.f.table[DFS.f.table$DFS_BH.adj < 0.2 & DFS.f.table$DFS_p.value < 0.05,]


### Figure 2A
#plot FDR-corrected p-values
DFS <- DFS.f.table[order(DFS.f.table$DFS_BH.adj),]
DFS$Variable <- as.character(DFS$Variable)
DFS$Variable <- gsub(".ASV.\\d+","", DFS$Variable)
DFS$Variable <- gsub("shannon.bin","Shannon diversity", DFS$Variable)
DFS$Variable <- factor(DFS$Variable, levels = rev(as.character(DFS$Variable)))

ggplot(DFS, aes(x=-1*log10(DFS_BH.adj),y=Variable, fill=Phylum)) + geom_bar(stat="identity") + theme_classic() +
  geom_vline(xintercept = -1*log10(0.05), colour="red", lty=2) + 
  geom_vline(xintercept = -1*log10(0.1), colour="grey", lty=2) +
  geom_vline(xintercept = -1*log10(0.2), colour="grey") +
  xlab("-log10(FDR)") + ylab("Genus") + scale_fill_manual(values = c("darkorange","lightblue"), na.value="grey47") +
  ggtitle("Log-rank test disease-free survival - fungal genera") + xlim(c(0,1.9))

### Figure 2B
#KM
surv.f.bin <- binary_taxa_surv(ps.f.start, agglom = "Species", read_depth = 200, prev = 0.9, surv_cols = c("DFS_months", "DFS_1_relapse_or_death"))
fit<- survfit(Surv(DFS_months, DFS_1_relapse_or_death, type = "right") ~ s__albicans.ASV.4, data = surv.f.bin$otus)
ggsurvplot(fit, pval = TRUE, conf.int = TRUE, risk.table = TRUE, palette = c("red", "blue"))



#GVHD-free, RFS
surv.f.list <- binary_taxa_surv(ps.f.start, agglom = "Genus", read_depth = 200, prev = 0.9, surv_cols = c("GVHD_free_relapse_free_time", "GVHD_free_relapse_free_survival"))
surv.f.bin <- surv.f.list$otus
#univariate
GRFS.pvals <- c()
for(i in 3:ncol(surv.f.bin)){
  if(length(unique(surv.f.bin[,i])) > 1){
    x <- survdiff(Surv(surv.f.bin[,1], surv.f.bin[,2], type = "right") ~ surv.f.bin[,i], data = surv.f.bin)
    p_value <- 1 - pchisq(x$chisq, length(x$n) - 1)
    GRFS.pvals <- c(GRFS.pvals, p_value)
  }
}
#table
GRFS.f.table <- data.frame(Variable = colnames(surv.f.bin)[3:ncol(surv.f.bin)], Phylum=c(surv.f.list$phyla,NA) , 
                         GRFS_p.value = as.numeric(GRFS.pvals), 
                         GRFS_BH.adj = p.adjust(as.numeric(GRFS.pvals), method = "fdr"))
#view
GRFS.f.table[GRFS.f.table$GRFS_BH.adj < 0.2 & GRFS.f.table$GRFS_p.value < 0.05,]

### Figure 2C
#plot FDR-corrected p-values

GRFS <- GRFS.f.table[order(GRFS.f.table$GRFS_BH.adj),]
GRFS$Variable <- as.character(GRFS$Variable)
GRFS$Variable <- gsub(".ASV.\\d+","", GRFS$Variable)
GRFS$Variable <- gsub("shannon.bin","Shannon diversity", GRFS$Variable)
GRFS$Variable <- factor(GRFS$Variable, levels = rev(as.character(GRFS$Variable)))

ggplot(GRFS, aes(x=-1*log10(GRFS_BH.adj),y=Variable, fill=Phylum)) + geom_bar(stat="identity") + theme_classic() +
  geom_vline(xintercept = -1*log10(0.05), colour="red", lty=2) + 
  geom_vline(xintercept = -1*log10(0.1), colour="grey", lty=2) +
  geom_vline(xintercept = -1*log10(0.2), colour="grey") +
  xlab("-log10(FDR)") + ylab("Genus") + scale_fill_manual(values = c("darkorange","lightblue"), na.value="grey47") +
  ggtitle("Log-rank test GVHD-free, relapse-free survival - fungal genera") + xlim(c(0,1.9))

### Figure 2D
#GVHD-free,DRS
surv.f.bin <- binary_taxa_surv(ps.f.start, agglom = "Species", read_depth = 200, prev = 0.9, surv_cols = c("GVHD_free_relapse_free_time", "GVHD_free_relapse_free_survival"))
#note: divide by 30.436875 as time provided in days not months
fit<- survfit(Surv(round(GVHD_free_relapse_free_time/30.436875, 1), GVHD_free_relapse_free_survival, type = "right") ~ s__albicans.ASV.4, data = surv.f.bin$otus)
ggsurvplot(fit, pval = TRUE, conf.int = TRUE, risk.table = TRUE, palette = c("red", "blue"))













### Figure S3
#include only samples with read counts of 10,000 and above

#OS
surv.f.list <- binary_taxa_surv(ps.f.start, agglom = "Species", read_depth = 10000, prev = 0.9, surv_cols = c("OS_months", "OS_1_death"))
surv.f.bin <- surv.f.list$otus


### Figure S3A
fit<- survfit(Surv(OS_months, OS_1_death, type = "right") ~ s__albicans.ASV.4, data = surv.f.bin)
ggsurvplot(fit, pval = TRUE, conf.int = TRUE, risk.table = TRUE, palette = c("red", "blue"))
dim(surv.f.bin)




#DFS
surv.f.list <- binary_taxa_surv(ps.f.start, agglom = "Species", read_depth = 10000, prev = 0.9, surv_cols = c("DFS_months", "DFS_1_relapse_or_death"))
surv.f.bin <- surv.f.list$otus



### Figure S3B
fit<- survfit(Surv(DFS_months, DFS_1_relapse_or_death, type = "right") ~ s__albicans.ASV.4, data = surv.f.bin)
ggsurvplot(fit, pval = TRUE, conf.int = TRUE, risk.table = TRUE, palette = c("red", "blue"))
dim(surv.f.bin)



#gvhd-free,rfs
surv.f.list <- binary_taxa_surv(ps.f.start, agglom = "Species", read_depth = 10000, prev = 0.9, surv_cols = c("GVHD_free_relapse_free_time", "GVHD_free_relapse_free_survival"))
surv.f.bin <- surv.f.list$otus


### Figure S3C
fit<- survfit(Surv(round(GVHD_free_relapse_free_time/30.436875,1), GVHD_free_relapse_free_survival, type = "right") ~ s__albicans.ASV.4, data = surv.f.bin)
ggsurvplot(fit, pval = TRUE, conf.int = TRUE, risk.table = TRUE, palette = c("red", "blue"))
dim(surv.f.bin)














###Bacterial analysis
#OS

surv.b.list <- binary_taxa_surv(ps.b.start, agglom = "Genus", read_depth = 200, prev = 0.9, 
                                surv_cols = c("OS_months", "OS_1_death"), tax_level = "Phylum")
surv.b.bin <- surv.b.list$otus


os.pvals <- c()
for(i in 3:ncol(surv.b.bin)){
  if(length(unique(surv.b.bin[,i])) > 1){
    x <- survdiff(Surv(surv.b.bin[,1], surv.b.bin[,2], type = "right") ~ surv.b.bin[,i], data = surv.b.bin)
    p_value <- 1 - pchisq(x$chisq, length(x$n) - 1)
    os.pvals <- c(os.pvals, p_value)
  }
}

OS.b.table <- data.frame(Variable = colnames(surv.b.bin)[3:ncol(surv.b.bin)], Phylum=c(surv.b.list$phyla,NA), 
                         OS_p.value = as.numeric(os.pvals), 
                         OS_BH.adj = p.adjust(as.numeric(os.pvals), method = "fdr"))

#reduce fdr p-value to 0.1 due to large number of genera
OS.b.table[OS.b.table$OS_BH.adj < 0.1 & OS.b.table$OS_p.value < 0.05,]



### Figure S2A
prac <- OS.b.table[order(OS.b.table$OS_BH.adj),]
prac$Variable <- gsub(".ASV.","", prac$Variable)
prac$Variable <- gsub("shannon.bin","Shannon diversity", prac$Variable)
prac$Variable <- factor(prac$Variable, levels = rev(as.character(prac$Variable)))

ggplot(prac, aes(x=-1*log10(OS_BH.adj),y=Variable, fill=Phylum)) + geom_bar(stat="identity") + theme_classic() +
  geom_vline(xintercept = -1*log10(0.05), colour="red", lty=2) + 
  geom_vline(xintercept = -1*log10(0.1), colour="grey", lty=2) +
  xlab("-log10(FDR)") + ylab("") + theme(legend.position = "none") +
  ggtitle("Log-rank test overall survival - bacterial genera") + xlim(c(0,1.8))


### Figure S2B
fit<- survfit(Surv(OS_months, OS_1_death, type = "right") ~ Lachnoclostridium.ASV.7, data = surv.b.bin)
ggsurvplot(fit, pval = TRUE, conf.int = TRUE, risk.table = TRUE, palette = c("red", "blue"))

### Figure S2C
fit<- survfit(Surv(OS_months, OS_1_death, type = "right") ~ shannon.bin, data = surv.b.bin)
ggsurvplot(fit, pval = TRUE, conf.int = TRUE, risk.table = TRUE, palette = c("red", "blue"))




#DFS
surv.b.list <- binary_taxa_surv(ps.b.start, agglom = "Genus", read_depth = 200, prev = 0.9, surv_cols = c("DFS_months", "DFS_1_relapse_or_death"))
surv.b.bin <- surv.b.list$otus

DFS.pvals <- c()
for(i in 3:ncol(surv.b.bin)){
  if(length(unique(surv.b.bin[,i])) > 1){
    x <- survdiff(Surv(surv.b.bin[,1], surv.b.bin[,2], type = "right") ~ surv.b.bin[,i], data = surv.b.bin)
    p_value <- 1 - pchisq(x$chisq, length(x$n) - 1)
    DFS.pvals <- c(DFS.pvals, p_value)
  }
}

DFS.b.table <- data.frame(Variable = colnames(surv.b.bin)[3:ncol(surv.b.bin)], Phylum=c(surv.b.list$phyla,NA) , 
                         DFS_p.value = as.numeric(DFS.pvals), 
                         DFS_BH.adj = p.adjust(as.numeric(DFS.pvals), method = "fdr"))

DFS.b.table[DFS.b.table$DFS_BH.adj < 0.1 & DFS.b.table$DFS_p.value < 0.05,]


### Figure S2D
fit<- survfit(Surv(DFS_months, DFS_1_relapse_or_death, type = "right") ~ shannon.bin, data = surv.b.bin)
ggsurvplot(fit, pval = TRUE, conf.int = TRUE, risk.table = TRUE, palette = c("red", "blue"))




#GVHD-free, RFS
surv.b.list <- binary_taxa_surv(ps.b.start, agglom = "Genus", read_depth = 200, prev = 0.9, surv_cols = c("GVHD_free_relapse_free_time", "GVHD_free_relapse_free_survival"))
surv.b.bin <- surv.b.list$otus

GRFS.pvals <- c()
for(i in 3:ncol(surv.b.bin)){
  if(length(unique(surv.b.bin[,i])) > 1){
    x <- survdiff(Surv(surv.b.bin[,1], surv.b.bin[,2], type = "right") ~ surv.b.bin[,i], data = surv.b.bin)
    p_value <- 1 - pchisq(x$chisq, length(x$n) - 1)
    GRFS.pvals <- c(GRFS.pvals, p_value)
  }
}

GRFS.b.table <- data.frame(Variable = colnames(surv.b.bin)[3:ncol(surv.b.bin)], Phylum=c(surv.b.list$phyla,NA) , 
                         GRFS_p.value = as.numeric(GRFS.pvals), 
                         GRFS_BH.adj = p.adjust(as.numeric(GRFS.pvals), method = "fdr"))

GRFS.b.table[GRFS.b.table$GRFS_BH.adj < 0.1 & GRFS.b.table$GRFS_p.value < 0.05,]

###Shannon diversity not significant on testing








### Figure S4 - samples with low read counts


prac_counts <- cbind(sample_data(ps.f.start), counts = rowSums(otu_table(ps.f.start)))


### Figure S24A - split on 200
prac_counts$counts.bin <- unlist(lapply(prac_counts$counts, function(x) if(x > 200) {"High"} else {"Low"}))
fit<- survfit(Surv(OS_months, OS_1_death, type = "right") ~ counts.bin, data = prac_counts)
ggsurvplot(fit, pval = TRUE, conf.int = TRUE, risk.table = TRUE, palette = c("red", "blue"))


### Figure S24B - split on median
prac_counts$counts.bin <- unlist(lapply(prac_counts$counts, function(x) if(x > median(prac_counts$counts)) {"High"} else {"Low"}))
fit<- survfit(Surv(OS_months, OS_1_death, type = "right") ~ counts.bin, data = prac_counts)
ggsurvplot(fit, pval = TRUE, conf.int = TRUE, risk.table = TRUE, palette = c("red", "blue"))



### Figure S24C - Candida albicans with <200 reads samples included as a separate group
surv.f.list <- binary_taxa_surv(ps.f.start, agglom = "Species", read_depth = 200, prev = 0.9, surv_cols = c("OS_months", "OS_1_death"))
surv.f.bin <- surv.f.list$otus


det.samp <- surv.f.bin[,c("OS_months", "OS_1_death", "s__albicans.ASV.4")]
undet.samp <- sample_data(ps.f.start)[!(rownames(sample_data(ps.f.start)) %in% rownames(surv.f.bin)),]
undet.samp$s__albicans.ASV.4 <- rep("lowCounts", nrow(undet.samp))
undet.samp <- undet.samp[,c("OS_months", "OS_1_death", "s__albicans.ASV.4")]
all.samp <- rbind(det.samp, undet.samp)

#plot KM
fit<- survfit(Surv(OS_months, OS_1_death, type = "right") ~ s__albicans.ASV.4, data = all.samp)
ggsurvplot(fit, pval = TRUE, conf.int = TRUE, risk.table = TRUE, palette = c("red", "blue", "grey"))




#########################################################

            ######## Multivariate ########

#########################################################

#### OS ####
library(survival)
library(survminer)



#subset samples to include those with > 200 reads and present having both bacterial and fungal sequencing available post filtering
ps.f.x <- ps.f.start
ps.f.x <- subset_samples(ps.f.x, rowSums(otu_table(ps.f.x)) > 200)
ps.f.x <- prune_taxa(colSums(otu_table(ps.f.x)) > 0, ps.f.x) # remove any zero columns

ps.b.x <- subset_samples(ps.b.start, rowSums(otu_table(ps.b.start)) > 200)
ps.b.x <- prune_taxa(colSums(otu_table(ps.b.x)) > 0, ps.b.x) # remove any zero columns
ps.b.x <- subset_samples(ps.b.x, rownames(sample_data(ps.b.x)) %in% rownames(sample_data(ps.f.x)))

ps.f.x <- subset_samples(ps.f.x, rownames(sample_data(ps.f.x)) %in% rownames(sample_data(ps.b.x)))
ps.b.x <- subset_samples(ps.b.x, rownames(sample_data(ps.b.x)) %in% rownames(sample_data(ps.f.x)))
all(rownames(sample_data(ps.b.x)) == rownames(sample_data(ps.f.x)))






surv.f.list <- binary_taxa_surv(ps.f.x, agglom = "Species", read_depth = 200, prev = 0.9, 
                                surv_cols = c("OS_months", "OS_1_death", "Antibiotics_at_commencement", "Age.bin",
                                              "DRI_2cat", "RIC_or_MAC", "Cell_source", "HLA_match"))

surv.f.bin <- surv.f.list$otus

surv.b.list <- binary_taxa_surv(ps.b.x, agglom = "Genus", read_depth = 200, prev = 0.9, 
                                surv_cols = c("OS_months", "OS_1_death", "Antibiotics_at_commencement", "Age.bin",
                                              "DRI_2cat", "RIC_or_MAC", "Cell_source", "HLA_match"))

surv.b.bin <- surv.b.list$otus

#simplify factors
surv.b.bin$Cell_source2 <- as.character(surv.b.bin$Cell_source)
surv.b.bin$Cell_source2[surv.b.bin$Cell_source == "CB"] <- "BMSC/CB"
surv.b.bin$Cell_source2[surv.b.bin$Cell_source == "BMSC"] <- "BMSC/CB"
surv.b.bin$DRI_2cat <- as.character(surv.b.bin$DRI_2cat)
surv.b.bin$DRI_2cat[is.na(surv.b.bin$DRI_2cat)] <- "low_int" # samples from AA with no DRI cat included as low_int
surv.b.bin$DRI_2cat <- factor(surv.b.bin$DRI_2cat, levels = c("low_int", "high_very_high"))
surv.b.bin$HLA_match <- as.character(surv.b.bin$HLA_match)
surv.b.bin$HLA_match[surv.b.bin$HLA_match == "MisUD"] <- "UD"
surv.b.bin$HLA_match[surv.b.bin$HLA_match == "MUD"] <- "UD"
surv.b.bin$HLA_match <- as.factor(surv.b.bin$HLA_match)
surv.b.bin$Antibiotics_at_commencement <- as.character(surv.b.bin$Antibiotics_at_commencement)
surv.b.bin$Antibiotics_at_commencement[surv.b.bin$Antibiotics_at_commencement == "Other"] <- "Antibiotics"
surv.b.bin$Antibiotics_at_commencement[surv.b.bin$Antibiotics_at_commencement == "ITM"] <- "Antibiotics"
#combine data-frames with OS, clinical variables and best bacterial (Lachnoclostridium) and fungal (Candida albicans) taxa 
surv.comb <- data.frame(surv.b.bin[,c("OS_months", "OS_1_death", "Antibiotics_at_commencement", "Age.bin",
                                      "DRI_2cat", "RIC_or_MAC", "Cell_source2", "HLA_match", 
                                      "Lachnoclostridium.ASV.7")], surv.f.bin[,c("s__albicans.ASV.4")])

#clean up names for plot
colnames(surv.comb)[3] <- "Antibiotics"
colnames(surv.comb)[4] <- "Age"
colnames(surv.comb)[5] <- "DRI"
colnames(surv.comb)[6] <- "Conditioning"
colnames(surv.comb)[7] <- "Cell source"
colnames(surv.comb)[8] <- "HLA match"
colnames(surv.comb)[9] <- "Lachnoclostridium"
colnames(surv.comb)[10] <- "Candida albicans"

#set reference levels
surv.comb$Lachnoclostridium <- relevel(surv.comb$Lachnoclostridium, ref = "Low")
surv.comb$'Candida albicans' <- relevel(surv.comb$'Candida albicans', ref = "Low")
surv.comb$Age <- relevel(surv.comb$Age, ref = "Under50")
surv.comb$Antibiotics <- as.factor(surv.comb$Antibiotics)
surv.comb$Antibiotics <- relevel(surv.comb$Antibiotics, ref = "None")
surv.comb$'HLA match' <- as.factor(surv.comb$'HLA match')
surv.comb$'HLA match' <- relevel(surv.comb$'HLA match', ref = "MSD")

#run coxph
fit.coxph <- coxph(Surv(surv.comb$OS_months,surv.comb$OS_1_death) ~
                     .,data=surv.comb)
summary(fit.coxph)

### Figure 3
#plot
ggforest(fit.coxph, data = surv.comb)

#check PH assumptions
cox.zph(fit.coxph)





#################################################################

                    ### Taxonomy ###

#################################################################


######compositional analysis - top X OTUs at genus level######

library(RColorBrewer)

### Fungal ###

ntaxa=20
ps2.prop <- ps.f.start
ps2.prop <- subset_samples(ps2.prop, rowSums(otu_table(ps2.prop)) > 200)
TopNOTUs <- names(sort(taxa_sums(ps2.prop), TRUE))[1:ntaxa]
ps2.prop <- transform_sample_counts(ps2.prop, function(x) x/sum(x))
ps2.prop   <- prune_taxa(TopNOTUs, ps2.prop)

remain.taxa <- 1 - rowSums(otu_table(ps2.prop))
taxa2 <- rbind(tax_table(ps2.prop), rep("xOther", 7))
rownames(taxa2)[ntaxa+1] <- "xOther"
taxa2 <- tax_table(taxa2)
otus2 <- otu_table(cbind(otu_table(ps2.prop), xOther = remain.taxa), taxa_are_rows = FALSE)

ps2.prop <- phyloseq(sample_data(ps2.prop), otu_table(otus2), tax_table(taxa2))

otu.hsct.h.genus.0 <- otu_table(ps2.prop)[sample_data(ps2.prop)$OS_1_death == "0",]
otu.hsct.h.genus.1 <- otu_table(ps2.prop)[sample_data(ps2.prop)$OS_1_death == "1",]


xorder.0 <- names(rev(rowSums(otu.hsct.h.genus.0[,grepl("albicans", tax_table(ps2.prop)[,7])])[order(rowSums(otu.hsct.h.genus.0[,grepl("albicans", tax_table(ps2.prop)[,7])]))]))
xorder.1 <- names(rev(rowSums(otu.hsct.h.genus.1[,grepl("albicans", tax_table(ps2.prop)[,7])])[order(rowSums(otu.hsct.h.genus.1[,grepl("albicans", tax_table(ps2.prop)[,7])]))]))


taxa.vec <- paste(gsub("g__","",tax_table(ps2.prop)[,6]), tax_table(ps2.prop)[,7], sep = "_")

taxa.vec <- unlist(lapply(taxa.vec, function(x) if(grepl("g__", x)) {paste0(gsub("g__", "", x), "(g)")} 
                          else if(grepl("f__", x)) {paste0(gsub("f__", "", x), "(f)")} 
                          else if(grepl("c__", x)) {paste0(gsub("c__", "", x), "(c)")} 
                          else if(grepl("o__", x)) {paste0(gsub("o__", "", x), "(o)")} 
                          else if(grepl("p__", x)) {paste0(gsub("p__", "", x), "(p)")} 
                          else if(grepl("k__", x)) {paste0(gsub("k__", "", x), "(k)")} 
                          else if(grepl("s__", x)) {paste0(gsub("s__", "", x), "(s)")} else {x}))

tax_table(ps2.prop) <- cbind(tax_table(ps2.prop), Taxa=taxa.vec)
colnames(tax_table(ps2.prop))[8] <- "Taxa"



p.f <- plot_bar(ps2.prop, x="Sample", y="Abundance", fill="Taxa") + theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), text = element_text(size=18)) + 
  guides(fill=guide_legend(ncol = 1)) + facet_grid(~OS_1_death, scales = "free") + 
  scale_x_discrete(limits=c(xorder.0, xorder.1))

g2 <- ggplot_build(p.f)
cols2 <- unique(g2$data[[1]]["fill"])$fill

dd <- unique(tax_table(ps2.prop)[,8])
names(cols2) <- levels(dd)


cols2[1] <- brewer.pal(6, "Set1")[1]
cols2[2] <- "red4"
cols2[3] <- "#FF689E"
cols2[4] <- "#00BC51"
cols2[5] <- brewer.pal(8, "Set1")[7]
cols2[6] <- "orange"
cols2[7] <- "navy"
#cols2[8] <- "lightgreen"
#cols2[9] <- "maroon"
cols2[10] <- brewer.pal(6, "Set1")[4]
cols2[11] <- "yellow"
cols2[12] <- "dodgerblue"
cols2[13] <- "tan"
cols2[15] <- "grey"

### Figure 1A
plot_bar(ps2.prop, x="Sample", y="Abundance", fill="Taxa") + theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), text = element_text(size=18)) + 
  guides(fill=guide_legend(ncol = 1)) + facet_grid(~OS_1_death, scales = "free") + 
  scale_x_discrete(limits=c(xorder.0, xorder.1)) + scale_fill_manual(values = cols2)




### Bacterial ###

#agglomerate at Genus
taxa <- "Genus"
ntaxa=20
ps2.prop <- tax_glom(ps.b.start, taxa, NArm = FALSE)
ps2.prop <- subset_samples(ps2.prop, rowSums(otu_table(ps2.prop)) > 100)
TopNOTUs <- names(sort(taxa_sums(ps2.prop), TRUE))[1:ntaxa]
ps2.prop <- transform_sample_counts(ps2.prop, function(x) x/sum(x))
ps2.prop   <- prune_taxa(TopNOTUs, ps2.prop)

remain.taxa <- 1 - rowSums(otu_table(ps2.prop))
taxa2 <- rbind(tax_table(ps2.prop), rep("xOther", 7))
rownames(taxa2)[ntaxa+1] <- "xOther"
taxa2 <- tax_table(taxa2)
otus2 <- otu_table(cbind(otu_table(ps2.prop), xOther = remain.taxa), taxa_are_rows = FALSE)

ps2.prop <- phyloseq(sample_data(ps2.prop), otu_table(otus2), tax_table(taxa2))

otu.hsct.h.genus.0 <- otu_table(ps2.prop)[sample_data(ps2.prop)$OS_1_death == "0",]
otu.hsct.h.genus.1 <- otu_table(ps2.prop)[sample_data(ps2.prop)$OS_1_death == "1",]


xorder.0 <- names(rev(rowSums(otu.hsct.h.genus.0[,grepl("Enterococcus", tax_table(ps2.prop)[,6])])[order(rowSums(otu.hsct.h.genus.0[,grepl("Enterococcus", tax_table(ps2.prop)[,6])]))]))
xorder.1 <- names(rev(rowSums(otu.hsct.h.genus.1[,grepl("Enterococcus", tax_table(ps2.prop)[,6])])[order(rowSums(otu.hsct.h.genus.1[,grepl("Enterococcus", tax_table(ps2.prop)[,6])]))]))


taxa.vec <- paste(gsub("g__","",tax_table(ps2.prop)[,6]), tax_table(ps2.prop)[,7], sep = "_")
temp.vec <- is.na(taxa.vec) | grepl("unidentified|unclassified|ambiguous|metagenome|uncultured|NA_NA", ignore.case = TRUE, taxa.vec)
taxa.vec[temp.vec] <- tax_table(ps2.prop)[,6][temp.vec]
temp.vec <- is.na(taxa.vec) | grepl("unidentified|unclassified|ambiguous|metagenome|uncultured|NA_NA", ignore.case = TRUE, taxa.vec)
taxa.vec[temp.vec] <- tax_table(ps2.prop)[,5][temp.vec]

taxa.vec <- unlist(lapply(taxa.vec, function(x) if(grepl("g__", x)) {paste0(gsub("g__", "", x), "(g)")} 
                          else if(grepl("f__", x)) {paste0(gsub("f__", "", x), "(f)")} 
                          else if(grepl("c__", x)) {paste0(gsub("c__", "", x), "(c)")} 
                          else if(grepl("o__", x)) {paste0(gsub("o__", "", x), "(o)")} 
                          else if(grepl("p__", x)) {paste0(gsub("p__", "", x), "(p)")} 
                          else if(grepl("k__", x)) {paste0(gsub("k__", "", x), "(k)")} 
                          else if(grepl("s__", x)) {paste0(gsub("s__", "", x), "(s)")} else {x}))

tax_table(ps2.prop) <- cbind(tax_table(ps2.prop), Taxa=taxa.vec)
colnames(tax_table(ps2.prop))[8] <- "Taxa"



p.f <- plot_bar(ps2.prop, x="Sample", y="Abundance", fill="Taxa") + theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), text = element_text(size=18)) + 
  guides(fill=guide_legend(ncol = 1)) + facet_grid(~OS_1_death, scales = "free") + 
  scale_x_discrete(limits=c(xorder.0, xorder.1))

g2 <- ggplot_build(p.f)
cols2 <- unique(g2$data[[1]]["fill"])$fill

dd <- unique(tax_table(ps2.prop)[,8])
names(cols2) <- levels(dd)


cols2[6] <- "lightblue"
cols2[11] <- "yellow"
cols2[12] <- "red"
cols2[13] <- "orange"
cols2[14] <- "lightgreen"
cols2[16] <- "red4"
cols2[18] <- "darkgreen"
cols2[20] <- "pink"
cols2[21] <- "grey"

### Figure S1
plot_bar(ps2.prop, x="Sample", y="Abundance", fill="Taxa") + theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), text = element_text(size=18)) + 
  guides(fill=guide_legend(ncol = 1)) + facet_grid(~OS_1_death, scales = "free") + 
  scale_x_discrete(limits=c(xorder.0, xorder.1)) + scale_fill_manual(values = cols2)





#######################################################################################################################################