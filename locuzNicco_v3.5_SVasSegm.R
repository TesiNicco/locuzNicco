#!/usr/bin/env Rscript

##################################################
#THIS IS LOCUZNICCO V3.5
##################################################

#LIBRARIES
args = commandArgs(trailingOnly=TRUE)
library(stringr)
library(data.table)
library(gridExtra)
##########

#FUNCTIONS
#function to extract chromosome and position from locus
function.extractChromosomeAndPos <- function(inp.variant){
  chr.pos <- as.data.frame(str_split_fixed(inp.variant, ":", 2))
  colnames(chr.pos) <- c("chr", "pos")
  chr.pos$chr <- as.numeric(as.character(chr.pos$chr))
  chr.pos$pos <- as.numeric(as.character(chr.pos$pos))
  
  return(chr.pos)
}

#function to take GWAS data depending on group chosen and reference genome
function.takeGWASofInterest <- function(groups, reference, chromosome){0
  if (reference == "--hg19"){
    if (groups == "--chc-lasa"){
      inputFile.path <- paste("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/data/assoc_hg19/assoc_chc_lasa/chr", chromosome, "_survival_chc__lasa.association.txt", sep="")
    } else if (groups == "--chc-lasa-scd"){
      inputFile.path <- paste("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/data/assoc_hg19/assoc_chc_lasa_scd/chr", chromosome, "_survival_chc__lasa_scd.association.txt", sep="")
    } else if (groups == "--igap"){
      inputFile.path <- paste("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/data/assoc_hg19/IGAP_2013/chr", chromosome, "_IGAP_S1.txt", sep="")
    } else if (groups == "--ad"){
      inputFile.path <- paste("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/data/assoc_hg19/assoc_ad_chc/chr", chromosome, "_assocAD_annot.txt", sep="")
    } else if (groups == "--grace"){
      inputFile.path <- paste("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/data/assoc_hg19/GRACE/GRACE/chr", chromosome, "_GRACE.txt", sep="")
    } else if (groups == "--igap-2019"){
      inputFile.path <- paste("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/data/assoc_hg19/IGAP_2019/chr", chromosome, "_IGAP_2k19.txt", sep="")
    } else if (groups == "--spigap"){
      inputFile.path <- paste("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/data/assoc_hg19/GRACE/SPIGAP/chr", chromosome, "_SPIGAP.txt", sep="")
    } else if (groups == "--spigap-ukb"){
      inputFile.path <- paste("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/data/assoc_hg19/GRACE/SPIGAP_UKB/chr", chromosome, "_SPIGAP_UKB.txt", sep="")
    }
  } else if (reference == "--hg38"){
    if (groups == "--chc-lasa"){
      inputFile.path <- paste("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/data/assoc_hg38/assoc_chc_lasa/chr", chromosome, "_chc_lasa_hg38_annot.bed", sep="")
    } else if (groups == "--chc-lasa-scd"){
      inputFile.path <- paste("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/data/assoc_hg38/assoc_chc_lasa_scd/chr", chromosome, "_chc_lasa_scd_hg38_annot.bed", sep="")
    } else if (groups == "--igap"){
      inputFile.path <- paste("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/data/assoc_hg38/IGAP_2013/chr", chromosome, "_IGAP_S1_hg38_annot.table", sep="")
    } else if (groups == "--ad"){
      inputFile.path <- paste("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/data/assoc_hg38/assoc_ad_chc/chr", chromosome, "_assocAD_annot_hg38.table", sep="")
    }
  }
  return(inputFile.path)
}

#function to manage association file in case is IGAP
function.manageIGAP <- function(groups, reference, inputFile.path, chromosome){
  if (reference == "--hg19"){
    #open association file
    assoc <- fread(inputFile.path, h=F)
    
    if (groups == "--igap"){
      #clean association file
      colnames(assoc) <- c("chr", "pos", "rsID", "a1", "a0", "beta", "se", "pvalue")
      #create supplementary column
      assoc$id <- paste(assoc$chr, assoc$pos, sep=":")
      
    } else if (groups == "--igap-2019"){
      #clean association file
      colnames(assoc) <- c("locID", "id", "chr", "pos", "a1", "a0", "beta", "se", "pvalue", "a1_frq", "or", "or_CI", "rsID", "stage")
      assoc$pos <- as.numeric(as.character(assoc$pos))
      #restrict to chromosome of interest
      assoc <- assoc[which(assoc$chr == chromosome),]
      assoc$beta <- as.numeric(as.character(assoc$beta))
    }
  } else if (reference == "--hg38"){
    #open association file
    assoc <- fread(inputFile.path, h=F)
    #clean association file
    colnames(assoc) <- c("chr", "pos.hg38", "pos.2.hg38", "id", "rsID", "a1", "a0", "beta", "se", "pvalue")
  }  
  return(assoc)
}

#function to select variant window and compute LD there, the output is a data frame merged with association results and ld informations
function.variantLD <- function(assoc, chr.pos, inp.variant, window.user){
  #define interval in which to compute LD
  max.value <- chr.pos$pos + window.user
  min.value <- chr.pos$pos - window.user
  beta.limit <- 4
  #subset association file to defined interval, and exclude variants with unreliable effect-size
  assoc.res <- assoc[which((assoc$pos >= min.value) & (assoc$pos <= max.value) & (abs(assoc$beta) <= beta.limit)),]
  
  #output list of variants to compute LD
  var.list <- as.data.frame(assoc.res$id)
  write.table(var.list, "list_variants.tmp", quote=F, row.names=F)
  
  #now compute LD
  print("Compute LD between variants in +- 150Kb window...")
  ##step 1 -- make bed file out of pfile
  pinputFile.path <- paste("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/data/genotypes_forLD/chr", chr.pos$chr, ".dose", sep="")
  indv.kept <- "/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/data/genotypes_forLD/all_indvs_kept.txt"
  cmd.conversion <- paste("plink2 --pfile ", pinputFile.path, " --extract list_variants.tmp --keep ", indv.kept, " --make-bed --out tmp", sep="")
  system(cmd.conversion)
  ##step 2 -- compute ld
  cmd.ld <- "plink --bfile tmp --r2 yes-really inter-chr --out out.ld"
  system(cmd.ld)
  
  #read LD info
  print("Reading LD info...")
  ld <- fread("out.ld.ld", h=T)
  
  #find all occurrences of the input variant
  occurrences <- c(grep(inp.variant, ld$SNP_A), grep(inp.variant, ld$SNP_B))
  
  #take set of variants that include the input variant
  ld.set <- ld[occurrences,]
  
  #re-format, i.e input variant always in one column, other variant in another
  ld.set.format <- as.data.frame(matrix(data = NA, nrow = nrow(ld.set), ncol = 3))
  colnames(ld.set.format) <- c("input.var", "other.var", "r2")
  for (i in 1:nrow(ld.set)){
    if (ld.set$SNP_A[i] == inp.variant){
      ld.set.format[i, ] <- c(ld.set$SNP_A[i], ld.set$SNP_B[i], ld.set$R2[i])
    } else {
      ld.set.format[i, ] <- c(ld.set$SNP_B[i], ld.set$SNP_A[i], ld.set$R2[i])
    }
  }
  
  #now merge with association results
  trial <- merge(assoc.res, ld.set.format, by.x="id", by.y="other.var", all.x=T)
  
  return(trial)
}

#function to select variant window and compute LD there, the output is a data frame merged with association results and ld informations -- this is for hg38
function.variantLD.hg38 <- function(assoc, chr.pos, inp.variant, groups, window.user){
  #define interval in which to compute LD
  max.value <- chr.pos$pos + window.user
  min.value <- chr.pos$pos - window.user
  beta.limit <- 4
  #subset association file to defined interval, and exclude variants with unreliable effect-size
  assoc.res <- assoc[which((assoc$pos.hg38 >= min.value) & (assoc$pos.hg38 <= max.value) & (abs(assoc$beta) <= beta.limit)),]
  
  #add locus id for hg38
  if (groups == "--igap"){
    assoc.res$id.hg38 <- paste(chr.pos$chr, assoc.res$pos.hg38, sep=":")
  } else {
    assoc.res$id.hg38 <- paste(assoc.res$chrom, assoc.res$pos.hg38, sep=":")
  }
  
  #identify same variant in hg19
  old.id <- assoc.res$id[which(assoc.res$id.hg38 == inp.variant)]
  
  #output list of variants to compute LD
  var.list <- as.data.frame(assoc.res$id)
  write.table(var.list, "list_variants.tmp", quote=F, row.names=F)
  
  #now compute LD
  print("Compute LD between variants in +- 150Kb window...")
  ##step 1 -- make bed file out of pfile
  pinputFile.path <- paste("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/data/genotypes_forLD/chr", chr.pos$chr, ".dose", sep="")
  indv.kept <- "/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/data/genotypes_forLD/all_indvs_kept.txt"
  cmd.conversion <- paste("plink2 --pfile ", pinputFile.path, " --extract list_variants.tmp --keep ", indv.kept, " --make-bed --out tmp", sep="")
  system(cmd.conversion)
  ##step 2 -- compute ld
  cmd.ld <- "plink --bfile tmp --r2 yes-really inter-chr --out out.ld"
  system(cmd.ld)
  
  #read LD info
  print("Reading LD info...")
  ld <- fread("out.ld.ld", h=T)
  
  #find all occurrences of the input variant
  occurrences <- c(grep(old.id, ld$SNP_A), grep(old.id, ld$SNP_B))
  
  #take set of variants that include the input variant
  ld.set <- ld[occurrences,]
  
  #re-format, i.e input variant always in one column, other variant in another
  ld.set.format <- as.data.frame(matrix(data = NA, nrow = nrow(ld.set), ncol = 3))
  colnames(ld.set.format) <- c("input.var", "other.var", "r2")
  for (i in 1:nrow(ld.set)){
    if (ld.set$SNP_A[i] == old.id){
      ld.set.format[i, ] <- c(ld.set$SNP_A[i], ld.set$SNP_B[i], ld.set$R2[i])
    } else {
      ld.set.format[i, ] <- c(ld.set$SNP_B[i], ld.set$SNP_A[i], ld.set$R2[i])
    }
  }
  
  #now merge with association results
  trial <- merge(assoc.res, ld.set.format, by.x="id", by.y="other.var", all.x=T)
  
  return(trial)
}

#function to compute LD, then select leading snp -- in case type of analysis is --gen
function.variantLD.Gene <- function(chromosome, assoc, start, end, reference){
  if (reference == "--hg19"){
    max.value <- end
    min.value <- start
    beta.limit <- 4
    #subset association file to defined interval, and exclude variants with unreliable effect-size
    assoc.res <- assoc[which((assoc$pos >= min.value) & (assoc$pos <= max.value) & (abs(assoc$beta) <= beta.limit)),]
    
    #output list of variants to compute LD
    var.list <- as.data.frame(assoc.res$id)
    write.table(var.list, "list_variants.tmp", quote=F, row.names=F)
    
    #now compute LD
    print("Compute LD between variants in +- 150Kb window...")
    ##step 1 -- make bed file out of pfile
    pinputFile.path <- paste("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/data/genotypes_forLD/chr", chromosome.n, ".dose", sep="")
    indv.kept <- "/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/data/genotypes_forLD/all_indvs_kept.txt"
    cmd.conversion <- paste("plink2 --pfile ", pinputFile.path, " --extract list_variants.tmp --keep ", indv.kept, " --make-bed --out tmp", sep="")
    system(cmd.conversion)
    ##step 2 -- compute ld
    cmd.ld <- "plink --bfile tmp --r2 yes-really inter-chr --out out.ld"
    system(cmd.ld)
    
    #select leading variant
    assoc.res$pvalue <- as.numeric(as.character(assoc.res$pvalue))
    assoc.res <- assoc.res[order(assoc.res$pvalue),]
    leading.snp <- assoc.res$id[1]
    
    #read LD info
    print("Reading LD info...")
    ld <- fread("out.ld.ld", h=T)
    
    #find all occurrences of the input variant
    occurrences <- c(grep(leading.snp, ld$SNP_A), grep(leading.snp, ld$SNP_B))
    
    #if there are no occurrences -- variant is too rare or not in LD with any other variant, take second as leading snp -- do this 3 times, then give a message
    if (length(occurrences) == 0){
      #take second leading snp
      leading.snp <- assoc.res$id[2]
      occurrences <- c(grep(leading.snp, ld$SNP_A), grep(leading.snp, ld$SNP_B))
      
      if (length(occurrences) == 0){
        #take third leading snp
        leading.snp <- assoc.res$id[3]
        occurrences <- c(grep(leading.snp, ld$SNP_A), grep(leading.snp, ld$SNP_B))
        
      }
    }
    
    #take set of variants that include the input variant
    ld.set <- ld[occurrences,]
    
    #re-format, i.e input variant always in one column, other variant in another
    ld.set.format <- as.data.frame(matrix(data = NA, nrow = nrow(ld.set), ncol = 3))
    colnames(ld.set.format) <- c("input.var", "other.var", "r2")
    for (i in 1:nrow(ld.set)){
      if (ld.set$SNP_A[i] == leading.snp){
        ld.set.format[i, ] <- c(ld.set$SNP_A[i], ld.set$SNP_B[i], ld.set$R2[i])
      } else {
        ld.set.format[i, ] <- c(ld.set$SNP_B[i], ld.set$SNP_A[i], ld.set$R2[i])
      }
    }
    
    #now merge with association results
    trial <- merge(assoc.res, ld.set.format, by.x="id", by.y="other.var", all.x=T)
  } else if (reference == "--hg38"){
    max.value <- end
    min.value <- start
    beta.limit <- 4
    #subset association file to defined interval, and exclude variants with unreliable effect-size
    assoc.res <- assoc[which((assoc$pos.hg38 >= min.value) & (assoc$pos.hg38 <= max.value) & (abs(assoc$beta) <= beta.limit)),]
    
    #output list of variants to compute LD
    var.list <- as.data.frame(assoc.res$id)
    write.table(var.list, "list_variants.tmp", quote=F, row.names=F)
    
    #now compute LD
    print("Compute LD between variants in +- 150Kb window...")
    ##step 1 -- make bed file out of pfile
    pinputFile.path <- paste("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/data/genotypes_forLD/chr", chromosome.n, ".dose", sep="")
    indv.kept <- "/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/data/genotypes_forLD/all_indvs_kept.txt"
    cmd.conversion <- paste("plink2 --pfile ", pinputFile.path, " --extract list_variants.tmp --keep ", indv.kept, " --make-bed --out tmp", sep="")
    system(cmd.conversion)
    ##step 2 -- compute ld
    cmd.ld <- "plink --bfile tmp --r2 yes-really inter-chr --out out.ld"
    system(cmd.ld)
    
    #select leading variant
    assoc.res$pvalue <- as.numeric(as.character(assoc.res$pvalue))
    assoc.res <- assoc.res[order(assoc.res$pvalue),]
    leading.snp <- assoc.res$id[1]
    print(paste("*****", leading.snp, sep=''))
    
    #read LD info
    print("Reading LD info...")
    ld <- fread("out.ld.ld", h=T)
    
    #find all occurrences of the input variant
    occurrences <- c(grep(leading.snp, ld$SNP_A), grep(leading.snp, ld$SNP_B))
    
    #if there are no occurrences -- variant is too rare or not in LD with any other variant, take second as leading snp -- do this 3 times, then give a message
    if (length(occurrences) == 0){
      #take second leading snp
      leading.snp <- assoc.res$id[2]
      occurrences <- c(grep(leading.snp, ld$SNP_A), grep(leading.snp, ld$SNP_B))
      
      if (length(occurrences) == 0){
        #take third leading snp
        leading.snp <- assoc.res$id[3]
        occurrences <- c(grep(leading.snp, ld$SNP_A), grep(leading.snp, ld$SNP_B))
        
      }
    }
    print(paste("*****", leading.snp, sep=''))
    
    #take set of variants that include the input variant
    ld.set <- ld[occurrences,]
    
    #re-format, i.e input variant always in one column, other variant in another
    ld.set.format <- as.data.frame(matrix(data = NA, nrow = nrow(ld.set), ncol = 3))
    colnames(ld.set.format) <- c("input.var", "other.var", "r2")
    for (i in 1:nrow(ld.set)){
      if (ld.set$SNP_A[i] == leading.snp){
        ld.set.format[i, ] <- c(ld.set$SNP_A[i], ld.set$SNP_B[i], ld.set$R2[i])
      } else {
        ld.set.format[i, ] <- c(ld.set$SNP_B[i], ld.set$SNP_A[i], ld.set$R2[i])
      }
    }
    
    #now merge with association results
    trial <- merge(assoc.res, ld.set.format, by.x="id", by.y="other.var", all.x=T)
    
  }  
  return(list(leading.snp, trial))
}

#function to load frequencies from file and merge file with association-ld file -- only for non-IGAP as IGAP is only common variants
function.combineFrequencies <- function(trial, inp.variant, chromosome){
  frq.path <- paste("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/data/frequencies/chr", chromosome, ".frq.afreq", sep="")
  frq <- fread(frq.path, h=T, check.names = F)
  
  #adjust frequencies
  frq$maf <- frq$ALT_FREQS
  frq$minor.allele <- frq$ALT
  frq$maf[which(frq$ALT_FREQS > 0.50)] <- 1 - frq$maf[which(frq$ALT_FREQS > 0.50)]
  frq$minor.allele[which(frq$ALT_FREQS > 0.50)] <- frq$REF[which(frq$ALT_FREQS > 0.50)]
  frq$locid <- paste(frq$ID, frq$ALT, sep=":")
  
  #merge with association and LD results
  trial.frq <- merge(trial, frq, by.x='locus.test', by.y='locid')
  
  #subset of minor allele frequency
  trial.frq <- trial.frq[which(trial.frq$maf >= 0.005),]
  
  return(trial.frq)
}

#function to prepare for plot -- assign colors, pch and sizes
function.plotPreparation <- function(trial, inp.variant, reference){
  #assign colors to snps depending on ld (r2)
  trial$col <- "navy"
  trial$col[which((trial$r2 >= 0.20) & (trial$r2 < 0.40))] <- "deepskyblue2"
  trial$col[which((trial$r2 >= 0.40) & (trial$r2 < 0.60))] <- "khaki"
  trial$col[which((trial$r2 >= 0.60) & (trial$r2 < 0.80))] <- "orange"
  trial$col[which((trial$r2 >= 0.80) & (trial$r2 <= 1))] <- "coral"
  trial$col[which(trial$id == inp.variant)] <- "black"
  
  #assign pch info
  trial$pch <- 21
  trial$pch[which(trial$id == inp.variant)] <- 23
  
  #assign pch size
  trial$size <- 1
  trial$size[which(trial$col == "deepskyblue2")] <- 1.20
  trial$size[which(trial$col == "khaki")] <- 1.40
  trial$size[which(trial$col == "orange")] <- 1.60
  trial$size[which(trial$col == "coral")] <- 1.80
  trial$size[which(trial$id == inp.variant)] <- 1.80
  
  #add log10P
  trial$pvalue <- as.numeric(as.character(trial$pvalue))
  trial$log10P <- -log10(trial$pvalue)
  
  #limits for plotting
  if (reference == "--hg19"){
    min.x <- min(trial$pos)
    max.x <- max(trial$pos)
    max.y <- max(trial$log10P[is.finite(trial$log10P)])
    interval <- (max.x - min.x)/3
  } else if (reference == "--hg38"){
    trial$col[which(trial$id.hg38 == inp.variant)] <- "black"
    trial$pch[which(trial$id.hg38 == inp.variant)] <- 23
    trial$size[which(trial$id.hg38 == inp.variant)] <- 1.80
    
    min.x <- min(trial$"pos.hg38")
    max.x <- max(trial$"pos.hg38")
    max.y <- max(trial$log10P[is.finite(trial$log10P)])
    interval <- (max.x - min.x)/3
  }
  #check if max.y is lower than 5 --> use 5 as maximum
  if (max.y < 5){
    max.y <- 5
  }
  
  return(list(trial, min.x, max.x, max.y, interval))
}

#function to parse gene location file and grep only region of interest
function.parseGeneLocation <- function(reference, chromosome, min.x, max.x){
  if (reference == "--hg19"){
    gene.loc <- fread("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/databases/hg19_geneListandPos", h=T)
    
    gene.loc.res <- subset(gene.loc, gene.loc$chrom == paste("chr", chromosome, sep=""))
    gene.loc.res$in.int <- 0
    gene.loc.res$in.int[which((gene.loc.res$txStart >= min.x) & (gene.loc.res$txEnd <= max.x))] <- 1
    gene.loc.res$in.int[which((gene.loc.res$txStart <= min.x) & (gene.loc.res$txEnd >= min.x))] <- 1
    gene.loc.res$in.int[which((gene.loc.res$txStart <= max.x) & (gene.loc.res$txEnd >= max.x))] <- 1
    genes <- gene.loc.res[which(gene.loc.res$in.int == 1),]
    genes <- genes[order(-genes$txEnd),]
    genes <- genes[!duplicated(genes$"#geneName"),]
    if ((nrow(genes) > 3) & (nrow(genes) <= 10)){
      genes$y <- seq(-1, -4, -1)
    } else if (nrow(genes) <= 3) {
      genes$y <- seq(-1, -3, -1)
    } else if (nrow(genes) > 10){
      genes$y <- seq(-1, -5, -1)
    }
  } else if (reference == "--hg38"){
    gene.loc <- fread("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/databases/trial_hg38_genes.txt", h=T)
    
    gene.loc.res <- subset(gene.loc, gene.loc$chrom == paste("chr", chromosome, sep=""))
    gene.loc.res$in.int <- 0
    gene.loc.res$in.int[which((gene.loc.res$txStart >= min.x) & (gene.loc.res$txEnd <= max.x))] <- 1
    gene.loc.res$in.int[which((gene.loc.res$txStart <= min.x) & (gene.loc.res$txEnd >= min.x))] <- 1
    gene.loc.res$in.int[which((gene.loc.res$txStart <= max.x) & (gene.loc.res$txEnd >= max.x))] <- 1
    genes <- gene.loc.res[which(gene.loc.res$in.int == 1),]
    genes <- genes[order(-genes$txEnd),]
    genes <- genes[!duplicated(genes$name2),]
    if ((nrow(genes) > 3) & (nrow(genes) < 10)){
      genes$y <- seq(-1, -4, -1)
    } else if (nrow(genes) <= 3){
      genes$y <- seq(-1, -3, -1)
    } else if (nrow(genes) > 10){
      genes$y <- seq(-1, -5, -1)
    }
  }
  
  return(genes)
}

#function to parse recombination peaks file and grep only region of interest
function.parseGeneticMap <- function(reference, chromosome, min.x, max.x, max.y){
  if (reference == "--hg19"){
    genmap.path <- paste("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/databases/recomb_rates/hg19/genetic_map_chr", chromosome, "_combined_b37.txt", sep="")
    genmap <- fread(genmap.path, h=T)
    colnames(genmap) <- c("position", "combined.rate", "cM")
    
    #normalize genmap between 0 and 8
    genmap$norm.rate <- (max.y - 0) * ((genmap$combined.rate - min(genmap$combined.rate)) / (max(genmap$combined.rate) - min(genmap$combined.rate)) + 0 )
    trial.norm <- data.frame(x = c(0, 20, 40, 60, 80, 100), y = c(0, 0, 0, 0, 0, 0))
    trial.norm$y <- (max.y - 0) * ((trial.norm$x - min(trial.norm$x)) / (max(trial.norm$x) - min(trial.norm$x)) + 0)
    
    #restrict to interval
    genmap <- genmap[which((genmap$position >= min.x) & (genmap$position <= max.x)),]
    genmap <- genmap[order(genmap$position),]
  } else if (reference == "--hg38"){
    genmap.path <- paste("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/databases/recomb_rates/hg38/recombinationRates_hg38_chr", chromosome, ".txt", sep="")
    genmap <- fread(genmap.path, h=T)
    
    #normalize genmap between 0 and 8
    genmap$norm.rate <- (max.y - 0) * ((genmap$rate - min(genmap$rate)) / (max(genmap$rate) - min(genmap$rate)) + 0 )
    trial.norm <- data.frame(x = c(0, 20, 40, 60, 80, 100), y = c(0, 0, 0, 0, 0, 0))
    trial.norm$y <- (max.y - 0) * ((trial.norm$x - min(trial.norm$x)) / (max(trial.norm$x) - min(trial.norm$x)) + 0)
    
    #restrict to interval
    genmap <- genmap[which((genmap$bp >= min.x) & (genmap$bp <= max.x)),]
    genmap <- genmap[order(genmap$bp),]
  } 
  return(list(genmap, trial.norm))
}

#function to plot png
function.plotPNG <- function(trial, trial.norm, chromosome, genes, genmap, min.x, max.x, max.y, inp.variant, interval, chromatin){
  par(mar=c(6, 6, 6, 6))
  
  #check whether there are genes in the plotted region
  min.y <- -1
  if (nrow(genes) > 0){
    min.y <- min(genes$y)
  }
  
  #set the minimun y-axis for p-value at 8 -- genome-wide significance
  max.y <- as.numeric(as.character(max.y))
  if (max.y < 8){
    max.y <- 8
  }
  
  if (nrow(genes) == 0){
    min.genes.y <- -3
  } else {
    min.genes.y <- min(genes$y)
  }
  
  #use the length as proportion for object sizes
  length.y <- abs(min.genes.y - max.y/4/4*6) + max.y 
  
  ##PLOT 1
  #title was --> main=paste("Chromosome ", chromosome, " ~ ", inp.variant, sep = "")
  
  #plot recombination map and its relative axis and label
  plot(genmap$position, genmap$norm.rate, type='l', col='darkolivegreen3', xaxt='none', ylim=c(min.genes.y - max.y/4/4*7, max.y), 
       xlim=c(min.x, max.x), xlab='Chromosomal Position (Mb)', ylab='-log10(P-value)', yaxt='none', cex.lab=1.40, lwd=2,
       main = inp.variant)
  axis(side = 4, at = trial.norm$y, labels=trial.norm$x, col='darkolivegreen3', col.axis = "darkolivegreen3", cex.axis=1.25)
  mtext("Recombination rate (cM/Mb)", side=4, line=3, cex=1.40, col='darkolivegreen3')
  
  #line at genome-wide significance
  abline(h=-log10(0.00000005), lty=2, lwd=2, col='dark grey')
  
  #plot point with low LD
  points(trial$pos[which(trial$col == "navy")], trial$log10P[which(trial$col == "navy")], 
         pch=trial$pch[which(trial$col == "navy")], bg=trial$col[which(trial$col == "navy")], 
         col=trial$col[which(trial$col == "navy")], cex=trial$size[which(trial$col == "navy")])
  #plot point in LD with my variant
  points(trial$pos[which(trial$col != "navy")], trial$log10P[which(trial$col != "navy")], 
         pch=trial$pch[which(trial$col != "navy")], bg=trial$col[which(trial$col != "navy")], 
         col=trial$col[which(trial$col != "navy")], cex=trial$size[which(trial$col != "navy")])
  #plot my variant
  points(trial$pos[which(trial$id == inp.variant)], trial$log10P[which(trial$id == inp.variant)], 
         pch=trial$pch[which(trial$id == inp.variant)], bg = trial$col[which(trial$id == inp.variant)], 
         cex=trial$size[which(trial$id == inp.variant)])
  #plot axes and legend
  axis(side = 1, at=c(min.x, min.x+(interval), min.x+(interval*2), max.x), cex.axis=1.25, 
       labels=c(round(min.x/1000000, 2), round((min.x+(interval))/1000000, 2), round((min.x+(interval*2))/1000000, 2), 
                round((max.x/1000000), 2)))
  axis(side = 2, at=seq(0, round(max.y, 0), round(max.y/4, 0)), labels=seq(0, round(max.y, 0), round(max.y/4, 0)), cex.axis=1.25)
  legend('topright', legend=c(0, 0.2, 0.4, 0.6, 0.8), pch=22, pt.bg=c("navy", "deepskyblue2", "khaki", "orange", "coral"), 
         cex=1.40, ncol=1, title=as.expression(bquote(r^2)) )
  
  #plot genes
  if (nrow(genes) > 0){
    for (g in 1:nrow(genes)){
      segments(x0 = genes$txStart[g], y0 = genes$y[g] - (length.y*1.12/100), x1 = genes$txEnd[g], y1 = genes$y[g] - (length.y*1.12/100), lwd=2)
      start <- as.data.frame(t(str_split_fixed(genes$exonStarts[g], ",", genes$exonCount[g]+1)))
      end <- as.data.frame(t(str_split_fixed(genes$exonEnds[g], ",", genes$exonCount[g]+1)))
      exons <- cbind(start, end)
      colnames(exons) <- c("start", "end")
      exons$start <- as.numeric(as.character(exons$start))
      exons$end <- as.numeric(as.character(exons$end))
      for (j in 1:nrow(exons)){
        rect(xleft=exons$start[j], ybottom=genes$y[g] - (length.y*1.12/100) - 0.15, xright=exons$end[j], ytop = genes$y[g] - (length.y*1.12/100) + 0.15, col='grey80')
      }
      text(x = genes$txStart[g] + (genes$txEnd[g] - genes$txStart[g])/2, y = genes$y[g] + 0.25, labels=genes$"#geneName"[g], font=4, cex=1.15)
      if (genes$strand[g] == "+"){
        arrows(x0 = genes$txEnd[g], y0 = genes$y[g] - (length.y*1.12/100), x1 = genes$txEnd[g] + 3000, y1 = genes$y[g] - (length.y*1.12/100), length=0.1, lwd=2, col='coral')
      } else if (genes$strand[g] == "-"){
        arrows(x0 = genes$txStart[g], y0 = genes$y[g] - (length.y*1.12/100), x1 = genes$txStart[g] - 3000, y1 = genes$y[g] - (length.y*1.12/100), length=0.1, lwd=2, col='coral')
      }
    }
  }  
  #new line before chromatin states
  abline(h=min.genes.y - (max.y/4/4*3), lty=1, lwd=1.50)
  
  #add plotting position
  chromatin$pos.plt <- min.genes.y-max.y/4/4*4
  
  #now add chromatin states
  chrom.colors <- str_split_fixed(chromatin$V5, ",", 3)
  chromatin$red <- as.numeric(chrom.colors[, 1])
  chromatin$green <- as.numeric(chrom.colors[, 2])
  chromatin$blue <- as.numeric(chrom.colors[, 3])
  for (row in 1:nrow(chromatin)){
    rect(xleft = chromatin$V2[row], ybottom = chromatin$pos.plt[row], xright = chromatin$V3[row], ytop = chromatin$pos.plt[row] - (length.y*4/100), 
         col=rgb(red = chromatin$red[row], green = chromatin$green[row], blue = chromatin$blue[row], maxColorValue = 255), border = "black",
         lwd = 0.50)
  }
  
  #read chromatin legend file
  chrom.legend <- read.table("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/databases/chromatin_states/chromLegend.txt", h=T, sep='\t', check.names = F)
  chrom.legend$colors <- paste("#", chrom.legend$colo, sep="")
  mtext("Chromatin\nstates", side=2, line=1, cex=1.10, col='black', at = min.genes.y-max.y/4/4*4.5, adj = 0.5)
  legend("bottom", legend = chrom.legend$state, pch=22, pt.bg = chrom.legend$colors, bty='n', ncol=nrow(chrom.legend)/2, 
         cex=0.90, pt.cex = 1.25, x.intersp = 0.40)
  
  
  return("PNG produced")
}

#function to plot pdf
function.plotPDF <- function(trial, trial.norm, chromosome, genes, genmap, min.x, max.x, max.y, inp.variant, interval, groups){
  par(mar=c(6, 6, 6, 6))
  
  ##PLOT 1
  #plot recombination map and its relative axis and label
  plot(genmap$position, genmap$norm.rate, type='l', col='darkolivegreen3', xaxt='none', ylim=c(min(genes$y), max.y), xlim=c(min.x, max.x), xlab='Chromosomal Position (Mb)', ylab='-log10(P-value)', yaxt='none', cex.lab=1.40, main=paste("Chromosome ", chromosome, " ~ ", inp.variant, sep = ""), lwd=2)
  axis(side = 4, at = trial.norm$y, labels=trial.norm$x, col='darkolivegreen3', col.axis = "darkolivegreen3", cex.axis=1.25)
  mtext("Recombination rate (cM/Mb)", side=4, line=3, cex=1.40, col='darkolivegreen3')
  #plot point with low LD
  points(trial$pos[which(trial$col == "navy")], trial$log10P[which(trial$col == "navy")], pch=trial$pch[which(trial$col == "navy")], bg=trial$col[which(trial$col == "navy")], col=trial$col[which(trial$col == "navy")], cex=trial$size[which(trial$col == "navy")])
  #plot point in LD with my variant
  points(trial$pos[which(trial$col != "navy")], trial$log10P[which(trial$col != "navy")], pch=trial$pch[which(trial$col != "navy")], bg=trial$col[which(trial$col != "navy")], col=trial$col[which(trial$col != "navy")], cex=trial$size[which(trial$col != "navy")])
  #plot my variant
  points(trial$pos[which(trial$id == inp.variant)], trial$log10P[which(trial$id == inp.variant)], pch=trial$pch[which(trial$id == inp.variant)], bg = trial$col[which(trial$id == inp.variant)], cex=trial$size[which(trial$id == inp.variant)])
  #plot axes and legend
  axis(side = 1, at=c(min.x, min.x+(interval), min.x+(interval*2), max.x), cex.axis=1.25, labels=c(round(min.x/1000000, 2), round((min.x+(interval))/1000000, 2), round((min.x+(interval*2))/1000000, 2), round((max.x/1000000), 2)))
  axis(side = 2, at=seq(0, 8, 1), labels=seq(0, 8, 1), cex.axis=1.25)
  legend('topright', legend=c(0, 0.2, 0.4, 0.6, 0.8), pch=22, pt.bg=c("navy", "deepskyblue2", "khaki", "orange", "coral"), cex=1.40, ncol=1, title=as.expression(bquote(r^2)) )
  
  #plot genes
  for (g in 1:nrow(genes)){
    segments(x0 = genes$txStart[g], y0 = genes$y[g], x1 = genes$txEnd[g], y1 = genes$y[g], lwd=3)
    start <- as.data.frame(t(str_split_fixed(genes$exonStarts[g], ",", genes$exonCount[g]+1)))
    end <- as.data.frame(t(str_split_fixed(genes$exonEnds[g], ",", genes$exonCount[g]+1)))
    exons <- cbind(start, end)
    colnames(exons) <- c("start", "end")
    exons$start <- as.numeric(as.character(exons$start))
    exons$end <- as.numeric(as.character(exons$end))
    for (j in 1:nrow(exons)){
      rect(xleft=exons$start[j], ybottom=genes$y[g]-0.15, xright=exons$end[j], ytop = genes$y[g]+0.15, col='grey80')
    }
    text(x = genes$txStart[g] + (genes$txEnd[g] - genes$txStart[g])/2, y = genes$y[g] + 0.50, labels=genes$"#geneName"[g], font=3, cex=0.85)
    if (genes$strand[g] == "+"){
      arrows(x0 = genes$txEnd[g], y0 = genes$y[g], x1 = genes$txEnd[g] + 3000, y1 = genes$y[g], length=0.1, lwd=2, col='coral')
    } else if (genes$strand[g] == "-"){
      arrows(x0 = genes$txStart[g], y0 = genes$y[g], x1 = genes$txStart[g] - 3000, y1 = genes$y[g], length=0.1, lwd=2, col='coral')
    }
  }
  
  #if group is igap, table should be slightly different
  if (groups == "--igap"){
    #prepare table
    toplot <- trial[, c("id", "a1", "beta", "se", "pvalue", "rsID")]
    toplot.sign <- toplot[which(toplot$pvalue <= 0.05),]
    toplot.sign <- toplot.sign[order(toplot.sign$pvalue),]
    toplot.sign <- as.data.frame(toplot.sign)
    toplot.sign <- na.omit(toplot.sign)
    toplot.sign <- toplot.sign[!duplicated(toplot.sign$id),]
    colnames(toplot.sign) <- c("LOCUS", "A1", "BETA", "SE", "P", "RSID")
    d <- seq(1, nrow(toplot.sign))
    a <- split(d, ceiling(seq_along(d)/22))
    for (page in 1:length(a)){
      df <- as.data.frame(a[page])
      #empty thing to be covered -- go to new page
      plot(1, bty='none', xaxt='none', yaxt='none', col='white', xlab='', ylab='')
      grid.table(toplot.sign[df[1, 1]:df[nrow(df), 1], ], rows = NULL)
    }
    
  } else {
    
    #prepare table
    toplot <- trial[, c("id", "a1", "beta", "odds.ratio", "se", "pvalue", "genotype", "rsID", "maf")]
    toplot.sign <- toplot[which(toplot$pvalue <= 0.05),]
    toplot.sign <- toplot.sign[order(toplot.sign$pvalue),]
    toplot.sign <- as.data.frame(toplot.sign)
    toplot.sign <- na.omit(toplot.sign)
    toplot.sign <- toplot.sign[!duplicated(toplot.sign$id),]
    colnames(toplot.sign) <- c("LOCUS", "A1", "BETA", "OR", "SE", "P", "GENOTYPE", "RSID", "MAF")
    d <- seq(1, nrow(toplot.sign))
    a <- split(d, ceiling(seq_along(d)/22))
    for (page in 1:length(a)){
      df <- as.data.frame(a[page])
      #empty thing to be covered -- go to new page
      plot(1, bty='none', xaxt='none', yaxt='none', col='white', xlab='', ylab='')
      grid.table(toplot.sign[df[1, 1]:df[nrow(df), 1], ], rows = NULL)
    }
  }
  return("PDF produced")
}

#function to clean temporary files
function.cleaning <- function(){
  cmd.clean.1 <- "rm list_variants.tmp"
  cmd.clean.2 <- "rm tmp.*"
  cmd.clean.3 <- "rm out.ld.*"
  system(cmd.clean.1)
  system(cmd.clean.2)
  system(cmd.clean.3)
  return("Done!")
}

#function to find genomic coordinates if analysis type is --gene
function.findGeneCoordinates <- function(inp.gene, reference, window.user){
  if (reference == "--hg19"){
    gene.loc <- fread("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/databases/hg19_geneListandPos", h=T, check.names = F)
    
    #find gene of interest, then order according to length
    gene.interest <- gene.loc[which(gene.loc$"#geneName" == inp.gene),]
    
    #check if results were found, otherwise return message to change gene
    if (nrow(gene.interest) == 0){
      print("***NO GENES FOUND, PLEASE CHANGE GENE NAME***")
    } else {
      #order according to length
      gene.interest <- gene.interest[order(-gene.interest$txEnd),]
      
      #selecte chromosome, start and end
      chromosome <- gene.interest$chrom[1]
      start <- gene.interest$txStart[1] - window.user
      end <- gene.interest$txEnd[1] + window.user
      chromosome.splt <- as.data.frame(str_split_fixed(chromosome, "chr", 2))
      chromosome.n <- as.numeric(as.character(chromosome.splt[1, 2]))
    }
  } else if (reference == "--hg38"){
    gene.loc <- fread("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/databases/trial_hg38_genes.txt", h=T, check.names = F)
    
    #find gene of interest, then order according to length
    gene.interest <- gene.loc[which(gene.loc$name2 == inp.gene),]
    
    #check if results were found, otherwise return message to change gene
    if (nrow(gene.interest) == 0){
      print("***NO GENES FOUND, PLEASE CHANGE GENE NAME***")
    } else {
      #order according to length
      gene.interest <- gene.interest[order(-gene.interest$txEnd),]
      
      #selecte chromosome, start and end
      chromosome <- gene.interest$chrom[1]
      start <- gene.interest$txStart[1] - window.user
      end <- gene.interest$txEnd[1] + window.user
      chromosome.splt <- as.data.frame(str_split_fixed(chromosome, "chr", 2))
      chromosome.n <- as.numeric(as.character(chromosome.splt[1, 2]))
    }
    
  }    
  return(list(chromosome.n, start, end))
}

#function to read jasper structural variations -- smaller variations
function.StructuralVar <- function(chromosome, min.value, max.value){
  #read input
  str.var <- fread("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/data/structual_variations/merged_w115_chm1_chm13.table", h=T)
  
  #take region of interest
  str.var.sbs <- str.var[which(str.var$"#reference" == paste("chr", chromosome, sep="")),]
  str.var.pos <- str.var.sbs[which((str.var.sbs$pos_start >= min.value) & (str.var.sbs$pos_end <= max.value)),]
  
  return(str.var.pos)  
}

#function to read jasper structural variations -- smaller variations
function.LargeStructuralVar <- function(chromosome, min.value, max.value){
  #read input
  large.var <- fread("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/data/structual_variations/5hap.all.bed", h=F)
  
  #take region of interest
  large.var.sbs <- large.var[which(large.var$V1 == paste("chr", chromosome, sep="")),]
  large.var.pos <- large.var.sbs[which((large.var.sbs$V2 >= min.value) & (large.var.sbs$V3 <= max.value)),]
  
  return(large.var.pos)  
}

#function to plot png -- hg38
function.plotPNG.hg38 <- function(trial, trial.norm, chromosome, genes, genmap, min.x, max.x, max.y, inp.variant, interval, str.var.pos, 
                                  large.var.pos, window.user, chromatin){
  par(mar=c(6, 9, 6, 9))
  
  #PLOT 1
  #plot recombination map and its relative axis and label
  plot(genmap$bp, genmap$norm.rate, type='l', col='darkolivegreen3', xaxt='none', ylim=c(min(genes$y)-max.y/4/4*9, max.y), xlim=c(min.x, max.x), xlab='Chromosomal Position (Mb)', ylab='-log10(P-value)', yaxt='none', cex.lab=1.40, main=paste("Chromosome ", chromosome, " ~ ", inp.variant, sep = ""), lwd=2)
  axis(side = 4, at = trial.norm$y, labels=trial.norm$x, col='darkolivegreen3', col.axis = "darkolivegreen3", cex.axis=1.25)
  mtext("Recombination rate (cM/Mb)", side=4, line=3, cex=1.40, col='darkolivegreen3', at = ((trial.norm$y[which(trial.norm$x == 60)] + trial.norm$y[which(trial.norm$x == 40)])/2))
  #plot point with low LD
  points(trial$pos.hg38[which(trial$col == "navy")], trial$log10P[which(trial$col == "navy")], pch=trial$pch[which(trial$col == "navy")], bg=trial$col[which(trial$col == "navy")], col=trial$col[which(trial$col == "navy")], cex=trial$size[which(trial$col == "navy")])
  #plot point in LD with my variant
  points(trial$pos.hg38[which(trial$col != "navy")], trial$log10P[which(trial$col != "navy")], pch=trial$pch[which(trial$col != "navy")], bg=trial$col[which(trial$col != "navy")], col=trial$col[which(trial$col != "navy")], cex=trial$size[which(trial$col != "navy")])
  #plot my variant
  points(trial$pos.hg38[which(trial$id == inp.variant)], trial$log10P[which(trial$id == inp.variant)], pch=trial$pch[which(trial$id == inp.variant)], bg = trial$col[which(trial$id == inp.variant)], cex=trial$size[which(trial$id == inp.variant)])
  #plot axes and legend
  axis(side = 1, at=c(min.x, min.x+(interval), min.x+(interval*2), max.x), cex.axis=1.25, labels=c(round(min.x/1000000, 2), round((min.x+(interval))/1000000, 2), round((min.x+(interval*2))/1000000, 2), round((max.x/1000000), 2)))
  axis(side = 2, at=seq(0, round(max.y, 0), round(max.y/4, 0)), labels=seq(0, round(max.y, 0), round(max.y/4, 0)), cex.axis=1.25)
  legend('topright', legend=c(0, 0.2, 0.4, 0.6, 0.8), pch=22, pt.bg=c("navy", "deepskyblue2", "khaki", "orange", "coral"), cex=1.40, ncol=1, title=as.expression(bquote(r^2)) )
  
  #plot genes
  for (g in 1:nrow(genes)){
    segments(x0 = genes$txStart[g], y0 = genes$y[g], x1 = genes$txEnd[g], y1 = genes$y[g], lwd=3)
    start <- as.data.frame(t(str_split_fixed(genes$exonStarts[g], ",", genes$exonCount[g]+1)))
    end <- as.data.frame(t(str_split_fixed(genes$exonEnds[g], ",", genes$exonCount[g]+1)))
    exons <- cbind(start, end)
    colnames(exons) <- c("start", "end")
    exons$start <- as.numeric(as.character(exons$start))
    exons$end <- as.numeric(as.character(exons$end))
    for (j in 1:nrow(exons)){
      rect(xleft=exons$start[j], ybottom=genes$y[g]-(max.y/4/4*0.15), xright=exons$end[j], ytop = genes$y[g]+(max.y/4/4*0.15), col='grey80')
    }
    text(x = genes$txStart[g] + (genes$txEnd[g] - genes$txStart[g])/2, y = genes$y[g] + (max.y/4/4*0.5), labels=genes$name2[g], font=3, cex=0.85)
    if (genes$strand[g] == "+"){
      arrows(x0 = genes$txEnd[g], y0 = genes$y[g], x1 = genes$txEnd[g] + (window.user*2/100), y1 = genes$y[g], length=0.1, lwd=2, col='coral')
    } else if (genes$strand[g] == "-"){
      arrows(x0 = genes$txStart[g], y0 = genes$y[g], x1 = genes$txStart[g] - (window.user*2/100), y1 = genes$y[g], length=0.1, lwd=2, col='coral')
    }
  }

  abline(h=min(genes$y) - (max.y/4/4*0.5), lty=1, lwd=1.50)
  
  #manage structural variations (new-way)
  large.var.pos$pos.plt <- min(genes$y)-max.y/4/4*2.5

  if (nrow(large.var.pos) > 0){
    for (v in 1:nrow(large.var.pos)){
      segments(x0 = large.var.pos$V2[v], y0 = large.var.pos$pos.plt[v], x1 = large.var.pos$V3[v], y1 = large.var.pos$pos.plt[v], lwd=4, col="springgreen4")
      segments(x0 = large.var.pos$V2[v], y0 = large.var.pos$pos.plt[v] + 0.10, x1 = large.var.pos$V2[v], y1 = large.var.pos$pos.plt[v] - 0.10, lwd=4, col='springgreen4')
      segments(x0 = large.var.pos$V3[v], y0 = large.var.pos$pos.plt[v] + 0.10, x1 = large.var.pos$V3[v], y1 = large.var.pos$pos.plt[v] - 0.10, lwd=4, col='springgreen4')
      text(x = ((large.var.pos$V3[v] - large.var.pos$V2[v])/2) + large.var.pos$V2[v], y = large.var.pos$pos.plt[v] + 0.30, labels = paste(large.var.pos$V4[v]), xpd=T, cex=0.80, adj = 0.5)
    }
  }

  mtext("Structural\nvariants", side=2, line=1, cex=1.10, col='black', at = min(genes$y)-max.y/4/4*2.5, adj = 0.5)
  #legend('topleft', legend=c("indels", "region"), lty = 1, lwd=5, col=c("blue", "red"), cex=1.40, ncol=1, title="SVs")
  
  #new line before chromatin states
  abline(h=min(genes$y) - (max.y/4/4*4.5), lty=1, lwd=1.50)
  
  #add plotting position
  chromatin$pos.plt <- min(genes$y)-max.y/4/4*6.5

  #now add chromatin states
  chrom.colors <- str_split_fixed(chromatin$V5, ",", 3)
  chromatin$red <- as.numeric(chrom.colors[, 1])
  chromatin$green <- as.numeric(chrom.colors[, 2])
  chromatin$blue <- as.numeric(chrom.colors[, 3])
  for (row in 1:nrow(chromatin)){
    rect(xleft = chromatin$V2[row], ybottom = chromatin$pos.plt[row], xright = chromatin$V3[row], ytop = chromatin$pos.plt[row] - 0.25, 
         col=rgb(red = chromatin$red[row], green = chromatin$green[row], blue = chromatin$blue[row], maxColorValue = 255), border = NA)
  }
  
  #read chromatin legend file
  chrom.legend <- read.table("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/databases/chromatin_states/chromLegend.txt", h=T, sep='\t', check.names = F)
  chrom.legend$colors <- paste("#", chrom.legend$colo, sep="")
  mtext("Chromatin\nstates", side=2, line=1, cex=1.10, col='black', at = min(genes$y)-max.y/4/4*7.5, adj = 0.5)
  legend("bottom", legend = chrom.legend$state, pch=22, pt.bg = chrom.legend$colors, bty='n', ncol=nrow(chrom.legend)/2, cex=0.80, pt.cex = 1.25)
  
  return("PNG produced")
}

#function to plot pdf -- hg38
function.plotPDF.hg38 <- function(trial, trial.norm, chromosome, genes, genmap, min.x, max.x, max.y, inp.variant, interval, str.var.pos, 
                                  large.var.pos, window.user){
  par(mar=c(6, 9, 6, 9))
  
  ##PLOT 1 -- pdf
  #plot recombination map and its relative axis and label
  plot(genmap$bp, genmap$norm.rate, type='l', col='darkolivegreen3', xaxt='none', ylim=c(min(genes$y)-max.y/4/4*8, max.y), xlim=c(min.x, max.x), xlab='Chromosomal Position (Mb)', ylab='-log10(P-value)', yaxt='none', cex.lab=1.40, main=paste("Chromosome ", chromosome, " ~ ", inp.variant, sep = ""), lwd=2)
  axis(side = 4, at = trial.norm$y, labels=trial.norm$x, col='darkolivegreen3', col.axis = "darkolivegreen3", cex.axis=1.25)
  mtext("Recombination rate (cM/Mb)", side=4, line=3, cex=1.40, col='darkolivegreen3', at = ((trial.norm$y[which(trial.norm$x == 60)] + trial.norm$y[which(trial.norm$x == 40)])/2))
  #plot point with low LD
  points(trial$pos.hg38[which(trial$col == "navy")], trial$log10P[which(trial$col == "navy")], pch=trial$pch[which(trial$col == "navy")], bg=trial$col[which(trial$col == "navy")], col=trial$col[which(trial$col == "navy")], cex=trial$size[which(trial$col == "navy")])
  #plot point in LD with my variant
  points(trial$pos.hg38[which(trial$col != "navy")], trial$log10P[which(trial$col != "navy")], pch=trial$pch[which(trial$col != "navy")], bg=trial$col[which(trial$col != "navy")], col=trial$col[which(trial$col != "navy")], cex=trial$size[which(trial$col != "navy")])
  #plot my variant
  points(trial$pos.hg38[which(trial$id == inp.variant)], trial$log10P[which(trial$id == inp.variant)], pch=trial$pch[which(trial$id == inp.variant)], bg = trial$col[which(trial$id == inp.variant)], cex=trial$size[which(trial$id == inp.variant)])
  #plot axes and legend
  axis(side = 1, at=c(min.x, min.x+(interval), min.x+(interval*2), max.x), cex.axis=1.25, labels=c(round(min.x/1000000, 2), round((min.x+(interval))/1000000, 2), round((min.x+(interval*2))/1000000, 2), round((max.x/1000000), 2)))
  axis(side = 2, at=seq(0, round(max.y, 0), round(max.y/4, 0)), labels=seq(0, round(max.y, 0), round(max.y/4, 0)), cex.axis=1.25)
  legend('topright', legend=c(0, 0.2, 0.4, 0.6, 0.8), pch=22, pt.bg=c("navy", "deepskyblue2", "khaki", "orange", "coral"), cex=1.40, ncol=1, title=as.expression(bquote(r^2)) )
  
  #plot genes
  for (g in 1:nrow(genes)){
    segments(x0 = genes$txStart[g], y0 = genes$y[g], x1 = genes$txEnd[g], y1 = genes$y[g], lwd=3)
    start <- as.data.frame(t(str_split_fixed(genes$exonStarts[g], ",", genes$exonCount[g]+1)))
    end <- as.data.frame(t(str_split_fixed(genes$exonEnds[g], ",", genes$exonCount[g]+1)))
    exons <- cbind(start, end)
    colnames(exons) <- c("start", "end")
    exons$start <- as.numeric(as.character(exons$start))
    exons$end <- as.numeric(as.character(exons$end))
    for (j in 1:nrow(exons)){
      rect(xleft=exons$start[j], ybottom=genes$y[g]-(max.y/4/4*0.15), xright=exons$end[j], ytop = genes$y[g]+(max.y/4/4*0.15), col='grey80')
    }
    text(x = genes$txStart[g] + (genes$txEnd[g] - genes$txStart[g])/2, y = genes$y[g] + (max.y/4/4*0.5), labels=genes$name2[g], font=3, cex=0.85)
    if (genes$strand[g] == "+"){
      arrows(x0 = genes$txEnd[g], y0 = genes$y[g], x1 = genes$txEnd[g] + (window.user*2/100), y1 = genes$y[g], length=0.1, lwd=2, col='coral')
    } else if (genes$strand[g] == "-"){
      arrows(x0 = genes$txStart[g], y0 = genes$y[g], x1 = genes$txStart[g] - (window.user*2/100), y1 = genes$y[g], length=0.1, lwd=2, col='coral')
    }
  }
  
  #now put structural variants in there
  #assign position where to plot
  str.var.pos$pos.plt <- 0
  str.var.pos$col <- "blue"
  str.var.pos$col[which(str.var.pos$type == "region")] <- "red"
  str.var.pos$pos.plt[which(str.var.pos$pheno == "w115.h1")] <- (min(genes$y) - (max.y/4/4*1.75))
  str.var.pos$pos.plt[which(str.var.pos$pheno == "w115.h2")] <- (min(genes$y) - (max.y/4/4*3))
  str.var.pos$pos.plt[which(str.var.pos$pheno == "chm1")] <- (min(genes$y) - (max.y/4/4*4.5))
  str.var.pos$pos.plt[which(str.var.pos$pheno == "cmh13")] <- (min(genes$y) - (max.y/4/4*6))
  if (nrow(str.var.pos) > 0){
    for (k in 1:nrow(str.var.pos)){
      segments(x0 = str.var.pos$pos_start[k], y0 = str.var.pos$pos.plt[k], x1 = str.var.pos$pos_end[k], y1 = str.var.pos$pos.plt[k], lwd=3, col=str.var.pos$col[k])
      #text(x = str.var.pos$pos_start[k], y = str.var.pos$pos.plt[k]+max.y/4/4*0.50, labels = paste(str.var.pos$diff_allele_size[k]), xpd=T, cex=0.50)
    }
  }
  #put labels and lines to divide
  abline(h=min(genes$y) - (max.y/4/4*0.5), lty=1, lwd=0.8)
  abline(h=min(genes$y) - (max.y/4/4*2), lty=2, lwd=0.4)
  abline(h=min(genes$y) - (max.y/4/4*3.25), lty=2, lwd=0.4)
  abline(h=min(genes$y) - (max.y/4/4*4.75), lty=2, lwd=0.4)
  abline(h=min(genes$y) - (max.y/4/4*6.25), lty=2, lwd=0.4)
  axis(side = 2, at = c(min(genes$y)-(max.y/4/4*1.75), min(genes$y)-(max.y/4/4*3), min(genes$y)-(max.y/4/4*4.5), min(genes$y)-(max.y/4/4*6)), labels = c("w115 haplo 1", "w115 haplo 2", "chm 1", "chm 13"), tick=F, las=1, cex.axis=0.90, xpd=T)
  
  #other structual variations
  large.var.pos$pos.plt <- min(genes$y)-max.y/4/4*7.5
  if (nrow(large.var.pos) > 0){
    for (v in 1:nrow(large.var.pos)){
      segments(x0 = large.var.pos$V2[v], y0 = large.var.pos$pos.plt[v], x1 = large.var.pos$V3[v], y1 = large.var.pos$pos.plt[v], lwd=3, col="darkolivegreen3")
      #text(x = large.var.pos$V2[v], y = large.var.pos$pos.plt[v]+max.y/4/4*0.5, labels = paste(large.var.pos$V4[v]), xpd=T, cex=0.50)
    }
  }
  abline(h=min(genes$y) - (max.y/4/4*7.75), lty=1, lwd=0.8)
  
  mtext("Structural\nvariants (SVs)", side=4, line=1, cex=1.10, col='black', at = -((abs(min(genes$y) - 2.5) + abs(min(genes$y) - 1))/2))      
  legend('topleft', legend=c("indels", "region"), lty = 1, lwd=5, col=c("blue", "red"), cex=1.40, ncol=1, title="SVs")
  
  ##PLOT 2
  #if group is igap, table should be slightly different
  if (groups == "--igap"){
    #prepare table
    toplot <- trial[, c("pos.hg38", "a1", "beta", "se", "pvalue", "rsID")]
    toplot.sign <- toplot[which(toplot$pvalue <= 0.05),]
    toplot.sign <- toplot.sign[order(toplot.sign$pvalue),]
    toplot.sign <- as.data.frame(toplot.sign)
    toplot.sign <- na.omit(toplot.sign)
    toplot.sign <- toplot.sign[!duplicated(toplot.sign$pos.hg38),]
    colnames(toplot.sign) <- c("POS", "A1", "BETA", "SE", "P", "RSID")
    d <- seq(1, nrow(toplot.sign))
    a <- split(d, ceiling(seq_along(d)/22))
    for (page in 1:length(a)){
      df <- as.data.frame(a[page])
      #empty thing to be covered -- go to new page
      plot(1, bty='none', xaxt='none', yaxt='none', col='white', xlab='', ylab='')
      grid.table(toplot.sign[df[1, 1]:df[nrow(df), 1], ], rows = NULL)
    }
    
  } else {
    
    #prepare table
    toplot <- trial[, c("pos.hg38", "a1", "beta", "odds.ratio", "se", "pvalue", "genotype", "rsID", "maf")]
    toplot.sign <- toplot[which(toplot$pvalue <= 0.05),]
    toplot.sign <- toplot.sign[order(toplot.sign$pvalue),]
    toplot.sign <- as.data.frame(toplot.sign)
    toplot.sign <- na.omit(toplot.sign)
    toplot.sign <- toplot.sign[!duplicated(toplot.sign$pos.hg38),]
    colnames(toplot.sign) <- c("POS", "A1", "BETA", "OR", "SE", "P", "GENOTYPE", "RSID", "MAF")
    d <- seq(1, nrow(toplot.sign))
    a <- split(d, ceiling(seq_along(d)/22))
    for (page in 1:length(a)){
      df <- as.data.frame(a[page])
      #empty thing to be covered -- go to new page
      plot(1, bty='none', xaxt='none', yaxt='none', col='white', xlab='', ylab='')
      grid.table(toplot.sign[df[1, 1]:df[nrow(df), 1], ], rows = NULL)
    }
  }
  return("PDF produced")
}

#function to read and save chromatin states information
function.parseChromatinStates <- function(chromosome, min.value, max.value, genome, window.user){
  #read input
  if (genome == "--hg19"){
    chromatin <- fread("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/databases/chromatin_states/wgEncodeBroadHmmGm12878HMM_hg19.bed", h=F, sep='\t')
  } else if (genome == "--hg38"){
    chromatin <- fread("/tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/script/tools/locuzNicco/databases/chromatin_states/wgEncodeBroadHmmGm12878HMM_hg38.bed", h=F, sep='\t')
  }
  
  #take region of interest
  chromatin.chr <- chromatin[which(chromatin$V1 == paste("chr", chromosome, sep="")),]
  chromatin.pos <- chromatin.chr[which((as.numeric(chromatin.chr$V2) >= as.numeric(min.value - window.user)) & (as.numeric(chromatin.chr$V3) <= as.numeric(max.value + window.user))),]
  
  return(chromatin.pos)
}

##########

#READ FIRST ARGUMENT -- ALWAYS AS FIRST THING
groups <- args[1]
##########

#DISPLAY HELP MESSAGE IN CASE FIRST ARGUMENT IS --help, OTHERWISE JUST LOAD ALL ARGUMENTS
if (groups == "--help"){
  print("Arg[1]: groups: --chc-lasa / --chc-lasa-scd / --igap / --ad / --grace / --igap-2019 / --spigap / --spigap-ukb")
  print("Arg[2]: genome version (hg19 or hg38) --> --hg19 / --hg38 / --hla")
  print("Arg[3]: analysis mode: variant (var) or gene (gene) --> --var --gene")
  print("Arg[4]: length of window to use (in Kb) upstream OR downstream --> 150 is 150Kb upstream and downstream")
  print("Arg[5]: input variant (chr:pos, hg19 coordinates) or gene (gene name) --> chr:pos / gene")
  print("Arg[6]: name of pdf output --> output name")
  print("**NB: --hla genome uses by default GRCh38, linkage analysis is not implemented yet")
  run.script = FALSE
} else {
  genome <- args[2]
  type <- args[3]
  window.user <- as.numeric(args[4])*1000
  target <- args[5]
  output <- args[6]
  run.script = TRUE
}
##########

#MAIN SCRIPT
#if there are all arguments, run the script
if (run.script == TRUE){
  print("*****ANALYSIS STARTED")
  print(paste("*****GROUPS: ", groups, sep=""))
  print(paste("*****GENOME: ", genome, sep=""))
  print(paste("*****TYPE: ", type, sep=""))
  print(paste("*****WINDOW: ", window.user, sep=""))
  print(paste("*****TARGET: ", target, sep=""))
  print(paste("*****OUTPUT: ", output, sep=""))
  print("************")
  print("************")
  
  #reference genome hg19
  if (genome == "--hg19"){
    
    #type of analysis is variant - based
    if (type == "--var"){
      print("*****EXTRACTING VARIANT INFORMATION, LOADING GWAS AND PERFORMING LD ANALYSIS ")
      
      #extract chromosome and position
      chr.pos <- function.extractChromosomeAndPos(target)
      
      #take GWAS of interest
      inputFile.path <- function.takeGWASofInterest(groups, genome, chr.pos$chr)
      
      #check whether user wants IGAP data, in case call the relative function
      if (groups %in% c("--igap", "--igap-2019")){
        assoc <- function.manageIGAP(groups, genome, inputFile.path, chr.pos$chr)
        
        #compute LD for the variants in 150 kb up/down-stream and merge association results with LD results
        trial <- function.variantLD(assoc, chr.pos, target, window.user)
        
      } else {
        assoc <- fread(inputFile.path, h=T)
        
        #if the group is GRACE, need to change column names as otherwise there will be problems
        if (groups %in% c("--grace", "--spigap", "spigap-ukb")){
          colnames(assoc) <- c("LOCID", "id", "chrom", "pos", "a1", "a2", "beta", "se", "pvalue", "freq", "or", "ci", "source")
          
          #select chromosome of interest
          assoc <- assoc[which(assoc$chrom == chr.pos$chr),]
        }

        #compute LD for the variants in 150 kb up/down-stream and merge association results with LD results
        trial <- function.variantLD(assoc, chr.pos, target, window.user)
        
        if (!(groups %in% c("--grace", "--spigap", "--spigap-ukb"))){
          #combine with frquency data -- subset to take loci where maf > 0.005
          trial <- function.combineFrequencies(trial, inp.variant, chr.pos$chr)
        }
      }
      
      print("*****PREPARING FOR PLOTTING")
      
      #run function to prepare for plotting
      results.function <- function.plotPreparation(trial, target, genome)
      trial <- as.data.frame(results.function[1])
      min.x <- as.numeric(results.function[2])
      max.x <- as.numeric(results.function[3])
      max.y <- as.numeric(results.function[4])
      interval <- as.numeric(results.function[5])
      
      #run function to grep region of interest in gene location file
      genes <- function.parseGeneLocation(genome, chr.pos$chr, min.x, max.x)
      
      #run function to grep region of interest in recombination peaks file
      results.function <- function.parseGeneticMap(genome, chr.pos$chr, min.x, max.x, max.y)
      genmap <- as.data.frame(results.function[1])
      trial.norm <- as.data.frame(results.function[2])
      
      #take chromatin information
      chromatin <- function.parseChromatinStates(chr.pos$chr, min.x, max.x, genome, window.user)
      
      print("*****PLOTTING AND CLEANING")
      #plot png
      out.name <- paste(output, ".png", sep="")
      png(out.name, height = 10, width=14, units='in', res=400)
      function.plotPNG(trial, trial.norm, chr.pos, genes, genmap, min.x, max.x, max.y, target, interval, chromatin)
      dev.off()
      
      #plot pdf
      #out.name <- paste(output, ".pdf", sep="")
      #pdf(out.name, height = 10, width=14)
      #function.plotPDF(trial, trial.norm, chr.pos, genes, genmap, min.x, max.x, max.y, target, interval, groups)
      #dev.off()
      
      #clean temporary files
      function.cleaning()
    } else if (type == "--gene"){
      
      print("*****EXTRACTING GENE INFORMATION, LOADING GWAS AND PERFORMING LD ANALYSIS ")
      
      #find gene coordinates
      results.function <- function.findGeneCoordinates(target, genome, window.user)
      chromosome.n <- as.numeric(results.function[1])
      start <- as.numeric(results.function[2])
      end <- as.numeric(results.function[3])
      
      #take GWAS of interest
      inputFile.path <- function.takeGWASofInterest(groups, genome, chromosome.n)
      print(inputFile.path)
      #check whether user wants IGAP data, in case call the relative function
      if (groups %in% c("--igap", "--igap-2019")){
        assoc <- function.manageIGAP(groups, genome, inputFile.path, chromosome.n)
        
        #compute LD for the variants in 150 kb up/down-stream and merge association results with LD results
        function.results <- function.variantLD.Gene(chromosome.n, assoc, start, end, genome)
        inp.variant <- as.character(function.results[1])
        trial <- as.data.frame(function.results[2])
        
      } else {
        assoc <- fread(inputFile.path, h=T)
        
        #compute LD for the variants in 150 kb up/down-stream and merge association results with LD results
        function.results <- function.variantLD.Gene(chromosome.n, assoc, start, end, genome)
        inp.variant <- as.character(function.results[1])
        trial <- as.data.frame(function.results[2])
        
        #combine with frquency data -- subset to take loci where maf > 0.005
        trial <- function.combineFrequencies(trial, inp.variant, chromosome.n)
      }
      
      print("*****PREPARING FOR PLOTTING")
      
      #run function to prepare for plotting
      results.function <- function.plotPreparation(trial, inp.variant, genome)
      trial <- as.data.frame(results.function[1])
      min.x <- as.numeric(results.function[2])
      max.x <- as.numeric(results.function[3])
      max.y <- as.numeric(results.function[4])
      interval <- as.numeric(results.function[5])
      
      #run function to grep region of interest in gene location file
      genes <- function.parseGeneLocation(genome, chromosome.n, min.x, max.x)
      
      #run function to grep region of interest in recombination peaks file
      results.function <- function.parseGeneticMap(genome, chromosome.n, min.x, max.x, max.y)
      genmap <- as.data.frame(results.function[1])
      trial.norm <- as.data.frame(results.function[2])
      
      #take chromatin information
      chromatin <- function.parseChromatinStates(chromosome.n, min.x, max.x, genome, window.user)
      
      print("*****PLOTTING AND CLEANING")
      #plot png
      out.name <- paste(output, ".png", sep="")
      png(out.name, height = 10, width=14, units='in', res=400)
      function.plotPNG(trial, trial.norm, chromosome.n, genes, genmap, min.x, max.x, max.y, target, interval, chromatin)
      dev.off()
      
      #plot pdf
      #out.name <- paste(output, ".pdf", sep="")
      #pdf(out.name, height = 10, width=14)
      #function.plotPDF(trial, trial.norm, chromosome.n, genes, genmap, min.x, max.x, max.y, target, interval, groups)
      #dev.off()
      
      #clean temporary files
      function.cleaning()
    }
    
  } else if (genome == "--hg38"){
    
    #type of analysis is variant - based
    if (type == "--var"){
      
      print("*****EXTRACTING VARIANT INFORMATION, LOADING GWAS AND PERFORMING LD ANALYSIS ")
      
      #extract chromosome and position
      chr.pos <- function.extractChromosomeAndPos(target)
      
      #take GWAS of interest
      inputFile.path <- function.takeGWASofInterest(groups, genome, chr.pos$chr)
      
      #check whether user wants IGAP data, in case call the relative function
      if (groups == "--igap"){
        assoc <- function.manageIGAP(genome, inputFile.path)
        
        #compute LD for the variants in 150 kb up/down-stream and merge association results with LD results
        trial <- function.variantLD.hg38(assoc, chr.pos, target, groups, window.user)
        
      } else {
        assoc <- fread(inputFile.path, h=T)
        
        #compute LD for the variants in 150 kb up/down-stream and merge association results with LD results
        trial <- function.variantLD.hg38(assoc, chr.pos, target, groups, window.user)
        
        #combine with frquency data -- subset to take loci where maf > 0.005
        trial <- function.combineFrequencies(trial, inp.variant, chr.pos$chr)
      }
      
      print("*****PREPARING FOR PLOTTING")
      
      #run function to prepare for plotting
      results.function <- function.plotPreparation(trial, target, genome)
      trial <- as.data.frame(results.function[1])
      min.x <- as.numeric(results.function[2])
      max.x <- as.numeric(results.function[3])
      max.y <- as.numeric(results.function[4])
      interval <- as.numeric(results.function[5])
      
      #run function to grep region of interest in gene location file
      genes <- function.parseGeneLocation(genome, chr.pos$chr, min.x, max.x)
      
      #run function to grep region of interest in recombination peaks file
      results.function <- function.parseGeneticMap(genome, chr.pos$chr, min.x, max.x, max.y)
      genmap <- as.data.frame(results.function[1])
      trial.norm <- as.data.frame(results.function[2])
      
      #take information from structural variants and large structural variants
      sv <- function.StructuralVar(chr.pos$chr, min.x, max.x)
      large.sv <- function.LargeStructuralVar(chr.pos$chr, min.x, max.x)
      
      #take chromatin information
      chromatin <- function.parseChromatinStates(chr.pos$chr, min.x, max.x, genome, window.user)

      print("*****PLOTTING AND CLEANING")
      
      #plot png
      out.name <- paste(output, ".png", sep="")
      png(out.name, height = 10, width=14, units='in', res=400)
      function.plotPNG.hg38(trial, trial.norm, chr.pos$chr, genes, genmap, min.x, max.x, max.y, target, interval, sv, large.sv, window.user, chromatin)
      dev.off()
      
      #plot pdf
      out.name <- paste(output, ".pdf", sep="")
      pdf(out.name, height = 10, width=14)
      function.plotPDF.hg38(trial, trial.norm, chr.pos$chr, genes, genmap, min.x, max.x, max.y, target, interval, sv, large.sv, window.user)
      dev.off()
      
      #clean temporary files
      function.cleaning()
    } else if (type == "--gene"){
      
      print("*****EXTRACTING GENE INFORMATION, LOADING GWAS AND PERFORMING LD ANALYSIS ")
      
      #find gene coordinates
      results.function <- function.findGeneCoordinates(target, genome, window.user)
      chromosome.n <- as.numeric(results.function[1])
      start <- as.numeric(results.function[2])
      end <- as.numeric(results.function[3])
      
      #take GWAS of interest
      inputFile.path <- function.takeGWASofInterest(groups, genome, chromosome.n)
      
      #check whether user wants IGAP data, in case call the relative function
      if (groups == "--igap"){
        assoc <- function.manageIGAP(genome, inputFile.path)
        
        #compute LD for the variants in 150 kb up/down-stream and merge association results with LD results
        function.results <- function.variantLD.Gene(chromosome.n, assoc, start, end, genome)
        inp.variant <- as.character(function.results[1])
        trial <- as.data.frame(function.results[2])
        
      } else {
        assoc <- fread(inputFile.path, h=T)
        
        #compute LD for the variants in 150 kb up/down-stream and merge association results with LD results
        function.results <- function.variantLD.Gene(chromosome.n, assoc, start, end, genome)
        inp.variant <- as.character(function.results[1])
        trial <- as.data.frame(function.results[2])
        
        #combine with frquency data -- subset to take loci where maf > 0.005
        trial <- function.combineFrequencies(trial, inp.variant, chromosome.n)
      }

      print("*****PREPARING FOR PLOTTING")
      
      #run function to prepare for plotting
      results.function <- function.plotPreparation(trial, inp.variant, genome)
      trial <- as.data.frame(results.function[1])
      min.x <- as.numeric(results.function[2])
      max.x <- as.numeric(results.function[3])
      max.y <- as.numeric(results.function[4])
      interval <- as.numeric(results.function[5])
      
      #run function to grep region of interest in gene location file
      genes <- function.parseGeneLocation(genome, chromosome.n, min.x, max.x)
      
      #run function to grep region of interest in recombination peaks file
      results.function <- function.parseGeneticMap(genome, chromosome.n, min.x, max.x, max.y)
      genmap <- as.data.frame(results.function[1])
      trial.norm <- as.data.frame(results.function[2])
      
      #take information from structural variants and large structural variants
      sv <- function.StructuralVar(chromosome.n, min.x, max.x)
      large.sv <- function.LargeStructuralVar(chromosome.n, min.x, max.x)

      #take chromatin information
      chromatin <- function.parseChromatinStates(chromosome.n, min.x, max.x, genome, window.user)
      print(max.x)
      print(tail(chromatin))
      print("*****PLOTTING AND CLEANING")
      
      #plot png
      out.name <- paste(output, ".png", sep="")
      png(out.name, height = 10, width=14, units='in', res=400)
      function.plotPNG.hg38(trial, trial.norm, chromosome.n, genes, genmap, min.x, max.x, max.y, inp.variant, interval, sv, large.sv, window.user, chromatin)
      dev.off()
      
      #plot pdf
      out.name <- paste(output, ".pdf", sep="")
      pdf(out.name, height = 10, width=14)
      function.plotPDF.hg38(trial, trial.norm, chromosome.n, genes, genmap, min.x, max.x, max.y, target, interval, sv, large.sv, window.user)
      dev.off()
      
      #clean temporary files
      function.cleaning()
    }
  }
}
###########



