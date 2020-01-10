###
#   File name : MakeRDA.R
#   Author    : Hyunjin Kim
#   Date      : Jan 8, 2020
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Make analysis-ready RDA file from Suciu-Foca Lab's raw data
#
#   Instruction
#               1. Source("MakeRDA.R")
#               2. Run the function "makeRDA" - specify the input files (raw counts and sample info)
#               3. The RDA file will be generated under the raw count files directory
#
#   Example
#               > source("The_directory_of_MakeRDA.R/MakeRDA.R")
#               > makeRDA(rCntDir="C:/Research/CUMC/Suciu-Foca_RNASeq/data/",
#                         sampleInfoPath="C:/Research/CUMC/Suciu-Foca_RNASeq/data/samples_matrix.xlsx")
###

# * All he raw count txt files for each sample should be located under the "rCntDir"

makeRDA <- function(rCntDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/suciu_foca/data/raw_counts_only/",
                    sampleInfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/suciu_foca/data/samples_matrix.xlsx") {
  
  ### load library
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(DESeq2, quietly = TRUE)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("DESeq2")
    require(DESeq2, quietly = TRUE)
  }
  
  ### get raw count files
  f <- list.files(path = rCntDir, pattern = "\\.txt$")
  
  ### extract sample names from the file names
  sample_names <- sapply(f, function(x) strsplit(x, split = ".", fixed = TRUE)[[1]][1])
  
  ### load the files
  rCnts <- vector("list", length = length(f))
  names(rCnts) <- sample_names
  for(file in f) {
    ### load the file
    rCnts[[sample_names[file]]] <- read.table(file = paste0(rCntDir, file),
                                              header = TRUE, sep = "\t", row.names = 1,
                                              stringsAsFactors = FALSE, check.names = FALSE)
  }
  
  ### check every entry has the same number of genes
  ### check every entry has the genes in same order
  for(i in 1:(length(rCnts)-1)) {
    if(!identical(rownames(rCnts[[i]]), rownames(rCnts[[i+1]]))) {
      writeLines(paste(names(rCnts)[i], "and", names(rCnts)[i+1], "do not have same genes."))
    }
  }
  
  ### combine the raw counts into one matrix
  raw_counts <- matrix(0, nrow(rCnts[[1]]), length(rCnts))
  rownames(raw_counts) <- rownames(rCnts[[1]])
  colnames(raw_counts) <- sample_names
  for(sample in sample_names) {
    for(gene in rownames(rCnts[[sample]])) {
      raw_counts[,sample] <- rCnts[[sample]][,sample]
    }
  }
  
  ### load the sample info
  sample_info_file <- read.xlsx2(file = sampleInfoPath, sheetIndex = 1, startRow = 7, row.names = 1,
                                 stringsAsFactors = FALSE, check.names = FALSE)
  
  ### make sample info with the info file
  sample_info <- matrix("", ncol(raw_counts), 4)
  rownames(sample_info) <- colnames(raw_counts)
  colnames(sample_info) <- c("Condition", "Time_Point", "Batch", "Technical_Replicates")
  sample_info[,"Condition"] <- rep(c("Untreated", "Low_Concentration", "High_Concentration"), 6)
  sample_info[,"Time_Point"] <- c(rep("24", 6), rep("6", 6), rep("12", 6))
  sample_info[,"Batch"] <- c(rep("Batch1", 6), rep("Batch2", 12))
  sample_info[,"Technical_Replicates"] <- c(rep(c("T1", "T2", "T3"), 2), rep(c("T4", "T5", "T6"), 2), rep(c("T7", "T8", "T9"), 2))
  
  ### sum up the counts for the same technical replicates
  deSeqData <- DESeqDataSetFromMatrix(countData=raw_counts, colData=sample_info, design= ~1)
  deSeqData <- collapseReplicates(deSeqData, factor(sample_info[,"Technical_Replicates"]))
  collapsed_raw_counts <- counts(deSeqData)
  
  ### make new sample info for the collapsed raw counts
  collapsed_sample_info <- matrix("", ncol(collapsed_raw_counts), 3)
  rownames(collapsed_sample_info) <- colnames(collapsed_raw_counts)
  colnames(collapsed_sample_info) <- c("Condition", "Time_Point", "Batch")
  collapsed_sample_info[,"Condition"] <- rep(c("Untreated", "Low_Concentration", "High_Concentration"), 3)
  collapsed_sample_info[,"Time_Point"] <- c(rep("24", 3), rep("6", 3), rep("12", 3))
  collapsed_sample_info[,"Batch"] <- c(rep("Batch1", 3), rep("Batch2", 6))
  
  ### set README function
  README <- function() {
    writeLines(paste(rep("#", 100), collapse = ""))
    writeLines("The \"raw_counts\" and \"sample_info\" are for original raw counts,")
    writeLines("and the \"collapsed_raw_counts\" and \"collapsed_sample_info\" are")
    writeLines("for summed-up version of the original ones.")
    writeLines("Because there are technical replicates, raw counts of the same replicates are")
    writeLines("summed up by using collapseReplicates() function of DESeq2 package.")
    writeLines(paste(rep("#", 100), collapse = ""))
  }
  
  ### save the raw counts and sample info as one RDA file
  save(list = c("raw_counts", "sample_info", "collapsed_raw_counts", "collapsed_sample_info", "README"),
       file = paste0(rCntDir, "raw_counts_suciu-foca.rda"))
  
}
