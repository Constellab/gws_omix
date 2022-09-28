#!/usr/bin/env Rscript

# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

# bash script salmon_index_cmd.sh :
#
# This script run salmon index program in salmonIndex method build_command()
# Arguments are given inside the build_command() in the cmd array
#
library("DESeq2")
library("tximport")

# Files importing
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) { # test if there is at least one argument: if not, return an error
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else {
  # Default output file
  input_mapping_file = args[1] # raw count table
  metadata_file = args[2] # metadata file
  #metadata_col = args[3] # metadata column name to be used for DE
  additional_design = paste(args[3])
  ##### DESIGN #####
  #Design = paste(" ~ ",metadata_col,additional_design) # " + " # see: https://cran.r-project.org/doc/manuals/R-intro.html#Formulae-for-statistical-models
  Design = paste(" ~ ",additional_design) # " + " # see: https://cran.r-project.org/doc/manuals/R-intro.html#Formulae-for-statistical-models
  output_file = args[4] # output file(s)
}

# Create tables
exp_design_initi = read.csv(metadata_file,header=TRUE,stringsAsFactors=F, sep="\t",fileEncoding="latin1")
salmon_raw_count_matrix_initi <- read.delim(input_mapping_file ,header=TRUE,stringsAsFactors=F, check.names = FALSE, row.names=1 , sep="\t",fileEncoding="latin1")
exp_design = exp_design_initi
rownames(exp_design) <- exp_design_initi[,1]
salmon_raw_count_matrix = salmon_raw_count_matrix_initi

# Get all existing groups in the choosen metadata column
uniq <- unique(exp_design[,metadata_col])

# Pairwaise comparison DESeq2 analysis: all groups versus all groups 
mylist <- c()
print("performing differential analysis")
for (i in 1:length(uniq)){
  tmp_list <- c(i)
  mylist <- append(mylist,tmp_list)
  for (j in 1:length(uniq)){
    if( i != j ){ 
      if( j %in% mylist == TRUE){
        next
      }
      else{ # create the sub-matrix to compare the group
        tmp_exp_matrix <- subset(exp_design, exp_design[,metadata_col]==uniq[i])
        tmp_exp_matrix2 <- subset(exp_design, exp_design[,metadata_col]==uniq[j])
        matrix_exp_combine <- rbind(tmp_exp_matrix,tmp_exp_matrix2)
        tmp_count_matrix<-c()
        list_row <- rownames(matrix_exp_combine)
        for (sample in list_row){
          tmp_count_matrix<- cbind(tmp_count_matrix, salmon_raw_count_matrix[[sample]])         
        }
        rownames(tmp_count_matrix) <- rownames(salmon_raw_count_matrix)
        colnames(tmp_count_matrix) <- list_row
        # DESeq2: differential expression analysis
        DESeq2_matrix <- DESeqDataSetFromMatrix(countData = round(tmp_count_matrix), colData = matrix_exp_combine, design = as.formula(Design))
        dds <- DESeq(DESeq2_matrix)
        res <- results(dds)
        write.table(res, file = paste(output_file,uniq[i],"vs",uniq[j],metadata_col,"and",additional_design,"complete","txt", sep='.'), quote=F, sep="\t")
        #res <- subset(res, res$padj < 0.05) # adjust p-value threshold to only keep differentaily expressed genes in a given condition
        # Output files: writing results according to the following filter: adjusted p-value from the Wald test
        #write.table(res, file = paste(output_file,uniq[i],"vs",uniq[j],metadata_col,"and",additional_design,"padj_0.05","txt", sep='.'), quote=F, sep="\t")
      }  
    }
  }
}
