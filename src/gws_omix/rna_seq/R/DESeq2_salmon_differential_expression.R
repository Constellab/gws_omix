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

args = commandArgs(trailingOnly=TRUE)

# Files importing

if (length(args)==0) { # test if there is at least one argument: if not, return an error
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else {
  # default output file
  input_mapping_file = args[1] 
  metadata_file = args[2]
  metadata_col = args[3]
  Design = paste(" ~ ",args[3])
  output_file = args[4]
}

# create tables

exp_design = read.csv(metadata_file,header=TRUE,stringsAsFactors=F,row.names=1, sep="\t")
salmon_raw_count_matrix <- read.delim(input_mapping_file ,header=TRUE,stringsAsFactors=F,row.names=1, sep="\t")

#get all existing groups in the choosen metadata column

uniq <- unique(exp_design[,metadata_col])


# Pairwaise comparison DESeq2 analysis: all groups versus all groups 
mylist <- c()

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
        res1 <- subset(res, res$padj < 0.1)
        res05 <- subset(res, res$padj < 0.05)
        res01 <- subset(res, res$padj < 0.01)

        # Output files: writing results according to the following filter: adjusted p-value from the Wald test
        write.table(res, file = paste(output_file,uniq[i],"versus",uniq[j],"all","txt", sep='.'), quote=F, sep="\t")
        write.table(res1, file = paste(output_file,uniq[i],"versus",uniq[j],"padj_0.1","txt", sep='.'), quote=F, sep="\t")
        write.table(res05, file = paste(output_file,uniq[i],"versus",uniq[j],"padj_0.05","txt", sep='.'), quote=F, sep="\t")
        write.table(res01, file = paste(output_file,uniq[i],"versus",uniq[j],"padj_0.01","txt", sep='.'), quote=F, sep="\t")       
      }  
    }
  }
}
