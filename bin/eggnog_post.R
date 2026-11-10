#!/usr/bin/env Rscript
############################
#######    EGGNOG   ########
############################

library(stringr)

args <- commandArgs(trailingOnly = TRUE)

file_path <- args[1]

if (is.na(file_path) || file_path == "-help" || file_path == "-h"){
  cat(" COMMAND : Rscript eggnog_post.R <PATH/TO/TARGET/ORTHOLOG/FILES>
      ")
  quit(save="no")
}

file_name <- list.files(file_path)

# READ in ortholog results from eggnog run
if (substr(file_path, nchar(file_path), nchar(file_path)) == "/"){
  raw_table <- read.csv(paste0(file_path,file_name), header = TRUE, sep = "\t", skip = 5)
}else{
  raw_table <- read.csv(paste0(file_path, "/",file_name), header = TRUE, sep = "\t", skip = 5)
}

# Table Formatting
new_table <- raw_table[-c(((nrow(raw_table)-2):(nrow(raw_table)))),-c(3:ncol(raw_table))]

colnames(new_table) <- c("QUERY", "MATCH")

write.csv(new_table, file = "./eggnog_std.csv", row.names = FALSE)



