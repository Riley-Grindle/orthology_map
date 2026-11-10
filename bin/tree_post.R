#!/usr/bin/env Rscript
############################
####### TREEGRAFTER ########
############################

library(stringr)

args <- commandArgs(trailingOnly = TRUE)

file_path <- args[1]

if (is.na(file_path) || file_path == "-help" || file_path == "-h"){
  cat(" COMMAND : Rscript tree_post.R <PATH/TO/TARGET/ORTHOLOG/FILES>
      ")
  quit(save="no")
}

file_name <- list.files(file_path)

# READ in ortholog results from Treegraft run
if (substr(file_path, nchar(file_path), nchar(file_path)) == "/"){
  raw_table <- read.csv(paste0(file_path,file_name), header = FALSE)
}else{
  raw_table <- read.csv(paste0(file_path, "/",file_name), header = FALSE)
}

raw_table <- as.data.frame(raw_table[,-2])

# Delete Spacer Rows
delete_rows <- which(substr(raw_table[,1], 1, 3) == "___")

if (length(delete_rows) > 0){
  raw_table <- raw_table[-delete_rows,]
}
raw_table <- as.data.frame(raw_table)

# Format in QUERY - MATCH data frame
#
# NOTE: this used to build new_table via `new_table <- rbind(new_table, addition)`
# inside the loop below. rbind() on a data.frame copies the entire (growing)
# table on every call, so for a ~170k-row input that loop is O(n^2) and can
# take many hours (or longer). Since raw_table always has exactly 1 column
# here (column 2 was dropped above), every iteration of this loop produces
# exactly one output row, so we can pre-allocate QUERY/MATCH vectors of
# length nrow(raw_table) and fill them by index instead -- same output,
# O(n) instead of O(n^2).
n <- nrow(raw_table)
query_vec <- character(n)
match_vec <- character(n)

for (i in 1:n){
  temp_str <- raw_table[i,1]
  if (substr(temp_str, 1, 7) != "gene_id"){
    if (i < n && substr(raw_table[i+1,1], 1, 7) == "gene_id"){
      query_vec[i] <- str_trim(gsub("_", "|", temp_str))
      match_vec[i] <- substr(raw_table[i+1, 1], 9, nchar(raw_table[i+1, 1]))
    }else{
      query_vec[i] <- str_trim(gsub("_", "|", temp_str))
      match_vec[i] <- "NO MATCH"
    }
  }else{
    for (j in 1:ncol(raw_table)){
      query_vec[i] <- "     "
      match_vec[i] <- substr(raw_table[i, j], 9, nchar(raw_table[i, j]))
    }
  }
}

new_table <- data.frame(QUERY = query_vec, MATCH = match_vec, stringsAsFactors = FALSE)

write.csv(new_table, file = "./tree_std.csv", row.names = FALSE)
