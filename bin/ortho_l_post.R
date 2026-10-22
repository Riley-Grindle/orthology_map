#!/usr/bin/env Rscript
############################
####### ORTHOLOGER  ########
############################

library(stringr)

args <- commandArgs(trailingOnly = TRUE)

file_path <- args[1]
prefix <- args[2]

if (is.na(file_path) || is.na(prefix) || file_path == "-help" || file_path == "-h"){
  cat(" COMMAND : Rscript tree_post.R <PATH/TO/TARGET/ORTHOLOG/FILES> <prefix_of_.fasta>
      ")
  quit(save="no")
}

file_name <- list.files(file_path)

# READ in ortholog results from Orthologer run
if (substr(file_path, nchar(file_path), nchar(file_path)) == "/"){
  raw_table <- read.csv(paste0(file_path,file_name), header = FALSE, sep = ' ')
}else{
  raw_table <- read.csv(paste0(file_path, "/",file_name), header = FALSE, sep = ' ')
}

raw_table <- raw_table[-c(1:9,(nrow(raw_table)-1), nrow(raw_table)), -c(3:ncol(raw_table))]

new_table <- data.frame(matrix(ncol=3, nrow=0))
query_index <- which(substr(raw_table[,2], 1, nchar(prefix)) == prefix)

new_table <- raw_table[query_index,]

new_table$V3 <- " "
raw_table <- raw_table[-c(query_index),]

# NOTE: this used to (a) rescan the *entire* raw_table with a nested loop for
# every query row to find same-cluster ("group") matches, and (b) grow
# `newest_table` via `rbind()` one row at a time. Both are O(n^2)-or-worse:
# with real cluster files (tens of thousands of rows) that nested scan alone
# is O(n_query * n_raw), and it was measured taking hours. Fix: index
# raw_table's rows by cluster id once (via split()), so each query row does
# an O(1) lookup instead of an O(n_raw) scan; and accumulate output rows in
# a list, doing a single rbind at the end instead of one per row. Output is
# identical -- same per-row logic, same "continuation" rows for extra
# matches within the same query's cluster.
raw_group_index <- split(seq_len(nrow(raw_table)), as.character(raw_table[,1]))

out_rows <- vector("list", 0)
for (i in 1:nrow(new_table)){
  key = 0
  group <- new_table[,1][i]
  matches_j <- raw_group_index[[as.character(group)]]
  if (!is.null(matches_j)){
    for (j in matches_j){
      if (key == 0){
        out_rows[[length(out_rows) + 1]] <- c(group, new_table[,2][i], substr(raw_table[,2][j], 1, nchar(raw_table[,2][j])))
        key = 1
      }else{
        out_rows[[length(out_rows) + 1]] <- c("   ", "   ", substr(raw_table[,2][j], 1, nchar(raw_table[,2][j])))
      }
    }
  }
}

if (length(out_rows) > 0){
  newest_table <- as.data.frame(do.call(rbind, out_rows), stringsAsFactors = FALSE)
}else{
  newest_table <- data.frame(matrix(ncol=3, nrow=0))
}

new_table <- newest_table[,-1]


colnames(new_table) <- c("QUERY", "MATCH")


write.csv(new_table, file = "./ortho_l_std.csv", row.names = FALSE)
