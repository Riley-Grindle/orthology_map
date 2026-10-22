#!/usr/bin/env Rscript
library(stringr)


args <- commandArgs(trailingOnly = TRUE)

file_path <- args[1]

if (is.na(file_path) || file_path == "-help" || file_path == "-h"){
  cat(" COMMAND : Rscript ortho_f_post.R <PATH/TO/TARGET/ORTHOLOG/FILES>
      ")
  quit(save="no")
}

# Change path name to ORTHO_f docker mounted results name ex: ORTHO_f_results (should only contain target prefix orthologs)
file_list <- as.array(list.files(file_path))

merge_tbl <- data.frame(matrix(ncol=3, nrow = 0))
colnames(merge_tbl) <- c("QUERY", "MATCH")
counter <- 0
for (i in 1:length(file_list)){
  
  # READ in ortholog results from orthofinder run 
  file_name <- file_list[i]
  if (substr(file_path, nchar(file_path), nchar(file_path)) == "/"){
    raw_table <- read.csv(paste0(file_path,file_name), header = TRUE, sep = "\t")
  }else{
    raw_table <- read.csv(paste0(file_path, "/",file_name), header = TRUE, sep = "\t")
  }
  
  raw_table <- raw_table[,-1]
  colnames(raw_table) <- c("QUERY", "MATCH")
  if (counter == 0){
    merge_tbl <- raw_table
  } else {
    merge_tbl <- merge(merge_tbl, raw_table, by = "QUERY", all = TRUE )
  }
  counter <- counter +1 
  
}

final_tbl <- data.frame(matrix(ncol=2, nrow=0))

for (i in 1:nrow(merge_tbl)){
  temp_str <- ''
  for (j in 2:ncol(merge_tbl)){
    if (j < ncol(merge_tbl)){
      temp_str <- paste0(temp_str, merge_tbl[i,j], ",")
    }else{
      temp_str <- paste0(temp_str, merge_tbl[i,j])
    }
  }
  final_tbl <- rbind(final_tbl, c(merge_tbl[i,1],temp_str))
}

colnames(final_tbl) <- c("QUERY", "MATCH")


# Loop through table and separate all orthologs of query sequences
new_table <- data.frame(matrix(ncol = 2, nrow = 0))
for (i in 1:nrow(final_tbl)){
  temp_list <- unlist(strsplit(final_tbl[,2][i], ","))
  key <- 0
  for (j in 1:length(temp_list)){
    if (key == 0){
      new_table <- rbind(new_table, c(final_tbl[,1][i], str_trim(temp_list[j])))
      key <- 1
    }else{
      new_table <- rbind(new_table, c('    ', str_trim(temp_list[j])))
    }
  }
}
colnames(new_table) <- c("QUERY", "MATCH")


newest_table <- data.frame(matrix(nrow = 0, ncol = 2))
non_na_list <- which(nchar(new_table[,1]) > 5)
non_na_list <- append(non_na_list, nrow(new_table)+1)
start <- 1
stop <- 2
key <- 0
for (i in 1:nrow(new_table)) {

  str_list <- unlist(str_split(new_table[i,1], ", "))

  if (str_list[1] != "    "){
    for (j in 1:length(str_list)){
      temp_df <- data.frame(str_list[j], new_table[non_na_list[start]:non_na_list[stop]-1, 2])
      if (nrow(temp_df) > 1 ) {
          temp_df[2:nrow(temp_df),1] <- "  "
      }
      colnames(temp_df) <- c("QUERY", "MATCH")
      if (key){
        temp_df[2,1] <- temp_df[1,1]
        temp_df <- temp_df[-1,]
      }
      newest_table <- rbind(newest_table, temp_df)
    }
    key <- 1
    start <- start + 1
    stop <- stop + 1
  }

}


# Extract Gene name from sequence headers
get_gene_info <- function(id) {
  result <- sub(".*GN=", "", id)
  
  return(result)
}

# Extract Taxa id from sequence headers
get_taxa_info <- function(id) {
  result <- sub(".*OX=(.*?)_.*", "\\1", id)
  
  return(result)
}

gene_list <- lapply(newest_table$MATCH, get_gene_info)
taxa_list <- lapply(newest_table$MATCH, get_taxa_info)

newest_table <- cbind(newest_table, unlist(gene_list), unlist(taxa_list))
colnames(newest_table) <- c("QUERY", "MATCH", "GENE", "SPECIES")



write.csv(newest_table, file = "./ortho_f_std.csv", row.names = FALSE)
