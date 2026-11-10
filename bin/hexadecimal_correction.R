#!/usr/bin/env Rscript
############################
## HEXADECIMAL Correction ##
############################

library(stringr)

args <- commandArgs(trailingOnly = TRUE)

fasta_path <- args[1]
gene_path <- args[2]
out_path <- args[3]

if (is.na(fasta_path) || fasta_path == "-help" || fasta_path == "-h"){
  cat(" COMMAND : Rscript tree_post.R <PATH/TO/TARGET/FASTA/FILES> <PATH/TO/HEADER_DIR> <PATH/TO/OUTS_TABLE>
      ")
  quit(save="no")
}

file_list <- list.files(fasta_path)
genes_list <- list.files(gene_path)

loger_out <- read.csv(paste0(out_path, "ortho_l_std.csv"))

# ORTHOLOGER sanitizes its species labels (replaces any char outside
# [A-Za-z0-9_] with '_') before writing them into ortho_l_std.csv's
# QUERY/MATCH columns (e.g. "ze-mu" -> "ze_mu"), but file_list here is read
# straight from the (unsanitized) odbdata/ filenames. search_ids must match
# the sanitized form or the match() lookup below silently finds nothing.
search_ids <- c()
for (file in file_list){
  temp_str <- unlist(strsplit(file, ".", fixed = TRUE))
  search_ids <- append(search_ids, gsub("[^A-Za-z0-9_]", "_", temp_str[1]))
}

replacement_list <- c()
for (file in genes_list){
  replacement_list <- append(replacement_list, read.csv(paste0(gene_path, file), header = FALSE))
}

# Pre-unlist each species' gene-header list once. The loops below used to call
# unlist(replacement_list[replacement_idx]) fresh on every one of the ~n rows,
# redoing the same work n times; doing it once per species up front turns that
# into a handful of calls instead of tens of thousands.
replacement_vecs <- lapply(seq_along(replacement_list), function(k) unlist(replacement_list[k]))

loger_out$FASTA <- unlist(strsplit(loger_out$MATCH, ":"))[c(TRUE, FALSE)]
loger_out$HEXI_M <- unlist(strsplit(loger_out$MATCH, ":"))[c(FALSE, TRUE)]

# NOTE: this whole section used to build `revised_names` via
# `revised_names <- append(revised_names, ...)` inside a for loop over every
# row (n ~ tens of thousands), on top of the redundant unlist() calls above --
# together those made this O(n^2)-ish and slow. The per-row species lookup
# (match()) and hex-string cleanup/decoding are fully vectorizable, so do
# those across the whole column at once; only the final "look up this row's
# value in that row's species vector" step still needs a loop, and it's now a
# cheap O(1)-per-row index into an already-unlisted vector.
fasta_idx_vec <- match(loger_out$FASTA, search_ids)

hexi_vec <- loger_out$HEXI_M
has_trailing_space <- substr(hexi_vec, nchar(hexi_vec), nchar(hexi_vec)) == " "
hexi_vec[has_trailing_space] <- substr(hexi_vec[has_trailing_space], 1, nchar(hexi_vec[has_trailing_space]) - 1)
hex_idx_vec <- as.numeric(strtoi(hexi_vec, base = 16)) + 1

n <- nrow(loger_out)
revised_names <- character(n)
for (i in 1:n){
  vec <- replacement_vecs[[fasta_idx_vec[i]]]
  val <- vec[hex_idx_vec[i]]
  revised_names[i] <- substr(val, 2, nchar(val))
}

revised_names <- as.data.frame(revised_names)


# Same fix as above, applied to QUERY: vectorize the split/lookup/hex-decode,
# keep a loop only for the final per-row vector index. Rows whose QUERY is
# the "   " placeholder (synteny continuation rows) are left untouched, same
# as the original `if (temp_list[1] != "   ")` guard.
query_split <- strsplit(loger_out$QUERY, ":")
query_prefix <- vapply(query_split, function(x) x[1], character(1))
query_hex_raw <- vapply(query_split, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))

query_idx_vec <- match(query_prefix, search_ids)
query_hex_idx_vec <- as.numeric(strtoi(query_hex_raw, base = 16)) + 1

new_query <- loger_out$QUERY
real_rows <- which(query_prefix != "   ")
for (i in real_rows){
  vec <- replacement_vecs[[query_idx_vec[i]]]
  val <- vec[query_hex_idx_vec[i]]
  new_query[i] <- substr(val, 2, nchar(val))
}
loger_out$QUERY <- new_query

new_table <- cbind(loger_out$QUERY, revised_names[,1])

colnames(new_table) <- c("QUERY", "MATCH")

new_table <- as.data.frame(new_table)
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

gene_list <- lapply(new_table$MATCH, get_gene_info)
taxa_list <- lapply(new_table$MATCH, get_taxa_info)

newest_table <- cbind(new_table, unlist(gene_list), unlist(taxa_list))
colnames(newest_table) <- c("QUERY", "MATCH", "GENE", "SPECIES")


write.csv(newest_table, file = "./ortho_l_std.csv", row.names = FALSE)
