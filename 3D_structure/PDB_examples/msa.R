library(muscle)
library(Biostrings)

# read all sequences
seqAll_path <- "/scicore/home/zavolan/zhu0006/3D_structure/PDB/seqAll/concatAll.fasta"
seqAll <- readAAStringSet(seqAll_path, format = "fasta")
# extract protein IDs from sequence names
proteinIDs <- sapply(strsplit(names(seqAll), "\\|"), `[`, 1)

# read clusters
clusters_path <- "/scicore/home/zavolan/zhu0006/3D_structure/PDB/clusters/clusters-by-entity-90.txt"
clusters_file <- file(clusters_path, open = "rt")
clusters <- readLines(clusters_file)
close(clusters_file)

# define msa folder
msa_folder <- "/scicore/home/zavolan/zhu0006/3D_structure/PDB/msa_90"


# run msa for clusters
cluster_id <- 0
for (cluster in clusters) {
  cluster_id <- cluster_id + 1 # count cluster ID
  if (cluster_id < 0) {next} # start from the stopped point
  identifiers <- unlist(strsplit(cluster, " "))
  if (length(identifiers) > 500) {next} # skip clusters with more than 500 sequences
  if (length(identifiers) < 80) {break} # break the loop if less than 80 sequences
  current_cluster <- seqAll[which(proteinIDs %in% identifiers)]
  current_msa <- muscle(stringset = current_cluster)
  
  # save the file
  current_msa_file <- paste0(msa_folder, "/cluster_", as.character(cluster_id), ".aln")
  writeXStringSet(as(current_msa, "AAStringSet"), filepath = current_msa_file)
}


# write a function to extract insertions out of MultipleAlignment object
### insertions should meet the following conditions:
### (1) they should not at the start or the end of a sequence
### (2) they should be longer than 8 residues
findInsertions <- function(aln, cluster_id) { # aln <- as(AAMultipleAlignment, "AAStringSet")
  results <- data.frame(Cluster = c(0, 0),
                        Sequence = c(0, 0),
                        Start = c(0, 0),
                        End = c(0, 0))
  
  for (i in seq_len(length(aln))) {
    seq <- as.character(aln[[i]]) # extract current sequence
    gapPos <- gregexpr("-", seq)[[1]]
    gapRegions <- unname(tapply(gapPos, cumsum(c(1, diff(gapPos)) != 1), range))
    for (j in seq_len(length(gapRegions))) { # loop in different regions
      current_region <- gapRegions[[j]]
      current_start <- current_region[1]
      current_end <- current_region[2]
      regionLength <- current_end - current_start + 1
      if ((regionLength > 8) & (current_start > 9) & (nchar(seq)-current_end > 9)) {
        current_insertions <- c(cluster_id, i, current_start, current_end)
        results[nrow(results)+1, ] <- current_insertions
      }
    }
  }
  
  results <- results[-c(1:2), ]
  return(results)
}



# find insertions in MSA files
msa_folder <- "/scicore/home/zavolan/zhu0006/3D_structure/PDB/msa_90"
res_folder <- "/scicore/home/zavolan/zhu0006/3D_structure/PDB/insertions_90"
msa_files <- list.files(msa_folder, full.names = TRUE)
for (msa_file in msa_files) {
  cluster_name <- unlist(strsplit(basename(msa_file), "\\."))[1]
  print(cluster_name)
  cluster_id <- unlist(strsplit(cluster_name, "_"))[2]
  current_msa <- readAAStringSet(msa_file, format = "fasta")
  current_res <- findInsertions(current_msa, cluster_id)
  if (length(row.names(current_res)) > 0) {
    csv_path <- paste0(res_folder, "/", cluster_name, ".csv")
    write.csv(current_res, csv_path, row.names = FALSE)
  }
}


