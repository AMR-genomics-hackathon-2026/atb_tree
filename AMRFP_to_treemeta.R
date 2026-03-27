#packages
library(data.table)
library(Matrix)

#import AMRfinderplus data for allthebacteria https://osf.io/ck7st 
amr_full <- fread("AMRFP_results.tsv", 
                  select = c("Name", c("Class","Subclass","Gene symbol","Element type")))
amr_full <- amr_full[`Element type` == "AMR"]

#target_cols <- c("Gene symbol", "Hierarchy node", "Subclass", "Class")
target_cols <- c("Class","Subclass","Gene symbol")

#Create presence absence csvs for different columns of AMRFP_results
for (column in target_cols){
  amr_data <- amr_full[, .(Name, get(column))]
  setnames(amr_data, 2, "Target")
  # Keep only unique combinations (binary presence-absence) observed more than once (e.g. for class)
  amr_unique <- unique(amr_data)
  name_factor <- as.factor(amr_unique$Name)
  col_factor <- as.factor(amr_unique$Target)
  
  # Build the sparse matrix
  sparse_pa <- sparseMatrix(i = name_factor, 
                            j = col_factor, 
                            x = 1)

  # Add row and column names
  rownames(sparse_pa) <- levels(name_factor)
  colnames(sparse_pa) <- levels(col_factor)
  
  ##prioritise common AMR determinants
  #Get the total occurrences per gene)
  gene_prevalence <- colSums(sparse_pa)
  
  #Get order from highest to lowest prevalence
  ordered_indices <- order(gene_prevalence, decreasing = TRUE)
  #Reorder the matrix columns
  sparse_pa_ordered <- sparse_pa[, ordered_indices]
  #reduce the size to the top 500 AMRvariants
  sparse_sub <- sparse_pa_ordered[, 1:min(500, ncol(sparse_pa_ordered))]
  # Convert the sparse subset to a dense matrix, then to a data.table
  sparse_dt <- as.data.table(as.matrix(sparse_sub))
  
  #Join reduced atb tree metadata to the unique list
  metadata <- fread("ATB_tree_taxonomic_assignments_GTDB.csv",
                  select = c("Sample", c("Family","Genus","Species")))
  metadata[, Sample := as.character(Sample)]
  metadata[, is_representative := NULL]
  metadata[, Root := NULL]
  #Add the 'Name' column back to join metadata
  sparse_dt[, Name := rownames(sparse_sub)]
  sparse_dt[, Name := as.character(Name)]
  #Merge matrix with ALL metadata columns
  final_output <- merge(metadata, sparse_dt, by.x = "Sample", by.y = "Name", all.y = TRUE)
  final_output[final_output == 0] <- NA
  
  ##Write out the CSV
  fwrite(final_output, paste0("AMRFP_res_for_atb_tree_", gsub(" ", "_", column), ".csv"))
  
  # Clear memory before the next iteration of the loop
  rm(sparse_pa, sparse_pa_ordered, sparse_sub, sparse_dt, amr_data, amr_unique)
  gc() 

}
