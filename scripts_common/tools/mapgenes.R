## Contains some functions for handling ARR data and gene IDs.
## They all require pathways/arr_hgnc.zip to be extracted to pathways/.
##
## readARRs(): reads pathways/arr_focal_pathways,7496x7496.csv,
## puts it into a tidy matrix, returns the matrix.
## pathways/arr_focal_pathways,7496x7496.csv is in pathways/arr_hgnc.zip.
##
## mapGeneIDs(): creates a mapping between Ensembl gene IDs and gene symbols,
## as well as labels those genes that are present in the ARR matrix.
## Uses the other file from pathways/arr_hgnc.zip.
## Other HGNC tables can be used instead, but this is discouraged,
## because the provided file was used when annotating the pathways used here.


readARRs <- function(out.dir) {
  library('data.table')
  ARR <- fread(file = file.path('pathways', 'arr_focal_pathways,7496x7496.csv'),
               sep = ',')
  arr.genes <- ARR$V1
  ARR <- as.matrix(ARR[,-1])
  # There was a pw names 'centr_pathway_SLC38A5, SN2.xlsx',
  # and the comma is gonna mess up CSV files. Replacing commas with _s.
  colnames(ARR) <- gsub(', ', '_', colnames(ARR))
  # The corresponding gene entry has a comma as well, actually...
  arr.genes <- gsub(',', '', arr.genes)
  rownames(ARR) <- arr.genes
  # Need to remove the pathways with abs(ARR) summing up to zero.
  ARR <- ARR[, colSums(abs(ARR)) != 0]
  # Also, need to return arr.genes somehow.
  write.csv(arr.genes, file.path(out.dir, 'arr_genes.csv'),
            quote = F, row.names = F)
  return(ARR)
}



mapGeneIDs <- function(ARR, expr.genes, project_id, genedata = NULL,
                       write.missing = FALSE, write.mapping = FALSE,
                       out.dir = NULL) {
  # Make a full table of all genes present both in expr.genes and hgnc.
  # With an additional column indicating whether it's in arr.genes.
  library('data.table')
  arr.genes <- rownames(ARR)
  gf <- file.path('pathways',
                  'hgnc_complete_set_Maxim_For_OB_samples_marking_13.07.17.txt')
  hgnc <- fread(gf)
  #hgnc <- hgnc[,c('symbol', 'ensembl_gene_id', 'locus_group')]
  hgnc <- hgnc[,c('symbol', 'ensembl_gene_id')]
  colnames(hgnc)[1] <- 'gene_symbol'
  hgnc <- subset(hgnc, ensembl_gene_id != '')
  #
  ### NEW BIT
  hgnc <- subset(hgnc, ensembl_gene_id %in% expr.genes)
  hgnc.1 <- hgnc
  # ^ why does this line exist? for more consistency between script versions
  hgnc.2 <- subset(hgnc, gene_symbol %in% arr.genes)
  # General plan:
  # 1. take hgnc
  # 2. filter to ensembl ids from expr
  # 3. label those that are present in the ARR data
  ###
  #
  # Interestingly, some of the missing gene entries are only missing because
  # they're actually several genes lumped together! E.g. 'AKR1D1 SRD5B1'.
  # I could split those, but since ARR is supposed to be a square matrix,
  # idk if I can handle that correctly.
  # So I'm just gonna exclude them.
  # There are only, like, 17 of those entries, after all.
  #
  # No, we do NOT need to use genedata here. Only the HGNC table.
  missing.genes <- arr.genes[!(arr.genes %in% hgnc.2$gene_symbol)]
  mapping <- hgnc.1
  inarr <- as.integer(mapping$gene_symbol %in% arr.genes)
  mapping$in_ARR <- inarr
  #
  # Sort
  mapping <- mapping[order(mapping$gene_symbol),]
  #
  ## Handling the ambiguity in mapping
  # If there's ambiguity in gene symbols, idk how to handle it
  if(nrow(mapping) != length(unique(mapping$gene_symbol))) {
    message('Error in mapGeneIDs(): ambiguous gene symbols')
    return(NULL)
  }
  # Now, need to exclude duplicated Ensembl IDs
  dupes <- unique(mapping$ensembl_gene_id[duplicated(mapping$ensembl_gene_id)])
  if(any(subset(mapping, ensembl_gene_id %in% dupes)$in_ARR == 1)) {
    message('Error in mapGeneIDs(): ambiguous Ensembl IDs in pathway genes')
    return(NULL)
  }
  # Since mapping is already sorted,
  # just take the 1st gene symbol for every dupe and discard the rest
  symbols.to.discard <- character()
  for(dupe in dupes) {
    m <- mapping[mapping$ensembl_gene_id == dupe,]
    symbols.to.discard <- c(symbols.to.discard, m[2:nrow(m)]$gene_symbol)
  }
  mapping <- subset(mapping, !(gene_symbol %in% symbols.to.discard))
  #
  #
  # Recording the missing genes and all pathways that contain them
  if(write.missing & is.null(out.dir)) {
    message('Warning from mapGeneIDs(): write.missing is TRUE but ',
            'output directory is not specified. Not writing anything.')
    write.missing <- F
  }
  if(write.missing) {
    missing.genes <- arr.genes[!(arr.genes %in% mapping$gene_symbol)]
    missing.genes.arr <- ARR[missing.genes,]
    missing.genes.pws <- c()
    for(pw.name in colnames(missing.genes.arr)) {
      if(any(missing.genes.arr[,pw.name] != 0)) {
        missing.genes.pws <- c(missing.genes.pws, pw.name)
      }
    }
    write.table(data.frame(V1 = missing.genes),
                file = file.path(out.dir,
                                 paste0('missing_genes_',
                                        project_id, '.txt')),
                quote = F, row.names = F, col.names = F)
    write.table(data.frame(V1 = missing.genes.pws),
                file = file.path(out.dir,
                                 paste0('missing_pathways_',
                                        project_id, '.txt')),
                quote = F, row.names = F, col.names = F)
  }
  #
  # Also recording the mapping itself
  if(write.mapping & is.null(out.dir)) {
    message('Warning from mapGeneIDs(): write.mapping is TRUE but ',
            'output directory is not specified. Not writing anything.')
    write.mapping <- F
  }
  if(write.mapping) {
    mapping.file <- file.path(out.dir,
                              paste0('gene_id_mapping_', project_id, '.csv'))
    write.csv(mapping, file = mapping.file,
              quote = F, row.names = F)
    message('Gene ID mapping data written to file ', mapping.file)
  }
  
  return(mapping)
}
