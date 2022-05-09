## Contains the functions for loading and normalizing expression data,
## as well as calculating PALs.
##
## normalizeExpr(): self-explanatory.
## calcPAL(): calculates PALs for gene-centric pathways.
## Both are called from dataprepTCGA().
## dataprepCPTAC() is not used (yet).
##
## Before running dataprepTCGA():
## > make sure that the two files from pathways/arr_hgnc.zip
##   are extracted to pathways/;
## > make sure that the TCGA data are downloaded to data/
##   with a proper directory structure, and the .tar.gz archives
##   in data/TCGA-*/Biospecimen/ and data/TCGA-*/Clinical/
##   are extracted to these respective directories.
## 
## dataprep.trad(): prepares expr data for oncoboxlib_calculate_scores,
## and generates the command itself.
## Running the command requires oncoboxlib installed:
## https://gitlab.com/oncobox/oncoboxlib
##

source('scripts_common/tools/write-session-info.R')
source('scripts_common/tools/mapgenes.R')


normalizeExpr <- function(E) {
  library('DESeq2')
  library('psych')
  library('preprocessCore')
  size.factors <- estimateSizeFactorsForMatrix(as.matrix(E+1))
  m3 <- t(t(as.matrix(E+1)) / size.factors)
  dimnames(m3) <- dimnames(E)
  return(m3)
}


calcPAL <- function(m3, ARR, missing_genes, use_normal = FALSE) {
  # Case-to-normal ratios
  library('psych')
  if(use_normal) {
    controls <- grep('^Norm_', colnames(m3), value = T)
    if(length(controls) == 0) {
      message('Warning: no normal samples, using mean tumor sample as control')
      use_normal <- FALSE
    }
  }
  if(use_normal) {
    m3.norm <- m3[,controls]
    m3.case <- m3[,!(colnames(m3) %in% colnames(m3.norm))]
    if(ncol(m3.case) == 0) {
      message('Error in calcPAL(): no tumor samples')
      return(NULL)
    }
    CNR <- m3.case / geometric.mean(t(as.matrix(m3.norm)))
  } else {
    CNR <- m3 / geometric.mean(t(as.matrix(m3)))
  }
  # Setting the CNRs of missing genes to 1 (so that their log10 is zero)
  CNR.missing.genes <- matrix(data = 1, nrow = length(missing_genes),
                              ncol = ncol(CNR))
  rownames(CNR.missing.genes) <- missing_genes
  colnames(CNR.missing.genes) <- colnames(CNR)
  CNR <- rbind(CNR, CNR.missing.genes)
  rm(CNR.missing.genes)
  # Aligning the matrix by gene symbol
  CNR <- CNR[rownames(ARR),]
  # Calculating PALs
  PAL <- t(t(as.matrix(log10(CNR))) %*% ARR) / colSums(abs(ARR))
  return(PAL)
}



dataprepTCGA <- function(projects, data.root, results.dir,
                         biotab.files, biotab.fields,
                         write.expr = FALSE, do.pal = FALSE,
                         main.histypes = NULL, use.normal = FALSE) {
  message('')
  program <- 'TCGA'
  if(use.normal) {
    sample.types <- c('Primary Tumor', 'Solid Tissue Normal')
  } else {
    sample.types <- c('Primary Tumor')
  }
  
  ## LOADING DATA
  library('data.table')
  library('jsonlite')
  message('Loading clinical data...')
  data.dir <- list()
  meta <- list()
  clin <- list()
  for(p in projects) {
    data.dir[[p]] <- file.path(data.root, p)
    meta.query <- '^metadata(.)*.json$'
    meta.file <- file.path(data.dir[[p]],
                           grep(meta.query,
                                list.files(data.dir[[p]]), value = T)[1])
    clin.file <- file.path(data.dir[[p]], 'Clinical', 'clinical.tsv')
    meta[[p]] <- fromJSON(meta.file)
    entities.id <- unlist(lapply(meta[[p]]$associated_entities,
                                 function(y) {return(y$entity_id)}))
    clin[[p]] <- fread(clin.file)
  }
  clin.df <- rbindlist(clin)
  # Cleaning up these two
  clin.cols <- c('case_submitter_id', 'project_id', 'primary_diagnosis',
                 'vital_status', 'days_to_death', 'days_to_last_follow_up')
  # 'tissue_or_organ_of_origin' for the loco
  clin.df <- unique(clin.df[,..clin.cols])
  colnames(clin.df)[2] <- 'project'
  colnames(clin.df)[3] <- 'histological_type'
  clin.cols.num <- grep('^days_', clin.cols, value = T)
  for(c in clin.cols.num) {
    clin.df[[c]] <- as.integer(clin.df[[c]])
  }
  clin.df[clin.df == "'--"] <- NA
  clin.df <- clin.df[order(clin.df$case_submitter_id),]
  
  message('Loading expression data...')
  first <- T
  for(p in projects) {
    for(i in 1:nrow(meta[[p]])) {
      if(length(meta[[p]][i,]$associated_entities) != 1) {
        message('Unexpected error in ', p, ' metadata, line ', i)
        next
      }
      eid <- meta[[p]][i,]$associated_entities[[1]]$entity_id
      esid <- meta[[p]][i,]$associated_entities[[1]]$entity_submitter_id
      expr1 <- fread(file.path(data.dir[[p]], 'STAR-Counts',
                               meta[[p]][i,]$file_id, meta[[p]][i,]$file_name))
      gene.data1 <- expr1[,c('gene_id','gene_name','gene_type')]
      expr1 <- expr1[,c('gene_id','unstranded')]
      if(substr(esid, 14, 15) == '11') {
        if(use.normal) {
          colnames(expr1)[2] <- paste0('Norm_', esid)
        } else { next }
      } else {
        colnames(expr1)[2] <- esid
      }
      if(first) {
        expr <- expr1
        gene.data <- gene.data1
        first <- F
      } else {
        if(colnames(expr1)[2] %in% colnames(expr)) {
          message('Warning: sample ', colnames(expr1)[2], ' already exists ',
                  'in the expression matrix!')
          next
        }
        expr <- merge(expr, expr1, by = 'gene_id', all = T)
        if(any(gene.data1 != gene.data)) {
          message('Warning: new gene mappings found in ', esid)
          gene.data <- merge(gene.data, gene.data1, all = T)
        }
      }
    }
  }
  
  # Make sure that all data are consistent
  cases.expr <- substr(gsub('^Norm_', '', colnames(expr)[-1]), 1, 12)
  samples.expr <- substr(gsub('^Norm_', '', colnames(expr)[-1]), 1, 16)
  cond1 <- all(cases.expr %in% clin.df$case_submitter_id)
  if(!cond1) {
    message(paste('All cases have clinical data: '), cond1)
    cols <- c(T, (cases.expr %in% clin.df$case_submitter_id))
    expr <- expr[,..cols]
    message('Warning: excluding cases with no clinical data ',
            'from the expression matrix (', sum(!cols),
            ' excluded, ', ncol(expr)-1, ' left)')
  }
  clin.df <- subset(clin.df, case_submitter_id %in% cases.expr)
  # Save this preliminary, filtered, tidy version of the clinical data
  # Just in case, idk why
  clin.df0 <- clin.df
  clin.outfile <- file.path(results.dir,
                            paste0('clinical_', program, '.tsv.old'))
  write.table(clin.df0, file = clin.outfile,
              quote = F, row.names = F, sep = '\t')
  message(paste('Initial clinical data written to file', clin.outfile))
  # Don't save as csv! there may be commas (in hist. type, for example)
  
  ## Cleaning up histological types, excluding the extra ones
  clin.df$histological_type <- gsub(',', '', clin.df$histological_type)
  if(!is.null(main.histypes) & !(all(is.na(clin.df$histological_type)))) {
    cases.to.discard <- subset(clin.df,
                               !(histological_type %in% main.histypes))
    cases.to.discard <- cases.to.discard$case_submitter_id
    clin.df <- subset(clin.df, !(case_submitter_id %in% cases.to.discard))
    barcodes.to.discard <- character()
    barcodes.expr <- gsub('^Norm_', '', colnames(expr)[-1])
    for(case in cases.to.discard) {
      barcodes.to.discard <- c(barcodes.to.discard,
                               grep(case, barcodes.expr, value = T))
    }
    cols <- !(gsub('^Norm_', '', colnames(expr)) %in% barcodes.to.discard)
    expr <- expr[,..cols]
  }
  
  # Exclude duplicates
  #
  # Let's only keep the FIRST aliquot (alphabetically).
  # A barcode goes like this:
  # TCGA-xx-PPPP-yyA-zzB-tttt-cc,
  # where xx is tissue source site,
  # PPPP is patient, yy is sample type
  # (01 for primary tumors, 11 for solid tissue normals),
  # A is vial, zz is portion, B is analyte, tttt is plate,
  # cc codes the center where the aliquot was analyzed.
  #
  # Gotta make sure that we keep ONE TUMOR SAMPLE from each patient
  # (multiple normal samples are okay).
  barcodes.to.discard <- character()
  barcodes.expr <- gsub('^Norm_', '', colnames(expr)[-1])
  cases.expr <- substr(barcodes.expr, 1, 12)
  cases.dupl <- unique(cases.expr[duplicated(cases.expr)])
  samples.expr <- substr(barcodes.expr, 1, 16)
  samples.dupl <- unique(samples.expr[duplicated(samples.expr)])
  for(case in cases.dupl) {
    dupl.group <- grep(case, barcodes.expr, value = T)
    dupl.group <- dupl.group[substr(dupl.group, 14, 15) == '01']
    # '01' = primary tumor, '11' = solid tissue normal
    dupl.group <- dupl.group[order(dupl.group)]
    barcodes.to.discard <- c(barcodes.to.discard, dupl.group[-1])
  }
  # Finally, excluding the recorded barcodes from expr
  cols <- !(gsub('^Norm_', '', colnames(expr)) %in% barcodes.to.discard)
  expr <- expr[,..cols]
  
  # Just printing some stats
  st.codes <- substr(gsub('^Norm_', '', colnames(expr)[-1]), 14, 15)
  message(paste('There are', sum(st.codes == '01'), 'tumor samples and',
                sum(st.codes == '11'), 'normal tissue samples.'))
  
  #
  ## COMPLETING CLINICAL DATA
  #
  message('Processing clinical data...')
  # Disabling the 'NAs introduced by coercion' warnings
  w <- getOption('warn')
  options(warn = -1)
  #
  # Initializing columns with days
  clin.df$days_fwup <- as.integer(NA)       # days to last follow-up
  clin.df$days_nte <- as.integer(NA)        # days to earliest new tumor event
  # Coercing these to integer is necessary.
  # Because by default, they were logical for some reason...
  clin.df$days_os <- clin.df$days_to_death  # days for overall survival analysis
  # Now we need to add NTE data and days to last known alive
  # (which have to be extracted from XMLs with clinical data)
  library('XML')
  clin.files <- character()
  for(p in projects) {
    clin.files <- c(clin.files, 
                    list.files(path = file.path(data.dir[[p]],
                                                'Clinical-XML'),
                               pattern = '*.xml$',
                               full.names = T, recursive = T))
  }
  # Shortcuts for XML element names
  el_fwups <- 'follow_ups'
  el_fwup_d <- 'days_to_last_followup'
  el_alive_d <- 'days_to_last_known_alive'
  el_nte <- 'new_tumor_events'
  el_nte_in <- 'new_tumor_event'
  el_nte_d <- 'days_to_new_tumor_event_after_initial_treatment'
  for(file in clin.files) {
    clin.doc <- xmlParse(file)
    clin.root <- xmlRoot(clin.doc)
    c <- xmlValue(clin.root[[2]][['bcr_patient_barcode']])
    if(!(c %in% clin.df$case_submitter_id)) { next }
    # Parsing follow-ups
    # There may be several follow-ups!
    value.fu <- NA   # the eventual days to follow-up value to be extracted
    # And the standalone el_fwup_d element should also be considered:
    value <- as.numeric(xmlValue(clin.root[[2]][[el_fwup_d]]))
    if(!is.na(value)) {
      fwups.days.values <- c(value)
    } else {
      fwups.days.values <- numeric()
    }
    # Looking through additional follow-ups
    n.fwups <- length(names(clin.root[[2]][[el_fwups]]))
    if(n.fwups > 0) {
      for(k in 1:n.fwups) {
        value <- xmlValue(clin.root[[2]][[el_fwups]][[k]][[el_fwup_d]])
        value <- as.numeric(value)
        if(!is.na(value)) {
          fwups.days.values <- c(fwups.days.values, value)
        }
      }
    }
    # Taking the latest follow-up date
    if(length(fwups.days.values) > 0) {
      value.fu <- max(fwups.days.values)
    }
    # Recording
    if(!is.na(value.fu)) {
      clin.df[clin.df$case_submitter_id == c, 'days_fwup'] <- value.fu
    }
    # Parsing days to last known alive, if applicable,
    # and assigning overall survival (os) days.
    # If no days to death, assign days to last known alive if available,
    # else assign days to last follow-up
    # (Might also make sense to check if all Dead patients have days to death)
    if(FALSE) {
      if(is.na(clin.df[clin.df$case_submitter_id == c,]$days_os)) {
        value <- as.numeric(xmlValue(clin.root[[2]][[el_alive_d]]))
        if(!is.na(value)) {
          clin.df[clin.df$case_submitter_id == c, 'days_os'] <- value
        } else {
          if(!is.na(value.fu)) {
            clin.df[clin.df$case_submitter_id == c, 'days_os'] <- value.fu
          }
        }
      }
    }
    # NOPE. Here's another version of that bit.
    # Let's still use days to last known alive
    # if no days to follow up or if that date is later.
    # Note that for days_pfs we should still use days to last fwup only,
    # because a patient might be 'known alive' but have an NTE.
    if(is.na(clin.df[clin.df$case_submitter_id == c,]$days_os)) {
      value.al <- as.numeric(xmlValue(clin.root[[2]][[el_alive_d]]))
      if(!is.na(value.fu)) {
        clin.df[clin.df$case_submitter_id == c, 'days_os'] <- value.fu
        if(!is.na(value.al) & value.al > value.fu) {
          clin.df[clin.df$case_submitter_id == c, 'days_os'] <- value.al
        }
      } else {
        if(!is.na(value.al)) {
          clin.df[clin.df$case_submitter_id == c, 'days_os'] <- value.al
        }
      }
    }
    #
    #
    # Parsing new tumor event data
    # if the NTE data are missing, value is gonna be NA anyway
    # There may be multiple NTEs!
    value.nte <- NA
    n.nte <- length(names(clin.root[[2]][[el_nte]]))
    if(n.nte > 0) {
      nte.days.values <- numeric()
      for(k in 1:n.nte) {
        if(names(clin.root[[2]][[el_nte]])[k] != el_nte_in) { next }
        value <- as.numeric(xmlValue(clin.root[[2]][[el_nte]][[k]][[el_nte_d]]))
        if(!is.na(value)) {
          nte.days.values <- c(nte.days.values, value)
        }
      }
      if(length(nte.days.values) > 0) {
        value.nte <- min(nte.days.values)
      }
    }
    if(!is.na(value.nte)) {
      clin.df[clin.df$case_submitter_id == c, 'days_nte'] <- value.nte
    }
  }
  
  #
  # Trying to extract missing NTE data from bcr biotab
  # This is there the biotab params come in.
  # biotab.files should be a list named after projects,
  # a list of LISTS of filenames.
  # biotab.fields should likewise be a named list,
  # but of LISTS, a list for each file, a list of character vectors
  # of all fields that should be extracted from the corresponding file.
  #
  biotab.data <- list()
  for(p in names(biotab.files)) {
    biotab.data[[p]] <- list()
    for(i in 1:length(biotab.files[[p]])) {
      filepath <- list.files(path = file.path(data.dir[[p]],
                                              'Clinical-Biotab'),
                             pattern = biotab.files[[p]][[i]],
                             full.names = T, recursive = T)
      if(length(filepath) > 1) {
        message(paste0('Warning: ', length(filepath), ' files supplied for ',
                       p, '-', i, ' (instead of 1). Taking the first one: ',
                       filepath[1]))
        filepath <- filepath[1]
      }
      biotab.data[[p]][[i]] <- fread(filepath)
    }
  }
  for(c in clin.df$case_submitter_id) {
    value.nte <- NA
    values <- numeric()
    p <- unlist(clin.df[clin.df$case_submitter_id == c, 'project'])
    if(!(p %in% names(biotab.files))) { next }
    for(i in 1:length(biotab.data[[p]])) {
      if(!(c %in% biotab.data[[p]][[i]]$bcr_patient_barcode)) { next }
      for(j in 1:length(biotab.fields[[p]][[i]])) {
        field <- biotab.fields[[p]][[i]][j]
        values1 <- subset(biotab.data[[p]][[i]],
                          bcr_patient_barcode == c)[[field]]
        values1 <- as.numeric(values1)
        values1 <- values1[!is.na(values1)]
        if(length(values1) > 0) {
          values <- c(values, values1)
        }
      }
    }
    if(length(values) > 0) {
      value.nte <- min(values)
    }
    if(is.na(value.nte)) { next }
    nte.missing <- is.na(clin.df[clin.df$case_submitter_id == c,
                                 'days_nte'])
    if(nte.missing) {
      clin.df[clin.df$case_submitter_id == c, 'days_nte'] <- value.nte
    } else {
      found.smaller.nte <- (value.nte<clin.df[clin.df$case_submitter_id == c,
                                              'days_nte'])
      if(found.smaller.nte) {
        clin.df[clin.df$case_submitter_id == c, 'days_nte'] <- value.nte
      }
    }
  }
  
  # Finally, making the days_pfs column (pfs = progression-free survival).
  # This one uses days_nte, if available;
  # if no days_nte and vital status is Alive, it uses days to last follow-up;
  # if no days_nte and vital status is Dead, it used days to death.
  # Event is either NTE of death.
  clin.df$days_pfs <- clin.df$days_nte
  for(c in clin.df$case_submitter_id) {
    if(is.na(clin.df[clin.df$case_submitter_id == c, 'days_pfs'])) {
      if(is.na(clin.df[clin.df$case_submitter_id == c, 'vital_status'])) {
        next
      }
      if(clin.df[clin.df$case_submitter_id == c, 'vital_status'] == 'Alive') {
        value.fu <- clin.df[clin.df$case_submitter_id == c, 'days_fwup']
        clin.df[clin.df$case_submitter_id == c, 'days_pfs'] <- value.fu
      }
      if(clin.df[clin.df$case_submitter_id == c, 'vital_status'] == 'Dead') {
        value <- clin.df[clin.df$case_submitter_id == c, 'days_to_death']
        clin.df[clin.df$case_submitter_id == c, 'days_pfs'] <- value
      }
    }
  }
  
  # Enabling all warnings again:
  options(warn = w)
  # Event columns
  clin.df$event_os <- 0
  clin.df$event_pfs <- 0
  filter1 <- clin.df$vital_status == 'Alive'
  filter2 <- clin.df$vital_status == 'Dead'
  clin.df$event_os[filter2] <- 1
  clin.df$event_pfs[filter2 | !is.na(clin.df$days_nte)] <- 1
  # Sorting
  clin.df <- clin.df[order(clin.df$case_submitter_id),]
  
  # One last check
  if(length(unique(clin.df$case_submitter_id)) != nrow(clin.df)) {
    message('Warning: duplicate patients left!')
  }
  
  ## Saving clinical data
  clin.outfile <- file.path(results.dir, paste0('clinical_', program, '.csv'))
  write.csv(clin.df, clin.outfile, quote = F, row.names = F)
  message(paste('Clinical data written to file', clin.outfile))
  
  ## Cleaning up expression data
  # Sort expr data by cols
  cols <- c(colnames(expr)[1], colnames(expr)[-1][order(colnames(expr)[-1])])
  expr <- expr[,..cols]
  # Exclude the N_ gene entries!
  ns <- grep('^N_', expr$gene_id, value = T)
  expr <- subset(expr, !(gene_id %in% ns))
  gene.data <- subset(gene.data, !(gene_id %in% ns))
  # Also exclude the Y chromosome genes
  gene.ids.Y <- grep('_PAR_Y$', expr$gene_id, value = T)
  expr <- subset(expr, !(gene_id %in% gene.ids.Y))
  gene.data <- subset(gene.data, !(gene_id %in% gene.ids.Y))
  # Sort by rows
  expr <- expr[order(expr$gene_id),]
  
  # Save the expr data
  if(write.expr) {
    exprn.outfile <- file.path(results.dir,
                               paste0('expr_', program, '.csv'))
    fwrite(expr, file = exprn.outfile, quote = F)
    message(paste('Expression data written to file', exprn.outfile))
  }
  
  if(!do.pal) {
    writeSessionInfo(file.path(results.dir,
                               paste0('sessionInfo_dataprep', program, '.txt')))
    return()
  }
  
  ## Converting expr data to a matrix
  expr.dt <- expr
  gene.ids <- expr.dt$gene_id
  expr <- as.matrix(expr[,-1])
  rownames(expr) <- gene.ids
  # Removing Ensembl ID versions, to be able to map them to gene symbols later
  # (e.g. ENSG00000000003.13 -> ENSG00000000003):
  rownames(expr) <- gsub('\\.(.)*', '', rownames(expr))

  
  #
  ## GENE ID MAPPING
  #
  # No, do this differently now.
  # Map BEFORE normalization.
  #
  message('Mapping gene IDs to gene symbols...')
  # Save the current gene data first
  gene.data0 <- gene.data
  gene.data.file <- file.path(results.dir,
                              paste0('gene_info_', program, '.csv'))
  write.csv(gene.data0, gene.data.file,
            quote = F, row.names = F)
  message('Gene info written to file ', gene.data.file)
  # Reading ARRs, generating the mapping
  ARR <- readARRs(results.dir)
  arr.genes <- read.csv(file.path(results.dir, 'arr_genes.csv'))[,1]
  ARR <- ARR[arr.genes,]
  message('ARR row name consistency: ', all(rownames(ARR) == arr.genes))
  gene.data$ensembl_gene_id <- gsub('\\.(.)*', '', gene.data$gene_id)
  mapping <- mapGeneIDs(ARR, rownames(expr), program, genedata = gene.data,
                        write.missing = T, write.mapping = T,
                        out.dir = results.dir)
  message('Pathway gene symbols mapped: ',
          sum(arr.genes %in% mapping$gene_symbol),
          ' / ', length(arr.genes))
  
  #
  ## NORMALIZING EXPRESSION DATA
  ## (exclude 1st and 2nd filtering if all genes need to be used)
  #
  # First filtering
  expr <- expr[mapping$ensembl_gene_id,]
  message('Total genes kept: ', nrow(expr))
  # The normalization
  message('Normalizing expression data (DESeq2)...')
  m3 <- normalizeExpr(expr)
  # Saving normalized expression data
  if(write.expr) {
    exprn.outfile <- file.path(results.dir,
                               paste0('expr_norm_', program, '.csv'))
    fwrite(cbind(data.table(gene_id = rownames(m3)), as.data.table(m3)),
           file = exprn.outfile, quote = F)
    message(paste('Normalized expression data written to file', exprn.outfile))
  }
  #
  # Second filtering, pre-PAL calculation
  mapping <- mapping[mapping$in_ARR == 1,]
  
  
  #
  ## CALCULATING PALS
  #
  # Renaming Ensembl IDs to gene symbols,
  # aligning m3 by gene symbols.
  m3 <- m3[mapping$ensembl_gene_id,]
  rownames(m3) <- mapping$gene_symbol
  message('Calculating PALs...')
  missing.genes <- read.csv(file.path(results.dir,
                                      paste0('missing_genes_',
                                             program, '.txt')),
                            header = F)[,1]
  PAL <- calcPAL(m3, ARR, missing.genes, use_normal = use.normal)
  if(is.null(PAL)) { return() }
  pal.file <- file.path(results.dir, paste0('PAL_', program, '.csv'))
  fwrite(cbind(data.table(V1 = rownames(PAL)), as.data.table(PAL)),
         file = pal.file, quote = F)
  message(paste('General PALs written to file', pal.file))
  
  #
  # Now we should also normalize the hist. types separately
  # and calculate PALs for them
  # (if there are multiple types, of course)
  tumor_types <- gsub('^TCGA-', '', unique(clin.df$project))
  K <- length(tumor_types)
  if(K < 2) {
    writeSessionInfo(file.path(results.dir,
                               paste0('sessionInfo_dataprep', program, '.txt')))
    return()
  }
  #
  ## Calculating PALs for survival analysis by type
  PAL.surv.list <- list()
  for(i in 1:K) {
    g <- tumor_types[i]      # 'g' for 'group'
    clin.df.g <- subset(clin.df, project == paste0('TCGA-', g))
    message(paste0('Group ', g, ': ', nrow(clin.df.g), ' samples'))
    # Filtering the expression matrix by this group
    barcodes.expr <- gsub('^Norm_', '', colnames(expr))
    cases.g <- clin.df.g$case_submitter_id
    barcodes.g <- character()
    for(case in cases.g) {
      barcodes.g <- c(barcodes.g, grep(case, barcodes.expr, value = T))
    }
    barcodes.g <- unique(barcodes.g)
    expr.g <- expr[,(gsub('^Norm_', '', colnames(expr)) %in% barcodes.g)]
    #
    #cases.g <- cases.g[order(cases.g)]
    #expr.g <- expr[,cases.g]
    if(use.normal) {
      n.norm <- sum(grepl('^Norm_', colnames(expr.g)))
      message(paste0('Group ', g, ': ', ncol(expr.g)-n.norm,
                     ' tumor samples and ', n.norm, ' normal samples used'))
    } else {
      message(paste0('Group ', g, ': ', ncol(expr.g), ' samples used'))
    }
    # Normalizing
    m3 <- normalizeExpr(expr.g)
    if(write.expr) {
      exprn.g.outfile <- file.path(results.dir,
                                   paste0('expr_norm_', program, '_',
                                          g, '.csv'))
      fwrite(cbind(data.table(gene_id = rownames(m3)), as.data.table(m3)),
             file = exprn.g.outfile, quote = F)
      message(paste('Normalized expression data for', g,
                    'written to file', exprn.g.outfile))
    }
    # Renaming Ensembl IDs to gene symbols,
    # aligning m3 by gene symbols.
    m3 <- m3[mapping$ensembl_gene_id,]
    rownames(m3) <- mapping$gene_symbol
    
    #
    # Calculating PALs
    PAL.surv.list[[g]] <- calcPAL(m3, ARR, missing.genes,
                                  use_normal = use.normal)
    if(is.null(PAL.surv.list[[g]])) { return() }
    pal.file <- file.path(results.dir,
                          paste0('PAL_', program, '_surv_', g, '.csv'))
    fwrite(cbind(data.table(V1 = rownames(PAL.surv.list[[g]])),
                 as.data.table(PAL.surv.list[[g]])),
           file = pal.file, quote = F)
    message(paste('PALs for group', g, 'written to file', pal.file))
  }
  
  # Printing sessionInfo
  writeSessionInfo(file.path(results.dir,
                             paste0('sessionInfo_dataprep', program, '.txt')))
  message('')
}



dataprepCPTAC <- function(project, data.root, loco, results.dir,
                          clinx.files, clinx.fields,
                          write.expr = FALSE, do.pal = FALSE,
                          use.normal = FALSE) {
  message('')
  program <- 'CPTAC'
  if(use.normal) {
    sample.types <- c('Primary Tumor', 'Solid Tissue Normal')
  } else {
    sample.types <- c('Primary Tumor')
  }
  origin <- list('lung' = c('Lung, NOS', 'Lower lobe, lung',
                            'Middle lobe, lung', 'Upper lobe, lung'),
                 'pancreas' = c('Body of pancreas', 'Head of pancreas',
                                'Pancreas, NOS', 'Tail of pancreas'),
                 'headneck' = c('Base of tongue, NOS', 'Cheek mucosa',
                                'Floor of mouth, NOS', 'Gum, NOS',
                                'Head, face or neck, NOS', 'Larynx, NOS',
                                'Lip, NOS', 'Oropharynx, NOS',
                                paste0('Overlapping lesion of lip, ',
                                       'oral cavity and pharynx'),
                                'Tongue, NOS', 'Tonsil, NOS'),
                 'kidney' = c('Kidney, NOS'),
                 'brain' = c('Brain, NOS', 'Frontal lobe', 'Occipital lobe',
                             'Parietal lobe', 'Temporal lobe'),
                 'corpusuteri' = c('Corpus uteri', 'Endometrium'))
  # CPTAC-3:
  # all head-n-neck cancers are squamous cell carcinomas
  # all lung cancers are either LUAD or LUSC (111 vs 108)
  # all kidney cancers are renal cell carcinomas
  # all pancreas cancers are infiltrating duct carcinomas
  # all uterus cancers are endometrioid adenocarcinomas
  # all brain cancers are glioblastomas
  tt <- list('pancreas' = 'IDCA', 'headneck' = 'HNSC',
             'kidney' = 'KIR', 'brain' = 'GBM', 'corpusuteri' = 'UCEAD')
  
  ## DATA PATHS
  data.dir <- file.path(data.root, project)
  if(use.normal) {
    meta.query <- '^metadata-with-normal(.)*.json$'
  } else {
    meta.query <- '^metadata(.)*.json$'
  }
  meta.file <- file.path(data.dir,
                         grep(meta.query,
                              list.files(data.dir), value = T)[1])
  clin.file <- file.path(data.dir,
                         'Clinical', 'clinical.tsv')
  spec.file <- file.path(data.dir,
                         'Biospecimen', 'sample.tsv')
  aliq.file <- file.path(data.dir,
                         'Biospecimen', 'aliquot.tsv')
  
  ## LOADING DATA
  library('jsonlite')
  library('data.table')
  message('Loading clinical and biospecimen data...')
  meta <- fromJSON(meta.file)
  cases <- unlist(lapply(meta$associated_entities,
                         function(y) {return(y$case_id)}))
  entities.id <- unlist(lapply(meta$associated_entities,
                               function(y) {return(y$entity_id)}))
  entities <- unlist(lapply(meta$associated_entities,
                            function(y) {return(y$entity_submitter_id)}))
  # all entity_type fields in meta$associated_entities are 'aliquot'.
  # looks like we need to use the aliquot file.
  clin.df <- fread(clin.file)
  aliq.df <- subset(fread(aliq.file), aliquot_id %in% entities.id)
  spec.df <- subset(fread(spec.file), sample_id %in% aliq.df$sample_id)
  # sooo, one aliquot may correspond to multiple samples,
  # all with different sample_submitter_id. wtf?
  #
  # for each case, keep only aliquots corresponding to solid primary tumors.
  # then pick the one with the latest sub id.
  #
  #
  # Filter by loco
  clin.df <- subset(clin.df, tissue_or_organ_of_origin %in% origin[[loco]])
  aliq.df <- subset(aliq.df, case_submitter_id %in% clin.df$case_submitter_id)
  # Merge (sa means spec+aliq)
  cols.aliq <- c('case_submitter_id', 'sample_submitter_id',
                 'aliquot_submitter_id')
  cols.spec <- c('case_submitter_id', 'sample_submitter_id',
                 'sample_type', 'composition')
  sa <- unique(merge(aliq.df[,..cols.aliq], spec.df[,..cols.spec], all.x = T))
  #
  # Keep only one aliquot per tumor case, make sure it's solid primary tumor;
  # all normal tissue aliquots are kept, if necessary
  aliqs.to.discard <- character()
  aliq.to.case <- list()
  aliq.to.sample <- list()
  cols <- c('case_submitter_id', 'sample_submitter_id', 'aliquot_submitter_id')
  for(case in clin.df$case_submitter_id) {
    case.df <- sa[sa$case_submitter_id==case & sa$composition=='Solid Tissue',]
    case.df <- unique(case.df[order(case.df$sample_submitter_id)])
    case.df.tum <- subset(case.df, sample_type == 'Primary Tumor')
    # Keep the 'oldest' tumor sample from the patient
    aliqs.to.keep <- c(case.df.tum$aliquot_submitter_id[1])
    #aliq.to.sample[[aliqs.to.keep[1]]] <- case.df.tum$sample_submitter_id[1]
    if(use.normal) {
      # Also keep all normal samples from that patient, if needed
      case.df.norm <- subset(case.df, sample_type == 'Solid Tissue Normal')
      aliqs.to.keep <- unique(c(aliqs.to.keep,
                                case.df.norm$aliquot_submitter_id))
    }
    # Discard the rest
    aliqs.to.discard <- c(aliqs.to.discard,
                          case.df$aliquot_submitter_id[
                            !(case.df$aliquot_submitter_id %in% aliqs.to.keep)
                          ])
    # Add to aliq-case-sample mapping
    for(a in aliqs.to.keep) {
      if(length(aliq.to.case) > 0) {
        if(a %in% names(aliq.to.case)) {
          message('Warning: ', a, ' already exists in aliq.to.case')
        }
      }
      aliq.to.case[[a]] <- case
      if(length(aliq.to.sample) > 0) {
        if(a %in% names(aliq.to.sample)) {
          message('Warning: ', a, ' already exists in aliq.to.sample')
        }
      }
      value <- subset(case.df, aliquot_submitter_id == a)$sample_submitter_id[1]
      aliq.to.sample[[a]] <- value
    }
  }
  aliq.df <- subset(aliq.df, !(aliquot_submitter_id %in% aliqs.to.discard))
  spec.df <- subset(spec.df,
                    sample_submitter_id %in% aliq.df$sample_submitter_id)
  # Determine which samples are normal tissue
  if(use.normal) {
    spec.df.norm <- subset(spec.df, sample_type == 'Solid Tissue Normal')
    aliq.df.norm <- subset(aliq.df, sample_id %in% spec.df.norm$sample_id)
    aliqs.norm <- unique(aliq.df.norm$aliquot_id)
  } else {
    aliqs.norm <- NULL
  }
  if(use.normal) {
    spec.df.tumor <- subset(spec.df, sample_type == 'Primary Tumor')
    spec.df.normal <- subset(spec.df, sample_type == 'Solid Tissue Normal')
    message(paste('There are', nrow(spec.df.tumor), 'tumor samples and',
                  nrow(spec.df.normal), 'normal tissue samples.'))
  }
  
  #
  # Load expression data
  message('Loading expression data...')
  first <- T
  for(i in 1:nrow(meta)) {
    if(length(meta[i,]$associated_entities) != 1) {
      message('Unexpected error in ', program, ' metadata, line ', i)
      next
    }
    eid <- meta[i,]$associated_entities[[1]]$entity_id
    esid <- meta[i,]$associated_entities[[1]]$entity_submitter_id
    if(!(eid %in% aliq.df$aliquot_id)) { next }
    case <- substr(aliq.to.sample[[esid]], 1, 9)
    if(length(case)==0) {next}
    if(!(case %in% clin.df$case_submitter_id)) { next }
    expr1 <- fread(file.path(data.dir, 'STAR-Counts',
                             meta[i,]$file_id, meta[i,]$file_name))
    gene.data1 <- expr1[,c('gene_id','gene_name','gene_type')]
    expr1 <- expr1[,c('gene_id','unstranded')]
    #colnames(expr1)[2] <- aliq.to.case[[esid]]
    if(eid %in% aliqs.norm) {
      colnames(expr1)[2] <- paste0('Norm_', aliq.to.sample[[esid]])
    } else {
      colnames(expr1)[2] <- aliq.to.sample[[esid]]
    }
    if(first) {
      expr <- expr1
      gene.data <- gene.data1
      first <- F
    } else {
      if(colnames(expr1)[2] %in% colnames(expr)) {
        message('Warning: sample ', colnames(expr1)[2], ' already exists ',
                'in the expression matrix!')
        next
      }
      expr <- merge(expr, expr1, by = 'gene_id', all = T)
      if(any(gene.data1 != gene.data)) {
        message('Warning: new gene mappings found in ', esid)
        gene.data <- merge(gene.data, gene.data1, all = T)
      }
    }
  }
  
  # Exclude the cases that didn't end up in expr data
  samples.expr <- gsub('^Norm_', '', colnames(expr)[-1])
  cases.expr <- subset(spec.df, sample_submitter_id %in% samples.expr)
  cases.expr <- unique(cases.expr$case_submitter_id)
  clin.df <- subset(clin.df, case_submitter_id %in% cases.expr)
  
  message('Processing clinical data...')
  #
  # Hey, should we consider cause of death?
  # If yes, make the 19 patients with 100% non-cancer-related deaths censored
  # (actually, let's disregard cause of death,
  # since we didn't consider it in TCGA)
  #
  
  clin.cols <- c('case_id', 'case_submitter_id',
                 'tissue_or_organ_of_origin', 'primary_diagnosis',
                 'vital_status', 'progression_or_recurrence', 'cause_of_death',
                 'days_to_death', 'days_to_recurrence',
                 'days_to_last_follow_up',
                 'last_known_disease_status')
  # just fyi, cause_of_death_source is not specified for any case
  # in all cases, days_to_last_known_disease_status is equal to days_fwup
  clin.df <- clin.df[,..clin.cols]
  clin.cols.num <- grep('^days_', clin.cols, value = T)
  for(c in clin.cols.num) {
    clin.df[[c]] <- as.integer(clin.df[[c]])
    # they all have fractional parts of 0.0, I checked
  }
  clin.df[clin.df == "'--"] <- NA
  spec.df[spec.df == "'--"] <- NA
  spec.cols <- colnames(spec.df)[!apply(spec.df, 2,
                                        function(x) {return(all(is.na(x)))})]
  spec.df <- unique(spec.df[,..spec.cols])
  
  # Save this preliminary, filtered, tidy version of the clinical data
  clin.df0 <- clin.df
  clin.outfile <- file.path(results.dir,
                            paste0('clinical_', program, '.tsv.old'))
  write.table(clin.df0, file = clin.outfile,
              quote = F, row.names = F, sep = '\t')
  message(paste('Initial clinical data written to file', clin.outfile))
  # Don't save as csv! there may be commas (in hist. types, for example)
  
  # Add localization (nah...)
  clin.site <- sapply(clin.df$tissue_or_organ_of_origin,
                      function(y) {
                        return(names(origin)[
                          unlist(lapply(origin,
                                        function(x) {
                                          return(y %in% x)
                                        }
                          ))])
                      })
  
  # There are some patients with vital_status 'Not Reported',
  # and no relevant values except for IDs, loco and diagnosis.
  # In CPTAC-3, there are 17 of them, including 9 with lung cancer.
  # Gotta exclude those (LATER, during surv analysis).
  # And btw, these cases are the only ones with missing days_fwup data.
  #clin.df <- subset(clin.df, vital_status %in% c('Alive', 'Dead'))
  #
  # Clean up
  clin.df1 <- clin.df
  clin.df$project <- project
  colnames(clin.df)[grep('^primary_diagnosis$',
                         colnames(clin.df))] <- 'histological_type'
  clin.df$histological_type <- gsub(', NOS$', '', clin.df$histological_type)
  colnames(clin.df)[grep('^days_to_recurrence$',
                         colnames(clin.df))] <- 'days_nte'
  colnames(clin.df)[grep('^days_to_last_follow_up$',
                         colnames(clin.df))] <- 'days_fwup'
  #clin.df$tumor_type <- NA
  if(loco == 'lung') {
    clin.df$tumor_type[
      clin.df$histological_type=='Adenocarcinoma'] <- 'LUAD'
    clin.df$tumor_type[
      clin.df$histological_type=='Squamous cell carcinoma'] <- 'LUSC'
  } else {
    clin.df$tumor_type <- tt[[loco]]
  }
  # Filling the important fields
  clin.df$days_os <- clin.df$days_to_death
  clin.df$days_os[is.na(clin.df$days_os)] <- clin.df$days_fwup[
    is.na(clin.df$days_os)]
  clin.df$days_pfs <- NA
  filter1 <- clin.df$vital_status == 'Alive'
  filter2 <- clin.df$vital_status == 'Dead'
  clin.df$days_pfs[filter1] <- clin.df$days_fwup[filter1]
  clin.df$days_pfs[filter2] <- clin.df$days_to_death[filter2]
  clin.df$days_pfs[!is.na(clin.df$days_nte)] <- clin.df$days_nte[
    !is.na(clin.df$days_nte)]
  # In the event cols, those with unrelated cause of death should be censored
  good.causes <- c(NA, 'Cancer Related', 'Not Reported', 'Unknown')
  # Filter for considering cause of death:
  #filter3 <- filter2 & (clin.df$cause_of_death %in% good.causes)
  # Filter disregarding cause of death:
  filter3 <- filter2
  clin.df$event_os <- 0
  clin.df$event_pfs <- 0
  clin.df$event_os[filter3] <- 1
  clin.df$event_pfs[filter3 | !is.na(clin.df$days_nte)] <- 1
  # Final
  clin.df <- clin.df[,c('case_submitter_id', 'project', 'tumor_type',
                        'histological_type', 'vital_status',
                        'days_to_death', 'days_fwup', 'days_nte',
                        'days_os', 'days_pfs', 'event_os', 'event_pfs')]
  clin.df <- clin.df[order(clin.df$case_submitter_id),]
  #
  # Save clinical data and specimen data
  clin.outfile <- file.path(results.dir,
                            paste0('clinical_', program, '.csv'))
  write.csv(clin.df, file = clin.outfile, quote = F, row.names = F)
  message(paste('Clinical data written to file', clin.outfile))
  spec.outfile <- file.path(results.dir, paste0('specimen_', program, '.csv'))
  write.csv(spec.df, file = spec.outfile, quote = F, row.names = F)
  message(paste('Specimen data written to file', spec.outfile))
  
  ## Cleaning up expression data
  # Sort expr data
  cols <- colnames(expr)[-1]
  cols <- c(colnames(expr)[1], cols[order(cols)])
  expr <- expr[,..cols]
  # Exclude the N_ gene entries!
  ns <- grep('^N_', expr$gene_id, value = T)
  expr <- subset(expr, !(gene_id %in% ns))
  gene.data <- subset(gene.data, !(gene_id %in% ns))
  # Also exclude the Y chromosome genes
  gene.ids.Y <- grep('_PAR_Y$', expr$gene_id, value = T)
  expr <- subset(expr, !(gene_id %in% gene.ids.Y))
  gene.data <- subset(gene.data, !(gene_id %in% gene.ids.Y))
  
  # Save the expr data
  if(write.expr) {
    exprn.outfile <- file.path(results.dir,
                               paste0('expr_', program, '.csv'))
    fwrite(expr, file = exprn.outfile, quote = F)
    message(paste('Expression data written to file', exprn.outfile))
  }
  
  if(!do.pal) {
    writeSessionInfo(file.path(results.dir,
                               paste0('sessionInfo_dataprep', program, '.txt')))
    return()
  }
  
  ## Converting expr data to a matrix
  expr.dt <- expr
  gene.ids <- expr.dt$gene_id
  expr <- as.matrix(expr[,-1])
  rownames(expr) <- gene.ids
  # Removing Ensembl ID versions, to be able to map them to gene symbols later
  # (e.g. ENSG00000000003.13 -> ENSG00000000003):
  rownames(expr) <- gsub('\\.(.)*', '', rownames(expr))
  
  #
  ## GENE ID MAPPING
  #
  message('Mapping gene IDs to gene symbols...')
  # Save the current gene data first
  gene.data0 <- gene.data
  gene.data.file <- file.path(results.dir,
                              paste0('gene_info_', program, '.csv'))
  write.csv(gene.data0, gene.data.file,
            quote = F, row.names = F)
  message('Gene info written to file ', gene.data.file)
  # Reading ARRs, generating the mapping
  ARR <- readARRs(results.dir)
  arr.genes <- read.csv(file.path(results.dir, 'arr_genes.csv'))[,1]
  ARR <- ARR[arr.genes,]
  message('ARR row name consistency: ', all(rownames(ARR) == arr.genes))
  gene.data$ensembl_gene_id <- gsub('\\.(.)*', '', gene.data$gene_id)
  mapping <- mapGeneIDs(ARR, rownames(expr), program, genedata = gene.data,
                        write.missing = T, write.mapping = T,
                        out.dir = results.dir)
  message('Pathway gene symbols mapped: ',
          sum(arr.genes %in% mapping$gene_symbol),
          ' / ', length(arr.genes))
  
  #
  ## NORMALIZING EXPRESSION DATA
  ## (disable 1st and 2nd filtering if all genes need to be used)
  #
  # First filtering
  expr <- expr[mapping$ensembl_gene_id,]
  message('Total genes kept: ', nrow(expr))
  # The normalization
  message('Normalizing expression data (DESeq2)...')
  m3 <- normalizeExpr(expr)
  # Saving normalized expression data
  if(write.expr) {
    exprn.outfile <- file.path(results.dir,
                               paste0('expr_norm_', program, '.csv'))
    fwrite(cbind(data.table(gene_id = rownames(m3)), as.data.table(m3)),
           file = exprn.outfile, quote = F)
    message(paste('Normalized expression data written to file', exprn.outfile))
  }
  #
  # Second filtering, pre-PAL calculation
  mapping <- mapping[mapping$in_ARR == 1,]
  
  
  #
  ## CALCULATING PALS
  #
  # Renaming Ensembl IDs to gene symbols,
  # aligning m3 by gene symbols.
  m3 <- m3[mapping$ensembl_gene_id,]
  rownames(m3) <- mapping$gene_symbol
  message('Calculating PALs...')
  missing.genes <- read.csv(file.path(results.dir,
                                      paste0('missing_genes_',
                                             program, '.txt')),
                            header = F)[,1]
  # Do the PAL calculation, write the result to file
  PAL <- calcPAL(m3, ARR, missing.genes, use_normal = use.normal)
  if(is.null(PAL)) { return() }
  pal.file <- file.path(results.dir, paste0('PAL_', program, '.csv'))
  fwrite(cbind(data.table(V1 = rownames(PAL)), as.data.table(PAL)),
         file = pal.file, quote = F)
  message(paste('General PALs written to file', pal.file))
  
  #
  # Now we should also normalize the hist. types separately
  # and calculate PALs for them
  # (if there are multiple types, of course)
  tumor_types <- unique(clin.df$tumor_type)
  K <- length(tumor_types)
  if(K < 2) {
    writeSessionInfo(file.path(results.dir,
                               paste0('sessionInfo_dataprep', program, '.txt')))
    return()
  }
  #
  ## Calculating PALs for survival analysis by type
  PAL.surv.list <- list()
  for(i in 1:K) {
    g <- tumor_types[i]      # 'g' for 'group'
    clin.df.g <- subset(clin.df, tumor_type == g)
    message(paste0('Group ', g, ': ', nrow(clin.df.g), ' samples'))
    # Filtering the expression matrix by this group
    spec.df.g <- subset(spec.df,
                        case_submitter_id %in% clin.df.g$case_submitter_id)
    cols <- (gsub('^Norm_', '',
                  colnames(expr)) %in% spec.df.g$sample_submitter_id)
    expr.g <- expr[,cols]
    if(use.normal) {
      n.norm <- sum(grepl('^Norm_', colnames(expr.g)))
      message(paste0('Group ', g, ': ', ncol(expr.g)-n.norm,
                     ' tumor samples and ', n.norm, ' normal samples used'))
    } else {
      message(paste0('Group ', g, ': ', ncol(expr.g), ' samples used'))
    }
    # Normalizing
    m3 <- normalizeExpr(expr.g)
    if(write.expr) {
      exprn.g.outfile <- file.path(results.dir,
                                   paste0('expr_norm_', program, '_',
                                          g, '.csv'))
      fwrite(cbind(data.table(gene_id = rownames(m3)), as.data.table(m3)),
             file = exprn.g.outfile, quote = F)
      message(paste('Normalized expression data for', g,
                    'written to file', exprn.g.outfile))
    }
    # Renaming Ensembl IDs to gene symbols,
    # aligning m3 by gene symbols.
    m3 <- m3[mapping$ensembl_gene_id,]
    rownames(m3) <- mapping$gene_symbol
    
    #
    # Calculating PALs
    PAL.surv.list[[g]] <- calcPAL(m3, ARR, missing.genes,
                                  use_normal = use.normal)
    if(is.null(PAL.surv.list[[g]])) { return() }
    pal.file <- file.path(results.dir,
                          paste0('PAL_', program, '_surv_', g, '.csv'))
    fwrite(cbind(data.table(V1 = rownames(PAL.surv.list[[g]])),
                 as.data.table(PAL.surv.list[[g]])),
           file = pal.file, quote = F)
    message(paste('PALs for group', g, 'written to file', pal.file))
  }
  
  # Printing sessionInfo
  writeSessionInfo(file.path(results.dir,
                             paste0('sessionInfo_dataprep', program, '.txt')))
  message('')
}



dataprep.trad <- function(results.dir, out.dir, tumor.types = NULL,
                          program = 'TCGA', use.normal = FALSE,
                          verbose = TRUE) {
  library('data.table')
  dir.create(out.dir, showWarnings = FALSE)
  bad.chars <- c('/', ',')

  # Make a list of all paths to relevant normalized expression files
  expr.files <- paste0('expr_norm_', program, '.csv')
  if(length(tumor.types) < 2) {
    tumor.types <- NULL
  }
  if(!is.null(tumor.types)) {
    for(tt in tumor.types) {
      tt0 <- tt
      for(c in bad.chars) {
        tt0 <- gsub(c, '-', tt0)
      }
      filename <- paste0('expr_norm_', program, '_', tt0, '.csv')
      expr.files <- c(expr.files, filename)
    }
  }
  
  # Process all those expression files
  library('psych')
  cmds <- character()
  for(expr.filename in expr.files) {
    expr.file <- file.path(results.dir, expr.filename)
    if(!file.exists(expr.file)) {
      message('Warning: file ', expr.file, ' does not exist, skipping it')
      next
    }
    if(verbose) {
      message('Processing ', expr.file, '...') 
    }
    
    # Read the expression file
    expr <- fread(expr.file)
    #expr.genes <- gsub('\\.(.)*', '', expr[[1]])
    expr.genes <- expr[[1]]
    expr <- as.matrix(expr[,-1])
    rownames(expr) <- expr.genes
    # Map gene Ensembl IDs to gene symbols
    mapping <- fread(file.path(results.dir,
                               paste0('gene_id_mapping_', program, '.csv')))
    expr <- expr[mapping$ensembl_gene_id,]
    rownames(expr) <- mapping$gene_symbol
    
    if(use.normal) {
      message('Error: normal tissue controls not implemented yet. Sorry.')
      return()
    } else {
      # Label the present samples as case
      samples <- colnames(expr)
      colnames(expr) <- paste0('Case_', samples)
      # Calculate controls
      expr.ctrl <- geometric.mean(t(expr))
      # Combine cases with controls
      if(FALSE) {
        N <- ncol(expr)
        for(sample in samples) {
          expr <- cbind(expr, expr.ctrl)
          N <- N + 1
          colnames(expr)[N] <- paste0('Control_', sample)
        }
      }
      expr <- cbind(expr, expr.ctrl)
      colnames(expr)[ncol(expr)] <- 'Norm_common'
    }
    
    # Write expr to file
    out.file <- file.path(out.dir, expr.filename)
    fwrite(cbind(data.table(SYMBOL = rownames(expr)), as.data.table(expr)),
           file = out.file, quote = FALSE)
    if(verbose) {
      message('New expression matrix written to file ', out.file)
    }
    
    # Generate command
    pal.file <- file.path(out.dir, gsub('^expr_norm', 'PAL', expr.filename))
    write.csv(data.frame(), pal.file)
    cmd <- paste0('oncoboxlib_calculate_scores',
                  ' --samples-file ', normalizePath(out.file),
                  ' --databases-dir ', normalizePath(file.path('pathways',
                                                               'databases')),
                  ' --databases-names all',
                  ' --results-file ', normalizePath(pal.file))
    unlink(pal.file)
    if(verbose) {
      message('Command to run in the terminal:')
      message(cmd)
    }
    cmds <- c(cmds, cmd)
    message('\n')
  }
  return(cmds)
}
