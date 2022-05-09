## The big function, doLoco() ("do localization"),
## that starts the data preparation function(s) (see scripts_common/dataprep.R),
## checks the data,
## and starts the analysis functions (see scripts_common/analysis/).
##
## Both data preparation and analysis can be skipped:
## see the description of the arguments under the function name.
## I usually do data preparation for several localizations,
## restart RStudio to free some memory, and then do the analysis.
##
## doLoco() for one tumor localization can be called from its manual.R script:
## locos/*/scripts/manual.R,
## where all args are specified and all fine tuning is done.
## Multiple localizations are processed by scripts_common/runs/runall_1,2.R.


source('scripts_common/dataprep.R')
source('scripts_common/analysis/analysis-diff.R')
source('scripts_common/analysis/analysis-surv.R')

doLoco <- function(project_ids, program = c('TCGA', 'CPTAC'),
                   data_root, loc, results_dir,
                   biotab_files, biotab_fields, doubt_patients,
                   do_prep = FALSE, write_expr = FALSE, do_pal = FALSE,
                   do_prep_trad = FALSE,
                   use_normal = FALSE, main_histypes = NULL,
                   days_cutoff = NULL,
                   run_diff = FALSE, run_surv = FALSE,
                   run_diff_trad = FALSE, run_surv_trad = FALSE,
                   run_diff_genes = FALSE, run_surv_genes = FALSE,
                   surv_mode = 'median',
                   surv_quant_probs = c(0.25, 0.75),
                   test = FALSE, surv.debug = FALSE) {
  ## Togglable parameters:
  # do_prep: whether to launch the dataprep function
  #          calculating normalized expression for genes
  #          and PALs for gene-centric pathways
  #          (see scripts_common/dataprep.R)
  # write_expr, do_pal: whether to write the expression files and calculate PALs
  #                     (passed to the dataprep function)
  # do_prep_trad: whether to calculate PALs for classical ("traditional")
  #               pathways using oncoboxlib.
  #               ATTENTION: if do_prep_trad is TRUE,
  #               then either set do_prep and write_expr to TRUE,
  #               or make sure that the expr_norm_* and gene_id_mapping_* files
  #               are already present in results_dir.
  # use_normal: if TRUE, normal controls will be used;
  #             otherwise, mean tumor controls.
  #             ATTENTION: normal controls are not implemented for classical
  #             pathways yet.
  # main_histypes: histological types to keep, commas should be removed
  #                (see primary_diagnosis from data/*/Clinical/clinical.tsv)
  # days_cutoff: if not NULL, then cases with longer follow-ups are censored
  # run_diff: run tumor type differentiation analysis (not used in the report)
  # run_surv: run survival analysis
  # (run_diff and run_surv: default for gene-centric pathways,
  # *_trad for classical pathways, *_genes for individual genes.)
  # surv_mode: how to split the cases into two groups.
  #            'median': variable >= or < than the median;
  #            'quantiles': variable <= than surv_quant_probs[1]
  #                         or >= than surv_quant_probs[2].
  # surv_quant_probs: see surv_mode.
  # test: if TRUE, then only the first 100 items (pathways, genes, etc.)
  #       will be analyzed.
  # surv.debug: if TRUE, then lengths of vectors with result values
  #             will be printed after survival analysis.
  #             It was useful during ""development"".
  
  library('data.table')
  bad.chars <- c('/', ',')
  flags <- as.logical(rep(NA, 100))
  message('')
  message(paste0(rep('_', 60), collapse = ''))
  message('')
  message('Processing ', loc)
  
  ## Prepare clinical data and calculate PALs, if needed
  if(do_prep) {
    message('')
    message('Preparing data for genes and gene-centric pathways...')
    dir.create(results_dir, showWarnings = FALSE)
    if(program == 'TCGA') {
      dataprepTCGA(project_ids, data_root, results_dir,
                   biotab_files, biotab_fields,
                   write.expr = write_expr, do.pal = do_pal,
                   main.histypes = main_histypes,
                   use.normal = use_normal)
    } else if(program == 'CPTAC') {
      dataprepCPTAC(project_ids, data_root, loc, results_dir,
                    write.expr = write_expr, do.pal = do_pal,
                    use.normal = use_normal)
    }
    gc()
  }
  if(do_prep_trad) {
    if(program == 'TCGA') {
      message('')
      message('Preparing data for classical pathways...')
      cmds <- dataprep.trad(results_dir, file.path(results_dir, 'tradpws'),
                            tumor.types = gsub('^TCGA-', '', project_ids),
                            use.normal = use_normal,
                            verbose = FALSE)
      for(cmd in cmds) {
        message('')
        message('Running:')
        message(cmd)
        system(cmd)
      }
    }
  }
  
  ## Read clinical data
  message('')
  message('Reading data...')
  clin.file <- file.path(results_dir, paste0('clinical_', program, '.csv'))
  clin.df <- read.csv(clin.file)
  if(program == 'TCGA') {
    clin.df$tumor_type <- gsub('^TCGA-', '', clin.df$project)
  }
  tumor_types <- unique(clin.df$tumor_type)
  if(!all(clin.df$days_nte <= clin.df$days_to_death, na.rm = TRUE)) {
    message('Warning: some NTE dates are later than death dates! ',
            'Consider setting doubt_patients.')
  }
  
  ## Read diff PAL data
  pal.prefix <- paste0('PAL_', program)
  pal.file <- file.path(results_dir, paste0(pal.prefix, '.csv'))
  if(!file.exists(pal.file)) {
    message('PAL file ', pal.file, ' does not exist, terminating.')
    return()
  }
  PAL.all <- as.data.frame(fread(pal.file))
  rownames(PAL.all) <- PAL.all[,1]
  PAL.all <- PAL.all[,-1]
  colnames(PAL.all) <- gsub('\\.', '-', colnames(PAL.all))
  if(program == 'TCGA') {
    colnames(PAL.all) <- substr(colnames(PAL.all), 1, 12)
  } else if(program == 'CPTAC') {
    colnames(PAL.all) <- substr(colnames(PAL.all), 1, 9)
  }
  PAL.all <- as.matrix(PAL.all)
  patients <- colnames(PAL.all)
  patients <- patients[order(patients)]
  PAL.all <- PAL.all[,patients]
  flags[1] <- !any(is.na(PAL.all))
  flags[2] <- !any(is.infinite(PAL.all))
  flags[3] <- all(colnames(PAL.all) == clin.df$case_submitter_id)
  message(paste('1. All general PALs are valid:', flags[1]))
  message(paste('2. All general PALs are finite:', flags[2]))
  message(paste('3. Sample names consistency for general PALs:', flags[3]))
  
  ## Read surv PAL data, split clinical data
  clin.df.g <- list()
  PAL <- list()
  j <- 3
  if(length(tumor_types) > 1) {
    for(g in tumor_types) {
      # Remove the illegal characters if present (for file names):
      g0 <- g
      for(c in bad.chars) {
        g0 <- gsub(c, '-', g0)
      }
      # Extract PALs for this group
      clin.df.g[[g]] <- subset(clin.df, tumor_type == g)
      PAL[[g]] <- as.data.frame(fread(file.path(results_dir,
                                                paste0(pal.prefix, '_surv_',
                                                       g0, '.csv'))))
      rownames(PAL[[g]]) <- PAL[[g]][,1]
      PAL[[g]] <- PAL[[g]][,-1]
      #print(dim(PAL[[g]]))
      colnames(PAL[[g]]) <- gsub('\\.', '-', colnames(PAL[[g]]))
      if(program == 'TCGA') {
        colnames(PAL[[g]]) <- substr(colnames(PAL[[g]]), 1, 12)
      } else if(program == 'CPTAC') {
        colnames(PAL[[g]]) <- substr(colnames(PAL[[g]]), 1, 9)
      }
      PAL[[g]] <- as.matrix(PAL[[g]])
      patients <- colnames(PAL[[g]])
      patients <- patients[order(patients)]
      PAL[[g]] <- PAL[[g]][,patients]
      flags[j+1] <- !any(is.na(PAL[[g]]))
      flags[j+2] <- !any(is.infinite(PAL[[g]]))
      flags[j+3] <- all(colnames(PAL[[g]]) == clin.df.g[[g]]$case_submitter_id)
      message(paste0(j+1, '. All ', g, ' PALs are valid: ', flags[j+1]))
      message(paste0(j+2, '. All ', g, ' PALs are finite: ', flags[j+2]))
      message(paste0(j+3, '. Sample names consistency for ', g, ' PALs: ',
                     flags[j+3]))
      j <- j+3
    }
    for(i in 2:length(tumor_types)) {
      g1 <- tumor_types[i-1]
      g2 <- tumor_types[i]
      j <- j+1
      flags[j] <- all(rownames(PAL[[g1]]) == rownames(PAL[[g2]]))
      message(paste0(j, '. Pathway names consistency (', g1, '/', g2, '): ',
                     flags[j]))
    }
  } else {
    clin.df.g[[tumor_types[1]]] <- clin.df
    PAL[[tumor_types[1]]] <- PAL.all
  }
  
  ## Read PALs for traditional pathways
  pal.file1 <- file.path(results_dir, 'tradpws', paste0(pal.prefix, '.csv'))
  if(run_diff_trad & !file.exists(pal.file1)) {
    message('Warning: PAL file ', pal.file1, ' does not exist, ',
            'cannot proceed with diff analysis of traditional pathways')
    run_diff_trad <- FALSE
  } else if(file.exists(pal.file1)) {
    PAL.all.trad <- as.data.frame(fread(pal.file1))
    #PAL.all.trad$database <- gsub(' ', '_', PAL.all.trad$database)
    pathways.trad <- PAL.all.trad$pathway
    for(i in 1:length(pathways.trad)) {
      new.name <- gsub(' ', '_', pathways.trad[i])
      new.name <- gsub(',', '.', new.name)
      pathways.trad[i] <- paste0(gsub(' ', '_', PAL.all.trad$database[i]),
                                 '_', new.name)
    }
    col.excl <- c(1, 2, grep('Norm', colnames(PAL.all.trad)))
    PAL.all.trad <- as.matrix(PAL.all.trad[,-col.excl])
    rownames(PAL.all.trad) <- pathways.trad
    colnames(PAL.all.trad) <- gsub('^Case_', '', colnames(PAL.all.trad))
    if(program == 'TCGA') {
      colnames(PAL.all.trad) <- substr(colnames(PAL.all.trad), 1, 12)
    } else if(program == 'CPTAC') {
      colnames(PAL.all.trad) <- substr(colnames(PAL.all.trad), 1, 9)
    }
    patients <- colnames(PAL.all.trad)
    patients <- patients[order(patients)]
    PAL.all.trad <- PAL.all.trad[,patients]
    flags[j+1] <- !any(is.na(PAL.all.trad))
    flags[j+2] <- !any(is.infinite(PAL.all.trad))
    flags[j+3] <- all(colnames(PAL.all.trad) == clin.df$case_submitter_id)
    message(j+1, '. All classical PALs are valid: ', flags[j+1])
    message(j+2, '. All classical PALs are finite: ', flags[j+2])
    message(j+3, '. Sample names consistency for classical PALs: ', flags[j+3])
    j <- j+3
  }
  
  ## Read PALs for traditional pathways by group
  PAL.trad <- list()
  if(length(tumor_types) > 1) {
    for(g in tumor_types) {
      # Remove the illegal characters if present (for file names):
      g0 <- g
      for(c in bad.chars) {
        g0 <- gsub(c, '-', g0)
      }
      # Extract PALs for this group
      pal.file2 <- file.path(results_dir, 'tradpws',
                             paste0(pal.prefix, '_', g0, '.csv'))
      if(run_surv_trad & !file.exists(pal.file2)) {
        message('Warning: PAL file ', pal.file2, ' does not exist, ',
                'cannot proceed with surv analysis of traditional pathways')
        run_surv_trad <- FALSE
      }
      if(!run_surv_trad) {
        next
      }
      PAL.trad[[g]] <- as.data.frame(fread(pal.file2))
      
      pathways.trad <- PAL.trad[[g]]$pathway
      for(i in 1:length(pathways.trad)) {
        new.name <- gsub(' ', '_', pathways.trad[i])
        new.name <- gsub(',', '.', new.name)
        pathways.trad[i] <- paste0(gsub(' ', '_', PAL.trad[[g]]$database[i]),
                                   '_', new.name)
      }
      col.excl <- c(1, 2, grep('Norm', colnames(PAL.trad[[g]])))
      PAL.trad[[g]] <- as.matrix(PAL.trad[[g]][,-col.excl])
      rownames(PAL.trad[[g]]) <- pathways.trad
      colnames(PAL.trad[[g]]) <- gsub('^Case_', '', colnames(PAL.trad[[g]]))
      if(program == 'TCGA') {
        colnames(PAL.trad[[g]]) <- substr(colnames(PAL.trad[[g]]), 1, 12)
      } else if(program == 'CPTAC') {
        colnames(PAL.trad[[g]]) <- substr(colnames(PAL.trad[[g]]), 1, 9)
      }
      patients <- colnames(PAL.trad[[g]])
      patients <- patients[order(patients)]
      PAL.trad[[g]] <- PAL.trad[[g]][,patients]
      flags[j+1] <- !any(is.na(PAL.trad[[g]]))
      flags[j+2] <- !any(is.infinite(PAL.trad[[g]]))
      flags[j+3] <- all(colnames(PAL.trad[[g]]) == clin.df.g[[g]]$case_submitter_id)
      message(j+1, '. All classical ', g, ' PALs are valid: ', flags[j+1])
      message(j+2, '. All classical ', g, ' PALs are finite: ', flags[j+2])
      message(j+3, '. Sample names consistency for classical ', g, 'PALs: ',
              flags[j+3])
      j <- j+3
    }
    if(run_surv_trad) {
      for(i in 2:length(tumor_types)) {
        g1 <- tumor_types[i-1]
        g2 <- tumor_types[i]
        j <- j+1
        flags[j] <- all(rownames(PAL.trad[[g1]]) == rownames(PAL.trad[[g2]]))
        message(paste0(j, '. Classical pathway names consistency (',
                       g1, '/', g2, '): ', flags[j]))
      }
    }
  } else {
    if(file.exists(pal.file1)) {
      PAL.trad[[tumor_types[1]]] <- PAL.all.trad
    }
  }
  
  
  
  ## Reading expression data
  # Total
  expr.total <- fread(file.path(results_dir,
                                paste0('expr_norm_', program, '.csv')))
  expr.genes.total <- expr.total[[1]]
  expr.total <- expr.total[,-1]
  expr.total <- as.matrix(expr.total)
  cols.norm <- grep('^Norm_', colnames(expr.total), value = T)
  expr.total <- expr.total[,!(colnames(expr.total) %in% cols.norm)]
  rownames(expr.total) <- expr.genes.total
  colnames(expr.total) <- gsub('\\.', '-', colnames(expr.total))
  if(program == 'TCGA') {
    colnames(expr.total) <- substr(colnames(expr.total), 1, 12)
  } else if(program == 'CPTAC') {
    colnames(expr.total) <- substr(colnames(expr.total), 1, 9)
  }
  patients <- colnames(expr.total)
  patients <- patients[order(patients)]
  expr.total <- expr.total[,patients]
  flags[j+1] <- !any(is.na(expr.total))
  flags[j+2] <- !any(is.infinite(expr.total))
  flags[j+3] <- all(colnames(expr.total) == clin.df$case_submitter_id)
  message(paste0(j+1, '. All total exprs are valid: ', flags[j+1]))
  message(paste0(j+2, '. All total exprs are finite: ', flags[j+2]))
  message(paste0(j+3, '. Sample names consistency for total exprs: ',
                 flags[j+3]))
  j <- j+3
  # Per tumor type
  expr <- list()
  if(length(tumor_types) > 1) {
    for(g in tumor_types) {
      # Remove the illegal characters if present (for file names):
      g0 <- g
      for(c in bad.chars) {
        g0 <- gsub(c, '-', g0)
      }
      # Read the expr
      expr[[g]] <- fread(file.path(results_dir,
                                   paste0('expr_norm_', program,
                                          '_', g0, '.csv')))
      expr.genes <- expr[[g]][[1]]
      expr[[g]] <- expr[[g]][,-1]
      expr[[g]] <- as.matrix(expr[[g]])
      cols.norm <- grep('^Norm_', colnames(expr[[g]]), value = T)
      expr[[g]] <- expr[[g]][,!(colnames(expr[[g]]) %in% cols.norm)]
      rownames(expr[[g]]) <- expr.genes
      colnames(expr[[g]]) <- gsub('\\.', '-', colnames(expr[[g]]))
      if(program == 'TCGA') {
        colnames(expr[[g]]) <- substr(colnames(expr[[g]]), 1, 12)
      } else if(program == 'CPTAC') {
        colnames(expr[[g]]) <- substr(colnames(expr[[g]]), 1, 9)
      }
      patients <- colnames(expr[[g]])
      patients <- patients[order(patients)]
      expr[[g]] <- expr[[g]][,patients]
      flags[j+1] <- !any(is.na(expr[[g]]))
      flags[j+2] <- !any(is.infinite(expr[[g]]))
      flags[j+3] <- all(colnames(expr[[g]]) == clin.df.g[[g]]$case_submitter_id)
      message(paste0(j+1, '. All ', g, ' exprs are valid: ', flags[j+1]))
      message(paste0(j+2, '. All ', g, ' exprs are finite: ', flags[j+2]))
      message(paste0(j+3, '. Sample names consistency for ', g, ' exprs: ',
                     flags[j+3]))
      j <- j+3
    }
  } else {
    expr[[tumor_types[1]]] <- expr.total
  }
  
  
  ## Checking if all the data are OK
  flags <- flags[!is.na(flags)]
  message(paste('All flags are OK:', all(flags)))
  if(!all(flags)) {
    message('Error: inconsistent data, aborting.')
    return()
  }
  message('')
  
  #
  #
  ## ANALYSIS
  #
  #
  
  ## Run differentiation analysis for gene-centric pathways
  if(run_diff & (length(tumor_types) > 1)) {
    message('')
    message(paste0('Running differentiation analysis for ', program,
                   ' (gene-centric pathways)...'))
    if(length(tumor_types) == 2) {
      analysis.diff(PAL.all, rownames(PAL.all), clin.df,
                    file.path(results_dir,
                              paste0('results_diff_', program, '.csv')),
                    mode = 'pal')
    } else {
      analysis.diff.multi(PAL.all, rownames(PAL.all), clin.df,
                          file.path(results_dir,
                                    paste0('results_diff_', program, '.csv')),
                          mode = 'pal')
    }
  }
  
  ## Run survival analysis for gene-centric pathways
  if(run_surv) {
    do_pfs <- (program %in% c('TCGA', 'CPTAC'))
    for(g in tumor_types) {
      # Remove the illegal characters if present (for file names)
      g0 <- g
      for(c in bad.chars) {
        g0 <- gsub(c, '-', g0)
      }
      message('')
      message(paste0('Running survival analysis for ', program, ' ', g,
                     ' (gene-centric pathways)...'))
      all.pathways <- rownames(PAL[[g]])
      if(test) {
        all.pathways <- all.pathways[1:100]
      }
      analysis.surv(PAL[[g]], all.pathways, clin.df.g[[g]], do_pfs,
                    doubt_patients[[g]],
                    file.path(results_dir,
                              paste0('results_surv_', program, '_',g0, '.csv')),
                    debug = surv.debug,
                    mode = surv_mode, quant.probs = surv_quant_probs,
                    days.cutoff = days_cutoff)
    }
  }
  
  ## Run differentiation analysis for classical pathways
  if(run_diff_trad & (length(tumor_types) > 1)) {
    message('')
    message('Running differentiation analysis for ', program,
            ' (classical pathways)...')
    if(length(tumor_types) == 2) {
      analysis.diff(PAL.all.trad, rownames(PAL.all.trad), clin.df,
                    file.path(results_dir, 'tradpws',
                              paste0('results_diff_', program, '.csv')),
                    mode = 'pal')
    } else {
      analysis.diff.multi(PAL.all.trad, rownames(PAL.all.trad), clin.df,
                          file.path(results_dir, 'tradpws',
                                    paste0('results_diff_', program, '.csv')),
                          mode = 'pal')
    }
  }
  
  ## Run survival analysis for classical pathways
  if(run_surv_trad) {
    do_pfs <- (program %in% c('TCGA', 'CPTAC'))
    for(g in tumor_types) {
      # Remove the illegal characters if present (for file names)
      g0 <- g
      for(c in bad.chars) {
        g0 <- gsub(c, '-', g0)
      }
      message('')
      message(paste0('Running survival analysis for ', program, ' ', g,
                     ' (classical pathways)...'))
      all.pathways <- rownames(PAL.trad[[g]])
      if(test) {
        all.pathways <- all.pathways[1:100]
      }
      analysis.surv(PAL.trad[[g]], all.pathways, clin.df.g[[g]], do_pfs,
                    doubt_patients[[g]],
                    file.path(results_dir, 'tradpws',
                              paste0('results_surv_', program, '_',g0, '.csv')),
                    debug = surv.debug,
                    mode = surv_mode, quant.probs = surv_quant_probs,
                    days.cutoff = days_cutoff)
    }
  }
  
  ## Analysis for genes
  # But first we might need to filter them out
  # ...actually, no
  
  ## Running control diff analysis for individual genes
  if(run_diff_genes & (length(tumor_types) > 1)) {
    genes.total <- expr.genes.total
    if(test) {
      genes.total <- genes.total[1:100]
    }
    outfile <- file.path(results_dir,
                         paste0('results_diff_', program, '_genes.csv'))
    if(nrow(expr.total) > 0 & length(genes.total) > 0) {
      message('')
      message(paste0('Running differentiation analysis for ', program,
                     ' (genes)...'))
      if(length(tumor_types) == 2) {
        analysis.diff(expr.total, genes.total, clin.df,
                      outfile, mode = 'expr')
      } else {
        analysis.diff.multi(expr.total, genes.total, clin.df,
                            outfile, mode = 'expr')
      }
    }
  }
  
  ## Running control survival analysis for individual genes
  if(program %in% c('TCGA', 'CPTAC')) {
    surv_types <- c('os', 'pfs')
  } else {
    surv_types <- c('os')
  }
  if(run_surv_genes) {
    for(g in tumor_types) {
      # Remove the illegal characters if present (for file names)
      g0 <- g
      for(c in bad.chars) {
        g0 <- gsub(c, '-', g0)
      }
      for(survtype in surv_types) {
        genes.g <- rownames(expr[[g]])
        if(test) {
          genes.g <- genes.g[1:100]
        }
        outfile <- file.path(results_dir,
                             paste0('results_surv_', program, '_',
                                    g0, '_genes_', survtype, '.csv'))
        if(nrow(expr[[g]]) > 0 & length(genes.g) > 0) {
          message('')
          message(paste0('Running ', toupper(survtype), ' analysis for ',
                         program, ' ', g, ' (genes)...'))
          analysis.surv.genes(expr[[g]], genes.g,
                              clin.df.g[[g]],
                              c(survtype),
                              doubt_patients[[g]],
                              outfile,
                              debug = surv.debug,
                              mode = surv_mode, quant.probs = surv_quant_probs,
                              days.cutoff = days_cutoff)
        }
      }
    }
  }
  
  gc()
  
}
