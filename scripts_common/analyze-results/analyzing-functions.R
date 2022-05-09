## Contains functions for summarizing the raw survival results;
## these are used by multiple scripts from scripts_common/analyze-results/.


overlapSignif <- function(num.to.check, n.group, m.group, N, M) {
  ## Permutation test returning p-val of an intersection
  m <- length(m.group)
  n <- length(n.group)
  iters <- 1e5
  intersections <- numeric(iters)
  for(i in 1:iters) {
    set.seed(i)
    intersections[i] <- length(intersect(sample(n, N), sample(m, M)))
  }
  p <- sum(intersections > num.to.check) / length(intersections)
  return(p)
}


printIfSmall <- function(v, fltr) {
  ## Prints the significant items if there aren't too many of them
  N.max <- 50
  v1 <- v[fltr]
  if(length(v1) > 0 & length(v1) <= N.max) {
    message(paste(c('#', v1), collapse = ' '))
  }
}


plotBoxes <- function(res.df, sample.size, log.y = FALSE, violins = TRUE) {
  ## Makes an untidy, exploratory violin plot with paired lines
  # Here's what helped:
  # https://stackoverflow.com/questions/49370705/implementing-paired-lines-into-boxplot-ggplot2
  library('ggplot2')
  plots <- list()
  tumor_types <- names(res.df)
  surv_types <- character()
  for(g in tumor_types) {
    surv_types <- c(surv_types, names(res.df[[g]]))
  }
  surv_types <- unique(surv_types)
  set.seed(42)
  for(g in tumor_types) {
    plots[[g]] <- list()
    for(st in surv_types) {
      if(is.null(nrow(res.df[[g]][[st]])) | nrow(res.df[[g]][[st]]) < 3) {
        plots[[g]][[st]] <- NULL
        next
      }
      # First, plot the boxes / violins themselves
      hrfield.g <- paste0('gene_', st, '_hr')
      hrfield.pw <- paste0('pw_', st, '_hr')
      lbl1 <- 'Genes'
      lbl2 <- 'Pathways'
      df <- rbind(data.frame(x = lbl1,
                             hr = res.df[[g]][[st]][[hrfield.g]]),
                  data.frame(x = lbl2,
                             hr = res.df[[g]][[st]][[hrfield.pw]]))
      df$x <- as.factor(df$x)
      p <- ggplot(df, aes(x = x, y = hr, fill = x))
      p <- p + scale_fill_brewer(palette = 'Accent')
      if(log.y) {
        p <- p + scale_y_continuous(trans = 'log10')
      }
      if(violins) {
        p <- p + geom_violin()
      }
      p <- p + geom_boxplot(width = 0.3, outlier.shape = NA)
      p <- p + labs(x = NULL, y = 'Cox Hazard Ratio')
      p <- p + scale_x_discrete(labels = c(lbl1 = 'Genes',
                                           lbl2 = 'Gene-centric Pathways'))
      p <- p + theme(legend.position = 'none')
      # Next, plot the cross-lines;
      # sample.size <= 0 would mean no cross-lines (and no points)
      if(sample.size > 0) {
        if(nrow(res.df[[g]][[st]]) > sample.size) {
          random.cases <- sample(1:nrow(res.df[[g]][[st]]), sample.size)
        } else {
          random.cases <- 1:nrow(res.df[[g]][[st]])
        }
        # Make jitter for the future points and lines
        jitt <- runif(nrow(df), -0.05, 0.05)
        jitt.g <- jitt[df$x == lbl1]
        jitt.pw <- jitt[df$x == lbl2]
        # Add the points and lines, pair by pair
        for(i in random.cases) {
          hrs <- c(subset(df, x == lbl1)$hr[i],
                   subset(df, x == lbl2)$hr[i])
          df1 <- data.frame(x = c(1+jitt.g[i], 2+jitt.pw[i]), hr = hrs)
          # how about highlighting the lines that cross the HR == 1 line
          hrs <- hrs - 1.0
          if((hrs[1]*hrs[2]) < 0) {
            clr <- 'red'
            opacity <- 0.2
          } else {
            clr <- 'black'
            opacity <- 0.2
          }
          p <- p + geom_line(aes(x = as.numeric(x), y = hr), df1,
                             inherit.aes = FALSE, color = clr, alpha = opacity)
          p <- p + geom_point(aes(x = as.numeric(x), y = hr), df1,
                              inherit.aes = FALSE, color = clr, size = 0.7)
        }
      }
      # Save the plot
      plots[[g]][[st]] <- p
    }
  }
  return(plots)
}


read.surv <- function(results.dir, program, tumor.types,
                      mode = c('gc', 'classic'),
                      keep.loci = c('protein-coding gene', 'pseudogene',
                                    'non-coding RNA', 'other')) {
  ## Reads and formats survival results for:
  ## genes and gene-centric pathways ('gc' mode) or classical pathways
  ## gc means gene-centric, classic is self-explanatory
  
  library('data.table')
  if(program == 'TCGA') {
    surv.types <- c('os', 'pfs')
  } else {
    surv.types <- c('os')
  }
  bad.chars <- c('/', ',')
  
  if(mode == 'classic') {
    res.surv.pws <- list()
    for(g in tumor.types) {
      g0 <- g
      for(c in bad.chars) {
        g0 <- gsub(c, '-', g0)
      }
      filename <- paste0('results_surv_', program, '_', g0, '.csv')
      res1 <- fread(file.path(results.dir, filename), sep = ',')
      res1 <- na.omit(res1)
      res.surv.pws[[g]] <- list()
      for(st in surv.types) {
        cols <- c(colnames(res1)[1],
                  grep(paste0('^', st, '_'), colnames(res1), value = T))
        res.surv.pws[[g]][[st]] <- res1[,..cols]
      }
    }
    return(res.surv.pws)
  }
  
  # Read gene ID-to-symbol mapping
  mapping <- fread(file.path(results.dir,
                             paste0('gene_id_mapping_', program, '.csv')))
  # Also read the original HGNC table used to create it,
  # because we may need to remove non-coding genes
  hgnc <- fread(file.path('pathways', paste0('hgnc_complete_set_Maxim_For_OB_',
                                             'samples_marking_13.07.17.txt')))
  hgnc <- hgnc[,c('symbol', 'locus_group')]
  hgnc <- hgnc[hgnc$locus_group != 'withdrawn',]
  
  # Load analysis results for genes
  res.surv.genes <- list()
  for(g in tumor.types) {
    g0 <- g
    for(c in bad.chars) {
      g0 <- gsub(c, '-', g0)
    }
    res.surv.genes[[g]] <- list()
    for(st in surv.types) {
      filename <- paste0('results_surv_', program, '_', g0,
                         '_genes_', st, '.csv')
      res1 <- fread(file.path(results.dir, filename))
      # Removing Ensembl ID versions, to be able to map them to gene symbols
      # (e.g. ENSG00000000003.13 -> ENSG00000000003)
      res1$gene <- gsub('\\.(.)*', '', res1$gene)
      # Attaching the symbols
      if(!all(res1$gene == mapping$ensembl_gene_id)) {
        message('Error: inconsistent gene IDs')
        return(NULL)
      }
      res1$symbol <- mapping$gene_symbol
      #
      # Removing the gene types we don't need
      res1 <- merge(res1, hgnc, by = 'symbol', all.x = T, all.y = F)
      res1 <- res1[(res1$locus_group %in% keep.loci),]
      cols <- colnames(res1)[colnames(res1) != 'locus_group']
      res1 <- res1[, ..cols]
      #
      # Sorting, excluding NAs
      id.cols <- c('symbol', 'gene')
      col.order <- c(id.cols, colnames(res1)[!(colnames(res1) %in% id.cols)])
      res1 <- res1[, ..col.order]
      res1 <- na.omit(res1[order(res1$symbol),])
      #
      # Readjusting p-values
      res1[[paste0(st,
                   '_pval.adj')]] <- p.adjust(res1[[paste0(st, '_pval')]],
                                              method = 'fdr')
      res1[[paste0(st,
                   '_hr.pval.adj')]] <- p.adjust(res1[[paste0(st, '_hr.pval')]],
                                                 method = 'fdr')
      #
      # Recording
      res.surv.genes[[g]][[st]] <- res1
    }
  }
  
  # Load pathway analysis results
  res.surv.pws <- list()
  for(g in tumor.types) {
    g0 <- g
    for(c in bad.chars) {
      g0 <- gsub(c, '-', g0)
    }
    filename <- paste0('results_surv_', program, '_', g0, '.csv')
    res1 <- fread(file.path(results.dir, filename))
    res1$symbol <- gsub('^centr_pathway_', '', res1$pathway)
    res1$symbol <- gsub('.xlsx$', '', res1$symbol)
    id.cols <- c('symbol', 'pathway')
    col.order <- c(id.cols, colnames(res1)[!(colnames(res1) %in% id.cols)])
    res1 <- res1[, ..col.order]
    res1 <- na.omit(res1[order(res1$symbol),])
    #res.surv.pws[[g]] <- res1
    res.surv.pws[[g]] <- list()
    for(st in surv.types) {
      cols <- c(colnames(res1)[1:length(id.cols)],
                grep(paste0('^', st, '_'), colnames(res1), value = T))
      res.surv.pws[[g]][[st]] <- res1[,..cols]
    }
  }
  
  # Just a quick check
  flags <- logical()
  for(g in tumor.types) {
    for(st in surv.types) {
      c1 <- colnames(res.surv.genes[[g]][[st]])[-c(1,2)]
      c2 <- colnames(res.surv.pws[[g]][[st]])[-c(1,2)]
      flags <- c(flags, all(c1 == c2))
    }
  }
  if(!all(flags)) {
    message('Warning: Column names consistency ', all(flags))
  }
  
  # Make a joint gene-pathway results table
  res.surv <- list()
  for(g in tumor.types) {
    res.surv[[g]] <- list()
    for(st in surv.types) {
      fields <- paste0(st, c('_pval.adj', '_hr',
                             '_hr.ci.lower', '_hr.ci.upper',
                             '_hr.pval', '_hr.pval.adj'))
      cols <- c('symbol', 'gene', fields)
      dt1 <- res.surv.genes[[g]][[st]][, ..cols]
      cols <- c('pathway', 'symbol', fields)
      dt2 <- res.surv.pws[[g]][[st]][, ..cols]
      colnames(dt1)[-c(1,2)] <- paste0('gene_', colnames(dt1)[-c(1,2)])
      colnames(dt2)[-c(1,2)] <- paste0('pw_', colnames(dt2)[-c(1,2)])
      res1 <- merge(dt1, dt2)
      id.cols <- c('symbol', 'gene', 'pathway')
      col.order <- c(id.cols,
                     colnames(res1)[!(colnames(res1) %in% id.cols)])
      res.surv[[g]][[st]] <- res1[,..col.order]
    }
  }
  # Return
  res <- list()
  res[['g']] <- res.surv.genes
  res[['pw']] <- res.surv.pws
  res[['gpw']] <- res.surv
  return(res)
}

read.diff <- function(results.dir, program, tumor.types,
                      mode = c('gc', 'classic'),
                      keep.loci = c('protein-coding gene', 'pseudogene',
                                    'non-coding RNA', 'other')) {
  ## Reads and formats diff. results for:
  ## genes and gene-centric pathways ('gc' mode), or classical pathways.
  ## gc means gene-centric, classical is self-explanatory
  
  library('data.table')
  
  if(mode == 'classic') {
    res <- fread(file.path(results.dir,
                           paste0('results_diff_', program, '.csv')), sep = ',')
    res <- na.omit(res)
    return(res)
  }
  
  # Read gene ID-to-symbol mapping
  mapping <- fread(file.path(results.dir,
                             paste0('gene_id_mapping_', program, '.csv')))
  # Also read the original HGNC table used to create it,
  # because we may need to remove non-coding genes
  hgnc <- fread(file.path('pathways', paste0('hgnc_complete_set_Maxim_For_OB_',
                                             'samples_marking_13.07.17.txt')))
  hgnc <- hgnc[,c('symbol', 'locus_group')]
  hgnc <- hgnc[hgnc$locus_group != 'withdrawn',]
  
  # Load analysis results
  filepath.g <- file.path(results.dir, paste0('results_diff_', program,
                                              '_genes.csv'))
  filepath.pw <- file.path(results.dir, paste0('results_diff_', program,
                                               '.csv'))
  if(!file.exists(filepath.g) | !file.exists(filepath.pw)) {
    message('Error in analyzeResults.diff(): results files are missing, ',
            'aborting.')
    return()
  } else {
    res.diff.genes <- fread(filepath.g)
    res.diff.pws <- fread(filepath.pw)
  }
  
  ## Add gene symbols
  colnames(res.diff.genes)[1] <- 'gene'
  # Removing Ensembl ID versions, to be able to map them to gene symbols
  # (e.g. ENSG00000000003.13 -> ENSG00000000003)
  res.diff.genes$gene <- gsub('\\.(.)*', '', res.diff.genes$gene)
  # Attaching the symbols
  # Attaching the symbols
  if(!all(res.diff.genes$gene == mapping$ensembl_gene_id)) {
    message('Error: inconsistent gene IDs')
    return()
  }
  res.diff.genes$symbol <- mapping$gene_symbol
  #
  # Removing the gene types we don't need
  res.diff.genes <- merge(res.diff.genes, hgnc, by = 'symbol',
                          all.x = T, all.y = F)
  res.diff.genes <- res.diff.genes[(res.diff.genes$locus_group %in% keep.loci),]
  cols <- colnames(res.diff.genes)[colnames(res.diff.genes) != 'locus_group']
  res.diff.genes <- res.diff.genes[, ..cols]
  #
  # Sorting, excluding NAs
  id.cols <- c('symbol', 'gene')
  col.order <- colnames(res.diff.genes)[!(colnames(res.diff.genes) %in% id.cols)]
  col.order <- c(id.cols, col.order)
  res.diff.genes <- res.diff.genes[, ..col.order]
  res.diff.genes <- res.diff.genes[order(res.diff.genes$symbol,
                                         res.diff.genes$gene),]
  res.diff.genes <- na.omit(res.diff.genes)
  #
  # Readjusting p-values
  res.diff.genes$w.pval.adj <- p.adjust(res.diff.genes$w.pval, method = 'fdr')
  #
  
  ## Load pathway analysis results
  res.diff.pws$symbol <- gsub('^centr_pathway_', '', res.diff.pws$pathway)
  res.diff.pws$symbol <- gsub('.xlsx$', '', res.diff.pws$symbol)
  id.cols <- c('symbol', 'pathway')
  col.order <- colnames(res.diff.pws)[!(colnames(res.diff.pws) %in% id.cols)]
  col.order <- c(id.cols, col.order)
  res.diff.pws <- na.omit(res.diff.pws[, ..col.order])
  res.diff.pws <- res.diff.pws[order(res.diff.pws$symbol),]
  
  # Combine
  dt1 <- res.diff.genes
  dt2 <- res.diff.pws
  colnames(dt1)[-c(1,2)] <- paste0('gene_', colnames(dt1)[-c(1,2)])
  colnames(dt2)[-c(1,2)] <- paste0('pw_', colnames(dt2)[-c(1,2)])
  res.diff <- merge(dt1, dt2)
  id.cols <- c('symbol', 'gene', 'pathway')
  col.order <- c(id.cols,
                 colnames(res.diff)[!(colnames(res.diff) %in% id.cols)])
  res.diff <- res.diff[,..col.order]
  # Return
  res <- list()
  res[['g']] <- res.diff.genes
  res[['pw']] <- res.diff.pws
  res[['gpw']] <- res.diff
  return(res)
}




analyzeResults.surv.gc <- function(results.dir, program, tumor.types,
                                   pval.max = 0.05,
                                   make.plots = TRUE, log_y = FALSE,
                                   sample_size = 200,
                                   keep_loci = c('protein-coding gene',
                                                 'pseudogene',
                                                 'non-coding RNA',
                                                 'other')) {
  ## Summarizes survival results for genes and gene-centric pathways
  ## for a single given localization
  
  library('data.table')
  library('cowplot')
  
  if(program %in% c('TCGA', 'CPTAC')) {
    surv.types <- c('os', 'pfs')
  } else {
    surv.types <- c('os')
  }
  bad.chars <- c('/', ',')
  
  # Read the raw analysis results
  res <- read.surv(results.dir, program, tumor.types,
                   mode = 'gc', keep.loci = keep_loci)
  if(is.null(res)) {
    message('Error: failed to read analysis results')
    return()
  }
  res.surv.genes <- res[['g']]
  res.surv.pws <- res[['pw']]
  res.surv <- res[['gpw']]
  
  ## SUMMARIZE
  summ <- list()
  res.surv.signif <- list()
  for(g in tumor.types) {
    summ[[g]] <- list()
    res.surv.signif[[g]] <- list()
    for(st in surv.types) {
      #
      # Printing a report, recording a summary
      summ[[g]][[st]] <- list()
      for(i in 1:2) { message('#') }
      message(paste('##', g, st))
      pfield1 <- paste0(st, '_pval.adj')
      #pfield2 <- paste0(st, '_hr.pval.adj')
      ci1 <- paste0(st, '_hr.ci.lower')
      ci2 <- paste0(st, '_hr.ci.upper')
      n.genes <- nrow(res.surv.genes[[g]][[st]])
      n.pws <- nrow(res.surv.pws[[g]][[st]])
      n.common <- nrow(res.surv[[g]][[st]])
      summ[[g]][[st]][['common.n.total']] <- n.common
      # gene count:
      cond1 <- (res.surv.genes[[g]][[st]][[pfield1]] < pval.max)
      cond2 <- (res.surv.genes[[g]][[st]][[ci1]] < 1.0 &
                  res.surv.genes[[g]][[st]][[ci2]] < 1.0) | (
                    res.surv.genes[[g]][[st]][[ci1]] > 1.0 &
                      res.surv.genes[[g]][[st]][[ci2]] > 1.0
                  )
      message('# Significant genes: ', sum(cond1 & cond2), ' / ', n.genes,
              ' (', round(100.0*sum(cond1 & cond2)/n.genes, 2), '%)')
      summ[[g]][[st]][['genes.n']] <- sum(cond1 & cond2)
      summ[[g]][[st]][['genes.n.total']] <- n.genes
      # gene top:
      printIfSmall(res.surv.genes[[g]][[st]]$symbol, cond1 & cond2)
      cols <- c('symbol', pfield1)
      tg <- res.surv.genes[[g]][[st]][,..cols]
      tg <- tg[cond1 & cond2,]
      tg <- tg[order(tg[[pfield1]]),]
      tg <- tg[['symbol']]
      tg <- tg[1:min(10, length(tg))]
      summ[[g]][[st]][['genes.top']] <- tg
      # genes with a pw:
      ci1.g <- paste0('gene_', ci1)
      ci2.g <- paste0('gene_', ci2)
      cond3 <- (res.surv[[g]][[st]][[paste0('gene_', pfield1)]] < pval.max)
      cond4 <- (res.surv[[g]][[st]][[ci1.g]] < 1.0 &
                  res.surv[[g]][[st]][[ci2.g]] < 1.0) | (
                    res.surv[[g]][[st]][[ci1.g]] > 1.0 &
                      res.surv[[g]][[st]][[ci2.g]] > 1.0
                  )
      message('# Including genes with a pathway: ', sum(cond3 & cond4),
              ' / ', n.common,
              ' (', round(100.0*sum(cond3 & cond4)/n.common, 2), '%)')
      summ[[g]][[st]][['genes.wpw.n']] <- sum(cond3 & cond4)
      printIfSmall(res.surv[[g]][[st]]$symbol, cond3 & cond4)
      # pathway count:
      cond1 <- (res.surv.pws[[g]][[st]][[pfield1]] < pval.max)
      cond2 <- (res.surv.pws[[g]][[st]][[ci1]] < 1.0 &
                  res.surv.pws[[g]][[st]][[ci2]] < 1.0) | (
                    res.surv.pws[[g]][[st]][[ci1]] > 1.0 &
                      res.surv.pws[[g]][[st]][[ci2]] > 1.0
                  )
      message('# Significant pathways: ', sum(cond1 & cond2), ' / ', n.pws,
              ' (', round(100.0*sum(cond1 & cond2)/n.pws, 2), '%)')
      summ[[g]][[st]][['pws.n']] <- sum(cond1 & cond2)
      summ[[g]][[st]][['pws.n.total']] <- n.pws
      # pathway top:
      printIfSmall(res.surv.pws[[g]][[st]]$symbol, cond1 & cond2)
      cols <- c('symbol', pfield1)
      tg <- res.surv.pws[[g]][[st]][,..cols]
      tg <- tg[cond1 & cond2,]
      tg <- tg[order(tg[[pfield1]]),]
      tg <- tg[['symbol']]
      tg <- tg[1:min(10, length(tg))]
      summ[[g]][[st]][['pws.top']] <- tg
      # pws with a gene:
      ci1.pw <- paste0('pw_', ci1)
      ci2.pw <- paste0('pw_', ci2)
      cond5 <- (res.surv[[g]][[st]][[paste0('pw_', pfield1)]] < pval.max)
      cond6 <- (res.surv[[g]][[st]][[ci1.pw]] < 1.0 &
                  res.surv[[g]][[st]][[ci2.pw]] < 1.0) | (
                    res.surv[[g]][[st]][[ci1.pw]] > 1.0 &
                      res.surv[[g]][[st]][[ci2.pw]] > 1.0
                  )
      message('# Including pathways with a gene: ', sum(cond5 & cond6),
              ' / ', n.common,
              ' (', round(100.0*sum(cond5 & cond6)/n.common, 2), '%)')
      summ[[g]][[st]][['pws.wg.n']] <- sum(cond5 & cond6)
      printIfSmall(res.surv[[g]][[st]]$symbol, cond5 & cond6)
      # Identifying the stable fraction
      hrfield.g <- paste0('gene_', st, '_hr')
      hrfield.pw <- paste0('pw_', st, '_hr')
      cond7 <- (res.surv[[g]][[st]][[hrfield.g]] < 1.0)
      cond8 <- (res.surv[[g]][[st]][[hrfield.pw]] < 1.0)
      cond9 <- (res.surv[[g]][[st]][[hrfield.g]] > 1.0)
      cond10 <- (res.surv[[g]][[st]][[hrfield.pw]] > 1.0)
      cond11 <- (cond7 & cond8) | (cond9 & cond10)
      cond12 <- cond3 & cond4 & cond5 & cond6
      cond13 <- cond11 & cond12
      # Significance
      if(sum(cond12) != 0) {
        p.overlap <- overlapSignif(sum(cond12),
                                   1:nrow(res.surv[[g]][[st]]),
                                   1:nrow(res.surv[[g]][[st]]),
                                   sum(cond3 & cond4), sum(cond5 & cond6))
        p.message <- paste0(' (p = ', p.overlap, ')')
        summ[[g]][[st]][['stable.pval']] <- p.overlap
      } else {
        p.message <- ''
      }
      message('# Stable genes: ', sum(cond12), ' / ', n.common, ' (',
              round(100.0*sum(cond12)/n.common, 2), '%) ', p.message)
      summ[[g]][[st]][['stable.n']] <- sum(cond12)
      printIfSmall(res.surv[[g]][[st]]$symbol, cond12)
      message('# Consistent stable genes: ', sum(cond13), ' / ', n.common,
              ' (', round(100.0*sum(cond13)/n.common, 2), '%)')
      summ[[g]][[st]][['stable.cst.n']] <- sum(cond13)
      printIfSmall(res.surv[[g]][[st]]$symbol, cond13)
      #
      #
      # Recording
      res.surv.signif[[g]][[st]] <- res.surv[[g]][[st]][cond12,]
    }
  }
  
  ## How about plotting
  if(make.plots) {
    out.dir <- file.path(results.dir, 'analysis')
    if(!dir.exists(out.dir)) {
      dir.create(out.dir)
    }
    # General
    plots.total <- plotBoxes(res.surv, 0, log.y = log_y)
    # Significant
    plots.signif <- tryCatch(
      plotBoxes(res.surv.signif, 99999, log.y = log_y, violins = FALSE),
      error = function(e) {
        message('Warning: stack overflow in plotBoxes(), ',
                'only ', sample_size, ' random lines will be shown')
        plotBoxes(res.surv.signif, sample_size, log.y = log_y, violins = FALSE)
      }
    )
    #
    # Writing
    for(g in tumor.types) {
      g0 <- g
      for(c in bad.chars) {
        g0 <- gsub(c, '-', g0)
      }
      for(st in surv.types) {
        exists.1 <- !is.null(plots.total[[g]][[st]])
        if(exists.1) {
          plot.file <- file.path(out.dir, paste0('plot_', g0, '_', st,
                                                 '_total.png'))
          ggsave(plot.file, plot = plots.total[[g]][[st]],
                 width = 1200, height = 1000, units = 'px',
                 dpi = 300)
          #message('Plot for ', g, ' ', st, ' (total)',
          #        ' written to file ', plot.file)
        }
        exists.2 <- !is.null(plots.signif[[g]][[st]])
        if(exists.2) {
          plot.file <- file.path(out.dir, paste0('plot_', g0, '_', st,
                                                 '_signif.png'))
          ggsave(plot.file, plot = plots.signif[[g]][[st]],
                 width = 1200, height = 1000, units = 'px',
                 dpi = 300)
          #message('Plot for ', g, ' ', st, ' (significant only)',
          #        ' written to file ', plot.file)
        }
        if(exists.1 & exists.2) {
          lims1 <- layer_scales(plots.signif[[g]][[st]])$y$get_limits()
          lims2 <- layer_scales(plots.total[[g]][[st]])$y$get_limits()
          lim.min <- min(lims1[1], lims2[1])
          lim.max <- max(lims1[2], lims2[2])
          plots.total[[g]][[st]] <- plots.total[[g]][[st]] +
            ylim(lim.min, lim.max) + ggtitle('All')
          plots.signif[[g]][[st]] <- plots.signif[[g]][[st]] +
            ylim(lim.min, lim.max) + ggtitle('Significant only') + ylab(NULL)
          pg <- plot_grid(plots.total[[g]][[st]],
                          plots.signif[[g]][[st]],
                          align = 'h', axis = 'tblr')
          plot.file <- file.path(out.dir, paste0('plot_', g0, '_', st,
                                                 '_all.png'))
          ggsave(plot.file, plot = pg,
                 width = 2500, height = 1000, units = 'px',
                 dpi = 300)
          #message('Two-panel plot for ', g, ' ', st,
          #        ' written to file ', plot.file)
        }
      }
    }
  }
  
  return(summ)
}



analyzeResults.surv.trad <- function(results.dir, program, tumor.types,
                                     pval.max = 0.05) {
  ## Summarizes survival results for classical pathways
  ## for a single given localization
  
  library('data.table')
  
  if(program == 'TCGA') {
    surv.types <- c('os', 'pfs')
  } else {
    surv.types <- c('os')
  }
  bad.chars <- c('/', ',')
  
  # Load pathway analysis results
  res.surv.pws <- read.surv(file.path(results.dir, 'tradpws'),
                            program, tumor.types,
                            mode = 'classic')
  
  # Extract significant results,
  # prepare a summary
  summ <- list()
  for(g in tumor.types) {
    summ[[g]] <- list()
    for(st in surv.types){
      n.pws <- nrow(res.surv.pws[[g]][[st]])
      #
      #
      # Printing a report, recording a summary
      summ[[g]][[st]] <- list()
      for(i in 1:2) { message('#') }
      message(paste('##', g, st))
      pfield1 <- paste0(st, '_pval.adj')
      #pfield2 <- paste0(st, '_hr.pval.adj')
      ci1 <- paste0(st, '_hr.ci.lower')
      ci2 <- paste0(st, '_hr.ci.upper')
      # pathway count:
      cond1 <- (res.surv.pws[[g]][[st]][[pfield1]] < pval.max)
      cond2 <- (res.surv.pws[[g]][[st]][[ci1]] < 1.0 &
                  res.surv.pws[[g]][[st]][[ci2]] < 1.0) | (
                    res.surv.pws[[g]][[st]][[ci1]] > 1.0 &
                      res.surv.pws[[g]][[st]][[ci2]] > 1.0
                  )
      message('# Significant pathways: ', sum(cond1 & cond2), ' / ', n.pws,
              ' (', round(100.0*sum(cond1 & cond2)/n.pws, 2), '%)')
      summ[[g]][[st]][['trad.pws.n']] <- sum(cond1 & cond2)
      summ[[g]][[st]][['trad.pws.n.total']] <- n.pws
      # pathway top:
      printIfSmall(res.surv.pws[[g]][[st]]$pathway, cond1 & cond2)
      cols <- c('pathway', pfield1)
      tg <- res.surv.pws[[g]][[st]][,..cols]
      tg <- tg[cond1 & cond2,]
      tg <- tg[order(tg[[pfield1]]),]
      tg <- tg[['pathway']]
      tg <- tg[1:min(10, length(tg))]
      summ[[g]][[st]][['trad.pws.top']] <- tg
    }
  }
  return(summ)
}



analyzeResults.diff <- function(results.dir, program, tumor.types,
                                pval.max = 0.05,
                                keep_loci = c('protein-coding gene',
                                              'pseudogene',
                                              'non-coding RNA',
                                              'other')) {
  library('data.table')
  
  # Read and format raw results
  res <- read.diff(results.dir, program, tumor.types,
                   mode = 'gc', keep.loci = keep_loci)
  res.diff.genes <- res[['g']]
  res.diff.pws <- res[['pw']]
  res.diff <- res[['gpw']]
  #
  #
  ## The report
  for(i in 1:2) { message('#') }
  n.genes <- nrow(res.diff.genes)
  n.common <- nrow(res.diff)
  n.pws <- nrow(res.diff.pws)
  if('w.pval.adj' %in% colnames(res.diff.genes)) {
    cond1 <- (res.diff.genes[['w.pval.adj']] < pval.max)
  } else {
    cond1 <- rep(TRUE, nrow(res.diff.genes))
  }
  #cond1 <- (res.diff.genes[['w.pval.adj']] < pval.max)
  cond2 <- (res.diff.genes[['auc']] >= 0.8)
  cond3 <- (res.diff.genes[['auc']] >= 0.9)
  message('# Significant genes with AUC >= 0.8: ', sum(cond1 & cond2), ' / ',
          n.genes, ' (', round(100.0*sum(cond1 & cond2)/n.genes, 2), '%)')
  printIfSmall(res.diff.genes[['symbol']], cond1 & cond2)
  cond4 <- (res.diff.genes[['symbol']][cond1 & cond2] %in% res.diff$symbol)
  message('# Including genes with a pathway: ', sum(cond4), ' / ',
          n.common, '( ', round(100.0*sum(cond4)/n.common, 2), '%)')
  printIfSmall(res.diff.genes[['symbol']][cond1 & cond2], cond4)
  message('# Significant genes with AUC >= 0.9: ', sum(cond1 & cond3), ' / ',
          n.genes, ' (', round(100.0*sum(cond1 & cond3)/n.genes, 2), '%)')
  printIfSmall(res.diff.genes[['symbol']], cond1 & cond3)
  cond5 <- (res.diff.genes[['symbol']][cond1 & cond3] %in% res.diff$symbol)
  message('# Including genes with a pathway: ', sum(cond5), ' / ',
          n.common, '( ', round(100.0*sum(cond5)/n.common, 2), '%)')
  printIfSmall(res.diff.genes[['symbol']][cond1 & cond3], cond5)
  message('#')
  #
  if('w.pval.adj' %in% colnames(res.diff.pws)) {
    cond6 <- (res.diff.pws[['w.pval.adj']] < pval.max)
  } else {
    cond6 <- rep(TRUE, nrow(res.diff.pws))
  }
  #cond6 <- (res.diff.pws[['w.pval.adj']] < pval.max)
  cond7 <- (res.diff.pws[['auc']] >= 0.8)
  cond8 <- (res.diff.pws[['auc']] >= 0.9)
  message('# Significant pathways with AUC >= 0.8: ', sum(cond6 & cond7), ' / ',
          n.pws, ' (', round(100.0*sum(cond6 & cond7)/n.pws, 2), '%)')
  printIfSmall(res.diff.pws[['symbol']], cond6 & cond7)
  cond9 <- (res.diff.pws[['symbol']][cond6 & cond7] %in% res.diff$symbol)
  message('# Including pathways with a gene: ', sum(cond9), ' / ',
          n.common, ' (', round(100.0*sum(cond9)/n.common, 2), '%)')
  printIfSmall(res.diff.pws[['symbol']][cond6 & cond7], cond9)
  message('# Significant pathways with AUC >= 0.9: ', sum(cond6 & cond8), ' / ',
          n.pws, ' (', round(100.0*sum(cond6 & cond8)/n.pws, 2), '%)')
  printIfSmall(res.diff.pws[['symbol']], cond6 & cond8)
  cond10 <- (res.diff.pws[['symbol']][cond6 & cond8] %in% res.diff$symbol)
  message('# Including pathways with a gene: ', sum(cond10), ' / ',
          n.common, ' (', round(100.0*sum(cond10)/n.common, 2), '%)')
  printIfSmall(res.diff.pws[['symbol']][cond6 & cond8], cond10)
  message('#')
  #
  if('gene_w.pval.adj' %in% colnames(res.diff)) {
    cond11 <- (res.diff[['gene_w.pval.adj']] < pval.max)
  } else {
    cond11 <- rep(TRUE, nrow(res.diff))
  }
  #cond11 <- (res.diff[['gene_w.pval.adj']] < pval.max)
  cond12 <- (res.diff[['gene_auc']] >= 0.8)
  cond13 <- (res.diff[['gene_auc']] >= 0.9)
  if('pw_w.pval.adj' %in% colnames(res.diff)) {
    cond14 <- (res.diff[['pw_w.pval.adj']] < pval.max)
  } else {
    cond14 <- rep(TRUE, nrow(res.diff))
  }
  #cond14 <- (res.diff[['pw_w.pval.adj']] < pval.max)
  cond15 <- (res.diff[['pw_auc']] >= 0.8)
  cond16 <- (res.diff[['pw_auc']] >= 0.9)
  #
  # Stability and significance for 0.8
  if(sum(cond11 & cond12 & cond14 & cond15) != 0) {
    p.overlap <- overlapSignif(sum(cond11 & cond12 & cond14 & cond15),
                               1:nrow(res.diff),
                               1:nrow(res.diff),
                               sum(cond11 & cond12), sum(cond14 & cond15))
    p.message <- paste0(' (p = ', p.overlap, ')')
  } else {
    p.message <- ''
  }
  message('# Stable genes with AUC >= 0.8: ',
          sum(cond11 & cond12 & cond14 & cond15), ' / ', n.common, ' (',
          round(100.0*sum(cond11 & cond12 & cond14 & cond15)/n.common, 2),
          ' %)', p.message)
  printIfSmall(res.diff[['symbol']], cond11 & cond12 & cond14 & cond15)
  #
  # Stability and significance for 0.9
  if(sum(cond11 & cond13 & cond14 & cond16) != 0) {
    p.overlap <- overlapSignif(sum(cond11 & cond13 & cond14 & cond16),
                               1:nrow(res.diff),
                               1:nrow(res.diff),
                               sum(cond11 & cond13), sum(cond14 & cond16))
    p.message <- paste0(' (p = ', p.overlap, ')')
  } else {
    p.message <- ''
  }
  message('# Stable genes with AUC >= 0.9: ',
          sum(cond11 & cond13 & cond14 & cond16), ' / ', n.common, ' (',
          round(100.0*sum(cond11 & cond13 & cond14 & cond16)/n.common, 2),
          ' %)', p.message)
  printIfSmall(res.diff[['symbol']], cond11 & cond13 & cond14 & cond16)
  message('#')
  #
  w1 <- wilcox.test(res.diff[['pw_auc']],
                    res.diff[['gene_auc']],
                    alternative = 'greater')
  message('# Wilcox test p-value for all: ', w1$p.value)
  if(sum(cond11 & cond14) >= 10) {
    w2 <- wilcox.test(res.diff[['pw_auc']][cond11 & cond14],
                      res.diff[['gene_auc']][cond11 & cond14],
                      alternative = 'greater')
    message('# Wilcox test p-value for stable signif: ', w2$p.value)
  }
  if(sum(cond11 & cond12 & cond14 & cond15) >= 10) {
    w3 <- wilcox.test(res.diff[['pw_auc']][cond11 & cond12 & cond14 & cond15],
                      res.diff[['gene_auc']][cond11 & cond12 & cond14 & cond15],
                      alternative = 'greater')
    message('# Wilcox test p-value for stable signif with AUC >= 0.8: ',
            w3$p.value)
  }
  if(sum(cond11 & cond13 & cond14 & cond16) >= 10) {
    w4 <- wilcox.test(res.diff[['pw_auc']][cond11 & cond13 & cond14 & cond16],
                      res.diff[['gene_auc']][cond11 & cond13 & cond14 & cond16],
                      alternative = 'greater')
    message('# Wilcox test p-value for stable signif with AUC >= 0.9: ',
            w4$p.value)
  }
  
}

