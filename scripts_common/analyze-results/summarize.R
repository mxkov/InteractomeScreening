## Makes a summary of survival results for all localizations;
## a lot of downstream analysis, plots, etc. will be based on this summary file


source('scripts_common/analyze-results/analyzing-functions.R')
source('scripts_common/locos-info.R')
source('scripts_common/input-locos.R')

summarize <- function(summ.file,
                      do.diff = FALSE, report.diff = FALSE,
                      do.surv.gc = TRUE, do.surv.trad = TRUE,
                      report.surv.gc = TRUE, report.surv.trad = TRUE,
                      make_plots = TRUE,
                      keep.loci = c('protein-coding gene', 'non-coding RNA'),
                      res.dirname = 'results_TCGA') {
  #locos.todo <- names(gs)    # gs is imported from locos-info.R
  locos.todo <- locos         # locos is imported from input-locos.R
  
  # Filling in the tumor types and result dir paths
  res.dirs <- list()
  for(loco in names(gs)) {
    res.dirs[[loco]] <- file.path('locos', loco, res.dirname)
  }
  
  s <- list()
  t <- list()
  d <- list()
  for(loco in locos.todo) {
    outdir <- file.path(res.dirs[[loco]], 'analysis')
    dir.create(outdir, showWarnings = FALSE)
    #
    # Diff analysis, if applicable
    if(do.diff & length(gs[[loco]]) > 1) {
      if(report.diff) {
        outfile <- file.path(outdir, 'report_diff.txt')
        capture.output((d1 <- analyzeResults.diff(res.dirs[[loco]],
                                                  'TCGA', gs[[loco]],
                                                  keep_loci = keep.loci)),
                       file = outfile, type = 'message')
        message('Report for ', loco, ' written to file ', outfile)
      } else {
        d1 <- analyzeResults.diff(res.dirs[[loco]], 'TCGA', gs[[loco]],
                                  keep_loci = keep.loci)
      }
      d[[loco]] <- d1
    }
    #
    # Surv analysis
    if(do.surv.gc) {
      message('Processing ', loco, '...')
      if(report.surv.gc) {
        outfile <- file.path(outdir, paste0('report_surv.txt'))
        capture.output((s1 <- analyzeResults.surv.gc(res.dirs[[loco]],
                                                     'TCGA', gs[[loco]],
                                                     make.plots = make_plots,
                                                     log_y = FALSE,
                                                     keep_loci = keep.loci)),
                       file = outfile, type = 'message')
        message('Report for ', loco, ' written to file ', outfile)
      } else {
        s1 <- analyzeResults.surv.gc(res.dirs[[loco]],
                                     'TCGA', gs[[loco]],
                                     make.plots = make_plots,
                                     log_y = FALSE,
                                     keep_loci = keep.loci)
      }
      s[[loco]] <- s1
    }
    #
    # Surv analysis for classical ("traditional") pathways
    if(do.surv.trad) {
      message('Processing ', loco, ' (classical pathways)...')
      outdir.trad <- file.path(res.dirs[[loco]], 'tradpws')
      if(!dir.exists(outdir.trad)) {
        message('Directory ', outdir.trad, ' does not exist, skipping ', loco)
        next
      }
      if(report.surv.trad) {
        outdir <- file.path(outdir.trad, 'analysis')
        dir.create(outdir, showWarnings = F)
        outfile <- file.path(outdir, paste0('report_surv.txt'))
        capture.output((t1 <- analyzeResults.surv.trad(res.dirs[[loco]],
                                                       'TCGA', gs[[loco]])),
                       file = outfile, type = 'message')
        message('Report for ', loco, ' written to file ', outfile)
      } else {
        t1 <- analyzeResults.surv.trad(res.dirs[[loco]],
                                       'TCGA', gs[[loco]])
      }
      t[[loco]] <- t1
    }
  }
  
  # Make the final summary
  surv.types <- c('os', 'pfs')
  cols.id <- c('localization', 'group', 'surv.type')
  cols.n <- c('genes.n', 'genes.wpw.n',
              'pws.n', 'pws.wg.n',
              'stable.n', 'stable.cst.n',
              'trad.pws.n')
  cols.n.total <- c('genes.n.total', 'pws.n.total', 'common.n.total',
                    'trad.pws.n.total')
  cols.top <- c('genes.top', 'pws.top', 'trad.pws.top')
  summ.df <- data.frame()
  for(col in cols.id) {
    summ.df[,col] <- character()
  }
  for(col in c(cols.n, cols.n.total)) {
    summ.df[,col] <- numeric()
  }
  for(col in cols.top) {
    summ.df[,col] <- character()
  }
  for(loco in locos.todo) {
    for(g in gs[[loco]]) {
      for(st in surv.types) {
        df <- data.frame(localization = loco, group = g, surv.type = st)
        entries <- list(s[[loco]][[g]][[st]], t[[loco]][[g]][[st]])
        for(entry in entries) {
          for(field in names(entry)) {
            if(field %in% cols.top) {
              df[,field] <- paste0(entry[[field]], collapse = '|')
            } else {
              df[,field] <- entry[[field]]
            }
          }
        }
        col.order <- colnames(summ.df)
        summ.df <- rbind(summ.df, df[,col.order])
      }
    }
  }
  write.csv(summ.df, summ.file, quote = F, row.names = F)
  message('')
  message('Final summary written to file ', summ.file)
}


summarize(file.path('locos', 'summary.csv'))
