## Wilcoxon's tests for distributions of hazard ratios

source('scripts_common/analyze-results/analyzing-functions.R')
source('scripts_common/locos-info.R')
library('data.table')


HRwilcox <- function(round.to = 4,
                     keep.loci = c('protein-coding gene', 'non-coding RNA'),
                     metric = c('mean', 'median')) {
  # round.to: how many decimal digits to keep in the output
  #
  surv.types <- c('os', 'pfs')
  locos.todo <- names(gs)    # gs is imported from locos-info.R
  res.dirs <- list()
  res.dirname <- 'results_TCGA'
  for(loco in names(gs)) {
    res.dirs[[loco]] <- file.path('locos', loco, res.dirname)
  }
  if(metric == 'mean') {
    mval <- function(x) { return(mean(x)) }
  } else if(metric == 'median') {
    mval <- function(x) { return(median(x)) }
  }
  
  # Read the summary to see which locos we need to plot
  message('Reading the summary file...')
  summ <- fread(file.path('locos', 'summary.csv'))
  summ.nonempty <- summ[apply(summ[,4:10], 1,
                              function(x) {return(any(unlist(x) != 0))}),]
  summ.nonempty <- summ.nonempty[summ.nonempty$group != 'LUAD',]
  
  # Read the survival results
  message('Reading raw survival results...')
  results <- list()
  for(loco in unique(summ.nonempty$localization)) {
    gs[[loco]] <- gs[[loco]][gs[[loco]] %in% summ.nonempty$group]
    res <- read.surv(res.dirs[[loco]], 'TCGA',
                     gs[[loco]], mode = 'gc',
                     keep.loci = keep.loci)
    res$trad <- read.surv(file.path(res.dirs[[loco]], 'tradpws'), 'TCGA',
                          gs[[loco]], mode = 'classic')
    results[[loco]] <- res
  }
  # Transform the data, extract significant
  message('Filtering and transforming data...')
  id.cols <- c('symbol', 'gene', 'pathway')
  first <- TRUE
  for(loco in names(results)) {
    for(cat in names(results[[loco]])) {
      if(cat == 'gpw') { next }
      if(cat %in% c('g', 'gpw')) {
        id.col <- 'symbol'
      } else {
        id.col <- 'pathway'
      }
      for(g in gs[[loco]]) {
        for(st in surv.types) {
          res <- results[[loco]][[cat]][[g]][[st]]
          # Clean up column names
          colnames(res) <- gsub(paste0('^', st, '_'), '', colnames(res))
          # Extract significant, exclude extra columns
          cond1 <- res$pval.adj < 0.05
          cond2 <- (res$hr.ci.lower-1)*(res$hr.ci.upper-1) > 0
          res <- res[cond1 & cond2,]
          id.cols.extra <- id.cols[!(id.cols == id.col)]
          cols <- colnames(res)[!(colnames(res) %in% id.cols.extra)]
          res <- res[,..cols]
          # Additional columns
          colnames(res)[1] <- 'item'
          res$item.type <- toupper(cat)
          res$group <- g
          res$surv.type <- toupper(st)
          # Sorting columns
          col.order <- c('group', 'surv.type', 'item', 'item.type',
                         'pval.adj', 'hr')
          res <- res[,..col.order]
          # Merging
          if(first) {
            results.all <- res
            first <- FALSE
          } else {
            results.all <- rbind(results.all, res)
          }
        }
      }
    }
  }
  # The joint gene-pathway tables need to be transformed differently
  first <- TRUE
  for(loco in names(results)) {
    for(g in gs[[loco]]) {
      for(st in surv.types) {
        res <- results[[loco]][['gpw']][[g]][[st]]
        filter <- rep(TRUE, nrow(res))
        for(item.lbl in c('gene_', 'pw_')) {
          filter <- filter & (res[[paste0(item.lbl, st, '_pval.adj')]] < 0.05)
          filter <- filter & ((
            res[[paste0(item.lbl, st, '_hr.ci.lower')]] < 1 &
              res[[paste0(item.lbl, st, '_hr.ci.upper')]] < 1) | (
                res[[paste0(item.lbl, st, '_hr.ci.lower')]] > 1 &
                  res[[paste0(item.lbl, st, '_hr.ci.upper')]] > 1)
          )
        }
        res <- res[filter,]
        if(nrow(res) == 0) { next }
        res1 <- rbind(data.table(group = g, surv.type = toupper(st),
                                 item = res$symbol, item.type = 'G',
                                 hr = res[[paste0('gene_', st, '_hr')]]),
                      data.table(group = g, surv.type = toupper(st),
                                 item = res$symbol, item.type = 'PW',
                                 hr = res[[paste0('pw_', st, '_hr')]]))
        if(first) {
          results.gpw <- res1
          first <- FALSE
        } else {
          results.gpw <- rbind(results.gpw, res1)
        }
      }
    }
  }
  
  # TEST
  # First, flip the ratios!
  results.all$hr[results.all$hr < 1] <- 1/results.all$hr[results.all$hr < 1]
  results.gpw$hr[results.gpw$hr < 1] <- 1/results.gpw$hr[results.gpw$hr < 1]
  # Now to the testing
  message('Testing...')
  item.pairs <- list(c('G', 'PW'), c('PW', 'TRAD'), c('G', 'TRAD'))
  item.types <- c('G', 'PW', 'TRAD')
  w <- list()
  first <- rep(TRUE, 2)
  for(g in unique(results.gpw$group)) {
    for(st in toupper(surv.types)) {
      #
      # Pairwise between three categories: w$all
      res0 <- subset(results.all, group == g & surv.type == st)
      for(i in 1:length(item.pairs)) {
        res1 <- subset(res0, item.type %in% item.pairs[[i]])
        if(!all(item.pairs[[i]] %in% res1$item.type)) { next }
        if(1 %in% table(res1$item.type)) { next }
        m <- c(0, 0)
        for(j in 1:2) {
          m[j] <- round(mval(res1$hr[res1$item.type == item.pairs[[i]][j]]),
                        round.to)
        }
        wtest <- wilcox.test(hr ~ item.type, data = res1)
        df <- data.frame(group = g, surv.type = st,
                         pair12 = paste0(item.pairs[[i]], collapse = '-'),
                         v1 = m[1], v2 = m[2],
                         v.diff = round(m[1]-m[2], round.to),
                         pval = signif(wtest$p.value, round.to))
        if(metric == 'mean') {
          colnames(df)[4:6] <- c('mean1', 'mean2', 'mean.diff')
        } else if(metric == 'median') {
          colnames(df)[4:6] <- c('median1', 'median2', 'median.diff')
        }
        if(first[1]) {
          w$all <- df
          first[1] <- FALSE
        } else {
          w$all <- rbind(w$all, df)
        }
      }
      #
      # Between genes and respective gene-centric pathways: w$gpw
      res2 <- subset(results.gpw, group == g & surv.type == st)
      if(!all(item.pairs[[1]] %in% res2$item.type)) { next }
      if(1 %in% table(res2$item.type)) { next }
      m <- c(0, 0)
      for(j in 1:2) {
        m[j] <- round(mval(res2$hr[res2$item.type == item.types[j]]),
                      round.to)
      }
      wtest <- wilcox.test(hr ~ item.type, data = res2)
      df <- data.frame(group = g, surv.type = st,
                       pair12 = paste0(item.types[1:2], collapse = '-'),
                       v1 = m[1], v2 = m[2],
                       v.diff = round(m[1]-m[2], round.to),
                       pval = signif(wtest$p.value, round.to))
      if(metric == 'mean') {
        colnames(df)[4:6] <- c('mean1', 'mean2', 'mean.diff')
      } else if(metric == 'median') {
        colnames(df)[4:6] <- c('median1', 'median2', 'median.diff')
      }
      if(first[2]) {
        w$gpw <- df
        first[2] <- FALSE
      } else {
        w$gpw <- rbind(w$gpw, df)
      }
    }
  }
  # Adjusting p-values
  for(id in names(w)) {
    w[[id]]$pval.adj <- signif(p.adjust(w[[id]]$pval, method = 'fdr'), round.to)
  }
  outdir <- file.path('report', 'HRboxes', 'wilcox')
  dir.create(outdir, showWarnings = F)
  # Writing to files
  for(id in names(w)) {
    outfile <- file.path(outdir, paste0('wilcox_', id, '_', metric, '.csv'))
    fwrite(w[[id]], file = outfile, quote = F, row.names = F)
    message('File written: ', outfile)
  }
  return(w)
}


run <- TRUE
metrics <- c('mean', 'median')
ids <- c('all', 'gpw')
w <- list()
for(mv in metrics) {
  if(run) {
    w[[mv]] <- HRwilcox(round.to = 4, metric = mv)
  } else {
    w[[mv]] <- list()
    for(id in ids) {
      w[[mv]][[id]] <- fread(file.path('report', 'HRboxes', 'wilcox',
                                       paste0('wilcox_', id, '_', mv, '.csv')))
    }
  }
}
#
mv <- 'median'
id <- 'all'
w1 <- w[[mv]][[id]]
diff.f <- paste0(mv, '.diff')
#
#w1[w1$pval.adj < 0.05,]
w1[w1$pval.adj < 0.05 & abs(w1[[diff.f]]) >= 0.1,]


