## Makes the violin plots with hazard ratio distributions

source('scripts_common/analyze-results/analyzing-functions.R')
source('scripts_common/locos-info.R')


# First, a function for generating labels with numbers of items:
makeNAnnot <- function(res, g, item.labs, survtypes, plot = c(1,2)) {
  gs.best.lbl.y1 <- list('HNSC' = c(2.2, 2.1, 1.95),
                         'KIRC' = c(3.5, 3.3, 2.9, 3.0, 3.0, 2.5),
                         'KIRP' = c(4.75, 4.6, 4.3, 4.3, 3.9, 3.8),
                         'LIHC' = c(3.35, 3.0, 3.0, 2.5, 2.5, 2.5),
                         'PRAD' = c(5.2, 3.6, 3.6))
  gs.best.lbl.y2 <- list('HNSC' = c(2.15, 2.1),
                         'KIRC' = c(3.3, 3.3, 3.0, 3.0),
                         'KIRP' = c(4.0, 4.0, 3.7, 3.2),
                         'LIHC' = c(3.0, 3.0, 2.5, 2.5),
                         'PRAD' = c(3.6, 3.6))
  n.annot <- data.frame(item.type = rep(item.labs, length(survtypes)),
                        surv.type = rep(survtypes, length(item.labs)),
                        n = '')
  n.annot$surv.type <- n.annot$surv.type[order(n.annot$surv.type)]
  if(plot == 1) {
    n.annot$y <- gs.best.lbl.y1[[g]]
  } else if(plot == 2) {
    n.annot$y <- gs.best.lbl.y2[[g]]
  }
  for(j in 1:nrow(n.annot)) {
    filter1 <- (res$item.type == n.annot$item.type[j])
    filter2 <- (res$surv.type == n.annot$surv.type[j])
    n.annot$n[j] <- paste0('n = ', nrow(res[filter1 & filter2,]))
  }
  return(n.annot)
}


plotHRs <- function(outdir, bigfig = TRUE, squished = FALSE,
                    keep.loci = c('protein-coding gene',
                                  'non-coding RNA')) {
  library('data.table')
  library('ggplot2')
  library('ggpubr')
  print(keep.loci)
  dir.create(outdir, recursive = T, showWarnings = F)
  
  surv.types <- c('os', 'pfs')
  locos.todo <- names(gs)    # gs is imported from locos-info.R
  
  res.dirname <- 'results_TCGA'
  res.dirs <- list()
  for(loco in names(gs)) {
    res.dirs[[loco]] <- file.path('locos', loco, res.dirname)
  }
  
  # Read the summary to see which locos we need to plot
  message('Reading the summary file...')
  summ <- fread(file.path('locos', 'summary.csv'))
  summ.nonempty <- summ[apply(summ[,4:10], 1,
                              function(x) {return(any(unlist(x) != 0))}),]
  
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
          filter1 <- res$pval.adj < 0.05
          filter2 <- (res$hr.ci.lower-1)*(res$hr.ci.upper-1) > 0
          res <- res[filter1 & filter2,]
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
  
  # Now plot
  # one plot for lung (2x2), one for kidney (2x2),
  # one for headneck and liver (2x2),
  # and mb another one.
  #
  # do this:
  # - bar plot for all summ.nonempty;
  # - selected plots for:
  # LUSC: pfs only, no cross
  # PRAD: pfs only, + cross
  # LIHC: os+pfs, both with cross
  # KIRC: os+pfs, both with cross
  # KIRP: os+pfs, pfs with cross
  # HNSC: os, + cross
  # no bars: BLCA os, HNSC pfs, LUAD os+pfs, LUSC pfs
  # discard: COAD, UCEC, SKCM, STAD
  
  # Okay. show the full plot (os+pfs, +cross) for LIHC and KIRC, the best ones.
  # Maybe also PRAD pfs, KIRP pfs, and HNSC os.
  # THe rest -- in suppl with no crosses: KIRP os, LUSC pfs, ...
  
  # Declaring titles/labels, some other plot parameters
  gs.best <- c('KIRC', 'KIRP', 'LIHC', 'HNSC', 'PRAD')    # order matters!
  title <- c('HNSC' = 'Head and Neck Squamous Cell Carcinoma (Overall Survival)',
             'LIHC' = 'Liver Hepatocellular Carcinoma',
             'KIRC' = 'Kidney Renal Clear Cell Carcinoma',
             'KIRP' = 'Kidney Renal Papillary Cell Carcinoma',
             'PRAD' = 'Prostate Adenocarcinoma (Progression-free Survival)')
  gs.best.survtypes <- list('LIHC' = c('OS', 'PFS'),
                            'KIRC' = c('OS', 'PFS'),
                            'KIRP' = c('OS', 'PFS'), 
                            'PRAD' = 'PFS', 'HNSC' = 'OS')
  cs <- 11    # common text size
  sample.size <- 50
  ylbl <- 'Hazard Ratio'
  survtype.labels <- c('Overall Survival', 'Progression-free Survival')
  names(survtype.labels) <- c('OS', 'PFS')
  item.label <- c('G' = 'Genes', 'PW' = 'Gene-centric\nPathways',
                  'TRAD' = 'Classical\nPathways')
  # Now, plotting
  p1.all <- list()
  p2.all <- list()
  p.all <- list()
  count <- 1
  for(g in gs.best) {
    n.st <- length(gs.best.survtypes[[g]])
    ## p1: all signif items, simple violins
    res1 <- subset(results.all, group == g)
    res1 <- subset(res1, surv.type %in% gs.best.survtypes[[g]])
    p1 <- ggplot(res1, aes(x = item.type, y = hr, fill = item.type)) +
      geom_violin() + geom_boxplot(width = 0.1, outlier.shape = NA) +
      geom_hline(yintercept = 1, linetype = 'dashed') +
      scale_fill_brewer(palette = 'Accent') +
      xlab(NULL) + ylab(ylbl) +
      scale_x_discrete(labels = item.label) +
      theme(axis.text.x = element_text(size = cs),
            axis.text.y = element_text(size = cs),
            legend.position = 'none',
            strip.text.y = element_text(size = ifelse(squished, cs-1, cs)))
    if(n.st > 1) {
      p1 <- p1 + facet_grid(surv.type ~ .,
                            labeller = labeller(surv.type = survtype.labels))
    }
    # text annotation
    n.annot1 <- makeNAnnot(res1, g, names(item.label),
                           gs.best.survtypes[[g]], plot = 1)
    p1 <- p1 + geom_text(aes(x = item.type, y = y, label = n),
                         n.annot1, size = 3.5)
    # save the y limits to align p2 later
    p1.ylims <- layer_scales(p1)$y$range$range
    ## p2: signif gpw, box+cross
    res2 <- subset(results.gpw, group == g)
    res2 <- subset(res2, surv.type %in% gs.best.survtypes[[g]])
    p2 <- ggplot(res2, aes(x = item.type, y = hr, fill = item.type)) +
      geom_violin() +
      geom_hline(yintercept = 1, linetype = 'dashed') +
      scale_fill_brewer(palette = 'Accent') +
      xlab(NULL) + ylab(ylbl) +
      scale_x_discrete(labels = item.label[1:2]) +
      theme(axis.text.x = element_text(size = cs),
            axis.text.y = element_text(size = cs),
            legend.position = 'none',
            strip.text.y = element_text(size = ifelse(squished, cs-1, cs)))
    if(n.st > 1) {
      p2 <- p2 + facet_grid(surv.type ~ .,
                            labeller = labeller(surv.type = survtype.labels))
    }
    # align the y axis with p1
    p2 <- p2 + ylim(p1.ylims[1], p1.ylims[2])
    # text annotation
    n.annot2 <- makeNAnnot(res2, g, names(item.label)[1:2],
                           gs.best.survtypes[[g]], plot = 2)
    p2 <- p2 + geom_text(aes(x = item.type, y = y, label = n),
                         n.annot2, size = 3.5)
    # points, jitter
    set.seed(42)
    for(st in gs.best.survtypes[[g]]) {
      res2.st <- subset(res2, surv.type == st)
      res2.st.g <- subset(res2.st, item.type == 'G')
      res2.st.pw <- subset(res2.st, item.type == 'PW')
      if(!all(res2.st.g$item == res2.st.pw$item)) {
        message('Warning: item names do not match in ', g, ' ', st, '!')
        # gotta make sure that we actually connect the paired points,
        # not just random points
      }
      jitt <- runif(nrow(res2.st), -0.05, 0.05)
      jitt.g <- jitt[res2.st$item.type == 'G']
      jitt.pw <- jitt[res2.st$item.type == 'PW']
      N <- nrow(res2.st)/2
      if(N < sample.size) {
        random.cases <- 1:N
      } else {
        random.cases <- sample(1:N, sample.size)
      }
      for(i in random.cases) {
        hrs <- c(res2.st.g$hr[i], res2.st.pw$hr[i])
        df <- data.frame(item.type = c(1+jitt.g[i], 2+jitt.pw[i]),
                         hr = hrs, surv.type = st)
        p2 <- p2 + geom_line(aes(x = as.numeric(item.type), y = hr), df,
                             inherit.aes = FALSE, alpha = 0.3, color = 'black',
                             lwd = 0.4)
        p2 <- p2 + geom_point(aes(x = as.numeric(item.type), y = hr), df,
                              inherit.aes = FALSE, size = 1.5, shape = 3,
                              alpha = 0.5)
      }
    }
    p1.all[[g]] <- p1
    p2.all[[g]] <- p2
    # Arrange p1 and p2 into one plot
    p <- ggarrange(p1, p2 + rremove('ylab'),
                   ncol = 2, nrow = 1, widths = c(3,2),
                   labels = LETTERS[count:(count+1)], legend = 'none')
    # bigfig controls the ABC labels,
    # in case we need to combine the plots into a big figure
    # (like in the report)
    if(bigfig) {
      count <- count + 2
      if(g == 'KIRP') { count <- 1 }
    }
    # Give it a title
    p <- annotate_figure(p, top = text_grob(title[g], color = 'black',
                                            face = 'bold', size = 12))
    # Save
    p.all[[g]] <- p
    if(n.st == 1) {
      w <- 2500
      h <- 1200
      h <- ifelse(squished, 900, 1200)
    } else if(n.st == 2) {
      w <- 2500
      h <- 2000
      h <- ifelse(squished, 1500, 2000)
    }
    plotfile <- file.path(outdir, paste0('hr_', g, '.png'))
    png(filename = plotfile, width = w, height = h,
        units = 'px', res = 350)
    print(p)
    dev.off()
    message('Plot saved: ', plotfile)
  }
}

plotHRs(file.path('report', 'HRboxes'), squished = F)

