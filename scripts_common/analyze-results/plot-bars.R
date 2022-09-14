## Makes a bar plot summarizing the results
## The actual function call is at the bottom

source('scripts_common/locos-info.R')

plotBars <- function(plotfile, plot.stable = FALSE,
                     lab.type = c('abbr', 'full')) {
  ## Parameters:
  # plotfile: path to the output file
  # plot.stable: whether to plot the numbers of significant genes
  #              paired with respective gene-centric pathways if the latter
  #              are also significant
  # lab.type: full for full x-axis tick labels (the cancer types),
  #           abbr for abbreviated
  #
  library('data.table')
  library('ggplot2')
  library('viridis')
  surv.types <- c('os', 'pfs')
  locos.todo <- names(gs)    # gs is imported from locos-info.R
  
  # Read the summary to see which locos we need to plot
  summ0 <- fread(file.path('locos', 'summary.csv'))[,2:14]
  summ <- summ0[apply(summ0[,3:9], 1,
                      function(x) {return(sum(unlist(x)) >= 10)}),]
  # Fill in some blanks
  # (i.e. if one survival type got discarded for a tumor type, add it back)
  for(g in unique(summ$group)) {
    g.st <- summ$surv.type[summ$group == g]
    g.st <- g.st[order(g.st)]
    if(!all(g.st == surv.types)) {
      summ <- rbind(summ, summ0[summ0$group==g,])
    }
  }
  summ <- unique(summ[order(summ$group, summ$surv.type),])
  
  # Add more columns:
  # percentages, labels
  round.to <- 1
  measure.cols <- c('genes.n', 'pws.n', 'stable.n', 'trad.pws.n')
  for(f in measure.cols) {
    f.pc <- gsub('\\.n$', '.pc', f)
    if(f == 'stable.n') {
      f.total <- 'common.n.total'
    } else {
      f.total <- paste0(f, '.total')
    }
    f.label <- paste0('label.', gsub('\\.n$', '', f))
    summ[[f.pc]] <- round(100.0*summ[[f]] / summ[[f.total]], round.to)
    pc <- as.character(summ[[f.pc]][summ[[f]] != 0])
    pc[pc == '0'] <- paste0('<', 10^(-round.to))
    summ[[f.label]] <- ''
    summ[[f.label]][summ[[f]] != 0] <- paste0(summ[[f]][summ[[f]] != 0],
                                              ' (', pc, '%)')
  }
  
  # Time to plot...
  # Declaring labels and stuff
  survtype.labs <- c('Overall Survival', 'Progression-free Survival')
  names(survtype.labs) <- c('os', 'pfs')
  legend.labs <- c('genes.n' = paste0('Signif. genes\n(',
                                      max(summ$genes.n.total), ' total)'),
                   'pws.n' = paste0('Signif. gene-centric pathways\n(',
                                    max(summ$pws.n.total), ' total)'),
                   'stable.n' = paste0('Signif. stable genes\n(',
                                       max(summ$common.n.total), ' total)'),
                   'trad.pws.n' = paste0('Signif. classical pathways\n(',
                                         max(summ$trad.pws.n.total), ' total)'))
  if(lab.type == 'abbr') {
    ls <- 11    # label size
    bar.labs <- c('BLCA' = 'Urothelial carc.',
                  'HNSC' = 'Head & Neck\nsq. c. carc.',
                  'KIRC' = 'Clear cell\nrenal cell carc.',
                  'KIRP' = 'Papillary\nrenal cell carc.',
                  'LIHC' = 'Hepatocellular\ncarc.',
                  'LUAD' = 'Lung\nadenocarc.',
                  'LUSC' = 'Lung sq. c. carc.',
                  'PRAD' = 'Prostate\nadenocarc.',
                  'SARC' = 'Sarcomas')
  } else if(lab.type == 'full') {
    ls <- 10
    bar.labs <- c('BLCA' = 'Urothelial\ncarcinoma',
                  'HNSC' = 'Head & Neck\nsquamous cell\ncarcinoma',
                  'KIRC' = 'Clear cell\nrenal cell carcinoma',
                  'KIRP' = 'Papillary renal\ncell carcinoma',
                  'LIHC' = 'Hepatocellular\ncarcinoma',
                  'LUAD' = 'Lung\nadenocarcinoma',
                  'LUSC' = 'Lung squamous\ncell carcinoma',
                  'PRAD' = 'Prostate\nadenocarcinoma',
                  'SARC' = 'Sarcomas')
  }
  cs = 11    # 'cs' means 'common size' (for text, in pt)
  # Need to reshape summ before plotting
  id.cols <- c('group', 'surv.type')
  first <- TRUE
  for(f in measure.cols) {
    f.label <- paste0('label.', gsub('\\.n$', '', f))
    summ1 <- summ[, ..id.cols]
    summ1$N <- summ[[f]]
    summ1$N[summ1$N == 0] <- 1    # for the log scale
    summ1$LAB <- summ[[f.label]]
    summ1$VAR <- f
    if(first) {
      # 'fp' = 'for plotting'
      summ.fp <- summ1
      first <- FALSE
    } else {
      summ.fp <- rbind(summ.fp, summ1)
    }
  }
  if(!plot.stable) {
    summ.fp <- summ.fp[summ.fp$VAR != 'stable.n']
    legend.labs <- legend.labs[names(legend.labs) != 'stable.n']
  }
  summ.fp$VAR <- factor(summ.fp$VAR,
                        levels = names(legend.labs), labels = legend.labs)
  #
  ## PLOTTING, for real this time
  p <- ggplot(data = summ.fp, aes(x = group, y = N, fill = VAR, width = .7)) +
    geom_bar(stat = 'identity', position = position_dodge()) +
    scale_fill_viridis(option = 'magma', discrete = T, begin = 0.1, end = 0.8) +
    geom_text(aes(label = LAB), angle = 90, hjust = -0.05, vjust = 0.45,
              size = 2.7, position = position_dodge(0.7)) +
    facet_grid(surv.type ~ ., labeller = labeller(surv.type = survtype.labs)) +
    labs(fill = NULL) + xlab(NULL) + ylab('Count') +
    scale_x_discrete(labels = bar.labs) +
    scale_y_continuous(trans = 'log10', limits = c(1, 10^5.4)) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = ls),
          axis.text.y = element_text(size = cs),
          axis.title.y = element_text(size = cs),
          axis.line = element_line(color = 'black', lineend = 'butt'),
          axis.ticks = element_line(color = 'black'),
          legend.position = 'top',
          legend.text = element_text(size = 10),
          strip.text.y = element_text(size = cs, color = 'black'),
          strip.background = element_rect(fill = 'gray'))
  
  # Saving the plot
  ggsave(plotfile, plot = p,
         width = 7.5, height = 6.6, units = 'in', dpi = 300)
  message('Plot saved: ', plotfile)
}

plotBars(file.path('report', 'summary_bars.tiff'), lab.type = 'full')
