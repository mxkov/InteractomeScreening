## Makes the tile plot with biomarker percentages
## (yes, 'pc' stands for 'percentage')
## The function call is at the very bottom

source('scripts_common/locos-info.R')

plotpcs <- function(plotfile,
                    size.scale = TRUE, round.to = 1,
                    keep.accur = TRUE, log.scale = TRUE) {
  ## The parameters:
  # plotfile: path to the output file
  # size.scale: max % in each case will have a larger label
  # round.to: number of decimal points for labels, >=0
  # keep.accur: disables '<10^(-round.to)%' labels
  # log.scale: log10 fill color scale
  
  library('data.table')
  library('ggplot2')
  library('ggnewscale')
  library('viridis')
  surv.types <- c('os', 'pfs')
  locos.todo <- names(gs)    # gs is imported from locos-info.R
  
  # Reading the summary to see which locos we need to plot
  summ0 <- fread(file.path('locos', 'summary.csv'))[,2:14]
  summ <- summ0[apply(summ0[,3:9], 1,
                      function(x) {return(any(unlist(x) != 0))}),]
  # Filling in some blanks
  # (i.e. if one survival type got discarded for a tumor type, add it back)
  for(g in unique(summ$group)) {
    g.st <- summ$surv.type[summ$group == g]
    g.st <- g.st[order(g.st)]
    if(!all(g.st == surv.types)) {
      summ <- rbind(summ, summ0[summ0$group==g,])
    }
  }
  summ <- unique(summ[order(summ$group, summ$surv.type),])
  
  # Removing unnecessary columns
  id.cols <- c('group', 'surv.type')
  measure.cols <- c('genes.n', 'pws.n', 'trad.pws.n')
  total.cols <- paste0(measure.cols, '.total')
  cols <- c(id.cols, measure.cols, total.cols)
  summ <- summ[, ..cols]
  
  # Adding percentages
  for(f in measure.cols) {
    f.pc <- gsub('\\.n$', '.pc', f)
    f.total <- paste0(f, '.total')
    summ[[f.pc]] <- 100.0*summ[[f]] / summ[[f.total]]
  }
  
  # Reshaping the data for plotting
  first <- TRUE
  for(f in measure.cols) {
    summ1 <- summ[,..id.cols]
    f.pc <- gsub('\\.n$', '.pc', f)
    summ1$PC <- summ[[f.pc]]
    summ1$VAR <- f
    if(first) {
      summ.fp <- summ1
      first <- FALSE
    } else {
      summ.fp <- rbind(summ.fp, summ1)
    }
  }
  group.levels <- unique(summ.fp$group)
  group.levels <- group.levels[order(group.levels, decreasing = TRUE)]
  summ.fp$group <- factor(summ.fp$group, levels = group.levels)
  
  # Label size? (larger for top performers)
  summ.fp$LAB <- ''
  if(size.scale) {
    summ.fp$LAB.sz <- 0
    for(g in unique(summ.fp$group)) {
      for(st in surv.types) {
        filter <- (summ.fp$group == g & summ.fp$surv.type == st)
        top.pc <- max(summ.fp$PC[filter])
        summ.fp$LAB.sz[filter] <- as.integer(summ.fp$PC[filter] == top.pc)
      }
    }
    summ.fp$LAB.sz <- as.factor(summ.fp$LAB.sz)
  }
  # Adding the labels themselves
  # (preserving accuracy wherever needed)
  summ.fp$LAB <- paste0(round(summ.fp$PC, round.to), '%')
  summ.fp$LAB[summ.fp$PC == 0.0] <- ''
  filter1 <- (summ.fp$LAB == '0%')
  if(keep.accur) {
    filter2 <- rep(TRUE, nrow(summ.fp))
  } else {
    filter2 <- (summ.fp$LAB.sz == 1)
  }
  r <- 0
  while(any(filter1 & filter2)) {
    summ.fp$LAB[filter1 & filter2] <- paste0(round(summ.fp$PC[filter1 & filter2],
                                                   round.to+r), '%')
    r <- r + 1
    filter1 <- (summ.fp$LAB == '0%')
  }
  summ.fp$LAB[filter1] <- paste0('<', 10^(-round.to), '%')
  # Additional special case:
  filter <- (summ.fp$group == 'SARC' & summ.fp$surv.type == 'os' &
               summ.fp$PC != 0.0)
  summ.fp$LAB[filter] <- paste0(round(summ.fp$PC[filter], round.to+r-1), '%')
  # Color for the labels: should be black on light tiles and white on dark tiles
  summ.fp$LAB.clr <- as.factor(as.integer(summ.fp$PC <= 5.0))
  
  # Log scale?
  if(log.scale) {
    summ.fp$PC <- summ.fp$PC + 1
  }
  
  # Plot the TILES
  # Then add the lines and labels
  item.label <- c('genes.n' = 'Genes', 'pws.n' = 'Gene-centric\nPathways',
                  'trad.pws.n' = 'Classical\nPathways')
  survtype.labs <- c('Overall Survival', 'Progression-free Survival')
  names(survtype.labs) <- c('os', 'pfs')
  cs = 11    # 'cs' means 'common size' (for text, in pt)
  
  p <- ggplot(summ.fp, aes(x = VAR, y = group)) +
    geom_tile(aes(fill = PC)) +
    facet_grid(~surv.type, labeller = labeller(surv.type = survtype.labs)) +
    labs(x = NULL, y = NULL) +
    scale_x_discrete(labels = item.label, expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    theme(axis.text.x = element_text(size = cs),
          axis.text.y = element_text(size = cs),
          legend.key.height = unit(1.0, 'cm'),
          legend.key.width = unit(0.8, 'cm'),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = cs),
          strip.text.x = element_text(size = cs))
  if(size.scale) {
    p <- p + geom_text(aes(x = VAR, y = group, label = LAB,
                           color = LAB.clr, size = LAB.sz), show.legend = F) +
      scale_color_manual(aesthetics = 'color', values = c('black', 'white')) +
      scale_size_discrete(range = c(3.5, 5), guide = 'none')
  } else {
    p <- p + geom_text(aes(x = VAR, y = group, label = LAB,
                           color = LAB.clr), show.legend = F) +
      scale_color_manual(aesthetics = 'color', values = c('black', 'white'))
  }
  if(log.scale) {
    log.breaks <- c(1, 2, 3, 6, (1:3)*10+1)
    p <- p + scale_fill_viridis(trans = 'log10', breaks = log.breaks,
                                labels = as.character(log.breaks-1),
                                name = '% of Markers')
  } else {
    p <- p + scale_fill_viridis(name = '% of Markers')
  }
  
  # Saving the plot
  ggsave(plotfile, plot = p,
         width = 7.5, height = 6.0, units = 'in', dpi = 300)
  message('Plot saved: ', plotfile)
}


plotpcs(file.path('report', 'summary_pc_tiles.tiff'))