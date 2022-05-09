## Just for making quick survival plots for presentations

survdraw <- function(clin, outfile, st = c('os', 'pfs'),
                     by = c('project', 'tumor_type'),
                     x.lim = NULL, days.cutoff = NULL, add.lbl = '') {
  library('ggplot2')
  library('survival')
  library('survminer')
  
  # Names, etc.
  days.col <- paste0('days_', st)
  event.col <- paste0('event_', st)
  survname <- list('os' = 'Overall Survival',
                   'pfs' = 'Progression-free Survival')
  ylbl <- survname[[st]]
  
  ## Data transforms
  # remove NAs
  clin <- clin[!is.na(clin[,days.col]) & !is.na(clin[,event.col]),]
  # remove negative days
  clin <- clin[clin[,days.col] >= 0,]
  # censor patients with long follow-ups, if requested
  if(!is.null(days.cutoff)) {
    filter <- (clin[,days.col] >= days.cutoff)
    clin[filter, event.col] <- 0
    clin[filter, days.col] <- days.cutoff
  }
  
  so <- Surv(clin[,days.col], clin[,event.col])
  sfit <- surv_fit(so ~ clin[,by], data = clin)
  pval <- surv_pvalue(sfit)$pval
  if(is.null(x.lim)) {
    p <- ggsurvplot(fit = sfit, data = clin, pval = TRUE,
                    risk.table = T, surv.median.line = 'hv',
                    ylab = ylbl)
  } else {
    p <- ggsurvplot(fit = sfit, data = clin, pval = TRUE,
                    risk.table = T, surv.median.line = 'hv',
                    ylab = ylbl, xlim = x.lim)
  }
  #print(pairwise_survdiff(so ~ clin[,by], data = clin))
  #print(layer_scales(p$plot)$x$range$range)
  #p$plot <- p$plot + geom_label(label = paste0('log-rank p-val = ',
  #                                             round(pval, 4)),
  #                              x = layer_scales(p$plot)$x$range$range[2]*0.7,
  #                              y = 0.8)
  png(filename = outfile, width = 2000, height = 2000, res = 300)
  print(p)
  dev.off()
}


## PARAMS
surv.types <- c('os', 'pfs')
locos <- c('lung', 'kidney', 'headneck')
#locos <- 'kidney'
split.kidney <- FALSE
gs <- list('lung' = c('LUAD', 'LUSC'), 'kidney' = c('KIRC', 'KIRP'),
           'headneck' = 'HNSC')
xlims <- list('lung' = c(0, 8000), 'kidney' = c(0, 6000),
              'headneck' = c(0, 6000))
projects <- c('TCGA', 'CPTAC')
out.dir <- 'seminar/2022-04-25'
#cut <- 1800
cut <- NULL

## RUN
by.what <- 'project'
if(by.what == 'tumor_type') {
  for(loco in locos) {
    loco.dir <- file.path('locos', loco)
    for(proj in projects) {
      clin.file <- file.path(loco.dir, paste0('results_', proj),
                             paste0('clinical_', proj, '.csv'))
      clin <- read.csv(clin.file)
      if(proj == 'TCGA') {
        clin$tumor_type <- gsub('^TCGA-', '', clin$project)
        clin <- clin[clin$tumor_type %in% gs[[loco]],]
      }
      for(stype in surv.types) {
        outfile <- file.path(out.dir, paste0('surv_', loco, '_',
                                             stype, '_', proj, '.png'))
        survdraw(clin, outfile, stype, by = 'tumor_type')
      }
    }
  }
}
if(by.what == 'project') {
  for(loco in locos) {
    loco.dir <- file.path('locos', loco)
    if(loco == 'kidney') {
      if(split.kidney) {
        lbl <- '_split'
      } else {
        lbl <- '_nosplit'
      }
    } else {
      lbl <- ''
    }
    for(g in gs[[loco]]) {
      first <- TRUE
      for(proj in projects) {
        clin.file <- file.path(loco.dir, paste0('results_', proj),
                               paste0('clinical_', proj, '.csv'))
        clin <- read.csv(clin.file)
        if(proj == 'TCGA') {
          clin$tumor_type <- gsub('^TCGA-', '', clin$project)
        }
        if(loco == 'kidney') {
          if(proj == 'TCGA') {
            if(split.kidney) {
              clin <- clin[clin$tumor_type %in% gs[[loco]],]
            } else {
              clin$project <- 'TCGA'
            }
          }
          g1 <- 'KIR'
        } else {
          g1 <- g
          clin <- clin[clin$tumor_type == g1,]
        }
        clin <- clin[, c('case_submitter_id', 'project', 'tumor_type',
                         paste0('days_', surv.types),
                         paste0('event_', surv.types))]
        if(first) {
          clin0 <- clin
          first <- FALSE
        } else {
          clin0 <- rbind(clin0, clin)
        }
      }
      for(stype in surv.types) {
        outfile <- file.path(out.dir, paste0('surv_', loco, '_', g1, lbl, '_',
                                             stype, '.png'))
        survdraw(clin0, outfile, stype, by = 'project')
      }
      if(g1 == 'KIR') { break }
    }
  }
}
