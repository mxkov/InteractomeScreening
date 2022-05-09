source('scripts_common/locos-info.R')
source('scripts_common/master.R')

#locos <- names(gs)
source('scripts_common/input-locos.R')

for(l in locos) {
  arg.script <- file.path('locos', l, 'scripts', 'get-args.R')
  if(!file.exists(arg.script)) {
    message('Warning: file ', arg.script, ' does not exist, skipping ', l)
    next
  }
  source(arg.script)
  argums <- getArgs()
  resdir <- file.path('locos', l, 'results_TCGA_test')
  doLoco(paste0('TCGA-',gs[[l]]), 'TCGA',
         'data', l, resdir,
         argums$bfil, argums$bfld, argums$dp,
         do_prep = TRUE, write_expr = TRUE, do_pal = TRUE,
         do_prep_trad = TRUE,
         use_normal = FALSE, main_histypes = argums$mht,
         run_surv = FALSE,
         run_surv_trad = FALSE,
         run_surv_genes = FALSE,
         test = FALSE)
  gc()
}