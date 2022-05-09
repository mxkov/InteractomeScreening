source('scripts_common/master.R')
source('scripts_common/tools/write-session-info.R')

source('locos/thymomas/scripts/get-args.R')
argums <- getArgs()

doLoco(c('TCGA-THYM'), 'TCGA',
       'data', 'thymomas', argums$resdir,
       argums$bfil, argums$bfld, argums$dp,
       do_prep = TRUE, write_expr = TRUE, do_pal = TRUE,
       do_prep_trad = TRUE,
       use_normal = FALSE, main_histypes = argums$mht,
       run_diff = TRUE, run_surv = TRUE,
       run_diff_trad = TRUE, run_surv_trad = TRUE,
       run_diff_genes = TRUE, run_surv_genes = TRUE,
       test = FALSE)
writeSessionInfo(file.path(resdir, 'sessionInfo_main.txt'))
