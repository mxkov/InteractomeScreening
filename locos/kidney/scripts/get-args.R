getArgs <- function() {
  biotab.files <- list('TCGA-KIRC' = list(), 'TCGA-KIRP' = list(),
                       'TCGA-KICH' = list())
  biotab.files[['TCGA-KIRC']][[1]] <- paste0('nationwidechildrens.org_',
                                             'clinical_nte_kirc.txt')
  biotab.files[['TCGA-KIRC']][[2]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v1.0_kirc.txt')
  biotab.files[['TCGA-KIRP']][[1]] <- paste0('nationwidechildrens.org_',
                                             'clinical_nte_kirp.txt')
  biotab.files[['TCGA-KIRP']][[2]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v1.0_kirp.txt')
  biotab.files[['TCGA-KICH']][[1]] <- paste0('nationwidechildrens.org_',
                                             'clinical_nte_kich.txt')
  biotab.files[['TCGA-KICH']][[2]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v4.4_',
                                             'nte_kich.txt')
  biotab.fields <- list('TCGA-KIRC' = list(), 'TCGA-KIRP' = list(), 'TCGA-KICH' = list())
  fld <- 'new_tumor_event_dx_days_to'
  for(p in names(biotab.files)) {
    for(i in 1:length(biotab.files[[p]])) {
      biotab.fields[[p]][[i]] <- fld
    }
  }
  doubt.patients <- list('KIRC' = character(), 'KIRP' = character(),
                         'KICH' = character())

  argums <- list()
  argums$resdir <- 'locos/kidney/results_TCGA'
  argums$bfil <- biotab.files
  argums$bfld <- biotab.fields
  argums$dp <- doubt.patients
  argums$mht <- NULL
  return(argums)
}
