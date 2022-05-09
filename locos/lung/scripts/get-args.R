getArgs <- function() {
  biotab.files <- list('TCGA-LUAD' = list(), 'TCGA-LUSC' = list())
  biotab.files[['TCGA-LUAD']][[1]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v1.0_luad.txt')
  biotab.files[['TCGA-LUAD']][[2]] <- paste0('nationwidechildrens.org_',
                                             'clinical_nte_luad.txt')
  biotab.files[['TCGA-LUSC']][[1]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v1.0_lusc.txt')
  biotab.files[['TCGA-LUSC']][[2]] <- paste0('nationwidechildrens.org_',
                                             'clinical_nte_lusc.txt')
  biotab.fields <- list('TCGA-LUAD' = list(), 'TCGA-LUSC' = list())
  fld <- 'new_tumor_event_dx_days_to'
  for(p in names(biotab.files)) {
    for(i in 1:length(biotab.files[[p]])) {
      biotab.fields[[p]][[i]] <- fld
    }
  }
  doubt.patients <- list('LUAD' = character(), 'LUSC' = character())
  
  argums <- list()
  argums$resdir <- 'locos/lung/results_TCGA'
  argums$bfil <- biotab.files
  argums$bfld <- biotab.fields
  argums$dp <- doubt.patients
  argums$mht <- NULL
  return(argums)
}
