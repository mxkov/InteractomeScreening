# BCR biotab exploration:
#
# 02d7, follow_up_v1.0_stad:
# new_tumor_event_dx_days_to
#
# ae85, nte_stad:
# new_tumor_event_dx_days_to

getArgs <- function() {
  biotab.files <- list('TCGA-STAD' = list())
  biotab.files[['TCGA-STAD']][[1]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v1.0_stad.txt')
  biotab.files[['TCGA-STAD']][[2]] <- paste0('nationwidechildrens.org_',
                                             'clinical_nte_stad.txt')
  biotab.fields <- list('TCGA-STAD' = list())
  fld <- 'new_tumor_event_dx_days_to'
  for(p in names(biotab.files)) {
    for(i in 1:length(biotab.files[[p]])) {
      biotab.fields[[p]][[i]] <- fld
    }
  }
  doubt.patients <- list('STAD' = character())

  argums <- list()
  argums$resdir <- 'locos/stomach/results_TCGA'
  argums$bfil <- biotab.files
  argums$bfld <- biotab.fields
  argums$dp <- doubt.patients
  argums$mht <- NULL
  return(argums)
}
