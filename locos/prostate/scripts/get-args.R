# BCR biotab exploration:
#
# 6b4a, nte_prad:
# new_tumor_event_dx_days_to
#
# e038, follow_up_v1.0_prad:
# new_tumor_event_dx_days_to

getArgs <- function() {
  biotab.files <- list('TCGA-PRAD' = list())
  biotab.files[['TCGA-PRAD']][[1]] <- paste0('nationwidechildrens.org_',
                                             'clinical_nte_prad.txt')
  biotab.files[['TCGA-PRAD']][[2]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v1.0_prad.txt')
  biotab.fields <- list('TCGA-PRAD' = list())
  fld <- 'new_tumor_event_dx_days_to'
  for(p in names(biotab.files)) {
    for(i in 1:length(biotab.files[[p]])) {
      biotab.fields[[p]][[i]] <- fld
    }
  }
  doubt.patients <- list('PRAD' = character())
  mht <- 'Adenocarcinoma NOS'

  argums <- list()
  argums$resdir <- 'locos/prostate/results_TCGA'
  argums$bfil <- biotab.files
  argums$bfld <- biotab.fields
  argums$dp <- doubt.patients
  argums$mht <- mht
  return(argums)
}
