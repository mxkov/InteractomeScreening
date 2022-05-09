# BCR biotab exploration:
#
# 0282, nte_skcm:
# new_tumor_event_dx_days_to, new_tumor_event_surgery_days_to
#
# c23a, follow_up_v2.0_skcm:
# new_tumor_event_dx_days_to, new_tumor_event_surgery_days_to

getArgs <- function() {
  biotab.files <- list('TCGA-SKCM' = list())
  biotab.files[['TCGA-SKCM']][[1]] <- paste0('nationwidechildrens.org_',
                                             'clinical_nte_skcm.txt')
  biotab.files[['TCGA-SKCM']][[2]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v2.0_skcm.txt')
  biotab.fields <- list('TCGA-SKCM' = list())
  fld1 <- 'new_tumor_event_dx_days_to'
  fld2 <- 'new_tumor_event_surgery_days_to'
  for(p in names(biotab.files)) {
    for(i in 1:length(biotab.files[[p]])) {
      biotab.fields[[p]][[i]] <- c(fld1, fld2)
    }
  }
  doubt.patients <- list('SKCM' = character())

  argums <- list()
  argums$resdir <- 'locos/skin/results_TCGA'
  argums$bfil <- biotab.files
  argums$bfld <- biotab.fields
  argums$dp <- doubt.patients
  argums$mht <- NULL
  return(argums)
}
