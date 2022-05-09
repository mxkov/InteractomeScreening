# BCR biotab exploration:
#
# ood8, follow_up_v4.4_nte_paad:
# new_tumor_event_dx_days_to,
# new_tumor_event_surgery_days_to
#
# 2cc6, nte_paad:
# new_tumor_event_dx_days_to,
# new_tumor_event_surgery_days_to
#
# f4af, follow_up_v4.4_paad:
# nothing

getArgs <- function() {
  biotab.files <- list('TCGA-PAAD' = list())
  biotab.files[['TCGA-PAAD']][[1]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v4.4_nte_paad.txt')
  biotab.files[['TCGA-PAAD']][[2]] <- paste0('nationwidechildrens.org_',
                                             'clinical_nte_paad.txt')
  biotab.fields <- list('TCGA-PAAD' = list())
  fld1 <- 'new_tumor_event_dx_days_to'
  fld2 <- 'new_tumor_event_surgery_days_to'
  for(p in names(biotab.files)) {
    for(i in 1:length(biotab.files[[p]])) {
      biotab.fields[[p]][[i]] <- c(fld1, fld2)
    }
  }
  doubt.patients <- list('PAAD' = character())
  mht <- 'Infiltrating duct carcinoma NOS'

  argums <- list()
  argums$resdir <- 'locos/pancreas/results_TCGA'
  argums$bfil <- biotab.files
  argums$bfld <- biotab.fields
  argums$dp <- doubt.patients
  argums$mht <- mht
  return(argums)
}
