# BCR biotab exploration (a whoooole lot of them):
#
# 0e25: follow_up_v4.0_nte_thca: new_tumor_event_dx_days_to
# (also new_tumor_event_dx_days_to exists but all values are equal or bigger)
# 
# 1aec: follow_up_v2.0_thca: nothing relevant
#
# 9c88: follow_up_v2.0_nte_thca: new_tumor_event_dx_days_to
# 
# 63b9: nte_thca: new_tumor_event_dx_days_to
# (also new_tumor_event_dx_days_to exists but all values are equal or bigger)
#
# 451b: follow_up_v4.0_thca: nothing relevant

getArgs <- function() {
  biotab.files <- list('TCGA-THCA' = list())
  biotab.files[['TCGA-THCA']][[1]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v4.0_nte_thca.txt')
  biotab.files[['TCGA-THCA']][[2]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v2.0_nte_thca.txt')
  biotab.files[['TCGA-THCA']][[3]] <- paste0('nationwidechildrens.org_',
                                             'clinical_nte_thca.txt')
  biotab.fields <- list('TCGA-THCA' = list())
  fld <- 'new_tumor_event_dx_days_to'
  for(p in names(biotab.files)) {
    for(i in 1:length(biotab.files[[p]])) {
      biotab.fields[[p]][[i]] <- fld
    }
  }
  doubt.patients <- list('THCA' = character())
  mht <- c('Papillary adenocarcinoma NOS',
           'Papillary carcinoma follicular variant',
           'Papillary carcinoma columnar cell')

  argums <- list()
  argums$resdir <- 'locos/thyroid/results_TCGA'
  argums$bfil <- biotab.files
  argums$bfld <- biotab.fields
  argums$dp <- doubt.patients
  argums$mht <- mht
  return(argums)
}
