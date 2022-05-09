# BCR biotab exploration:
#
# 99c2, nte_blca:
# new_tumor_event_dx_days_to;
# new_tumor_event_surgery_days_to has no new data and all values are >=
#
# c5fe, follow_up_v2.0_blca:
# new_tumor_event_dx_days_to
#
# d54d, follow_up_v4.0_nte_blca:
# new_tumor_event_dx_days_to;
# new_tumor_event_surgery_days_to has no new data and all values are >=
#
# f913, follow_up_v4.0_blca:
# nothing relevant

getArgs <- function() {
  biotab.files <- list('TCGA-BLCA' = list())
  biotab.files[['TCGA-BLCA']][[1]] <- paste0('nationwidechildrens.org_',
                                             'clinical_nte_blca.txt')
  biotab.files[['TCGA-BLCA']][[2]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v2.0_blca.txt')
  biotab.files[['TCGA-BLCA']][[3]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v4.0_nte_blca.txt')
  biotab.fields <- list('TCGA-BLCA' = list())
  fld <- 'new_tumor_event_dx_days_to'
  for(p in names(biotab.files)) {
    for(i in 1:length(biotab.files[[p]])) {
      biotab.fields[[p]][[i]] <- fld
    }
  }
  doubt.patients <- list('BLCA' = character())
  mht <- c('Transitional cell carcinoma',
           'Papillary transitional cell carcinoma')

  argums <- list()
  argums$resdir <- 'locos/bladder/results_TCGA'
  argums$bfil <- biotab.files
  argums$bfld <- biotab.fields
  argums$dp <- doubt.patients
  argums$mht <- mht
  return(argums)  
}
