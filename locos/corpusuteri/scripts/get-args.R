# BCR biotab exploration:
#
# 6b35, follow_up_v4.0_ucec: nothing relevant
# 
# 8cc6, nte_ucec:
# new_tumor_event_dx_days_to;
# also new_tumor_event_surgery_days_to, but it has no new data
#
# 60e7, follow_up_v1.7_ucec:
# new_tumor_event_dx_days_to
#
# 7710, follow_up_v2.0_ucec:
# new_tumor_event_dx_days_to
#
# bdde, follow_up_v4.0_nte_ucec:
# new_tumor_event_dx_days_to;
# also new_tumor_event_surgery_days_to, but it has no new data

getArgs <- function() {
  biotab.files <- list('TCGA-UCEC' = list())
  biotab.files[['TCGA-UCEC']][[1]] <- paste0('nationwidechildrens.org_',
                                             'clinical_nte_ucec.txt')
  biotab.files[['TCGA-UCEC']][[2]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v1.7_ucec.txt')
  biotab.files[['TCGA-UCEC']][[3]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v2.0_ucec.txt')
  biotab.files[['TCGA-UCEC']][[4]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v4.0_nte_ucec.txt')
  biotab.fields <- list('TCGA-UCEC' = list())
  fld <- 'new_tumor_event_dx_days_to'
  for(p in names(biotab.files)) {
    for(i in 1:length(biotab.files[[p]])) {
      biotab.fields[[p]][[i]] <- fld
    }
  }
  doubt.patients <- list('UCEC' = character())
  mht <- 'Endometrioid adenocarcinoma NOS'
  
  argums <- list()
  argums$resdir <- 'locos/corpusuteri/results_TCGA'
  argums$bfil <- biotab.files
  argums$bfld <- biotab.fields
  argums$dp <- doubt.patients
  argums$mht <- mht
  return(argums)
}
