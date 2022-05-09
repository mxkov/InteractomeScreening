# BCR biotab exploration:
#
# 5be2, follow_up_v4.0_nte_cesc:
# new_tumor_event_dx_days_to;
# new_tumor_event_surgery_days_to has no new data and all values are >=
#
# 143b, follow_up_v4.0_cesc:
# nothing relevant
#
# 7104, follow_up_v2.0_cesc:
# new_tumor_event_dx_days_to
#
# ac90, nte_cesc:
# new_tumor_event_dx_days_to;
# new_tumor_event_surgery_days_to has no new data and all values are >=

getArgs <- function() {
  biotab.files <- list('TCGA-CESC' = list())
  biotab.files[['TCGA-CESC']][[1]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v4.0_nte_cesc.txt')
  biotab.files[['TCGA-CESC']][[2]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v2.0_cesc.txt')
  biotab.files[['TCGA-CESC']][[3]] <- paste0('nationwidechildrens.org_',
                                             'clinical_nte_cesc.txt')
  biotab.fields <- list('TCGA-CESC' = list())
  fld <- 'new_tumor_event_dx_days_to'
  for(p in names(biotab.files)) {
    for(i in 1:length(biotab.files[[p]])) {
      biotab.fields[[p]][[i]] <- fld
    }
  }
  doubt.patients <- list('CESC' = character())

  mht <- c('Squamous cell carcinoma keratinizing NOS',
           'Squamous cell carcinoma large cell nonkeratinizing NOS',
           'Squamous cell carcinoma NOS',
           'Papillary squamous cell carcinoma',
           'Basaloid squamous cell carcinoma')
  
  argums <- list()
  argums$resdir <- 'locos/cervixuteri/results_TCGA'
  argums$bfil <- biotab.files
  argums$bfld <- biotab.fields
  argums$dp <- doubt.patients
  argums$mht <- mht
  return(argums)
}
