#c <- testClin('TCGA-ESCA', 'locos/esophagus')
# no histo type data

# BCR biotab exploration:
#
# 16e4, follow_up_v4.0_esca:
# nothing relevant
#
# 7837, follow_up_v4.0_nte_esca:
# new_tumor_event_dx_days_to;
# new_tumor_event_surgery_days_to has no new data and all values are >=
#
# b8ab, nte_esca:
# new_tumor_event_dx_days_to;
# new_tumor_event_surgery_days_to has no new data and all values are >=

getArgs <- function() {
  biotab.files <- list('TCGA-ESCA' = list())
  biotab.files[['TCGA-ESCA']][[1]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v4.0_nte_esca.txt')
  biotab.files[['TCGA-ESCA']][[2]] <- paste0('nationwidechildrens.org_',
                                             'clinical_nte_esca.txt')
  biotab.fields <- list('TCGA-ESCA' = list())
  fld <- 'new_tumor_event_dx_days_to'
  for(p in names(biotab.files)) {
    for(i in 1:length(biotab.files[[p]])) {
      biotab.fields[[p]][[i]] <- fld
    }
  }
  doubt.patients <- list('ESCA' = character())
  mht <- NULL
  
  argums <- list()
  argums$resdir <- 'locos/esophagus/results_TCGA'
  argums$bfil <- biotab.files
  argums$bfld <- biotab.fields
  argums$dp <- doubt.patients
  argums$mht <- mht
  return(argums)
}
