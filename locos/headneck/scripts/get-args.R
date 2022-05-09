# BCR biotab exploration:
#
# 6864, nte_hnsc:
# new_tumor_event_dx_days_to;
# days_to_new_tumor_event_additional_surgery_procedure has no new info
#
# 1799, follow_up_v4.8_hnsc:
# nothing relevant
#
# c44d, follow_up_v1.0_hnsc:
# new_tumor_event_dx_days_to
#
# e025, follow_up_v4.8_nte_hnsc:
# new_tumor_event_dx_days_to;
# new_tumor_event_surgery_days_to has no new data and all days are >=

getArgs <- function() {
  biotab.files <- list('TCGA-HNSC' = list())
  biotab.files[['TCGA-HNSC']][[1]] <- paste0('nationwidechildrens.org_',
                                             'clinical_nte_hnsc.txt')
  biotab.files[['TCGA-HNSC']][[2]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v1.0_hnsc.txt')
  biotab.files[['TCGA-HNSC']][[3]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v4.8_nte_hnsc.txt')
  biotab.fields <- list('TCGA-HNSC' = list())
  fld <- 'new_tumor_event_dx_days_to'
  for(p in names(biotab.files)) {
    for(i in 1:length(biotab.files[[p]])) {
      biotab.fields[[p]][[i]] <- fld
    }
  }
  doubt.patients <- list('HNSC' = character())
  mht <- NULL

  argums <- list()
  argums$resdir <- 'locos/headneck/results_TCGA_normal'
  argums$bfil <- biotab.files
  argums$bfld <- biotab.fields
  argums$dp <- doubt.patients
  argums$mht <- mht
  return(argums)
}
