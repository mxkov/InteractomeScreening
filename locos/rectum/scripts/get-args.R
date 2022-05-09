# BCR biotab exploration:
#
# 2e7b, nte_read:
# contains a single case entry; new_tumor_event_dx_days_to;
# days_to_new_tumor_event_additional_surgery_procedure is >=
#
# afbb, follow_up_v1.0_read:
# nothing relevant
#
# fe40, follow_up_v1.0_nte_read:
# new_tumor_event_dx_days_to;
# days_to_new_tumor_event_additional_surgery_procedure has no new data and all >=

getArgs <- function() {
  biotab.files <- list('TCGA-READ' = list())
  biotab.files[['TCGA-READ']][[1]] <- paste0('nationwidechildrens.org_',
                                             'clinical_nte_read.txt')
  biotab.files[['TCGA-READ']][[2]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v1.0_nte_read.txt')
  biotab.fields <- list('TCGA-READ' = list())
  fld <- 'new_tumor_event_dx_days_to'
  for(p in names(biotab.files)) {
    for(i in 1:length(biotab.files[[p]])) {
      biotab.fields[[p]][[i]] <- fld
    }
  }
  doubt.patients <- list('READ' = character())

  argums <- list()
  argums$resdir <- 'locos/rectum/results_TCGA'
  argums$bfil <- biotab.files
  argums$bfld <- biotab.fields
  argums$dp <- doubt.patients
  argums$mht <- NULL
  return(argums)
}
