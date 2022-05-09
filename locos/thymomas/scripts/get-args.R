# BCR biotab exploration:
#
# 6fe5, follow_up_v4.0_thym:
# nothing relevant
#
# 296c, nte_thym:
# new_tumor_event_dx_days_to
#
# d77b, follow_up_v4.0_nte_thym:
# new_tumor_event_dx_days_to

getArgs <- function() {
  biotab.files <- list('TCGA-THYM' = list())
  biotab.files[['TCGA-THYM']][[1]] <- paste0('nationwidechildrens.org_',
                                             'clinical_nte_thym.txt')
  biotab.files[['TCGA-THYM']][[2]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v4.0_nte_thym.txt')
  biotab.fields <- list('TCGA-THYM' = list())
  fld <- 'new_tumor_event_dx_days_to'
  for(p in names(biotab.files)) {
    for(i in 1:length(biotab.files[[p]])) {
      biotab.fields[[p]][[i]] <- fld
    }
  }
  doubt.patients <- list('THYM' = character())
  mht <- c('Thymoma type A malignant', 'Thymoma type A NOS',
           'Thymoma type AB malignant', 'Thymoma type AB NOS',
           'Thymoma type B1 malignant', 'Thymoma type B1 NOS',
           'Thymoma type B2 malignant', 'Thymoma type B2 NOS',
           'Thymoma type B3 malignant')

  argums <- list()
  argums$resdir <- 'locos/thymomas/results_TCGA'
  argums$bfil <- biotab.files
  argums$bfld <- biotab.fields
  argums$dp <- doubt.patients
  argums$mht <- mht
  return(argums)
}
