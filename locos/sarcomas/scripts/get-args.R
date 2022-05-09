# BCR biotab exploration:
#
# 53f2, nte_sarc:
# new_tumor_event_dx_days_to; 
# new_tumor_event_surgery_days_to has no new data and one negative date
# (TCGA-3B-A9HO)
#
# 168d, follow_up_v4.0_nte_sarc:
# new_tumor_event_dx_days_to, new_tumor_event_surgery_days_to;
# the latter has no new data but one smaller date (for TCGA-3B-A9HR)
#
# b221, follow_up_v4.0_sarc:
# nothing relevant

getArgs <- function() {
  biotab.files <- list('TCGA-SARC' = list())
  biotab.files[['TCGA-SARC']][[1]] <- paste0('nationwidechildrens.org_',
                                             'clinical_nte_sarc.txt')
  biotab.files[['TCGA-SARC']][[2]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v4.0_nte_sarc.txt')
  biotab.fields <- list('TCGA-SARC' = list())
  fld1 <- 'new_tumor_event_dx_days_to'
  fld2 <- 'new_tumor_event_surgery_days_to'
  biotab.fields[['TCGA-SARC']][[1]] <- fld1
  biotab.fields[['TCGA-SARC']][[2]] <- c(fld1, fld2)
  doubt.patients <- list('SARC' = character())

  argums <- list()
  argums$resdir <- 'locos/sarcomas/results_TCGA'
  argums$bfil <- biotab.files
  argums$bfld <- biotab.fields
  argums$dp <- doubt.patients
  argums$mht <- NULL
  return(argums)
}
