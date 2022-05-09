# BCR biotab exploration:
#
# 0abf, follow_up_v1.0_ov:
# nothing relevant
#
# 5ea3, follow_up_v1.0_nte_ov:
# new_tumor_event_dx_days_to,
# days_to_new_tumor_event_additional_surgery_procedure (has a lot of new data)
#
# d922, nte_ov:
# new_tumor_event_dx_days_to;
# days_to_new_tumor_event_additional_surgery_procedure are all NA

getArgs <- function() {
  biotab.files <- list('TCGA-OV' = list())
  biotab.files[['TCGA-OV']][[1]] <- paste0('nationwidechildrens.org_',
                                           'clinical_follow_up_v1.0_nte_ov.txt')
  biotab.files[['TCGA-OV']][[2]] <- paste0('nationwidechildrens.org_',
                                           'clinical_nte_ov.txt')
  biotab.fields <- list('TCGA-OV' = list())
  fld1 <- 'new_tumor_event_dx_days_to'
  fld2 <- 'days_to_new_tumor_event_additional_surgery_procedure'
  biotab.fields[['TCGA-OV']][[1]] <- c(fld1, fld2)
  biotab.fields[['TCGA-OV']][[2]] <- fld1
  doubt.patients <- list('OV' = character())
  mht <- 'Serous cystadenocarcinoma NOS'

  argums <- list()
  argums$resdir <- 'locos/ovary/results_TCGA'
  argums$bfil <- biotab.files
  argums$bfld <- biotab.fields
  argums$dp <- doubt.patients
  argums$mht <- mht
  return(argums)
}
