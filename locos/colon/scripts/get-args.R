# BCR biotab exploration:
#
# 807f, nte_coad: new_tumor_event_dx_days_to
# (also days_to_new_tumor_event_additional_surgery_procedure but it's all >=)
#
# c535, follow_up_v1.0_nte_coad: new_tumor_event_dx_days_to, AND
# days_to_new_tumor_event_additional_surgery_procedure!
# some barcodes are missing nte days but have surgery days.
#
# ce1b, follow_up_v1.0_coad: nothing relevant
#

getArgs <- function() {
  biotab.files <- list('TCGA-COAD' = list())
  biotab.files[['TCGA-COAD']][[1]] <- paste0('nationwidechildrens.org_',
                                             'clinical_nte_coad.txt')
  biotab.files[['TCGA-COAD']][[2]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v1.0_nte_coad.txt')
  biotab.fields <- list('TCGA-COAD' = list())
  fld1 <- 'new_tumor_event_dx_days_to'
  fld2 <- 'days_to_new_tumor_event_additional_surgery_procedure'
  biotab.fields[['TCGA-COAD']][[1]] <- fld1
  biotab.fields[['TCGA-COAD']][[2]] <- c(fld1, fld2)
  doubt.patients <- list('COAD' = character())
  mht <- c('Adenocarcinoma NOS', 'Mucinous adenocarcinoma',
           'Adenocarcinoma with mixed subtypes',
           'Adenocarcinoma with neuroendocrine differentiation',
           'Papillary adenocarcinoma NOS')
  
  argums <- list()
  argums$resdir <- 'locos/colon/results_TCGA'
  argums$bfil <- biotab.files
  argums$bfld <- biotab.fields
  argums$dp <- doubt.patients
  argums$mht <- mht
  return(argums)
}
