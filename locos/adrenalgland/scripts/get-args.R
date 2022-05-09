# BCR biotab exploration:
#
# 5a0d, nte_pcpg:
# days_to_new_tumor_event_after_initial_treatment
#
# 970f, follow_up_v4.0_pcpg:
# nothing relevant
#
# f1f9, follow_up_v4.0_nte_pcpg:
# days_to_new_tumor_event_after_initial_treatment

getArgs <- function() {
  biotab.files <- list('TCGA-PCPG' = list())
  biotab.files[['TCGA-PCPG']][[1]] <- paste0('nationwidechildrens.org_',
                                             'clinical_nte_pcpg.txt')
  biotab.files[['TCGA-PCPG']][[2]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v4.0_nte_pcpg.txt')
  biotab.fields <- list('TCGA-PCPG' = list())
  fld <- 'days_to_new_tumor_event_after_initial_treatment'
  for(p in names(biotab.files)) {
    for(i in 1:length(biotab.files[[p]])) {
      biotab.fields[[p]][[i]] <- fld
    }
  }
  doubt.patients <- list('PCPG' = character())
  mht <- c('Pheochromocytoma NOS', 'Pheochromocytoma malignant')
  
  argums <- list()
  argums$resdir <- 'locos/adrenalgland/results_TCGA'
  argums$bfil <- biotab.files
  argums$bfld <- biotab.fields
  argums$dp <- doubt.patients
  argums$mht <- mht
  return(argums)
}
