# BCR biotab exploration:
#
# a762, follow_up_v4.0_nte_lihc:
# new_tumor_event_dx_days_to;
# new_tumor_event_surgery_days_to has no new data and all values are >=,
# additionally, there are TWO cols with that name but they're complementary
# (i.e. the values that are NA in one are non-NA in the other, and vice versa)
#
# be89, follow_up_v4.0_lihc:
# nothing relevant
#
# ef5b, nte_lihc:
# new_tumor_event_dx_days_to;
# new_tumor_event_surgery_days_to: same as a762

getArgs <- function() {
  biotab.files <- list('TCGA-LIHC' = list())
  biotab.files[['TCGA-LIHC']][[1]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v4.0_nte_lihc.txt')
  biotab.files[['TCGA-LIHC']][[2]] <- paste0('nationwidechildrens.org_',
                                             'clinical_nte_lihc.txt')
  biotab.fields <- list('TCGA-LIHC' = list())
  fld <- 'new_tumor_event_dx_days_to'
  for(p in names(biotab.files)) {
    for(i in 1:length(biotab.files[[p]])) {
      biotab.fields[[p]][[i]] <- fld
    }
  }
  doubt.patients <- list('LIHC' = character())
  mht <- 'Hepatocellular carcinoma NOS'

  argums <- list()
  argums$resdir <- 'locos/liver/results_TCGA'
  argums$bfil <- biotab.files
  argums$bfld <- biotab.fields
  argums$dp <- doubt.patients
  argums$mht <- mht
  return(argums)
}
