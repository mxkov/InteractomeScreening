# BCR biotab exploration:
#
# fwup 1.5: has nte days
# fwup 4.0: no nte days
# fwup 2.1: has nte days
# But the actual data we need are in the two files with 'nte' in the name;
# they have all data from the fwup files (with consistency) and more.
# (and the overlap between the nte files and the fwup files is tiny anyway,
# it's like 1-3 samples)
# And yes, we need both nte files.
#
# About the fields:
# follow_up_v4.0_nte_brca (a9ba) has new_tumor_event_surgery_days_to,
# but it doesn't have any new data compared to new_tumor_event_dx_days_to,
# and all values are >= than new_tumor_event_dx_days_to.
# nte_brca (a88c) has new_tumor_event_surgery_days_to, with one new value,
# but the same patient also has a smaller new_tumor_event_dx_days_to value.

getArgs <- function() {
  biotab.files <- list('TCGA-BRCA' = list())
  biotab.files[['TCGA-BRCA']][[1]] <- paste0('nationwidechildrens.org_',
                                             'clinical_nte_brca.txt')
  biotab.files[['TCGA-BRCA']][[2]] <- paste0('nationwidechildrens.org_',
                                             'clinical_follow_up_v4.0_nte_brca.txt')
  biotab.fields <- list('TCGA-BRCA' = list())
  fld <- 'new_tumor_event_dx_days_to'
  for(i in 1:length(biotab.files[['TCGA-BRCA']])) {
    biotab.fields[['TCGA-BRCA']][[i]] <- fld
  }
  doubt.patients <- list('BRCA' = character())
  mht <- 'Infiltrating duct carcinoma NOS'
  
  argums <- list()
  argums$resdir <- 'locos/breast/results_TCGA'
  argums$bfil <- biotab.files
  argums$bfld <- biotab.fields
  argums$dp <- doubt.patients
  argums$mht <- mht
  return(argums) 
}
