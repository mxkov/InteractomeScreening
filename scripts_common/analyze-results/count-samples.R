## Just a quick tool for counting samples;
## mostly those used for survival analysis, with all the filters applied,
## but also those initially present in the project if tables.full is TRUE

locos.all <- list.dirs('locos', full.names = F, recursive = F)
source('scripts_common/locos-info.R')

for(loco in names(gs)) {
  if(loco == 'gliomas') { next }
  clin <- read.csv(file.path('locos', loco, 'results_TCGA',
                             'clinical_TCGA.csv'))
  clin <- clin[!is.na(clin$days_os) & !is.na(clin$days_pfs),]
  clin <- clin[!is.na(clin$event_os) & !is.na(clin$event_pfs),]
  clin <- clin[clin$days_os >= 0 & clin$days_pfs >= 0,]
  message('')
  message(loco, ' ', nrow(clin))
  if(loco %in% c('kidney', 'lung')) {
    print(table(clin$project))
  } else {
    print(table(clin$histological_type))
  }
  message('')
}

tables.full <- FALSE
if(tables.full) {
  source('scripts_common/locos-info.R')
  for(g in unlist(gs)) {
    message('')
    message(g)
    cf <- file.path('data', paste0('TCGA-', g), 'Clinical', 'clinical.tsv')
    clin <- fread(cf)
    clin <- unique(clin[,c('case_submitter_id', 'primary_diagnosis',
                           'vital_status')])
    print(table(clin$primary_diagnosis))
  }
}
