## Contains the mapping between localization names and cancer types

# 'gs' means 'groups' I guess
gs <- list('adrenalgland' = 'PCPG', 'bladder' = 'BLCA', 'breast' = 'BRCA',
           'cervixuteri' = 'CESC', 'corpusuteri' = 'UCEC',
           'colon' = 'COAD', 'esophagus' = 'ESCA',
           'headneck' = 'HNSC', 'kidney' = c('KIRC', 'KIRP', 'KICH'),
           'liver' = 'LIHC', 'lung' = c('LUAD', 'LUSC'),
           'ovary' = 'OV', 'pancreas' = 'PAAD',
           'prostate' = 'PRAD', 'rectum' = 'READ', 'sarcomas' = 'SARC',
           'skin' = 'SKCM', 'stomach' = 'STAD',
           'thymomas' = 'THYM', 'thyroid' = 'THCA')

# discarded:
# 'testis' = 'TGCT'
