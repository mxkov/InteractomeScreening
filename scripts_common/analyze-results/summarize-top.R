## Extracts a table with top items from the summary file

library('data.table')
library('openxlsx')
# Reading
summ <- fread('locos/summary.csv')
# Transforming
summ <- summ[apply(summ[,4:10], 1,
                   function(x) {return(any(unlist(x) != 0))}),]
cols.top <- grep('\\.top$', colnames(summ), value = T)
cols <- c('group', 'surv.type', cols.top)
summ <- summ[,..cols]
summ[is.na(summ)] <- ''
summ$surv.type <- toupper(summ$surv.type)
for(f in cols.top) {
  summ[[f]] <- gsub('\\|', ', ', summ[[f]])
  summ[[f]] <- gsub('_', ' ', summ[[f]])
}
summ <- summ[order(summ$group, summ$surv.type),]
# Putting into a workbook
wb <- createWorkbook()
addWorksheet(wb, 'Sheet1')
writeData(wb, 'Sheet1', summ, colNames = TRUE)
# Writing
outfile <- file.path('report', 'summary_tops.xlsx')
saveWorkbook(wb, outfile, overwrite = TRUE)
message('File written: ', outfile)
