source('scripts_common/locos-info.R')

# By default, all tumor localizations from locos-info.R are used:
locos <- names(gs)
# But on not-so-powerful machines, processing all at once can take very long.
# So you might want to specify a shorter vector. Like this:
locos <- c('lung', 'pancreas')
