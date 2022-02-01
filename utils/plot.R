
library(tidyverse)
library(cowplot)
library(extrafont)

sweeps_save <- function(file, fig, width=4.5, asp=4/3) {
    save_plot(file, fig, base_height=NULL, base_width=width, base_asp=asp)
    embed_fonts(file)
}

sweeps_theme <- theme_minimal() +
    theme(
        text=element_text(size=12, family="Arial Narrow"),
	axis.title.x=element_text(hjust=1),
        panel.grid.minor=element_blank(),
	strip.text=element_text(hjust=0, face='italic')
    )

sweeps_colour <- scale_colour_brewer(palette='Dark2')

target_factor <- function(target) {
    result <- str_replace_all(target, c(
        'log-sel-strength'='Sel. strength',
        'sweep-mode'='Sweep mode',
        'hard-vs-soft'='Hard vs. Soft',
        'rnm-vs-sgv'='RNM vs. SGV'
    ))
    result <- factor(result, levels=c('Sel. strength', 'Sweep mode', 'Hard vs. Soft', 'RNM vs. SGV'))
    return(result)
}

sweepmode_factor <- function(mode) {
    result <- str_replace_all(mode, c(
        'hard'='Hard sweeps',
        'rnm \\(true\\)'='RNM sweeps',
        'sgv \\(true\\)'='SGV sweeps'
    ))
    result <- factor(result, levels=c('Hard sweeps', 'RNM sweeps', 'SGV sweeps'))
    return(result)
}


sweepmode_factor_short <- function(mode) {
    result <- str_replace_all(mode, c(
        'hard'='Hard',
        'rnm \\(true\\)'='RNM',
        'sgv \\(true\\)'='SGV'
    ))
    result <- factor(result, levels=c('Hard', 'RNM', 'SGV'))
    return(result)
}