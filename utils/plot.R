
library(tidyverse)
library(cowplot)

sweeps_save <- function(file, fig, width=4.5, asp=4/3) {
    save_plot(file, fig, base_height=NULL, base_width=width, base_asp=asp)
}

sweeps_theme <- theme_minimal() +
    theme(
        text=element_text(size=11),
        panel.grid.minor=element_blank()
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