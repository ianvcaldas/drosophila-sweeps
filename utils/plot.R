
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
        panel.grid.minor=element_blank(),
        strip.text=element_text(hjust=0, face='italic'),
        panel.spacing.x = unit(1, "lines"),
        plot.title=element_text(size=12, family="Arial Narrow", hjust=0, face='italic')
    )

sweeps_colour <- scale_colour_brewer(palette='Dark2')
sweeps_fill <- scale_fill_brewer(palette='Dark2')

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
        'sgv \\(true\\)'='SGV sweeps',
        'soft'='Soft sweeps'
    ))
    result <- factor(result, levels=c('Hard sweeps', 'RNM sweeps', 'SGV sweeps', 'Soft sweeps'))
    return(result)
}


sweepmode_factor_short <- function(mode) {
    result <- str_replace_all(mode, c(
        'hard'='Hard',
        'rnm \\(true\\)'='RNM',
        'sgv \\(true\\)'='SGV',
        'soft'='Soft'
    ))
    result <- factor(result, levels=c('Hard', 'RNM', 'SGV', 'Soft'))
    return(result)
}

rmse <- function(true, pred) {
    result <- sqrt(mean((true - pred)^2))
    return(result)
}

mean_relative_error <- function(log_true, log_pred) {
    true <- 10^log_true
    pred <- 10^log_pred
    result <- mean(abs((true - pred)/true))
    return(result)
}

accuracy <- function(true, pred) {
    result <- mean(true == pred)
    return(result)
}