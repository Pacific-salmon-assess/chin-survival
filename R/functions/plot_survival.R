# function to plot cumulative survival w/ or w/out sims
plot_surv <- function(x, show_mcmc = TRUE) {
  #subset for plotting
  set.seed(123)
  trials <- sample(1:length(unique(x$iter)), 50, replace = F)
  dat <- x %>% 
    filter(iter %in% trials)
  mu_dat <- x %>% 
    dplyr::select(-iter, -est) %>% 
    distinct()
  
  fill_pal <- c("white", "red")
  names(fill_pal) <- c("phi", "beta")
  
  p <- ggplot() +
    geom_pointrange(data = mu_dat, 
                    aes(x = fct_reorder(segment_name, segment), 
                        y = median, ymin = low, ymax = up, fill = par),
                    shape = 21) +
    ggsidekick::theme_sleek() +
    scale_fill_manual(values = fill_pal) +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          legend.text=element_text(size = 9),
          legend.title = element_blank(),
          axis.text.x = element_text(size = rel(.8))) +
    guides(fill = "none") +
    lims(y = c(0, 1)) 
  
  if (show_mcmc == TRUE) {
    p <- p +
      geom_line(data = dat, 
                aes(x = segment_name, y = est, group = iter), alpha = 0.1)
  }
  
  if (!is.na(x$agg_name_f[1])) {
    title <- as.character(unique(x$agg_name_f))
    p <- p +
      labs(title = title)
  }
  
  return(p)
}
