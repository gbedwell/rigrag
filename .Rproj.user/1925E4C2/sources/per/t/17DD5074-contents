#' Quantifies and plots the relative proportions of each outlier type per sample
#'
#'
#'@param df The dataframe generated with the model calculation functions
#'@param output.type Defines whether the output values are \code{"relative"} or \code{"absolute"}.
#'\code{"relative"} values express the population fractions relative to total integration.
#'\code{"absolute"} values expression the population fraction relative to genic integration.
#'@param sample.levels Character string defining the factor levels of the samples. Used for ordering the plot.
#'@param sample.names Character string defining the sample names along the x-axis.
#'
#' @examples genic_subpop_plot()
#'
#' @export
genic_subpop_plot <- function(df, output.type=c("relative","absolute"), sample.levels, sample.names){
  if(output.type == "relative"){
    dat <- df %>%
      dplyr::mutate(genic.frac = genic.sites/total.sites) %>%
      dplyr::group_by(sample, outlier.type) %>%
      dplyr::mutate(int.sum = sum(int.counts)) %>%
      dplyr::mutate(pop.frac = (int.sum/genic.sites)*genic.frac) %>%
      dplyr::select(sample, genic.frac, outlier.type, pop.frac) %>%
      dplyr::distinct(sample, outlier.type, .keep_all=TRUE) %>%
      dplyr::ungroup()
    dat$sample <- factor(dat$sample, levels = sample.levels)

    dat$outlier.type <- factor(dat$outlier.type,
                               levels =c("UPPER","LOWER","NONE"))

    plot <- ggplot(dat, aes(x=sample, y=pop.frac, fill=outlier.type)) +
      geom_bar(stat="identity", color = "black", width = 0.75) +
      scale_fill_brewer("Genic Subpopulations",
                        palette = "Blues",
                        guide = guide_legend(reverse = FALSE,
                                             title.position = "top",
                                             nrow = 1,
                                             title.hjust=0.5),
                        labels = c("RIGs", "RAGs", "Non-outliers")) +
      theme_bw() + theme(panel.grid.major.x=element_blank(),
                         text=element_text(size=12),
                         legend.position="top",
                         legend.key.size = unit(0.5, "cm"),
                         legend.key.width = unit(0.5,"cm"),
                         legend.text = element_text(size=12),
                         legend.title=element_text(size=12),
                         legend.box.margin=margin(-10,-10,-10,-10),
                         axis.text.x=element_text(size=12),
                         axis.text.y=element_text(size=14),
                         axis.title=element_text(size=16)) +
      scale_x_discrete(labels = sample.names) +
      scale_y_continuous(limits = c(0,1)) +
      labs(x = "Sample", y = "Fraction of Integration into Genes")

    return(plot) }

  if(output.type == "absolute"){
    dat <- df %>%
      dplyr::mutate(genic.frac = genic.sites/total.sites) %>%
      dplyr::group_by(sample, outlier.type) %>%
      dplyr::mutate(int.sum = sum(int.counts)) %>%
      dplyr::mutate(pop.frac = (int.sum/genic.sites)) %>%
      dplyr::select(sample, genic.frac, outlier.type, pop.frac) %>%
      dplyr::distinct(sample, outlier.type, .keep_all=TRUE) %>%
      dplyr::ungroup()
    dat$sample <- factor(dat$sample, levels = sample.levels)

    dat$outlier.type <- factor(dat$outlier.type,
                               levels =c("UPPER","LOWER","NONE"))

    plot <- ggplot(dat, aes(x=sample, y=pop.frac, fill=outlier.type)) +
      geom_bar(stat="identity", color = "black", width = 0.75) +
      scale_fill_brewer("Genic Subpopulations",
                        palette = "Blues",
                        guide = guide_legend(reverse = FALSE,
                                             title.position = "top",
                                             nrow = 1,
                                             title.hjust=0.5),
                        labels = c("RIGs", "RAGs", "Non-outliers")) +
      theme_bw() + theme(panel.grid.major.x=element_blank(),
                         text=element_text(size=12),
                         legend.position="top",
                         legend.key.size = unit(0.5, "cm"),
                         legend.key.width = unit(0.5,"cm"),
                         legend.text = element_text(size=12),
                         legend.title=element_text(size=12),
                         legend.box.margin=margin(-10,-10,-10,-10),
                         axis.text.x=element_text(size=12),
                         axis.text.y=element_text(size=14),
                         axis.title=element_text(size=16)) +
      scale_x_discrete(labels = sample.names) +
      scale_y_continuous(limits = c(0,1)) +
      labs(x = "Sample", y = "Fraction of Genic Integration")

    return(plot) }

}
