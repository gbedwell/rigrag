#' Visualization of cross-sample occurrence
#'
#' Plots the cross-sample overlap between samples.
#'
#'@param df The dataframe generated with \code{random_frag_model()}
#'@param population.names Character string defining the sample names along the x-axis.
#'@param sample.num The number of samples being analyzed. This is used to extend the color palette in case there are more samples than colors.
#'
#' @examples total_cso_plot()
#'
#' @export
total_cso_plot <- function(df, sample.num, population.names){

  if(sample.num > 9){
    require(RColorBrewer)

    df <- df %>% dplyr::group_by(sample) %>%
      dplyr::arrange(cso)

    df$cso <- factor(df$cso, levels=as.character(1:sample.num))

    mycolors <- colorRampPalette(brewer.pal(9, "Blues"))(sample.num)


    dat <- ggplot(df, aes(x=sample, y=frac.total, fill=forcats::fct_rev(cso))) +
      geom_bar(stat="identity", color = "black", width = 0.75) +
      scale_fill_manual("Total CSO",
                        values = mycolors,
                        guide = guide_legend(reverse =  TRUE,
                                             title.position = "top",
                                             nrow = 1,
                                             title.hjust=0.5)) +
      theme_bw() + theme(panel.grid.major.x=element_blank(),
                         text=element_text(size=12),
                         legend.position="top",
                         legend.key.size = unit(0.5, "cm"),
                         legend.key.width = unit(0.5,"cm"),
                         legend.text = element_text(size=12),
                         legend.title=element_text(size=12),
                         legend.box.margin=margin(-10,-10,-10,-10),
                         legend.justification="center",
                         axis.text.x=element_text(size=16),
                         axis.text.y=element_text(size=16),
                         axis.title=element_text(size=18)) +
      scale_x_discrete(labels = population.names) +
      labs(x = "Sample", y = "Fraction of Population")

    return(dat) }

  else{
    df <- df %>% dplyr::group_by(sample) %>%
      dplyr::arrange(cso)

    df$cso <- factor(df$cso, levels=as.character(1:sample.num))

    dat <- ggplot(df, aes(x=sample, y=frac.total, fill=forcats::fct_rev(as.character(cso)))) +
      geom_bar(stat="identity", color = "black", width = 0.75) +
      scale_fill_brewer("Total CSO",
                        palette = "Blues",
                        guide = guide_legend(reverse = FALSE,
                                             title.position = "top",
                                             nrow = 1,
                                             title.hjust=0.5)) +
      theme_bw() + theme(panel.grid.major.x=element_blank(),
                         text=element_text(size=12),
                         legend.position="top",
                         legend.key.size = unit(0.5, "cm"),
                         legend.key.width = unit(0.5,"cm"),
                         legend.text = element_text(size=12),
                         legend.title=element_text(size=12),
                         legend.box.margin=margin(-10,-10,-10,-10),
                         axis.text.x=element_text(size=16),
                         axis.text.y=element_text(size=16),
                         axis.title=element_text(size=18)) +
      scale_x_discrete(labels = population.names) +
      labs(x = "Sample", y = "Fraction of Population")

    return(dat) }


}
