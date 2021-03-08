#' Visualization of cross-sample occurrence
#'
#' Plots the cross-sample overlap between samples.
#'
#'@param df The dataframe generated with \code{random_frag_model()}
#'@param sample.names Character string defining the sample names along the x-axis.
#'@param sample.num The number of samples being analyzed. This is used to extend the color palette in case there are more samples than colors.
#'@param gene.type Used for the y-axis title, denotes either \code{"RIGs"} or \code{"RAGs"}.
#'@param write.out Determines whether or not to save the generated plot. To save the plots, define as \code{"yes"}.
#'@param h The desired height of the generated plot
#'@param w The desired width of the generated plot
#'@param output.dir The directory in which the files will be output. If the working directory is the desired directory, enter \code{NULL}.
#'@param file.name The name of the output file. Enter as a character string.
#'
#' @examples cso_plot()
#'
#' @export
cso_plot <- function(df, sample.num, sample.names, gene.type=c("RIGs", "RAGs"), write.out=c("yes", "no"), h=NULL, w=NULL, output.dir=NULL, file.name=NULL){

  if(sample.num > 9){
    require(RColorBrewer)

    mycolors <- colorRampPalette(brewer.pal(9, "Blues"))(sample.num)

    df <- df %>% dplyr::group_by(sample) %>%
      dplyr::arrange(cso)

    df$cso <- factor(df$cso, levels=as.character(1:sample.num))


    dat <- ggplot(df, aes(x=sample, y=frac.total, fill=forcats::fct_rev(cso))) +
      geom_bar(stat="identity", color = "black", width = 0.75) +
      scale_fill_manual("CSO",
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
                         axis.text.x=element_text(size=12),
                         axis.text.y=element_text(size=14),
                         axis.title=element_text(size=16)) +
      scale_x_discrete(labels = sample.names) +
      labs(x = "Sample", y = paste0("Fraction of Total ", gene.type))

    if(write.out=="yes"){
      if(is.null(output.dir)){
        ggsave(filename=paste0(getwd(), "/", file.name), plot=dat, height=h, width=w, units="in", device="tiff") }
      else{
        ggsave(filename=paste0(getwd(), "/", output.dir, "/", file.name), plot=plots, height=h, width=w, units="in", device="tiff") } }
    if(write.out=="no"){
      return(dat) } }

  else{
    df <- df %>% dplyr::group_by(sample) %>%
      dplyr::arrange(cso)

    df$cso <- factor(df$cso, levels=as.character(1:sample.num))

    dat <- ggplot(df, aes(x=sample, y=frac.total, fill=forcats::fct_rev(cso))) +
      geom_bar(stat="identity", color = "black", width = 0.75) +
      scale_fill_brewer("CSO",
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
                         axis.text.x=element_text(size=12),
                         axis.text.y=element_text(size=14),
                         axis.title=element_text(size=16)) +
      scale_x_discrete(labels = sample.names) +
      labs(x = "Sample", y = paste0("Fraction of Total ", gene.type))

    if(write.out=="yes"){
      if(is.null(output.dir)){
        ggsave(filename=paste0(getwd(), "/", file.name), plot=dat, height=h, width=w, units="in", device="tiff") }
      else{
        ggsave(filename=paste0(getwd(), "/", output.dir, "/", file.name), plot=plots, height=h, width=w, units="in", device="tiff") } }
    if(write.out=="no"){
      return(dat) } }
}
