#' Outlier visualization
#'
#' Generates separate plots of genic integration frequency vs. gene length for each sample being analyzed.
#' Depicts model boundaries within gray shaded areas and colors upper and lower outliers blue and red, respectively.
#'
#'@param df The dataframe generated with \code{random_frag_model()}
#'@param labels The labels to be assigned to each panel. Use \code{"AUTO"} to automatically assign alphabetical labels starting with "A". For custom labels, enter them as a character vector.
#'@param label.size The size of the labels on the combined plot output.
#'@param output.type If a grid of plots is desired define as \code{"combo"}, if individual plots are desired define as \code{"single"}.
#'@param write.out Determines whether or not to save the generated plot. To save the plots, define as \code{"yes"}.
#'@param num.col The desired number of columns in the cowplot output
#'@param h The desired height of the generated plot
#'@param w The desired width of the generated plot
#'@param output.dir The directory in which the files will be output. If the working directory is the desired directory, enter \code{NULL}.
#'@param file.name The name of the output file. Enter as a character string.
#'
#'
#' @examples model_plot()
#'
#' @export
model_plot <- function(df, output.type=c("single", "combo"), labels=NULL, label.size=NULL, write.out=c("yes", "no"), num.col=NULL, h=NULL, w=NULL, output.dir=NULL, file.name=NULL){
  dat <- df %>%
    dplyr::group_by(sample) %>%
    tidyr::nest() %>%
    dplyr::mutate(plot = map2(data, sample, ~ggplot(data=.x) + aes(x=gene.length, y=expected.frac) +
                                geom_ribbon(aes(ymin = expected.frac, ymax = upper.limit), fill = "grey85") +
                                geom_ribbon(aes(ymin =lower.limit, ymax = expected.frac), fill = "grey85")  +
                                geom_line(aes(y=upper.limit), size = 0.3) +
                                geom_line(aes(y=lower.limit), size = 0.3) +
                                geom_line(size=0.3) +
                                geom_point(data=., aes(x=gene.length, y=frac.int), shape=21, size=2.5) +
                                geom_point(data=., aes(x=gene.length, y=frac.int, fill=outlier.type),
                                           color="black", alpha=0.7, shape=21, size=2.5) +
                                scale_fill_manual(values = c("LOWER"="darkred","NONE"="gray70","UPPER"="darkblue")) +
                                theme_bw() + theme(plot.title = element_text(face="bold", size=12),
                                                   axis.text=element_text(size=12),
                                                   axis.title=element_text(size=14),
                                                   legend.position = "none") +
                                scale_x_continuous(breaks=c(0,500000,1000000,1500000,2000000),
                                                   labels=c(0, 0.5, 1, 1.5, 2)) +
                                scale_y_continuous(limits=c(0,0.004)) +
                                labs(x = "Gene Length (Mbp)", y = "Fraction of Genic Integration",
                                     title=paste0(sample,
                                                  " - ", scales::comma_format()(unique(.x$total.sites))," total sites"))))
  if(output.type=="combo"){
    plots <- cowplot::plot_grid(plotlist=dat$plot, ncol=num.col, labels=labels, label_size=label.size)
    if(write.out=="yes"){
      if(is.null(output.dir)){
        ggsave(filename=paste0(getwd(), "/", file.name), plot=plots, height=h, width=w, units="in", device="tiff") }
      else{
        ggsave(filename=paste0(getwd(), "/", output.dir, "/", file.name), plot=plots, height=h, width=w, units="in", device="tiff") } }
    if(write.out=="no"){
      return(plots) } }

  if(output.type=="single") {
    if(write.out=="yes"){
      plots <- map2(dat$plot, dat$sample, ~ggsave(filename=paste0(getwd(), "/", output.dir, "/", mgsub(as.character(.y), c(" ","#","+"), c("_","","")), ".tiff"),
                                                  plot=.x, width = w, height = h, units = "in", device="tiff")) }
    if(write.out=="no"){
      plots <- dat$plot
      return(plots) } }

}

