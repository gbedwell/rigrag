#' Single file outliers output
#'
#' Outputs a single unique gene list for the combined upper and lower outliers for all samples
#'
#'@param df The dataframe generated with \code{random_frag_model()} that contains the outlier information.
#'@param output.dir The directory in which the file will be output. If the working directory is the desired directory, enter \code{NULL}.
#'@param file.prefix The desired prefix of the output files as a character string. The function will add "upper/lower_outliers.txt" to the prefix.
#'
#' @examples outliers_output()
#'
#' @export
outliers_output <- function(df, output.dir=NULL, file.prefix){
  if(is.null(output.dir)){
    dat <- df %>%
      dplyr::filter(outlier.type == "UPPER") %>%
      dplyr::select(gene.name) %>%
      dplyr::distinct(gene.name) %>%
      dplyr::arrange(gene.name) %>%
      dplyr::group_walk(~write.table(.x, file = paste0(getwd(), "/", file.prefix,"_upper_outliers.txt"),
                                     append = FALSE, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE))

    dat <- df %>%
      dplyr::filter(outlier.type == "LOWER") %>%
      dplyr::select(gene.name) %>%
      dplyr::distinct(gene.name) %>%
      dplyr::arrange(gene.name) %>%
      dplyr::group_walk(~write.table(.x, file = paste0(getwd(), "/", file.prefix,"_lower_outliers.txt"),
                                     append = FALSE, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)) }
  else{
    dat <- df %>%
      dplyr::filter(outlier.type == "UPPER") %>%
      dplyr::select(gene.name) %>%
      dplyr::distinct(gene.name) %>%
      dplyr::arrange(gene.name) %>%
      dplyr::group_walk(~write.table(.x, file = paste0(getwd(), "/", output.dir, "/", file.prefix,"_upper_outliers.txt"),
                                     append = FALSE, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE))

    dat <- df %>%
      dplyr::filter(outlier.type == "LOWER") %>%
      dplyr::select(gene.name) %>%
      dplyr::distinct(gene.name) %>%
      dplyr::arrange(gene.name) %>%
      dplyr::group_walk(~write.table(.x, file = paste0(getwd(), "/", output.dir, "/", file.prefix,"_lower_outliers.txt"),
                                     append = FALSE, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)) }
}
