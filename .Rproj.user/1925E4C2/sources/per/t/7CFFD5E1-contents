#' Sample-specific outliers output
#'
#' Outputs a gene list for each sample analyzed and the specified outlier types.
#'
#'@param df The dataframe generated with \code{random_frag_model()} that contains the outlier information.
#'@param output.dir The directory in which the files will be output. If the working directory is the desired directory, enter \code{NULL}.
#'@param outliers The type of outliers that you want output (\code{"both"}, \code{"upper"}, or \code{"lower"}).
#'@param data.depth \code{"list"} returns a list of genes classified as outliers for each sample. \code{complete} returns all variables for each sample.
#'
#' @examples separate_outliers_output()
#'
#' @export
separate_outliers_output <- function(df, outliers=c("both","upper","lower"), output.dir=NULL){
  if(is.null(output.dir)){
    if(outliers=="BOTH"){
      dat <- df %>% dplyr::group_by(sample) %>%
        dplyr::filter(outlier.type == "UPPER") %>%
        dplyr::select(gene.name) %>%
        dplyr::group_walk(~ write.table(.x, file = paste0(getwd(),"/", mgsub(as.character(.y$sample),c(" ", "#","+"),c("_","","")),"_upper_outliers.txt")),
                          append = FALSE, sep = "\t", row.names=FALSE,
                          col.names=TRUE, quote=FALSE)

      dat <- df %>% dplyr::group_by(sample) %>%
        dplyr::filter(outlier.type == "LOWER") %>%
        dplyr::select(gene.name) %>%
        dplyr::group_walk(~ write.table(.x, file = paste0(getwd(),"/", mgsub(as.character(.y$sample),c(" ", "#","+"),c("_","","")),"_lower_outliers.txt")),
                          append = FALSE, sep = "\t", row.names=FALSE,
                          col.names=TRUE, quote=FALSE) }
    if(outliers=="UPPER"){
      dat <- df %>% dplyr::group_by(sample) %>%
        dplyr::filter(outlier.type == "UPPER") %>%
        dplyr::select(gene.name) %>%
        dplyr::group_walk(~ write.table(.x, file = paste0(getwd(),"/", mgsub(as.character(.y$sample),c(" ", "#","+"),c("_","","")),"_upper_outliers.txt")),
                          append = FALSE, sep = "\t", row.names=FALSE,
                          col.names=TRUE, quote=FALSE) }
    if(outliers=="LOWER"){
      dat <- df %>% dplyr::group_by(sample) %>%
        dplyr::filter(outlier.type == "LOWER") %>%
        dplyr::select(gene.name) %>%
        dplyr::group_walk(~ write.table(.x, file = paste0(getwd(),"/", mgsub(as.character(.y$sample),c(" ", "#", "+"),c("_","","")),"_lower_outliers.txt")),
                          append = FALSE, sep = "\t", row.names=FALSE,
                          col.names=TRUE, quote=FALSE) } }

  else{
    if(outliers=="BOTH"){
      dat <- df %>% dplyr::group_by(sample) %>%
        dplyr::filter(outlier.type == "UPPER") %>%
        dplyr::select(gene.name) %>%
        dplyr::group_walk(~ write.table(.x, file = paste0(getwd(),"/", output.dir, "/", mgsub(as.character(.y$sample),c(" ", "#","+"),c("_","","")),"_upper_outliers.txt")),
                          append = FALSE, sep = "\t", row.names=FALSE,
                          col.names=TRUE, quote=FALSE)

      dat <- df %>% dplyr::group_by(sample) %>%
        dplyr::filter(outlier.type == "LOWER") %>%
        dplyr::select(gene.name) %>%
        dplyr::group_walk(~ write.table(.x, file = paste0(getwd(),"/", output.dir, "/", mgsub(as.character(.y$sample),c(" ", "#","+"),c("_","","")),"_lower_outliers.txt")),
                          append = FALSE, sep = "\t", row.names=FALSE,
                          col.names=TRUE, quote=FALSE) }

    if(outliers=="UPPER"){
      dat <- df %>% dplyr::group_by(sample) %>%
        dplyr::filter(outlier.type == "UPPER") %>%
        dplyr::select(gene.name) %>%
        dplyr::group_walk(~ write.table(.x, file = paste0(getwd(),"/", output.dir, "/", mgsub(as.character(.y$sample),c(" ", "#","+"),c("_","","")),"_upper_outliers.txt")),
                          append = FALSE, sep = "\t", row.names=FALSE,
                          col.names=TRUE, quote=FALSE) }
    if(outliers=="LOWER"){
      dat <- df %>% dplyr::group_by(sample) %>%
        dplyr::filter(outlier.type == "LOWER") %>%
        dplyr::select(gene.name) %>%
        dplyr::group_walk(~ write.table(.x, file = paste0(getwd(),"/", output.dir, "/", mgsub(as.character(.y$sample),c(" ", "#", "+"),c("_","","")),"_lower_outliers.txt")),
                          append = FALSE, sep = "\t", row.names=FALSE,
                          col.names=TRUE, quote=FALSE) } }
}
