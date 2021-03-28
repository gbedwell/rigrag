#' Integration site file import and summary statistics
#'
#' Imports integration site files.
#' Calculates the number of genic integration sites and the fraction of genic integration events per sample.
#' Calculates gene length per gene.
#'
#' @param files List of integration site files containing the data to be analyzed. Generated with \code{generate_file_list()}
#' @param sample.ids List of sample identifiers for each file in the file list.
#' Defined with \code{list()}.
#' The order of the sample identifiers must match the order of the files listed in the file list.
#' Spaces and special characters like '+' and '#' are OK in the sample identifiers.
#' Files written within package functions named using sample identifiers will have their names cleaned to remove spaces and select special characters.
#' @param collapsed If \code{"TRUE"}, the number of occurrences of each gene per sample are counted and the collapsed data is returned. If \code{"FALSE"}, every integration event into every gene is returned.
#'
#' @examples int_file_import()
#'
#' @export
int_file_import <- function(files, sample.ids, collapsed = c(TRUE, FALSE)){
  sample.list = lapply(files, data.table::fread)
  sample.list <- Map(as.data.frame, sample.list)
  sample.list <- Map(cbind, sample.list, sample = sample.ids)

  if(collapsed==FALSE){
    sample.df <- data.table::rbindlist(sample.list) %>%
      dplyr::group_by(sample) %>%
      magrittr::set_colnames(c("chromosome","int.start","int.end","int.strand","gene.name",
                               "gene.start","gene.end","gene.strand","total.sites","sample")) %>%
      dplyr::mutate(genic.sites=n_distinct(int.start)) }

  if(collapsed==TRUE){
    sample.df <- data.table::rbindlist(sample.list) %>%
      dplyr::group_by(sample) %>%
      magrittr::set_colnames(c("chromosome","int.start","int.end","int.strand","gene.name",
                               "gene.start","gene.end","gene.strand","total.sites","sample")) %>%
      dplyr::mutate(genic.sites=n_distinct(int.start))

    sample.df <- sample.df %>%
      dplyr::group_by(gene.name, sample) %>%
      dplyr::add_count(gene.name, name = "int.counts") %>%
      dplyr::distinct(gene.name, .keep_all=TRUE) %>%
      dplyr::ungroup() %>%
      dplyr::select(-c(chromosome, int.start, int.end, int.strand)) %>%
      dplyr::mutate(frac.int = int.counts/genic.sites) %>%
      dplyr::mutate(gene.length = gene.end - gene.start) %>%
      dplyr::arrange(sample) }

  return(sample.df)

}
