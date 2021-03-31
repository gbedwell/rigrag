#' Sample cross-sample occurrence
#'
#' Calculates the overlap in outlier genes across each sample.
#'
#'@param df The dataframe generated with \code{random_frag_model()} or \code{pattern_frag_model()} that contains the outlier information.
#'@param output.type Determines whether the output is the cross-sample occurrence across each gene in each sample (\code{"by_gene"}) or a summary of the cross-sample occurrences across samples (\code{"summary"}).
#'@param outlier.type Determines whether cross-sample occurrence is calculated for each outlier type (\code{"UPPER"}, \code{"LOWER"}, \code{"NONE"}).
#'
#'Note that in order to obtain unique non-outliers, the \code{"NONE"} output will have to be further filtered to remove any/all overlap with the \code{"UPPER"} and \code{"LOWER"} outputs.
#'
#' @examples cross_sample_occurrence()
#'
#' @export
cross_sample_occurrence <- function(df, output.type=c("by_gene", "summary"), outlier.type=c("UPPER", "LOWER", "NONE")){
  if(outlier.type=="UPPER"){
    if(output.type=="by_gene"){
      upper.outliers <- df %>%
        dplyr::filter(outlier.type == "UPPER") %>%
        dplyr::select(gene.name, sample) %>%
        dplyr::group_by(gene.name) %>%
        dplyr::add_count(gene.name, name="cso") %>%
        dplyr::ungroup()

      return(upper.outliers) }

    if(output.type=="summary"){
      upper.outliers <- df %>%
        dplyr::filter(outlier.type == "UPPER") %>%
        dplyr::select(gene.name, sample) %>%
        dplyr::group_by(gene.name) %>%
        dplyr::add_count(gene.name, name="cso") %>%
        dplyr::ungroup() %>%
        dplyr::group_by(sample, cso) %>%
        dplyr::add_count(cso, name="total.count") %>%
        dplyr::distinct(sample, total.count) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(sample) %>%
        dplyr::mutate(frac.total = total.count/sum(total.count)) %>%
        dplyr::arrange(sample, cso) %>%
        dplyr::ungroup()

      upper.outliers.total <- upper.outliers %>%
        dplyr::group_by(cso) %>%
        dplyr::summarize(total.count = sum(total.count)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(frac.total = total.count/sum(total.count)) %>%
        dplyr::mutate(sample = "Total") %>%
        dplyr::select(sample, cso, total.count, frac.total)

      upper.outliers <- rbind(upper.outliers, upper.outliers.total) }

    return(upper.outliers) }

  if(outlier.type=="LOWER"){
    if(output.type=="by_gene"){
      lower.outliers <- df %>%
        dplyr::filter(outlier.type == "LOWER") %>%
        dplyr::select(gene.name, sample) %>%
        dplyr::group_by(gene.name) %>%
        dplyr::add_count(gene.name, name="cso") %>%
        dplyr::ungroup()

      return(lower.outliers) }

    if(output.type=="summary"){
      lower.outliers <- df %>%
        dplyr::filter(outlier.type == "LOWER") %>%
        dplyr::select(gene.name, sample) %>%
        dplyr::group_by(gene.name) %>%
        dplyr::add_count(gene.name, name="cso") %>%
        dplyr::ungroup() %>%
        dplyr::group_by(sample, cso) %>%
        dplyr::add_count(cso, name="total.count") %>%
        dplyr::distinct(sample, total.count) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(sample) %>%
        dplyr::mutate(frac.total = total.count/sum(total.count)) %>%
        dplyr::arrange(sample, cso) %>%
        dplyr::ungroup()

      lower.outliers.total <- lower.outliers %>%
        dplyr::group_by(cso) %>%
        dplyr::summarize(total.count = sum(total.count)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(frac.total = total.count/sum(total.count)) %>%
        dplyr::mutate(sample = "Total") %>%
        dplyr::select(sample, cso, total.count, frac.total)

      lower.outliers <- rbind(lower.outliers, lower.outliers.total) }

    return(lower.outliers) }

  if(outlier.type=="NONE"){
    if(output.type=="by_gene"){
      print("To obtain unique outliers, this output should be filtered to remove genes overlapping upper and/or lower outliers!")
      non.outliers <- df %>%
        dplyr::filter(outlier.type == "NONE") %>%
        dplyr::select(gene.name, sample) %>%
        dplyr::group_by(gene.name) %>%
        dplyr::add_count(gene.name, name="cso") %>%
        dplyr::ungroup()

      return(non.outliers) }

    if(output.type=="summary"){
      print("This output might contain upper or lower outliers that were identified as such in other samples!")
      print("Use the by_gene option and filter against upper and lower outliers to obtain unique non-outlier counts!")
      non.outliers <- df %>%
        dplyr::filter(outlier.type == "NONE") %>%
        dplyr::select(gene.name, sample) %>%
        dplyr::group_by(gene.name) %>%
        dplyr::add_count(gene.name, name="cso") %>%
        dplyr::ungroup() %>%
        dplyr::group_by(sample, cso) %>%
        dplyr::add_count(cso, name="total.count") %>%
        dplyr::distinct(sample, total.count) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(sample) %>%
        dplyr::mutate(frac.total = total.count/sum(total.count)) %>%
        dplyr::arrange(sample, cso) %>%
        dplyr::ungroup()

      non.outliers.total <- non.outliers %>%
        dplyr::group_by(cso) %>%
        dplyr::summarize(total.count = sum(total.count)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(frac.total = total.count/sum(total.count)) %>%
        dplyr::mutate(sample = "Total") %>%
        dplyr::select(sample, cso, total.count, frac.total)

      non.outliers <- rbind(non.outliers, non.outliers.total) }

    return(non.outliers) }

}

