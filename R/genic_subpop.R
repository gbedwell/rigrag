#' Quantifies the relative proportions of each outlier type per sample
#'
#'
#'@param df The dataframe generated with the model calculation functions
#'@param output.type Defines whether the output values are \code{"relative"} or \code{"absolute"}.
#'\code{"relative"} values express the population fractions relative to total integration.
#'\code{"absolute"} values expression the population fraction relative to genic integration.
#'
#' @examples genic_subpop()
#'
#' @export
genic_subpop <- function(df, output.type=c("relative","absolute")){
  if(output.type == "relative"){
    dat <- df %>%
      dplyr::mutate(genic.frac = genic.sites/total.sites) %>%
      dplyr::group_by(sample, outlier.type) %>%
      dplyr::mutate(int.sum = sum(int.counts)) %>%
      dplyr::mutate(pop.frac = (int.sum/genic.sites)*genic.frac) %>%
      dplyr::select(sample, genic.frac, outlier.type, pop.frac) %>%
      dplyr::distinct(sample, outlier.type, .keep_all=TRUE) %>%
      dplyr::ungroup() }

  if(output.type == "absolute"){
    dat <- df %>%
      dplyr::mutate(genic.frac = genic.sites/total.sites) %>%
      dplyr::group_by(sample, outlier.type) %>%
      dplyr::mutate(int.sum = sum(int.counts)) %>%
      dplyr::mutate(pop.frac = (int.sum/genic.sites)) %>%
      dplyr::select(sample, genic.frac, outlier.type, pop.frac) %>%
      dplyr::distinct(sample, outlier.type, .keep_all=TRUE) %>%
      dplyr::ungroup() }

  return(dat)
}

