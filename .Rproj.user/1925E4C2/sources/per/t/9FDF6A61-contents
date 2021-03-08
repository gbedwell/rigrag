#' Sample cross-sample occurrence
#'
#' Calculates the overlap in outlier genes across each sample.
#'
#'@param df The dataframe generated with \code{random_frag_model()} that contains the outlier information.
#'@param uo_list The combined list of genes corresponding to upper outliers. Requires a "gene.name" column header.
#'@param lo_list The combined list of genes corresponding to lower outliers. Requires a "gene.name" column header.
#'
#'
#'
#' @examples total_cso()
#'
#' @export
total_cso <- function(df, uo_list, lo_list){
  dat1 <- df %>%
    dplyr::inner_join(., uo_list, by="gene.name") %>%
    dplyr::select(gene.name) %>%
    dplyr::group_by(gene.name) %>%
    dplyr::add_count(gene.name, name="cso") %>%
    dplyr::ungroup() %>%
    dplyr::distinct(gene.name, .keep_all=TRUE) %>%
    dplyr::group_by(cso) %>%
    dplyr::add_count(cso, name="total.count") %>%
    dplyr::select(cso, total.count) %>%
    dplyr::distinct() %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(cso)) %>%
    dplyr::mutate(frac.total = total.count/sum(total.count),
                  sample="rigs")

  dat2 <- df %>%
    dplyr::inner_join(., lo_list, by="gene.name") %>%
    dplyr::select(gene.name) %>%
    dplyr::group_by(gene.name) %>%
    dplyr::add_count(gene.name, name="cso") %>%
    dplyr::ungroup() %>%
    dplyr::distinct(gene.name, .keep_all=TRUE) %>%
    dplyr::group_by(cso) %>%
    dplyr::add_count(cso, name="total.count") %>%
    dplyr::select(cso, total.count) %>%
    dplyr::distinct() %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(cso)) %>%
    dplyr::mutate(frac.total = total.count/sum(total.count),
                  sample="rags")


  dat3 <- df %>%
    dplyr::anti_join(., uo_list, by="gene.name") %>%
    dplyr::anti_join(., lo_list, by="gene.name") %>%
    dplyr::select(gene.name) %>%
    dplyr::group_by(gene.name) %>%
    dplyr::add_count(gene.name, name="cso") %>%
    dplyr::ungroup() %>%
    dplyr::distinct(gene.name, .keep_all=TRUE) %>%
    dplyr::group_by(cso) %>%
    dplyr::add_count(cso, name="total.count") %>%
    dplyr::select(cso, total.count) %>%
    dplyr::distinct() %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(cso)) %>%
    dplyr::mutate(frac.total = total.count/sum(total.count),
                  sample="nos")

  dat.combo <- rbind(dat1, dat2, dat3)

  return(dat.combo)

}

