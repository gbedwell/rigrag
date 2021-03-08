#' Finds overrepresented genes in patient-derived integration datasets
#'
#' Compares patient integration frequencies to the corresponding frequencies in in vitro samples.
#'
#'@param df1 Dataframe containing patient-derived genes and the corresponding integration frequency values. This should be a "final" list, with all averaging already done (if necessary).
#'@param column Name of the column containing the patient integration frequency values
#'
#' @examples patient_overrep()
#'
#' @export
patient_overrep <- function(df, df2, column){
  data(mean_int_freq)
  column <- enquo(column)

  patient.z.scores <- dplyr::full_join(df, mean_int_freq, by="gene.name") %>%
    tidyr::drop_na() %>%
    dplyr::mutate(fold.change = !!column/mean.freq,
                  z.score = (!!column - mean.freq)/sd.freq,
                  p.from.z = pnorm(z.score, lower.tail=FALSE, log.p = FALSE),
                  q.value = p.adjust(p.from.z, method = "BH", n = length(p.from.z))) %>%
    dplyr::select(gene.name,  fold.change, z.score, p.from.z, q.value) %>%
    dplyr::arrange(gene.name) %>%
    dplyr::mutate(index = row_number())

  return(patient.z.scores)

  }
