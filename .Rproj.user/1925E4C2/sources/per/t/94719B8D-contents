#' Model generation based on random fragmentation
#'
#' Calculates the expected integration frequency, the upper boundary, and the lower boundary for each sample.
#' Defines outliers based on the relation to the calculated upper and lower boundaries.
#'
#'@param df The dataframe generated with \code{int_file_import}
#'@param total.sites The column in the \code{df} dataframe containing the total number of integration sites.
#'@param gene.length The column in the \code{df} dataframe containing the gene length of a given gene.
#'
#' @examples random_frag_model()
#'
#' @export
random_frag_model <- function(df){
  print("Pattern is random")
  dat <- df %>%
    dplyr::mutate(expected.frac = case_when(

      total.sites < exp(10.051250567) ~
        (exp(-21.125319132 - 10.051250567*0.036524692)*total.sites^0.036524692)*gene.length +
        ((1.188053e+01*total.sites^-1.520168e+00) + -3.238539e-07),

      total.sites >= exp(10.051250567) ~
        (exp(-21.125319132 - 10.051250567*0.001037773)*total.sites^0.001037773)*gene.length +
        ((1.188053e+01*total.sites^-1.520168e+00) + -3.238539e-07)),


      upper.spread = case_when(

        total.sites < exp(10.63798039) ~ ((2.124986e-01*total.sites^-6.817895e-01) + -3.264780e-06) +
          (((5.748894e+03*total.sites^-2.210068e+00) + 1.751955e-07)*gene.length^
             (exp(-0.81554604 - 10.63798039*0.22855264)*total.sites^0.22855264)),

        total.sites >= exp(10.63798039) ~ ((2.124986e-01*total.sites^-6.817895e-01) + -3.264780e-06) +
          (((5.748894e+03*total.sites^-2.210068e+00) + 1.751955e-07)*gene.length^
             (exp(-0.81554604 - 10.63798039*0.03961336)*total.sites^0.03961336))),


      lower.spread = upper.spread,


      upper.limit = expected.frac + upper.spread,

      lower.limit = expected.frac - lower.spread,

      lower.limit = replace(lower.limit, lower.limit < 0, 0),

      lower.limit = replace(lower.limit, gene.length < (1.543353e+07*total.sites^-3.866596e-01) + -6.214747e+04, 0),

      outlier.type = case_when(frac.int < upper.limit & frac.int > lower.limit ~ "NONE",
                               frac.int > upper.limit ~ "UPPER",
                               frac.int < lower.limit ~ "LOWER")) %>%
    dplyr::select(-c(upper.spread, lower.spread))

  return(dat) }
