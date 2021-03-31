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

  mean.slope.c = 1.035733e+01
  mean.slope.b0 = -2.107253e+01
  mean.slope.b1 = 2.122253e-02
  mean.slope.b2 = 4.996048e-04
  mean.int.N = 3.996978e+00
  mean.int.a = -1.404334e+00
  mean.int.k = -3.790485e-07
  spread.a.c = 10.40810900
  spread.a.b0 = -0.82714022
  spread.a.b1 = 0.31269120
  spread.a.b2 = 0.04003754
  spread.N.N = 1.472670e+04
  spread.N.a = -2.294790e+00
  spread.N.k = 1.802925e-07
  resid.uw.N = 2.332219e-01
  resid.uw.a = -6.870340e-01
  resid.uw.k = -2.549764e-06
  zeroes.N = 1.431823e+07
  zeroes.a = -3.774231e-01
  zeroes.k = -6.683311e+04

  dat <- df %>%
    dplyr::mutate(expected.frac = case_when(

      total.sites < exp(mean.slope.c) ~
        (exp(mean.slope.b0 - mean.slope.c*mean.slope.b1)*total.sites^mean.slope.b1)*gene.length +
        ((mean.int.N*total.sites^mean.int.a) + mean.int.k),

      total.sites >= exp(mean.slope.c) ~
        (exp(mean.slope.b0 - mean.slope.c*mean.slope.b2)*total.sites^mean.slope.b2)*gene.length +
        ((mean.int.N*total.sites^mean.int.a) + mean.int.k)),


      upper.spread = case_when(

        total.sites < exp(spread.a.c) ~ ((resid.uw.N*total.sites^resid.uw.a) + resid.uw.k) +
          (((spread.N.N*total.sites^spread.N.a) + spread.N.k)*gene.length^
             (exp(spread.a.b0 - spread.a.c*spread.a.b1)*total.sites^spread.a.b1)),

        total.sites >= exp(spread.a.c) ~ ((resid.uw.N*total.sites^resid.uw.a) + resid.uw.k) +
          (((spread.N.N*total.sites^spread.N.a) + spread.N.k)*gene.length^
             (exp(spread.a.b0 - spread.a.c*spread.a.b2)*total.sites^spread.a.b2))),

      lower.spread = upper.spread,

      upper.limit = expected.frac + upper.spread,

      lower.limit = expected.frac - lower.spread,

      lower.limit = replace(lower.limit, lower.limit < 0, 0),

      lower.limit = replace(lower.limit, gene.length < (zeroes.N*total.sites^zeroes.a) + zeroes.k, 0),

      outlier.type = case_when(frac.int < upper.limit & frac.int > lower.limit ~ "NONE",
                               frac.int > upper.limit ~ "UPPER",
                               frac.int < lower.limit ~ "LOWER")) %>%
    dplyr::select(-c(upper.spread, lower.spread))

  return(dat) }
