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

  mean.slope.c = 9.72733328
  mean.slope.b0 = -21.07957158
  mean.slope.b1 = 0.06159223
  mean.slope.b2 = 0.00190597
  mean.int.N = 9.859066e+00
  mean.int.a = -1.496102e+00
  mean.int.k = -3.301019e-07
  spread.a.c = 10.21677389
  spread.a.b0 = -0.85257513
  spread.a.b1 = 0.33256170
  spread.a.b2 = 0.04632329
  spread.N.N = 9.359475e+03
  spread.N.a = -2.255893e+00
  spread.N.k = 1.905592e-07
  resid.uw.N = 2.647977e-01
  resid.uw.a = -6.995866e-01
  resid.uw.k = -1.602546e-06
  zeroes.N = 1.527221e+07
  zeroes.a = -3.853157e-01
  zeroes.k = -6.294615e+04

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
