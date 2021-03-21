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

  mean.slope.c = 9.478231367
  mean.slope.b0 = -21.080243617
  mean.slope.b1 = 0.099308243
  mean.slope.b2 = 0.001392479
  mean.int.N = 1.372085e+01
  mean.int.a = -1.535794e+00
  mean.int.k = -1.076369e-06
  spread.a.c = 10.41544528
  spread.a.b0 = -0.81456301
  spread.a.b1 = 0.24888476
  spread.a.b2 = 0.04385694
  spread.N.N = 2.930537e+02
  spread.N.a = -1.906554e+00
  spread.N.k = 1.210055e-07
  spread.uw.N = 2.949172e-01
  spread.uw.a = -7.094517e-01
  spread.uw.k = 2.482763e-06
  zeroes.N = 1.633124e+07
  zeroes.a = -3.887341e-01
  zeroes.k = -6.077641e+04

  dat <- df %>%
    dplyr::mutate(expected.frac = case_when(

      total.sites < exp(mean.slope.c) ~
        (exp(mean.slope.b0 - mean.slope.c*mean.slope.b1)*total.sites^mean.slope.b1)*gene.length +
        ((mean.int.N*total.sites^mean.nt.a) + mean.int.k),

      total.sites >= exp(mean.slope.c) ~
        (exp(mean.slope.b0 - mean.slope.c*mean.slope.b2)*total.sites^mean.slope.b2)*gene.length +
        ((mean.int.N*total.sites^mean.int.a) + mean.int.k)),


      upper.spread = case_when(

        total.sites < exp(spread.a.c) ~ ((spread.uw.N*total.sites^spread.uw.a) + spread.uw.k) +
          (((spread.N.N*total.sites^spread.N.a) + spread.N.k)*gene.length^
             (exp(spread.a.b0 - spread.a.c*spread.a.b1)*total.sites^spread.a.b1)),

        total.sites >= exp(spread.a.c) ~ ((spread.uw.N*total.sites^spread.uw.a) + spread.uw.k) +
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
