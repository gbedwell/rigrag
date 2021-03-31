#' Model generation based on pattern fragmentation
#'
#' Calculates the expected integration frequency, the upper boundary, and the lower boundary for each sample.
#' Defines outliers based on the relation to the calculated upper and lower boundaries.
#'
#'@param df The dataframe generated with \code{int_file_import}
#'@param pattern The character string defining the enzymes used to digest the genomic DNA. Currenty available patterns are \code{"MB"} for MseI/BglII digestion and \code{"NASB"} for NheI/AvrII/SpeI/BamHI digestion.
#'
#' @examples pattern_frag_model()
#'
#' @export
pattern_frag_model <- function(df, pattern=c("NASB", "MB")){
  print(paste0("Pattern is ", pattern))
  if(pattern=="NASB"){

      mean.slope.c = 10.046065270
      mean.slope.b0 = -21.100594914
      mean.slope.b1 = 0.024715032
      mean.slope.b2 = 0.001236263
      mean.int.N = 1.680105e+00
      mean.int.a = -1.304389e+00
      mean.int.k = 7.883625e-07
      spread.a.c = 10.4352508
      spread.a.b0 = -0.8436968
      spread.a.b1 = 0.3005003
      spread.a.b2 = 0.0418972
      spread.N.N = 2.300702e+04
      spread.N.a = -2.337391e+00
      spread.N.k = 2.118615e-07
      resid.uw.N = 2.252665e-01
      resid.uw.a = -6.865222e-01
      resid.uw.k = 8.206164e-06
      zeroes.N = 1.601822e+07
      zeroes.a = -3.921552e-01
      zeroes.k = -5.664448e+04 }

  if(pattern=="MB"){

      mean.slope.c = 11.137546445
      mean.slope.b0 = -21.063221901
      mean.slope.b1 = 0.009444708
      mean.slope.b2 = 0.000148406
      mean.int.N = 1.144371e+00
      mean.int.a = -1.273171e+00
      mean.int.k = -7.486481e-07
      spread.a.c = 10.41440181
      spread.a.b0 = -0.82287674
      spread.a.b1 = 0.27678204
      spread.a.b2 = 0.04178149
      spread.N.N = 1.175781e+03
      spread.N.a = -2.039190e+00
      spread.N.k = 1.459199e-07
      resid.uw.N = 2.796058e-01
      resid.uw.a = -7.044471e-01
      resid.uw.k = -3.751581e-09
      zeroes.N = 1.600084e+07
      zeroes.a = -3.906015e-01
      zeroes.k = -5.920290e+04 }

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
