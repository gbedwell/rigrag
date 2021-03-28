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

      mean.slope.c = 10.035189537
      mean.slope.b0 = -21.100853269
      mean.slope.b1 = 0.024386533
      mean.slope.b2 = 0.001593322
      mean.int.N = 1.663138e+00
      mean.int.a = -1.303294e+00
      mean.int.k = 7.883465e-07
      spread.a.c = 10.43594848
      spread.a.b0 = -0.84366465
      spread.a.b1 = 0.30025068
      spread.a.b2 = 0.04189944
      spread.N.N = 1.462464e+04
      spread.N.a = -2.288880e+00
      spread.N.k = 2.004174e-07
      resid.uw.N = 2.256293e-01
      resid.uw.a = -6.867147e-01
      resid.uw.k = 8.262533e-06
      zeroes.N = 1.594914e+07
      zeroes.a = -3.915795e-01
      zeroes.k = -5.695653e+04 }

  if(pattern=="MB"){

      mean.slope.c = 1.094382e+01
      mean.slope.b0 = -2.106403e+01
      mean.slope.b1 = 9.524159e-03
      mean.slope.b2 = 9.048685e-04
      mean.int.N = 1.081564e+00
      mean.int.a = -1.267149e+00
      mean.int.k = -7.477234e-07
      spread.a.c = 10.41568101
      spread.a.b0 = -0.82246561
      spread.a.b1 = 0.27689327
      spread.a.b2 = 0.04165547
      spread.N.N = 1.177050e+03
      spread.N.a = -2.039308e+00
      spread.N.k = 1.456944e-07
      resid.uw.N = 2.802989e-01
      resid.uw.a = -7.047310e-01
      resid.uw.k = 5.283514e-08
      zeroes.N = 1.538530e+07
      zeroes.a = -3.858750e-01
      zeroes.k = -6.132808e+04 }

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
