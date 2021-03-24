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

      mean.slope.c = 1.062591e+01
      mean.slope.b0 = -2.109951e+01
      mean.slope.b1 = 1.419552e-02
      mean.slope.b2 = 7.373111e-04
      mean.int.N = 8.907227e-01
      mean.int.a = -1.236737e+00
      mean.int.k = 7.155503e-07
      spread.a.c = 10.41224339
      spread.a.b0 = -0.83053680
      spread.a.b1 = 0.29211443
      spread.a.b2 = 0.03722969
      spread.N.N = 1.216210e+04
      spread.N.a = -2.281313e+00
      spread.N.k = 1.906715e-07
      spread.uw.N = 2.843457e-01
      spread.uw.a = -7.115858e-01
      spread.uw.k = 1.260564e-05
      zeroes.N = 1.500378e+07
      zeroes.a = -3.843191e-01
      zeroes.k = -5.899074e+04 }

  if(pattern=="MB"){

      mean.slope.c = 10.300985961
      mean.slope.b0 = -21.068387171
      mean.slope.b1 = 0.028349258
      mean.slope.b2 = 0.001733579
      mean.int.N = 1.926378e+00
      mean.int.a = -1.316746e+00
      mean.int.k = -7.921794e-07
      spread.a.c = 10.60512491
      spread.a.b0 = -0.80427803
      spread.a.b1 = 0.22991466
      spread.a.b2 = 0.03864262
      spread.N.N = 2.096778e+03
      spread.N.a = -2.105591e+00
      spread.N.k = 1.569561e-07
      spread.uw.N = 2.590415e-01
      spread.uw.a = -6.984847e-01
      spread.uw.k = 8.456193e-07
      zeroes.N = 1.706085e+07
      zeroes.a = -3.981275e-01
      zeroes.k = -5.489989e+04 }

  dat <- df %>%
    dplyr::mutate(expected.frac = case_when(

      total.sites < exp(mean.slope.c) ~
        (exp(mean.slope.b0 - mean.slope.c*mean.slope.b1)*total.sites^mean.slope.b1)*gene.length +
        ((mean.int.N*total.sites^mean.int.a) + mean.int.k),

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
