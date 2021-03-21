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

      mean.slope.c = 10.012633458
      mean.slope.b0 = -21.108738592
      mean.slope.b1 = 0.019483919
      mean.slope.b2 = 0.002354043
      mean.int.N = 6.303895e-01
      mean.int.a = -1.196897e+00
      mean.int.k = 6.339807e-08
      spread.a.c = 10.64781989
      spread.a.b0 = -0.81298244
      spread.a.b1 = 0.21913184
      spread.a.b2 = 0.04130661
      spread.N.N = 9.572647e+02
      spread.N.a = -2.020024e+00
      spread.N.k = 1.522057e-07
      spread.uw.N = 2.846307e-01
      spread.uw.a = -7.115993e-01
      spread.uw.k = 1.518783e-05
      zeroes.N = 1.667800e+07
      zeroes.a = -3.916833e-01
      zeroes.k = -5.680890e+04 }

  if(pattern=="MB"){

      mean.slope.c = 1.012329e+01
      mean.slope.b0 = -2.106895e+01
      mean.slope.b1 = 1.570697e-02
      mean.slope.b2 = 3.726597e-04
      mean.int.N = 1.557550e+00
      mean.int.a = -1.309916e+00
      mean.int.k = -1.481097e-06
      spread.a.c = 10.3771405
      spread.a.b0 = -0.8089277
      spread.a.b1 = 0.2667871
      spread.a.b2 = 0.0487804
      spread.N.N = 1.356586e+04
      spread.N.a = -2.310493e+00
      spread.N.k = 1.527404e-07
      spread.uw.N = 2.421844e-01
      spread.uw.a = -6.887260e-01
      spread.uw.k = 1.389093e-06
      zeroes.N = 1.661793e+07
      zeroes.a = -3.895366e-01
      zeroes.k = -5.995046e+04 }

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
