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

      mean.slope.c = 10.337732037
      mean.slope.b0 = -21.102360621
      mean.slope.b1 = 0.019917463
      mean.slope.b2 = 0.001181648
      mean.int.N = 1.503708e+00
      mean.int.a = -1.289750e+00
      mean.int.k = 7.807597e-07
      spread.a.c = 10.15665914
      spread.a.b0 = -0.86566007
      spread.a.b1 = 0.42778151
      spread.a.b2 = 0.04614291
      spread.N.N = 5.146977e+04
      spread.N.a = -2.415631e+00
      spread.N.k = 2.059273e-07
      resid.uw.N = 2.256866e-01
      resid.uw.a = -6.884679e-01
      resid.uw.k = 9.166680e-06
      zeroes.N = 1.585392e+07
      zeroes.a = -3.908851e-01
      zeroes.k = -5.747605e+04 }

  if(pattern=="MB"){

      mean.slope.c = 1.062621e+01
      mean.slope.b0 = -2.106555e+01
      mean.slope.b1 = 1.223730e-02
      mean.slope.b2 = 2.815437e-04
      mean.int.N = 1.217188e+00
      mean.int.a = -1.281141e+00
      mean.int.k = -7.675881e-07
      spread.a.c = 10.69577187
      spread.a.b0 = -0.80094104
      spread.a.b1 = 0.20089455
      spread.a.b2 = 0.04068415
      spread.N.N = 2.811103e+02
      spread.N.a = -1.896742e+00
      spread.N.k = 1.292952e-07
      resid.uw.N = 2.926493e-01
      resid.uw.a = -7.095572e-01
      resid.uw.k = 4.385521e-07
      zeroes.N = 1.648939e+07
      zeroes.a = -3.939781e-01
      zeroes.k = -5.787765e+04 }

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
