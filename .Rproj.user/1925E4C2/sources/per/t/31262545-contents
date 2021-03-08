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
    dat <- df %>%
      dplyr::mutate(expected.frac = case_when(

        total.sites < exp(10.023112729) ~
          (exp(-21.150989669 - 10.023112729*0.042760180)*total.sites^0.042760180)*gene.length +
          ((4.593904e+00*total.sites^-1.409699e+00) + 7.066821e-07),

        total.sites >= exp(10.023112729) ~
          (exp(-21.150989669 - 10.023112729*0.001398569)*total.sites^0.001398569)*gene.length +
          ((4.593904e+00*total.sites^-1.409699e+00) + 7.066821e-07)),


        upper.spread = case_when(

          total.sites < exp(10.68573493) ~ ((2.428667e-01*total.sites^-6.999891e-01) + 9.016225e-06) +
            (((3.200194e+03*total.sites^-2.136769e+00) + 1.661411e-07)*gene.length^
               (exp(-0.81498460 - 10.68573493*0.24759905)*total.sites^0.24759905)),

          total.sites >= exp(10.68573493) ~ ((2.428667e-01*total.sites^-6.999891e-01) + 9.016225e-06) +
            (((3.200194e+03*total.sites^-2.136769e+00) + 1.661411e-07)*gene.length^
               (exp(-0.81498460 - 10.68573493*0.03551543)*total.sites^0.03551543))),


        lower.spread = upper.spread,


        upper.limit = expected.frac + upper.spread,

        lower.limit = expected.frac - lower.spread,

        lower.limit = replace(lower.limit, lower.limit < 0, 0),

        lower.limit = replace(lower.limit, gene.length < (1.588258e+07*total.sites^-3.910966e-01) + -5.742514e+04, 0),

        outlier.type = case_when(frac.int < upper.limit & frac.int > lower.limit ~ "NONE",
                                 frac.int > upper.limit ~ "UPPER",
                                 frac.int < lower.limit ~ "LOWER")) %>%
      dplyr::select(-c(upper.spread, lower.spread)) }

  if(pattern=="MB"){
    dat <- df %>%
      dplyr::mutate(expected.frac = case_when(

        total.sites < exp(10.039421304) ~
          (exp(-21.245236519 - 10.039421304*0.051181669)*total.sites^0.051181669)*gene.length +
          ((4.055190e-01*total.sites^-1.126081e+00) + 3.641413e-06),

        total.sites >= exp(10.039421304) ~
          (exp(-21.245236519 - 10.039421304*0.003898347)*total.sites^0.003898347)*gene.length +
          ((4.055190e-01*total.sites^-1.126081e+00) + 3.641413e-06)),


        upper.spread = case_when(

          total.sites < exp(10.37650916) ~ ((2.568801e-01*total.sites^-7.070787e-01) + 1.034859e-05) +
            (((1.456483e+04*total.sites^-2.280417e+00) + 2.454798e-07)*gene.length^
               (exp(-0.88266617 - 10.37650916*0.32259087)*total.sites^0.32259087)),

          total.sites >= exp(10.37650916) ~ ((2.568801e-01*total.sites^-7.070787e-01) + 1.034859e-05) +
            (((1.456483e+04*total.sites^-2.280417e+00) + 2.454798e-07)*gene.length^
               (exp(-0.88266617 - 10.37650916*0.03579225)*total.sites^0.03579225))),


        lower.spread = upper.spread,


        upper.limit = expected.frac + upper.spread,

        lower.limit = expected.frac - lower.spread,

        lower.limit = replace(lower.limit, lower.limit < 0, 0),

        lower.limit = replace(lower.limit, gene.length < (1.488892e+07*total.sites^-3.836516e-01) + -6.454891e+04, 0),

        outlier.type = case_when(frac.int < upper.limit & frac.int > lower.limit ~ "NONE",
                                 frac.int > upper.limit ~ "UPPER",
                                 frac.int < lower.limit ~ "LOWER")) %>%
      dplyr::select(-c(upper.spread, lower.spread)) }

  return(dat) }
