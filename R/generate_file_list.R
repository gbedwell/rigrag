#' File list construction
#'
#' Generates a list of integration site files in a given directory for import into rigrag.
#'
#' @param directory The directory to search. Entered as a character string. Can be called with \code{getwd()} if the directory of interest is the current working directory.
#' The function is written to search recursively, so all subdirectories within the named directory are searched.
#'
#' @param pattern The common file naming pattern to search for. Input as a character string. Typically something like \code{"*genic.bed"}.
#'
#' @examples generate_file_list()
#'
#' @export
generate_file_list <- function(directory, pattern){
  files <- list.files(directory,
                      pattern=pattern,
                      full.names = TRUE,
                      recursive = TRUE)
}
