#' runSH
#' run the bash command line
#' @param sh character value representing the bash command line
#'
#' @return NULL
#' @export
#'
#' @examples
runSH = function( sh )
{
  cat(sh,'\n\n')
  system(sh)
}
