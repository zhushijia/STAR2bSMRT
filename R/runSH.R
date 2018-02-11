#' runSH
#'
#' @param sh 
#'
#' @return
#' @export
#'
#' @examples
runSH = function( sh )
{
  cat(sh)
  system(sh)
}
