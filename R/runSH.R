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
  cat(sh,'\n\n')
  system(sh)
}
