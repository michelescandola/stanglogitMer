#' loo for stanglogitFit objects
#'
#' @param object a stanglogitFit object
#' @return an ""psis_loo" "loo" object
#' @export

loo.stanglogitFit = function(object,...){

  if (requireNamespace("loo", quietly = TRUE)) {
    ans = loo::loo(extract_log_lik(object[[2]]),...)
  }

  return(ans)

}


#' bridge_sampler for stanglogitFit objects
#'
#' @param object a stanglogitFit object
#' @return an ""psis_loo" "loo" object
#' @export

bridge_sampler.stanglogitFit = function(object, ...){

  if (requireNamespace("bridgesampling", quietly = TRUE)) {
    ans = bridgesampling::bridge_sampler(object[[2]], ...)
  }

  return(ans)

}

#' waic for stanglogitFit objects
#'
#' @param object a stanglogitFit object
#' @return an ""psis_loo" "waic" object
#' @export

waic.stanglogitFit = function(object, ...){

  if (requireNamespace("bridgesampling", quietly = TRUE)) {
    ans = bridgesampling::waic(object[[2]], ...)
  }

  return(ans)

}
