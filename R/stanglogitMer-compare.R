#' loo for stanglogitFit objects
#'
#' @param object a stanglogitFit object
#' @return an ""psis_loo" "loo" object
#' @export

loo.stanglogitFit = function(object,...){

  if (requireNamespace("loo", quietly = TRUE)) {
    rel_n_eff = loo::relative_eff(exp(loo::extract_log_lik(object[[2]])),
                             chain_id = rep(1:length(object[[2]]@stan_args),
                                            each=object[[2]]@stan_args[[1]]$iter-object[[2]]@stan_args[[1]]$warmup),...)
    ans = loo::loo(loo::extract_log_lik(object[[2]]),rel_n=rel_n_eff,...)
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
    ans = bridgesampling::bridge_sampler(object[[2]],...)
  }

  return(ans)

}

#' waic for stanglogitFit objects
#'
#' @param object a stanglogitFit object
#' @return an ""psis_loo" "waic" object
#' @export

waic.stanglogitFit = function(object, ...){

  if (requireNamespace("loo", quietly = TRUE)) {
    ans = loo::waic(loo::extract_log_lik(object[[2]]), ...)
  }

  return(ans)

}
