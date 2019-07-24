#' Generalized Logit Function.
#'
#' \code{glogit} returns the value on the generalized logit curve,
#' along the covariate \code{x}, between the \code{a}
#' and \code{k} parameters, with growth \code{g} and shift \code{s}.
#'
#' \deqn{a + \frac{k-a}{1+e^{-g \times (x-s)}}}
#'
#' @param x the independent continuous variable.
#' @param a the minimum value.
#' @param k the maximum value.
#' @param g the growth parameter.
#' @param s the shift (intercept) parameter.
#' @return The value from the generalized logit function.
#' @seealso \code{\link{invglogit}}
#' @examples
#' glogit(-100:100,-70,100,0.5,10)
#' glogit(-100:100,0,1,0.5,10)
#' plot(x=-100:100,y=glogit(-100:100,-70,100,0.2,-10),type="l")
#' @export

glogit = function(x, a, k, g, s)  return (a + ((k-a)/(1+exp(-g*(x-s)))))

#' Inverse Generalized Logit Function.
#'
#' \code{invglogit} returns the value from the inverse generalized logit function,
#' along the continuous \code{y}, coming from a Generalized Logit Function, between the \code{a}
#' and \code{k} parameters, with growth \code{g} and shift \code{s}.
#'
#' \deqn{a + \frac{k-a}{1+e^{-g \times (x-s)}}}
#'
#' @param y the dependent variable / a variable coming from a Generalized Logit Function.
#' @param a the minimum value.
#' @param k the maximum value.
#' @param g the growth parameter.
#' @param s the shift (intercept) parameter.
#' @return The value from the inverse generalized logit function.
#' @seealso \code{\link{glogit}}
#' @examples
#' y = glogit(-100:100,-70,100,0.5,10)
#' x = invglogit(y,-70,100,0.5,10)
#' plot(x=x,y=y,type="l")
#' @export

invglogit = function(y, a, k, g, s)  return ((-1/g)*log(((k-a)/(y-a))-1)+s)

#' Summary for stanglogitFit objects
#'
#' @param object a stanglogitFit object
#' @return an "summary.stanglogitFit" object
#' @export

summary.stanglogitFit = function(object){
  mypaste = function(k){
    out = paste(k[1],k[2])
    if(length(k)>2)
      for(i in 3:length(k))
        out = paste(out,"*",k[i])
      return(out)
  }

  testProblem = function(k,N=nrow(as.matrix(object[[2]]))){

    k$summary = as.data.frame(k$summary)
    k$summary$" " =ifelse(k$summary[,"Rhat"]>1.1,"Rhat too high (greater than 1.1).","")
    k$summary$" " =ifelse(k$summary[,"n_eff"]<(N/10),paste(k$summary$" ","ESS too low (less than 10% of total iterations)."),k$summary$" ")
    return(k)
  }

  if(class(object)[2]!="stanglogitFit") stop("Not a valid stanglogitFit object.")

  sum.gbeta = testProblem(summary(object[[2]],pars="growth"))
  sum.sbeta = testProblem(summary(object[[2]],pars="shift"))
  sum.gsd   = testProblem(summary(object[[2]],pars="sigma_g"))
  sum.ssd   = testProblem(summary(object[[2]],pars="sigma_s"))
  sum.esd   = testProblem(summary(object[[2]],pars="sigma_e"))
  sum.rand  = testProblem(summary(object[[2]],pars="sigma_u"))
  sum.rcor  = testProblem(summary(object[[2]],pars="Cor_1"))

  dep.name = object[[9]]

  out_waic = NULL

  if (requireNamespace("loo", quietly = TRUE)) {
    out_waic = loo::waic(extract_log_lik(object[[2]]))
  }

  ans =
    list(
      "Formula Growth FE" = paste(dep.name,mypaste(object[[4]])),
      "Formula Shift FE"  = paste(dep.name,mypaste(object[[5]])),
      "RE" = as.character(object[[6]]),
      "Growth parameters" = sum.gbeta,
      "Growth variances"  = sum.gsd,
      "Shift parameters"  = sum.sbeta,
      "Shift variances"   = sum.ssd,
      "Error parameter"   = sum.esd,
      "Random parameters" = sum.rand,
      "Random cor"        = sum.rcor,
      "WAIC"              = out_waic,
      "Asymptoms formula" = as.character(object[[10]])
    )

  class(ans) = append(class(ans),"summary.stanglogitFit")

  return(ans)
}

print.summary.stanglogitFit = function(object,digits=4){

  cat("\n\nBAYESIAN GENERALISED LOGIT MULTILEVEL MODEL\n\n")

  cat("\n\n    WAIC: \n")

  print(object$`WAIC`)

  cat("\n\n   Asymptoms found by: ", object$`Asymptoms formula`,"\n\n")

  cat("\n\nFIXED EFFECTS OR POPULATION-LEVEL EFFECTS\n\n")

  cat("   Sigma: \n\n")
  print(object$`Error parameter`$summary,digits=digits)

  cat("\n\n   Formula for growth parameters: ", object$`Formula Growth FE`,"\n\n")

  cat("   Growth parameters:\n\n")
  print(rbind(object$`Growth parameters`$summary,
              object$`Growth variances`$summary),digits=digits)

  cat("\n\n   Formula for shift parameters: ", object$`Formula Shift FE`,"\n\n")

  cat("   Shift parameters:\n\n")
  print(rbind(object$`Shift parameters`$summary,
                object$`Shift variances`$summary),digits=digits)

  cat("\n\nRANDOM EFFECTS OR GROUP-LEVEL EFFECTS\n\n")

  cat("   Formula for random parameters: ", object$`RE`,"\n\n")

  cat("   Estimates of variances\n\n")
  print(object$`Random parameters`$summary,digits=digits)
  cat("\n\n   Correlations\n\n")
  tmp = object$`Random cor`$summary
  tmp2= tmp[!is.nan(tmp[,"Rhat"]),]
  tmp3= strsplit(rownames(tmp2),",")
  tmp4= matrix(unlist(lapply(tmp3,function(x){gsub(pattern="Cor_1",replace="",x,ignore.case = TRUE)})),ncol=2,byrow=TRUE)
  tmp4= gsub(pattern="[^0-9///' ]",replace="",tmp4)

  tmp5 = t(apply(tmp4,1,function(x){x[order(x)]}))
  tmp6 = apply(tmp5,1,function(x){paste(x[1],x[2],sep="-")})
  sel = NULL
  for(i in 1:(length(tmp6)-1)){
    for(j in (i+1):length(tmp6)){
      if(tmp6[i]==tmp6[j]) sel = c(sel,j)
    }
  }


  sel = c(sel,which(tmp4[,1]==tmp4[,2]))

  print(tmp2[-sel,],digits = digits)
}

print.stanglogitFit = function(x){
  print(x[[2]])
}
