#' plot_sigmoid for stanglogitFit objects
#'
#' \code{plot_sigmoid()} plots the resulting generalised logit curves from the posterior distributions of the model.
#' For each combination of population-level factors for the Growth and Shift parameters, 9 curves will be plotted.
#' One curve will come from the median value of the posterior distributions of the Growth and Shift parameters,
#' the other eight parameters will be obtained by the combination of the boundaries of the Credible
#' Interval for the Growth and Shift parameters. The dimension of the Credible Interval can be
#' set by means of the \code{CI} parameter.
#'
#' A sigmoid function can be seen as a cumulative gaussian function. In order to plot the gaussian curves,
#' use \code{plot_gaussian()}.
#'
#' @param object a stanglogitFit object
#' @param CI the dimension of the credible interval taken into account. Default 0.95.
#' @return a ggplot2 object or a base-R plot
#' @seealso \code{\link{plot_gaussian}}
#' @export

plot_sigmoid = function(object,CI=0.95){

  gs = as.data.frame(extract(object[[2]],c("g","s")))

  gS = apply(gs,2,function(x){quantile(x,probs=c((1-CI)/2,0.5,1-((1-CI)/2)))})
  gS = as.data.frame(cbind(t(gS[,1:(ncol(gS)/2)]),
                           t(gS[,(ncol(gS)/2+1):(ncol(gS))]))
                     )
  colnames(gS) = c(paste0("g_",colnames(gS)[1:3]),paste0("s_",colnames(gS)[4:6]))

  dat = object[[1]]

  forG= gsub("~","",x=as.character(object[[4]]))
  forS= gsub("~","",x=as.character(object[[5]]))

  distance = object[[8]]

  maximum= unlist(as.data.frame(extract(object[[2]],c("maximum")))[1,])
  minimum= unlist(as.data.frame(extract(object[[2]],c("minimum")))[1,])

  ag = unlist(strsplit(forG,"\\+"))[-1]
  ag = trimws(ag[!grepl(":",ag)])

  as = unlist(strsplit(forS,"\\+"))[-1]
  as = trimws(as[!grepl(":",as)])

  selG = colnames(dat)%in%ag
  selS = colnames(dat)%in%as

  tmp0 = data.frame(
    g1 = gS[,1],
    g2 = gS[,2],
    g3 = gS[,3],
    s1 = gS[,4],
    s2 = gS[,5],
    s3 = gS[,6],
    factgrowth= interaction(dat[,selG]),
    factshift = interaction(dat[,selS]),
    distance=distance,
    maximum = maximum,
    minimum = minimum
  )

  tmp0$divide = paste("Growth:",as.character(tmp0$factgrowth),"\nShift:",as.character(tmp0$factshift))

  # tmp0$y0 = glogit(x=tmp0$distance*10,a=(tmp0$minimum),k=(tmp0$maximum),g=tmp0$g2, s=tmp0$s2)
  # tmp0$y1 = glogit(x=tmp0$distance*10,a=(tmp0$minimum),k=(tmp0$maximum),g=tmp0$g1, s=tmp0$s1)
  # tmp0$y2 = glogit(x=tmp0$distance*10,a=(tmp0$minimum),k=(tmp0$maximum),g=tmp0$g3, s=tmp0$s3)
  # tmp0$y3 = glogit(x=tmp0$distance*10,a=(tmp0$minimum),k=(tmp0$maximum),g=tmp0$g1, s=tmp0$s2)
  # tmp0$y4 = glogit(x=tmp0$distance*10,a=(tmp0$minimum),k=(tmp0$maximum),g=tmp0$g1, s=tmp0$s3)
  # tmp0$y5 = glogit(x=tmp0$distance*10,a=(tmp0$minimum),k=(tmp0$maximum),g=tmp0$g2, s=tmp0$s1)
  # tmp0$y6 = glogit(x=tmp0$distance*10,a=(tmp0$minimum),k=(tmp0$maximum),g=tmp0$g2, s=tmp0$s3)
  # tmp0$y7 = glogit(x=tmp0$distance*10,a=(tmp0$minimum),k=(tmp0$maximum),g=tmp0$g3, s=tmp0$s1)
  # tmp0$y8 = glogit(x=tmp0$distance*10,a=(tmp0$minimum),k=(tmp0$maximum),g=tmp0$g3, s=tmp0$s2)

  tmp0$y0 = glogit(x=tmp0$distance,a=(tmp0$minimum),k=(tmp0$maximum),g=tmp0$g2, s=tmp0$s2)
  tmp0$y1 = glogit(x=tmp0$distance,a=(tmp0$minimum),k=(tmp0$maximum),g=tmp0$g1, s=tmp0$s1)
  tmp0$y2 = glogit(x=tmp0$distance,a=(tmp0$minimum),k=(tmp0$maximum),g=tmp0$g3, s=tmp0$s3)
  tmp0$y3 = glogit(x=tmp0$distance,a=(tmp0$minimum),k=(tmp0$maximum),g=tmp0$g1, s=tmp0$s2)
  tmp0$y4 = glogit(x=tmp0$distance,a=(tmp0$minimum),k=(tmp0$maximum),g=tmp0$g1, s=tmp0$s3)
  tmp0$y5 = glogit(x=tmp0$distance,a=(tmp0$minimum),k=(tmp0$maximum),g=tmp0$g2, s=tmp0$s1)
  tmp0$y6 = glogit(x=tmp0$distance,a=(tmp0$minimum),k=(tmp0$maximum),g=tmp0$g2, s=tmp0$s3)
  tmp0$y7 = glogit(x=tmp0$distance,a=(tmp0$minimum),k=(tmp0$maximum),g=tmp0$g3, s=tmp0$s1)
  tmp0$y8 = glogit(x=tmp0$distance,a=(tmp0$minimum),k=(tmp0$maximum),g=tmp0$g3, s=tmp0$s2)

  N = round(sqrt(length(unique(tmp0$divide))))

  tmp0 = tmp0[order(tmp0$distance),]

  ans=NULL

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ans = ggplot(tmp0)+
      geom_smooth(aes(y=y1,x=distance),se=FALSE,col="grey81",method = 'loess')+
      geom_smooth(aes(y=y2,x=distance),se=FALSE,col="grey81",method = 'loess')+
      geom_smooth(aes(y=y3,x=distance),se=FALSE,col="grey81",method = 'loess')+
      geom_smooth(aes(y=y4,x=distance),se=FALSE,col="grey81",method = 'loess')+
      geom_smooth(aes(y=y5,x=distance),se=FALSE,col="grey81",method = 'loess')+
      geom_smooth(aes(y=y6,x=distance),se=FALSE,col="grey81",method = 'loess')+
      geom_smooth(aes(y=y7,x=distance),se=FALSE,col="grey81",method = 'loess')+
      geom_smooth(aes(y=y8,x=distance),se=FALSE,col="grey81",method = 'loess')+
      geom_smooth(aes(y=y0,x=distance),se=FALSE,col="steelblue1",method = 'loess')+
      facet_wrap(~divide,ncol=N)+xlab("")+ylab("")
  }else{
    i = 1
    old.par = par(no.readonly = TRUE)
    par(mfrow=c(N,N))
    afterFirst = FALSE
    for(ww in unique(tmp0$divide)){
      plot(x=tmp0$distance[tmp0$divide==ww],y=tmp0$y0[tmp0$divide==ww],col=i,type="l",xlab="",ylab="",main=ww)
      points(x=tmp0$distance[tmp0$divide==ww],y=tmp0$y1[tmp0$divide==ww],col=i,type="l",lty=3)
      points(x=tmp0$distance[tmp0$divide==ww],y=tmp0$y2[tmp0$divide==ww],col=i,type="l",lty=3)
      points(x=tmp0$distance[tmp0$divide==ww],y=tmp0$y3[tmp0$divide==ww],col=i,type="l",lty=3)
      points(x=tmp0$distance[tmp0$divide==ww],y=tmp0$y4[tmp0$divide==ww],col=i,type="l",lty=3)
      points(x=tmp0$distance[tmp0$divide==ww],y=tmp0$y5[tmp0$divide==ww],col=i,type="l",lty=3)
      points(x=tmp0$distance[tmp0$divide==ww],y=tmp0$y6[tmp0$divide==ww],col=i,type="l",lty=3)
      points(x=tmp0$distance[tmp0$divide==ww],y=tmp0$y7[tmp0$divide==ww],col=i,type="l",lty=3)
      points(x=tmp0$distance[tmp0$divide==ww],y=tmp0$y8[tmp0$divide==ww],col=i,type="l",lty=3)

      i = i+1
    }

    par(old.par)
  }

  return(ans)

}

#' plot_gaussian for stanglogitFit objects
#'
#' \code{plot_gaussian()} plots the resulting generalised logit curves from the posterior distributions of the model,
#' as gaussian curves, with mean the Shift parameter, and standard deviation the inverse of the
#' Growth parameter.
#' For each combination of population-level factors for the Growth and Shift parameters, 9 curves will be plotted.
#' One curve will come from the median value of the posterior distributions of the Growth and Shift parameters,
#' the other eight parameters will be obtained by the combination of the boundaries of the Credible
#' Interval for the Growth and Shift parameters. The dimension of the Credible Interval can be
#' set by means of the \code{CI} parameter.
#'
#' @param object a stanglogitFit object
#' @param CI the dimension of the credible interval taken into account. Default 0.95.
#' @return a ggplot2 object or a base-R plot
#' @seealso \code{\link{plot_sigmoid}}
#' @export

plot_gaussian = function(object,CI=0.95){

  gs = as.data.frame(extract(object[[2]],c("g","s")))

  gS = apply(gs,2,function(x){quantile(x,probs=c((1-CI)/2,0.5,1-((1-CI)/2)))})
  gS = as.data.frame(cbind(t(gS[,1:(ncol(gS)/2)]),
                           t(gS[,(ncol(gS)/2+1):(ncol(gS))]))
  )
  colnames(gS) = c(paste0("g_",colnames(gS)[1:3]),paste0("s_",colnames(gS)[4:6]))

  dat = object[[1]]

  forG= gsub("~","",x=as.character(object[[4]]))
  forS= gsub("~","",x=as.character(object[[5]]))

  distance = object[[8]]

  maximum= unlist(as.data.frame(extract(object[[2]],c("maximum")))[1,])
  minimum= unlist(as.data.frame(extract(object[[2]],c("minimum")))[1,])

  ag = unlist(strsplit(forG,"\\+"))[-1]
  ag = trimws(ag[!grepl(":",ag)])

  as = unlist(strsplit(forS,"\\+"))[-1]
  as = trimws(as[!grepl(":",as)])

  selG = colnames(dat)%in%ag
  selS = colnames(dat)%in%as

  tmp0 = data.frame(
    g1 = gS[,1],
    g2 = gS[,2],
    g3 = gS[,3],
    s1 = gS[,4],
    s2 = gS[,5],
    s3 = gS[,6],
    factgrowth= interaction(dat[,selG]),
    factshift = interaction(dat[,selS]),
    distance=distance,
    maximum = maximum,
    minimum = minimum
  )

  tmp0$divide = paste("Growth:",as.character(tmp0$factgrowth),"\nShift:",as.character(tmp0$factshift))

  # tmp0$y0 = dnorm(x=tmp0$distance*10,sd=1/tmp0$g2, mean=tmp0$s2)
  # tmp0$y1 = dnorm(x=tmp0$distance*10,sd=1/tmp0$g1, mean=tmp0$s1)
  # tmp0$y2 = dnorm(x=tmp0$distance*10,sd=1/tmp0$g3, mean=tmp0$s3)
  # tmp0$y3 = dnorm(x=tmp0$distance*10,sd=1/tmp0$g1, mean=tmp0$s2)
  # tmp0$y4 = dnorm(x=tmp0$distance*10,sd=1/tmp0$g1, mean=tmp0$s3)
  # tmp0$y5 = dnorm(x=tmp0$distance*10,sd=1/tmp0$g2, mean=tmp0$s1)
  # tmp0$y6 = dnorm(x=tmp0$distance*10,sd=1/tmp0$g2, mean=tmp0$s3)
  # tmp0$y7 = dnorm(x=tmp0$distance*10,sd=1/tmp0$g3, mean=tmp0$s1)
  # tmp0$y8 = dnorm(x=tmp0$distance*10,sd=1/tmp0$g3, mean=tmp0$s2)

  tmp0$y0 = dnorm(x=tmp0$distance,sd=1/tmp0$g2, mean=tmp0$s2)
  tmp0$y1 = dnorm(x=tmp0$distance,sd=1/tmp0$g1, mean=tmp0$s1)
  tmp0$y2 = dnorm(x=tmp0$distance,sd=1/tmp0$g3, mean=tmp0$s3)
  tmp0$y3 = dnorm(x=tmp0$distance,sd=1/tmp0$g1, mean=tmp0$s2)
  tmp0$y4 = dnorm(x=tmp0$distance,sd=1/tmp0$g1, mean=tmp0$s3)
  tmp0$y5 = dnorm(x=tmp0$distance,sd=1/tmp0$g2, mean=tmp0$s1)
  tmp0$y6 = dnorm(x=tmp0$distance,sd=1/tmp0$g2, mean=tmp0$s3)
  tmp0$y7 = dnorm(x=tmp0$distance,sd=1/tmp0$g3, mean=tmp0$s1)
  tmp0$y8 = dnorm(x=tmp0$distance,sd=1/tmp0$g3, mean=tmp0$s2)

  N = round(sqrt(length(unique(tmp0$divide))))

  tmp0 = tmp0[order(tmp0$distance),]

  ans=NULL

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ans = ggplot(tmp0)+
      geom_smooth(aes(y=y1,x=distance),se=FALSE,col="grey81",method = 'loess')+
      geom_smooth(aes(y=y2,x=distance),se=FALSE,col="grey81",method = 'loess')+
      geom_smooth(aes(y=y3,x=distance),se=FALSE,col="grey81",method = 'loess')+
      geom_smooth(aes(y=y4,x=distance),se=FALSE,col="grey81",method = 'loess')+
      geom_smooth(aes(y=y5,x=distance),se=FALSE,col="grey81",method = 'loess')+
      geom_smooth(aes(y=y6,x=distance),se=FALSE,col="grey81",method = 'loess')+
      geom_smooth(aes(y=y7,x=distance),se=FALSE,col="grey81",method = 'loess')+
      geom_smooth(aes(y=y8,x=distance),se=FALSE,col="grey81",method = 'loess')+
      geom_smooth(aes(y=y0,x=distance),se=FALSE,col="steelblue1",method = 'loess')+
      geom_vline(aes(xintercept=y1),col="steelblue1")+
      facet_wrap(~divide,ncol=N)+xlab("")+ylab("")
  }else{
    i = 1
    afterFirst = FALSE
    for(ww in unique(tmp2$label)){
      if(afterFirst){
        points(x=tmp2$x[tmp2$label==ww],y=tmp2$y[tmp2$label==ww],col=i,type="l")
        points(x=tmp2$x[tmp2$label==ww],y=tmp2$ym[tmp2$label==ww],col=i,type="l",lty=2)
        points(x=tmp2$x[tmp2$label==ww],y=tmp2$yM[tmp2$label==ww],col=i,type="l",lty=3)
        abline(v=tmp2$s[tmp2$label==ww],col=i)
        abline(v=tmp2$sm[tmp2$label==ww],col=i,lty=2)
        abline(v=tmp2$sM[tmp2$label==ww],col=i,lty=3)
      }else{
        plot(x=tmp2$x[tmp2$label==ww],y=tmp2$y[tmp2$label==ww],col=i,type="l",xlab="",ylab="",ylim=c(0,0.2))
        points(x=tmp2$x[tmp2$label==ww],y=tmp2$ym[tmp2$label==ww],col=i,type="l",lty=2)
        points(x=tmp2$x[tmp2$label==ww],y=tmp2$yM[tmp2$label==ww],col=i,type="l",lty=3)
        abline(v=tmp2$s[tmp2$label==ww],col=i)
        abline(v=tmp2$sm[tmp2$label==ww],col=i,lty=2)
        abline(v=tmp2$sM[tmp2$label==ww],col=i,lty=3)
        afterFirst = TRUE
      }

      i = i+1
    }

    legend(min(tmp2$x),0.19,col=1:(i-1),unique(tmp2$label),lty=1)
  }

  return(ans)

}


#' plot for stanglogitFit objects
#'
#' @param object a stanglogitFit object
#' @param pars an optional vector of strings, with possible values:
#'     \describe{
#'         \item{all} all the parameters
#'         \item{growth} only the growth parameters
#'         \item{shift} only the shift parameters
#'         \item{random} only the random variances
#'         \item{cor_random} only the correlations between random variables
#'         \item{sigma} only the overall sigma of the model
#'     }
#' @param type to select the type of plot:
#'         \item{combo} a density and a trace plot
#'         \item{area} a density plot
#'         \item{hist} a density histogram
#'         \item{intervals} an 95% intervals plot
#'         \item{pairs} a pairs plot (only with \code{bayesplot})
#'         \item{trace} a trace/catterpillar plot
#'         \item{violin} a violin plot (only with \code{bayesplot})
#'
#' @param ... other bayesplot or ggplot2 arguments
#'
#' @return a plot, a ggplot2 object, or a bayesplot object
#' @export

plot.stanglogitFit = function(object,pars=NULL,type="combo",...){

  ans=list()

  if (requireNamespace("bayesplot", quietly = TRUE)) {
    if(is.null(pars)||("growth"%in%pars)||("all"%in%pars)){
      if(type == "combo"){
        ans[["growth"]]     = bayesplot::mcmc_combo(As.mcmc.list(object[[2]]),regex_pars = c("growth","sigma_g"),...)
      } else if(type == "area"){
        ans[["growth"]]     = bayesplot::mcmc_areas(As.mcmc.list(object[[2]]),regex_pars = c("growth","sigma_g"),...)
      } else if(type == "intervals"){
        ans[["growth"]]     = bayesplot::mcmc_intervals(As.mcmc.list(object[[2]]),regex_pars = c("growth","sigma_g"),...)
      } else if(type == "hist"){
        ans[["growth"]]     = bayesplot::mcmc_hist(As.mcmc.list(object[[2]]),regex_pars = c("growth","sigma_g"),...)
      } else if(type == "pairs"){
        ans[["growth"]]     = bayesplot::mcmc_pairs(As.mcmc.list(object[[2]]),regex_pars = c("growth","sigma_g"),...)
      } else if(type == "trace"){
        ans[["growth"]]     = bayesplot::mcmc_trace(As.mcmc.list(object[[2]]),regex_pars = c("growth","sigma_g"),...)
      } else if(type == "violin"){
        ans[["growth"]]     = bayesplot::mcmc_violin(As.mcmc.list(object[[2]]),regex_pars = c("growth","sigma_g"),...)
      }
    }

    if(is.null(pars)||("shift"%in%pars)||("all"%in%pars)){
      if(type == "combo"){
        ans[["shift"]]     = bayesplot::mcmc_combo(As.mcmc.list(object[[2]]),regex_pars = c("shift","sigma_s"),...)
      } else if(type == "area"){
        ans[["shift"]]     = bayesplot::mcmc_areas(As.mcmc.list(object[[2]]),regex_pars = c("shift","sigma_s"),...)
      } else if(type == "intervals"){
        ans[["shift"]]     = bayesplot::mcmc_intervals(As.mcmc.list(object[[2]]),regex_pars = c("shift","sigma_s"),...)
      } else if(type == "hist"){
        ans[["shift"]]     = bayesplot::mcmc_hist(As.mcmc.list(object[[2]]),regex_pars = c("shift","sigma_s"),...)
      } else if(type == "pairs"){
        ans[["shift"]]     = bayesplot::mcmc_pairs(As.mcmc.list(object[[2]]),regex_pars = c("shift","sigma_s"),...)
      } else if(type == "trace"){
        ans[["shift"]]     = bayesplot::mcmc_trace(As.mcmc.list(object[[2]]),regex_pars = c("shift","sigma_s"),...)
      } else if(type == "violin"){
        ans[["shift"]]     = bayesplot::mcmc_violin(As.mcmc.list(object[[2]]),regex_pars = c("shift","sigma_s"),...)
      }
    }

    if(is.null(pars)||("random"%in%pars)||("all"%in%pars)){
      if(type == "combo"){
        ans[["random"]]     = bayesplot::mcmc_combo(As.mcmc.list(object[[2]]),regex_pars = "sigma_u",...)
      } else if(type == "area"){
        ans[["random"]]     = bayesplot::mcmc_areas(As.mcmc.list(object[[2]]),regex_pars = "sigma_u",...)
      } else if(type == "intervals"){
        ans[["random"]]     = bayesplot::mcmc_intervals(As.mcmc.list(object[[2]]),regex_pars = "sigma_u",...)
      } else if(type == "hist"){
        ans[["random"]]     = bayesplot::mcmc_hist(As.mcmc.list(object[[2]]),regex_pars = "sigma_u",...)
      } else if(type == "pairs"){
        ans[["random"]]     = bayesplot::mcmc_pairs(As.mcmc.list(object[[2]]),regex_pars = "sigma_u",...)
      } else if(type == "trace"){
        ans[["random"]]     = bayesplot::mcmc_trace(As.mcmc.list(object[[2]]),regex_pars = "sigma_u",...)
      } else if(type == "violin"){
        ans[["random"]]     = bayesplot::mcmc_violin(As.mcmc.list(object[[2]]),regex_pars = "sigma_u",...)
      }
    }

    if(is.null(pars)||("cor_random"%in%pars)||("all"%in%pars)){
      if(type == "combo"){
        ans[["cor_random"]]     = bayesplot::mcmc_combo(As.mcmc.list(object[[2]]),regex_pars = "Cor_1",...)
      } else if(type == "area"){
        ans[["cor_random"]]     = bayesplot::mcmc_areas(As.mcmc.list(object[[2]]),regex_pars = "Cor_1",...)
      } else if(type == "intervals"){
        ans[["cor_random"]]     = bayesplot::mcmc_intervals(As.mcmc.list(object[[2]]),regex_pars = "Cor_1",...)
      } else if(type == "hist"){
        ans[["cor_random"]]     = bayesplot::mcmc_hist(As.mcmc.list(object[[2]]),regex_pars = "Cor_1",...)
      } else if(type == "pairs"){
        ans[["cor_random"]]     = bayesplot::mcmc_pairs(As.mcmc.list(object[[2]]),regex_pars = "Cor_1",...)
      } else if(type == "trace"){
        ans[["cor_random"]]     = bayesplot::mcmc_trace(As.mcmc.list(object[[2]]),regex_pars = "Cor_1",...)
      } else if(type == "violin"){
        ans[["cor_random"]]     = bayesplot::mcmc_violin(As.mcmc.list(object[[2]]),regex_pars = "Cor_1",...)
      }
    }

    if(is.null(pars)||("sigma"%in%pars)||("all"%in%pars)){
      if(type == "combo"){
        ans[["sigma_e"]]     = bayesplot::mcmc_combo(As.mcmc.list(object[[2]]),regex_pars = "sigma_e",...)
      } else if(type == "area"){
        ans[["sigma_e"]]     = bayesplot::mcmc_areas(As.mcmc.list(object[[2]]),regex_pars = "sigma_e",...)
      } else if(type == "intervals"){
        ans[["sigma_e"]]     = bayesplot::mcmc_intervals(As.mcmc.list(object[[2]]),regex_pars = "sigma_e",...)
      } else if(type == "hist"){
        ans[["sigma_e"]]     = bayesplot::mcmc_hist(As.mcmc.list(object[[2]]),regex_pars = "sigma_e",...)
      } else if(type == "pairs"){
        ans[["sigma_e"]]     = bayesplot::mcmc_pairs(As.mcmc.list(object[[2]]),regex_pars = "sigma_e",...)
      } else if(type == "trace"){
        ans[["sigma_e"]]     = bayesplot::mcmc_trace(As.mcmc.list(object[[2]]),regex_pars = "sigma_e",...)
      } else if(type == "violin"){
        ans[["sigma_e"]]     = bayesplot::mcmc_violin(As.mcmc.list(object[[2]]),regex_pars = "sigma_e",...)
      }
    }

  }else if (requireNamespace("ggplot2", quietly = TRUE)&&requireNamespace("gridExtra", quietly = TRUE)) {
    tmp0 = As.mcmc.list(object[[2]])
    tmp1 = as.data.frame(as.matrix(object[[2]]))
    tmp1$ID = 1:nrow(tmp1)

    GimmeGraph = function(what,tmp0,tmp1,tmp2,type,...){
      what = c(what,"ID")
      sel = NULL
      N = length(tmp0)
      for(ww in what){
        sel = c(sel,which(grepl(ww,colnames(tmp1))))
      }

      tmp2 = reshape(tmp1[,sel],direction="long",idvar="ID",varying=list(1:(length(sel)-1)),times=colnames(tmp1[,sel])[-length(sel)],v.names="y")
      tmp2$chain = as.factor(rep(1:N,each=nrow(tmp1)/N))
      tmp2$n     = rep(1:(nrow(tmp1)/N),N)

      if(type=="combo"){
        out = gridExtra::arrangeGrob(ggplot2::ggplot(tmp2,aes(x=y),...)+
                                 geom_density(fill="blue")+
                                 facet_grid(time~., scales="free_y")+ylab("")+xlab(""),
                               ggplot2::ggplot(tmp2,aes(y=y,x=n,colour=chain),...)+
                                 geom_path()+
                                 facet_grid(time~., scales="free_y")+ylab("")+xlab(""))
      }else if(type=="area"){
        out = ggplot2::ggplot(tmp2,aes(x=y),...)+
                                       geom_density(fill="blue")+
                                       facet_grid(time~., scales="free_y")+ylab("")+xlab("")
      }else if(type=="hist"){
        out = ggplot2::ggplot(tmp2,aes(x=y),...)+
          geom_histogram(fill="blue")+
          facet_grid(time~., scales="free_y")+ylab("")+xlab("")
      }else if(type=="trace"){
        out = ggplot2::ggplot(tmp2,aes(y=y,x=n,colour=chain),...)+
          geom_path()+
          facet_grid(time~., scales="free_y")+ylab("")+xlab("")
      }else if(type=="intervals"){
        out = ggplot2::ggplot(tmp2,aes(y=y,x=0),...)+
          stat_summary(geom="pointrange",fun.y=mean,fun.ymin=function(x){quantile(x,probs=0.025)},fun.ymax=function(x){quantile(x,probs=0.975)})+
          facet_grid(time~., scales="free_y")+ylab("")+xlab("")+coord_flip()
      }


      return(out)
    }

    if(is.null(pars)||("growth"%in%pars)||("all"%in%pars))
      ans[["growth"]]     = GimmeGraph(what=c("growth","sigma_g"),tmp0,tmp1,tmp2,type=type)
    if(is.null(pars)||("shift"%in%pars)||("all"%in%pars))
      ans[["shift"]]      = GimmeGraph(what=c("shift","sigma_s"),tmp0,tmp1,tmp2,type=type)
    if(is.null(pars)||("random"%in%pars)||("all"%in%pars))
      ans[["random"]]     = GimmeGraph(what=c("sigma_u"),tmp0,tmp1,tmp2,type=type)
    if(is.null(pars)||("cor_random"%in%pars)||("all"%in%pars))
      ans[["cor_random"]] = GimmeGraph(what=c("Cor_1"),tmp0,tmp1,tmp2,type=type)
    if(is.null(pars)||("sigma"%in%pars)||("all"%in%pars))
      ans[["sigma_e"]] = GimmeGraph(what=c("sigma_e"),tmp0,tmp1,tmp2,type=type)
  }else{
    tmp0 = As.mcmc.list(object[[2]])
    tmp1 = as.data.frame(as.matrix(tmp0))
    tmp1$ID = 1:nrow(tmp1)

    GimmeGraph = function(what,tmp0,tmp1,tmp2){

      what = c(what,"ID")
      sel = NULL
      N = length(tmp0)
      for(ww in what){
        sel = c(sel,which(grepl(ww,colnames(tmp1))))
      }

      tmp2 = reshape(tmp1[,sel],direction="long",idvar="ID",varying=list(1:(length(sel)-1)),times=colnames(tmp1[,sel])[-length(sel)],v.names="y")
      tmp2$chain = as.factor(rep(1:N,each=nrow(tmp1)/N))
      tmp2$n     = rep(1:(nrow(tmp1)/N),N)

      old.par = par(no.readonly = TRUE)

      if(type=="combo"){
        par(mfrow=c(length(unique(tmp2$time)),2))

        for(ww in unique(tmp2$time)){
          plot(density(tmp2$y[tmp2$time==ww]),main=ww)

          plot(x=tmp2$n[tmp2$time==ww],y=tmp2$y[tmp2$time==ww],col=tmp2$chain[tmp2$time==ww],type="l",xlab="",ylab="")
        }
      }else if(type=="area"){
        par(mfrow=c(length(unique(tmp2$time)),1))

        for(ww in unique(tmp2$time)){
          plot(density(tmp2$y[tmp2$time==ww]),main=ww)
        }
      }else if(type=="hist"){
        par(mfrow=c(length(unique(tmp2$time)),1))

        for(ww in unique(tmp2$time)){
          plot(density(tmp2$y[tmp2$time==ww]),main=ww,type="h")
        }
      }else if(type=="trace"){
        par(mfrow=c(length(unique(tmp2$time)),1))

        for(ww in unique(tmp2$time)){
          plot(x=tmp2$n[tmp2$time==ww],y=tmp2$y[tmp2$time==ww],col=tmp2$chain[tmp2$time==ww],type="l",xlab="",ylab="",main=ww)
        }
      }else if(type=="intervals"){
        par(mfrow=c(length(unique(tmp2$time)),1))

        ll = range(tmp2$y)

        for(ww in unique(tmp2$time)){
          plot(x=quantile(tmp2$y[tmp2$time==ww],probs=c(0.025,0.975)),y=c(0,0),type="l",xlab="",ylab="",main=ww,xlim=ll)
          points(x=mean(tmp2$y[tmp2$time==ww]),y=0,type="p",xlab="",ylab="",main=ww)
        }
      }

      par(old.par)


      return()
    }

    if(is.null(pars)||("growth"%in%pars)||("all"%in%pars))
      ans[["growth"]]     = GimmeGraph(what=c("growth","sigma_g"),tmp0,tmp1,tmp2,type=type)
    if(is.null(pars)||("shift"%in%pars)||("all"%in%pars))
      ans[["shift"]]      = GimmeGraph(what=c("shift","sigma_s"),tmp0,tmp1,tmp2,type=type)
    if(is.null(pars)||("random"%in%pars)||("all"%in%pars))
      ans[["random"]]     = GimmeGraph(what=c("sigma_u"),tmp0,tmp1,tmp2,type=type)
    if(is.null(pars)||("cor_random"%in%pars)||("all"%in%pars))
      ans[["cor_random"]] = GimmeGraph(what=c("Cor_1"),tmp0,tmp1,tmp2,type=type)
    if(is.null(pars)||("sigma"%in%pars)||("all"%in%pars))
      ans[["sigma_e"]] = GimmeGraph(what=c("sigma_e"),tmp0,tmp1,tmp2,type=type)
  }

  return(ans)

}

#' pp_check for stanglogitFit objects
#'
#' \code{plot_sigmoid()} plots the resulting generalised logit curves from the posterior distributions of the model.
#' For each combination of population-level factors for the Growth and Shift parameters, 9 curves will be plotted.
#' One curve will come from the median value of the posterior distributions of the Growth and Shift parameters,
#' the other eight parameters will be obtained by the combination of the boundaries of the Credible
#' Interval for the Growth and Shift parameters. The dimension of the Credible Interval can be
#' set by means of the \code{CI} parameter.
#'
#' A sigmoid function can be seen as a cumulative gaussian function. In order to plot the gaussian curves,
#' use \code{plot_gaussian()}.
#'
#' @param object a stanglogitFit object
#' @param CI the dimension of the credible interval taken into account. Default 0.95.
#' @return a ggplot2 object or a base-R plot
#' @seealso \code{\link{plot_gaussian}}
#' @export

pp_check = function(object,type="dens"){

  y_rep = extract(object[[2]], pars = "y_rep")

  y     = object[[3]]

  ans=list()

  if (requireNamespace("bayesplot", quietly = TRUE)) {
    if(type=="hist"){
       ans = ppc_hist(y, y_rep[[1]][1:8, ])+xlim(c(-sd(y)*3,sd(y)*3))
    }else if(type=="mean"){
       ans = ppc_stat(y = y, yrep = y_rep[[1]], stat = "mean")+xlim(c(-sd(y)*3,sd(y)*3))
    }else{
       ans = ppc_dens_overlay(y, y_rep[[1]][1:200, ])+xlim(c(-sd(y)*3,sd(y)*3))
    }
  }else if (requireNamespace("ggplot2", quietly = TRUE)){
    if(type=="hist"){
      nn  = ncol(y_rep[[1]])
      tmp = data.frame(
        y = c(y,y_rep[[1]][1,],y_rep[[1]][2,],y_rep[[1]][3,],y_rep[[1]][4,],y_rep[[1]][5,],y_rep[[1]][6,],y_rep[[1]][7,],y_rep[[1]][8,]),
        f = rep(c("y","yrep"),c(length(y),nn*8)),
        c = rep(0:8,c(length(y),nn,nn,nn,nn,nn,nn,nn,nn))
      )
      nnn = length(unique(tmp$c))

      ans = ggplot(tmp,aes(x=y,fill=f))+
        geom_histogram()+
        facet_wrap(~c,ncol = round(sqrt(nnn)))+
        xlim(c(-sd(y)*3,sd(y)*3))+ylab("")+theme(legend.title = element_blank())
    }else if(type=="mean"){
      ym = mean(y)
      yr = data.frame(apply(y_rep[[1]],1,mean))

      ans = ggplot(tmp,aes(x=y))+
        geom_histogram()+
        xlim(c(-sd(y)*3,sd(y)*3))+ylab("")+theme(legend.title = element_blank())+
        geom_vline(xintercept = ym,size=2)+xlab("")
    }else{
      nn  = ncol(y_rep[[1]])
      tmp = data.frame(
        y = c(y,y_rep[[1]][1,],y_rep[[1]][2,],y_rep[[1]][3,],y_rep[[1]][4,],y_rep[[1]][5,],y_rep[[1]][6,],y_rep[[1]][7,],y_rep[[1]][8,]),
        f = rep(c("y","yrep"),c(length(y),nn*8)),
        c = rep(0:8,c(length(y),nn,nn,nn,nn,nn,nn,nn,nn))
      )

      ans = ggplot(tmp,aes(x=y,colour=f,group=c))+
        geom_density()+
        xlim(c(-sd(y)*3,sd(y)*3))+ylab("")+theme(legend.title = element_blank())
    }
  }else{
    if(type=="hist"){
      nn  = ncol(y_rep[[1]])
      tmp = data.frame(
        y = c(y,y_rep[[1]][1,],y_rep[[1]][2,],y_rep[[1]][3,],y_rep[[1]][4,],y_rep[[1]][5,],y_rep[[1]][6,],y_rep[[1]][7,],y_rep[[1]][8,]),
        f = rep(c("y","yrep"),c(length(y),nn*8)),
        c = rep(0:8,c(length(y),nn,nn,nn,nn,nn,nn,nn,nn))
      )
      nnn = length(unique(tmp$c))

      old.par = par(no.readonly = TRUE)

      par(mfrow=c(round(sqrt(nnn)),round(sqrt(nnn))))

      for(cc in unique(tmp$c)){
        hist(x=tmp$y[tmp$c==cc],main=tmp$f[tmp$c==cc][1],xlab="",xlim=c(-sd(y)*3,sd(y)*3),breaks=round(diff(range(tmp$y[tmp$c==cc]))/50))
      }

      par(old.par)

      ans = ggplot(tmp,aes(x=y,fill=f))+
        geom_histogram()+
        facet_wrap(~c,ncol = round(sqrt(nnn)))+
        xlim(c(-sd(y)*3,sd(y)*3))+ylab("")+theme(legend.title = element_blank())
    }else if(type=="mean"){
      ym = mean(y)
      yr = data.frame(apply(y_rep[[1]],1,mean))

      ans = ggplot(tmp,aes(x=y))+
        geom_histogram()+
        xlim(c(-sd(y)*3,sd(y)*3))+ylab("")+theme(legend.title = element_blank())+
        geom_vline(xintercept = ym,size=2)+xlab("")
    }else{
      nn  = ncol(y_rep[[1]])
      tmp = data.frame(
        y = c(y,y_rep[[1]][1,],y_rep[[1]][2,],y_rep[[1]][3,],y_rep[[1]][4,],y_rep[[1]][5,],y_rep[[1]][6,],y_rep[[1]][7,],y_rep[[1]][8,]),
        f = rep(c("y","yrep"),c(length(y),nn*8)),
        c = rep(0:8,c(length(y),nn,nn,nn,nn,nn,nn,nn,nn))
      )

      ans = ggplot(tmp,aes(x=y,colour=f,group=c))+
        geom_density()+
        xlim(c(-sd(y)*3,sd(y)*3))+ylab("")+theme(legend.title = element_blank())
    }

  }

  return(ans)
}
