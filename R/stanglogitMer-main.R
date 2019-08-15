#' stanglogitMer: Bayesian Generalized Logit Multilevels Models using 'Stan'
#'
#' The \strong{stanglogitMer} package provides an interface to fit Bayesian Generalized Logit Multilevel Models.
#' These models are useful when the dependent data has a sigmoid-like shape, with boundaries that are not
#' limited within the 0-1 range, as expressed by the following function:
#'
#' \deqn{a + \frac{k-a}{1+e^{-g \times (x-s)}}}
#'
#' By means of these functions, you can estimate the linear regression describing the growth parameter (\code{g}),
#' and the shift parameter (\code{s}), taking into account of a multilevel structure, also known as
#' random effects.
#'
#' @section Details:
#' The main function of \strong{stanglogitMer} is \code{\link{stanglogitMer}}, which uses formula syntax to
#' specify your model.
#'
#' @docType package
#' @name stanglogitMer
NULL


#' Fit Bayesian Generalized Logit  Multilevel Models
#'
#' \code{stanglogitMer} fits the Generalized Logit Function by means of Stan.
#'
#' @param growth.formula the formula for the fixed effects of the growth parameter.
#' @param shift.formula the formula for the fixed effects of the shift parameter.
#' @param random.formula the formula for the random effects (slopes).
#' @param grouping.random the categorial factor to group the random effects (the intercept).
#' @param covariate the x-axis covariate (distance, time, etc...) along which the dependent variable is distributed.
#' You must not include it in the fixed or random effects.
#' @param data the data frame
#' @param cores Number of cores to use when executing the chains in parallel, which
#' defaults to 1 but we recommend setting the mc.cores option to be as many processors
#' as the hardware and RAM allow (up to the number of chains).
#' For non-Windows OS in non-interactive R sessions, forking is used
#' instead of PSOCK clusters.
#' @param chains Number of Markov chains (defaults to 4).
#' @param warmup Number of total iterations per chain (including warmup; defaults to 2000).
#' @param iter A positive integer specifying number of warmup (aka burnin) iterations.
#' This also specifies the number of iterations used for stepsize adaptation, so warmup
#' samples should not be used for inference. The number of warmup should not
#' be larger than iter and the default is 4000
#' @param seed The seed for random number generation to make results reproducible.
#' If NA (the default), Stan will set the seed randomly.
#' @param control A named list of parameters to control the sampler's behavior.
#' It defaults to NULL so all the default values are
#' used. For a comprehensive overview see \strong{stan}.
#'
#' @details
#' The necessity to split the formula in three different pieces may be confusing.
#' In order to better explain how it works, let's say that we want to
#' fit the following formula, with the classic \code{lmer} syntax:
#' \code{y ~ Condition * Group + (Condition|Participant)}.
#'
#' In order to fit with \code{stanglogitMer} that model we will set:
#' \itemize{
#'     \item \code{growth.formula = ~ Condition * Group}
#'     \item \code{shift.formula = ~ Condition * Group}
#'     \item \code{random.formula = ~ Condition}
#'     \item \code{grouping = ~ Participant}
#' }
#'
#' @export

stanglogitMer = function(dependent,
                         growth.formula,
                         shift.formula,
                         random.formula,
                         grouping.random,
                         asymptoms.formula = NULL,
                         covariate,
                         data,
                         cores = 1, chains = 4, warmup = 2000, iter = 4000, seed = NA, control = list(adapt_delta=0.8)){

  if(grouping.random==~1){
    stanmodel = "
    functions {
      real glogit(real x, real a, real k, real g, real s){
        return a + ((k-a)/(1+exp(-g*(x-s))));
      }
    }
    data {
      int<lower=1> Nobs;                          //number of data points
      int<lower=1> NpredsG;                       //number of predictors for growth
      int<lower=1> NpredsS;                       //number of predictors for shift
      real y[Nobs];                               //dependent variable
      matrix[Nobs,NpredsG] XFg;                   //contrast matrix for fixed effects, independent variables
      matrix[Nobs,NpredsS] XFs;                   //contrast matrix for fixed effects, independent variables
      real a[Nobs];                               //minimum for subject
      real k[Nobs];                               //maximum for subject
      real dist[Nobs];                            //the distances of the experiment
    }
    transformed data{
      real standard_deviation = sd(y);
    }
    parameters {
      vector[NpredsG] growth;               // growth predictor parameters
      vector[NpredsS] shift;                // shift predictor parameters
      real<lower=0> sigma_e;
      real<lower=0> sigma_s;
      real<lower=0> sigma_g;
    }
    transformed parameters{
      real mu[Nobs];
      real g[Nobs];
      real s[Nobs];
      real minimum[Nobs];
      real maximum[Nobs];

      for(i in 1:Nobs){
        g[i] = dot_product(growth,XFg[i,]);
        s[i] = dot_product(shift,XFs[i,]);
        minimum[i] = a[i];
        maximum[i] = k[i];
        mu[i] = glogit(dist[i], a[i], k[i], g[i], s[i]);
      }
    }
    model {
      //priors
      target += cauchy_lpdf(sigma_g | 0,  standard_deviation);
      target += cauchy_lpdf(sigma_s | 0,  standard_deviation);
      target += cauchy_lpdf(sigma_e | 0,  sqrt(standard_deviation));

      target += cauchy_lpdf(growth | 0, sigma_g);
      target += cauchy_lpdf(shift  | 0, sigma_s);

      //likelihood
      target += cauchy_lpdf(y | mu, sigma_e);
    }
    generated quantities {
      vector[Nobs] log_lik;
      real y_rep[Nobs];

      for (n in 1:Nobs){
        log_lik[n] = cauchy_lpdf(y[n]| mu[n], sigma_e);
        y_rep[n] = cauchy_rng(mu[n], sigma_e);
      }
    }
    "

  }else if (random.formula==~1){##only random intercept
    stanmodel = "
    functions {
      real glogit(real x, real a, real k, real g, real s){
        return a + ((k-a)/(1+exp(-g*(x-s))));
      }
    }
    data {
      int<lower=1> Nobs;                          //number of data points
      int<lower=1> NpredsG;                       //number of predictors for growth
      int<lower=1> NpredsS;                       //number of predictors for shift
      int<lower=1> Nrands;                        //number of random coefficients
      int<lower=1> Ngroup;                        //number of levels of the grouping factor
      real y[Nobs];                               //dependent variable
      matrix[Nobs,NpredsG] XFg;                   //contrast matrix for fixed effects, independent variables
      matrix[Nobs,NpredsS] XFs;                   //contrast matrix for fixed effects, independent variables
      matrix[Nobs,Nrands] XR;                     //contrast matrix for random effects
      int<lower=1, upper=Ngroup> grouping[Nobs];  //grouping factor
      real a[Nobs];                               //minimum for subject
      real k[Nobs];                               //maximum for subject
      real dist[Nobs];                            //the distances of the experiment
    }
    transformed data{
      real standard_deviation = sd(y);
    }
    parameters {
      vector[NpredsG] growth;               // growth predictor parameters
      vector[NpredsS] shift;                // shift predictor parameters
      real<lower=0> sigma_u;                // random effects sd
      real<lower=0> sigma_e;
      real<lower=0> sigma_s;
      real<lower=0> sigma_g;
      vector[Ngroup] u;                     // random effects
    }
    transformed parameters{
      real mu[Nobs];
      real g[Nobs];
      real s[Nobs];
      real minimum[Nobs];
      real maximum[Nobs];

      for(i in 1:Nobs){
        g[i] = dot_product(growth,XFg[i,]);
        s[i] = dot_product(shift,XFs[i,]);
        minimum[i] = a[i];
        maximum[i] = k[i];
        mu[i] = glogit(dist[i], a[i], k[i], g[i], s[i]) + u[grouping[i]];
      }
    }
    model {
      //priors
      target += normal_lpdf(u | 0, sigma_u);

      target += cauchy_lpdf(sigma_g | 0,  standard_deviation);
      target += cauchy_lpdf(sigma_s | 0,  standard_deviation);
      target += cauchy_lpdf(sigma_e | 0,  sqrt(standard_deviation));

      target += cauchy_lpdf(growth | 0, sigma_g);
      target += cauchy_lpdf(shift  | 0, sigma_s);

      //likelihood
      target += cauchy_lpdf(y | mu, sigma_e);
    }
    generated quantities {
      vector[Nobs] log_lik;
      real y_rep[Nobs];

      for (n in 1:Nobs){
        log_lik[n] = cauchy_lpdf(y[n]| mu[n], sigma_e);
        y_rep[n] = cauchy_rng(mu[n], sigma_e);
      }
    }
    "

  }else{
    stanmodel = "
    functions {
      real glogit(real x, real a, real k, real g, real s){
      return a + ((k-a)/(1+exp(-g*(x-s))));
      }
    }
    data {
      int<lower=1> Nobs;                          //number of data points
      int<lower=1> NpredsG;                       //number of predictors for growth
      int<lower=1> NpredsS;                       //number of predictors for shift
      int<lower=1> Nrands;                        //number of random coefficients
      int<lower=1> Ngroup;                        //number of levels of the grouping factor
      real y[Nobs];                               //dependent variable
      matrix[Nobs,NpredsG] XFg;                   //contrast matrix for fixed effects, independent variables
      matrix[Nobs,NpredsS] XFs;                   //contrast matrix for fixed effects, independent variables
      matrix[Nobs,Nrands] XR;                     //contrast matrix for random effects
      int<lower=1, upper=Ngroup> grouping[Nobs];  //grouping factor
      real a[Nobs];                               //minimum for subject
      real k[Nobs];                               //maximum for subject
      real dist[Nobs];                            //the distances of the experiment
    }
    transformed data{
      real standard_deviation = sd(y);
    }
    parameters {
      vector[NpredsG] growth;               // growth predictor estimates
      vector[NpredsS] shift;                // shift predictor estimates
      vector<lower=0>[Nrands] sigma_u;      // random effects sd
      real<lower=0> sigma_e;
      real<lower=0> sigma_s;
      real<lower=0> sigma_g;
      cholesky_factor_corr[Nrands] L_Omega;
      matrix[Nrands,Ngroup] z_u;
    }
    transformed parameters{
      matrix[Nrands,Ngroup] u;
      real mu[Nobs];
      real g[Nobs];
      real s[Nobs];
      real minimum[Nobs];
      real maximum[Nobs];
      u = (diag_pre_multiply(sigma_u, L_Omega) * z_u); //random effects

      for(i in 1:Nobs){
        g[i] = dot_product(growth,XFg[i,]);
        s[i] = dot_product(shift,XFs[i,]);
        minimum[i] = a[i];
        maximum[i] = k[i];
        mu[i] = glogit(dist[i], a[i], k[i], g[i], s[i]) + dot_product(u[,grouping[i]],XR[i,]);
      }
    }
    model {
      //priors
      target += lkj_corr_cholesky_lpdf(L_Omega | 1);
      target += normal_lpdf(to_vector(z_u) | 0, 1);

      target += cauchy_lpdf(sigma_g | 0,  standard_deviation);
      target += cauchy_lpdf(sigma_s | 0,  standard_deviation);
      target += cauchy_lpdf(sigma_e | 0,  sqrt(standard_deviation));

      target += cauchy_lpdf(growth | 0,  sigma_g);
      target += cauchy_lpdf(shift  | 0,  sigma_s);

      //likelihood
      target += cauchy_lpdf(y | mu,  sigma_e);
    }
    generated quantities {
      vector[Nobs] log_lik;
      real y_rep[Nobs];
      corr_matrix[Nrands] Cor_1 = multiply_lower_tri_self_transpose(L_Omega);

      for (n in 1:Nobs){
        log_lik[n] = cauchy_lpdf(y[n]| mu[n],  sigma_e);
        y_rep[n] = cauchy_rng(mu[n],  sigma_e);
      }
    }
    "
  }

  if(is.null(asymptoms.formula)) asymptoms.formula = grouping.random

  grouping.string = as.character(grouping.random)[2]
  gg = 1
  if(grouping.string!="1") gg = as.numeric(data[,grouping.string])

  minmax = .findMaxMin(asymptoms.formula,data,dependent)

  y    = dependent
  XFg  = model.matrix(growth.formula,data=data)
  XFs  = model.matrix(shift.formula,data=data)
  XR   = model.matrix(random.formula,data=data)

  data.list = list(
    Nobs    = nrow(data),
    NpredsG = ncol(XFg),
    NpredsS = ncol(XFs),
    Nrands  = ncol(XR),
    Ngroup  = ifelse(grouping.string!="1",length(levels(data[,grouping.string])),1),

    y       = y,
    XFg     = XFg,
    XFs     = XFs,
    XR      = XR,

    grouping= gg,

    dist    = covariate,

    a       = minmax$min,
    k       = minmax$max
  )

  mdl = stan(model_code=stanmodel, data=data.list,
             warmup=warmup, iter=iter, control=control,
             chains=chains, cores = cores,
             seed=seed)

  out = list(data,mdl,dependent,
             growth.formula,
             shift.formula,
             random.formula,
             grouping.random,
             covariate,
             deparse(substitute(dependent)),
             asymptoms.formula,
             XFg,XFs,XR)

  class(out) = append(class(out),"stanglogitFit")

  return(out)
}
