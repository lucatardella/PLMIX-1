#' Gibbs sampling for a Bayesian mixture of Plackett-Luce models
#' 
#' Perform Gibbs sampling simulation for a Bayesian mixture of Plackett-Luce models fitted to partial orderings. Differently from the \code{gibbsPLMIX} function it contains a \strong{normalization step} for the support parameters \strong{which should be proven to be equivalent} to the case where the normalization step is skipped and the support parameter space is unidentified up to a proportionality constant.
#' 
#' The size \eqn{L} of the final MCMC sample is equal to \code{n_iter}-\code{n_burn}.
#' 
#' @param pi_inv An object of class \code{top_ordering}, collecting the numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings, or an object that can be coerced with \code{\link{as.top_ordering}}.
#' @param K Number of possible items. 
#' @param G Number of mixture components.
#' @param init List of named objects with initialization values: \code{z} is a numeric \eqn{N}\eqn{\times}{x}\eqn{G} matrix of binary mixture component memberships; \code{p} is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific support parameters. If starting values are not supplied (\code{NULL}), they are randomly generated with a uniform distribution. Default is \code{NULL}.
#' @param n_iter Total number of MCMC iterations.
#' @param n_burn Number of initial burn-in drawings removed from the returned MCMC sample.
#' @param hyper List of named objects with hyperparameter values for the conjugate prior specification: \code{shape0} is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of shape hyperparameters; \code{rate0} is a numeric vector of \eqn{G} rate hyperparameters; \code{alpha0} is a numeric vector of \eqn{G} Dirichlet hyperparameters. Default is vague prior setting.
#' @param centered_start Logical: whether a random start whose support parameters and weights should be centered around the observed relative frequency that each item has been ranked top. Default is \code{FALSE}. Ignored when \code{init} is not \code{NULL}.
#'
#' @return A list of S3 class \code{gsPLMIX} with named elements:
#' 
#'  \item{\code{W}}{ Numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with MCMC samples of the mixture weights.}
#'  \item{\code{P}}{ Numeric \eqn{L}\eqn{\times}{x}\eqn{(G*K)} matrix with MCMC samples of the component-specific support parameters.}
#'  \item{\code{log_lik}}{ Numeric vector of \eqn{L} posterior log-likelihood values.}
#'  \item{\code{deviance}}{ Numeric vector of \eqn{L} posterior deviance values (\eqn{-2 * }\code{log_lik}).}
#'  \item{\code{objective}}{ Numeric vector of \eqn{L} objective function values (that is the kernel of the log-posterior distribution).}
#'  \item{\code{call}}{ The matched call.}
#' 
#' @references 
#' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#' 
#' @author Cristina Mollica and Luca Tardella
#' 
#' @examples
#' 
#' data(d_carconf)
#' GIBBS <- gibbsPLMIX_with_norm(pi_inv=d_carconf, K=ncol(d_carconf), G=3, n_iter=30, n_burn=10)
#' str(GIBBS)
#' GIBBS$P
#' GIBBS$W
#' 
#' @export 


gibbsPLMIX_with_norm <- function(pi_inv,K,G,
                                 init=list(z=NULL,p=NULL),
                                 n_iter=1000,
                                 n_burn=500,
                                 hyper=list(shape0=matrix(1,nrow=G,ncol=K),rate0=rep(0.001,G),alpha0=rep(1,G)),
                                 centered_start=FALSE){
  
  cl=match.call()
  
  if(class(pi_inv)[1]!="top_ordering"){
    if(class(pi_inv)[1]=="RankData"){
      pi_inv=as.top_ordering(data=pi_inv)
    }
    if(class(pi_inv)[1]=="rankings"){
      pi_inv=as.top_ordering(data=pi_inv)
    }
    if(class(pi_inv)[1]=="matrix" | class(pi_inv)[1]=="data.frame"){
      pi_inv=as.top_ordering(data=pi_inv,format_input="ordering",aggr=FALSE)
    }
  }
  pi_inv <- fill_single_entries(data=pi_inv)
  N <- nrow(pi_inv)
  n_rank <-  howmanyranked(pi_inv)
  rho <- matrix(1:K,nrow=G,ncol=K,byrow=TRUE)
  
  if(is.null(init$z)){
    z <- binary_group_ind(class=sample(1:G,size=N,replace=TRUE),G=G)
  }else{
    z <- init$z
  }
  
  
  omega <- colMeans(z)
  
  if(is.null(init$p)){
    if(centered_start){
      
      print("CENTERED START !!")
      
      # omega <- rdirichlet(1,rep(1,G))
      mle1comp <- matrix(prop.table(table(factor(pi_inv[,1],levels=1:K))),nrow=1)
      p <- random_start(mlesupp=mle1comp, givenweights=omega)
      # p <- p/rowSums(p)
      
    }else{
      
      print("COMPLETELY RANDOM (uniform support, rescaled) START")                    
      
      p <- matrix(rgamma(n=G*K,shape=1,rate=1),nrow=G,ncol=K)
      
    }
    
  }else{
    
    p <- init$p
  }
  
  shape0 <- hyper$shape0
  rate0 <- hyper$rate0
  alpha0 <- hyper$alpha0
  
  u_bin <- umat(pi_inv=pi_inv)
  
  
  log_lik <- c(loglikPLMIX(p=p,ref_order=rho,weights=omega,pi_inv=pi_inv),
               rep(NA,n_iter))
  
  log_prior <- c(log(ddirichlet(omega,alpha0))+sum(dgamma(p,shape=shape0,rate=rate0,log=TRUE)),
                 rep(NA,n_iter))
  
  objective <- log_lik+log_prior
  
  Pi <- array(NA,dim=c(G,K,n_iter+1))
  Pi[,,1] <- p
  Zeta <- z
  Omega <- matrix(NA,nrow=n_iter+1,ncol=G)
  Omega[1,] <- omega
  
  
  for(l in 1:n_iter){
    #    print(l)
    
    if(l%%500==0){
      print(paste("GIBBS iteration",l))
    }
    
    Omega[l+1,] <- rdirichlet(n=1,alpha=alpha0+colSums(Zeta))
    
    temprate <- CompRateYpartial(p=adrop(Pi[,,l,drop=FALSE],3),pi_inv=pi_inv,ref_order=rho,z=Zeta,n_rank=n_rank)
    Ypsilon <- SimYpsilon(rate=temprate,n_rank=n_rank)
    
    Pi[,,l+1] <- matrix(rgamma(n=G*K,shape=shape0+gammamat(u_bin=u_bin,z_hat=Zeta),
                               rate <- CompRateP(pi_inv=pi_inv, Y=Ypsilon, z=Zeta, u_bin=u_bin, n_rank=n_rank, rate0=rate0)),nrow=G,ncol=K)
    tempPi <- adrop(Pi[,,l+1,drop=FALSE],drop=3)
    Pi[,,l+1] <- tempPi/rowSums(tempPi)
    
    Zeta <- binary_group_ind(apply(CompProbZpartial(p=adrop(Pi[,,l+1,drop=FALSE],3),pi_inv=pi_inv,Y=Ypsilon, u_bin=u_bin,n_rank,omega=Omega[l+1,]),1,FUN=sample,x=1:G,replace=TRUE,size=1),G=G)
    
    log_lik[l+1] <- loglikPLMIX(p=adrop(Pi[,,l+1,drop=FALSE],3),ref_order=rho,weights=Omega[l+1,],
                                pi_inv=pi_inv)
    
    log_prior[l+1] <- log(ddirichlet(Omega[l+1,],alpha0))+sum(dgamma(adrop(Pi[,,l+1,drop=FALSE],3),shape=shape0,rate=rate0,log=TRUE))
    
    objective[l+1] <- log_lik[l+1]+log_prior[l+1]
    
  }
  
  log_lik <- log_lik[-c(1:(n_burn+1))]
  
  objective <- objective[-c(1:(n_burn+1))]
  
  Omega <- Omega[-c(1:(n_burn+1)),,drop=FALSE]
  colnames(Omega) <- paste0("w_",1:G)
  
  
  Pi <- array(apply(Pi,3,FUN=function(x)x/rowSums(x)),c(G,K,n_iter+1))	
  
  Pi=t(apply(Pi,3,c))[-c(1:(n_burn+1)),]
  colnames(Pi) <- paste0("p_",rep(1:G,K),rep(1:K,each=G))
  
  out=list(W=Omega,P=Pi,log_lik=log_lik,deviance=-2*log_lik,objective=objective,call=cl)
  
  class(out)="gsPLMIX"
  
  return(out)
  
}
