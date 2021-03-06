% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misturafina.R
\name{em_multivar_normal_mix}
\alias{em_multivar_normal_mix}
\title{algoritmo EM para misturas de normais com variancia fixa E variável}
\usage{
em_multivar_normal_mix(
  parini = rep(0, 6),
  sigmas0 = list(diag(3), diag(3)),
  prop = c(0.5, 0.5),
  amostra,
  cc0 = 1e-04,
  rep0 = 10000,
  sigmavar = FALSE
)
}
\arguments{
\item{parini}{initial mean parameters 'mu'}

\item{sigmas0}{initial variance matrix parameters 'sigmas0'}

\item{prop}{proportion of each distribution mix}

\item{amostra}{sample data}

\item{cc0}{convergence criterion}

\item{rep0}{number of replications}

\item{sigmavar}{variable/unknown variance matrix. default=FALSE}
}
\value{
lista de saidas

p: parametros mu's otimizado

rep: numero de repeticoes/iteracoes realizadas pelo algoritmo

valor_loglike: valor final da funcao loglikelihood

sigmas0: parametros sigma otimizado

pi1: parametro pi final (probabilidade a observação pertencer a mistura 1)
}
\description{
Executa o algoritmo EM para otimização numerica de uma função de máxima-verossimilhança
para misturas de normais com variancia fixa E variável e obtém a probabilidade
de uma observação pertencer a mistura 1.
}
\examples{

sample = matrix(rnorm(40*3),nrow=40,ncol=3)
mu_ini <- c(3,.5,-2,0,2,-0.5)
em_multivar_normal_mix(parini = mu_ini,amostra = sample)

mu_ini <- c(3,.5,-2,0,2,-0.5)
mu <- list(mu1=c(0,0,0),mu2=c(3,3,3))
rvarcov1 <- matrix(rWishart(1,3,diag(3)),ncol=3,nrow=3)
rvarcov2 <- matrix(rWishart(1,3,diag(3)),ncol=3,nrow=3)
sigmas0 <- list(rvarcov1,rvarcov2)
set.seed(1234)
sample <- rmultivar_normal_mix(n=40,mu=mu,sigma = sigmas0)
em_multivar_normal_mix(parini = mu_ini,amostra = sample,sigmavar=TRUE)
}
