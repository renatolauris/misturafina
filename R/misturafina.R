# library(mixtools)

#' funcao geradora de dados aleatorios mistura de distribuicao normal multivariada
#'
#' geração de uma amostra aleatória para uma mistura de distribuição normal multivariada com médias igual a 'mu' e variancias igual a 'sigma'.
#'
#' @param n number of observations.
#' @param prop proportion of each distribution mix
#' @param mu list of vector of means for each distribution mix
#' @param sigma list of variance-covariance matrix for each distribution mix
#'
#' @return generates multivariate normal distribution random samples
#' @export
#'
#' @examples
#' rmultivar_normal_mix(n=40,mu=list(mu1=c(0,0,0),mu2=c(1,1,3)))
rmultivar_normal_mix <- function(n=40,prop=c(.5,.5),mu=list(rep(0,3),rep(0,3)),
                                 sigma=list(diag(3),diag(3))){
  # argumentos
  # n: tamanho da amostra
  # prop: proporcao de cada mistura de distribuicoes
  # mu: parametros de media das amostras de distr normal multivariada
  # sigma: parametros de matriz de covariancia das amostras de distr normal multivariada

  # amostra mistura 1
  amostra1 <- mixtools::rmvnorm(n=ceiling(n*prop[1]), mu = mu[[1]], sigma = sigma[[1]])
  # amostra mistura 2
  amostra2 <- mixtools::rmvnorm(n=ceiling(n*prop[2]), mu = mu[[2]],sigma = sigma[[2]])
  # uniao das amostras
  amostra <- rbind(amostra1,amostra2)

  # saida
  # amostra: amostra especificada
  return(amostra)
}

#' log-verossimilhanca de uma mistura de distribuição normal multivariada
#'
#' apresenta o resultado de uma função log-verossimilhanca para uma mistura de distribuição normal multivariada
#' com vetor de médias 'mu' e matriz de variacia 'sigma'.
#'
#' @param mu vector of means for each distribution mix
#' @param sigma list of variance-covariance matrix for each distribution mix
#' @param x multivariate normal distribution sample data
#' @param prop proportion of each distribution mix
#'
#' @return multivariate normal distribution log-likelihood result
#' @export
#'
#' @examples
#' sample = matrix(rnorm(40*3),nrow=40,ncol=3)
#' loglike_multivar_normal_mix(mu=c(0,0,0,3,3,3),sigma=list(diag(3),diag(3)),x=sample)
loglike_multivar_normal_mix <-  function(mu=rep(0,6),
                                         sigma=list(diag(3),
                                                    diag(3)), x, prop = c(0.5, 0.5)){
  # argumentos
  # mu: parametro de media da distr normal multivariada
  # sigma: parametro da matriz de covarianca e variancia da distr normal multivariada
  # x: dados amostrais
  # prop: proporcao de cada mistura de distribuicoes

  dim <- length(mu)
  mu=list(mu[1:(dim/2)],mu[(dim/2+1):(dim)])

  # saida
  # funcao de logverossimilhanca da distrib normal multivariada
  sum(log(prop[1]*mixtools::dmvnorm(x,mu = mu[[1]],sigma = sigma[[1]]) +
            prop[2]*mixtools::dmvnorm(x,mu = mu[[2]],sigma = sigma[[2]])))
}

#' algoritmo simulated annealing
#'
#' Executa o algoritmo simulated anneling para otimização numerica de uma função objetivo especificada
#'
#' @param parini initial parameters
#' @param x sample data
#' @param f_objt objective function to be optimized
#' @param n.iter number of iterations
#' @param temp cooling schedule function
#' @param d uniform distribution parameter for random variable 'e'
#'
#' @return
#' output list
#'
#' teta.ot: parametros teta otimizado
#'
#' tx_accept: taxa de aceitacao de novos valores de teta
#'
#' teta: matriz de valores percorridos do parametro teta
#'
#' f.eval: vetor de valores  percorridos da funcao objetivo
#'
#' @export
#'
#' @examples
#' sample = matrix(rnorm(40*3),nrow=40,ncol=3)
#' temp = function(i){1/log(1+i)}
#' d = function(i){0.3*sqrt(temp(i))}
#' mu_ini <- c(3,.5,-2,0,2,-0.5)
#' estima_anneling <- sim_annealing(parini = mu_ini,x = sample,
#' f_objt = loglike_multivar_normal_mix,n.iter = 500,temp=temp,
#' d = d)

sim_annealing <- function(parini=rep(0,6),x,f_objt,
                          n.iter=500,
                          temp = function(i){1/log(1+i)},d = function(i){0.1*sqrt(temp(i))}){
  # argumentos
  # parini: parametros iniciais
  # x: dados amostrais
  # f_objt: funcao objetivo a ser otimizada
  # n.iter: numero de iteracoes do algoritmo
  # temp: funcao cooling schedule - sequencia da temperatura (Tt)
  # d: parametro da va 'e' uniforme

  # matriz que armazena o valor dos parametros ao longo da iteracao do algoritmo
  teta = matrix(0,nrow=n.iter,ncol=length(parini))

  # armazena os parametros inicias na matriz de iteracao
  teta[1,] = parini

  # vetor dos valores da funcao objetivo avaliados no algoritmo
  # accept: vetor para verificar a taxa de aceitacao de um valor maior da funcao ao longo da iteracao
  f.eval = accept = rep(0,n.iter)
  # armazena o valor da funcao objetivo com parametros iniciais
  f.eval[1] = f_objt(mu=teta[1,],x=x)
  mu=teta[1,]

  for(i in 2:n.iter){
    # funcao cooling schedule - sequencia da temperatura (Tt)
    # temp = 1/log(1+i)
    # d: parametro da va 'e' uniforme
    # d = 0.5*sqrt(temp)
    # 'e': va com distrib uniforme (gt)
    e = runif(length(mu),-d(i),d(i))
    # 'te1': possivel novo de teta
    te1 = teta[(i-1),]+e
    # 'prob': probabilidade de percorrer um caminho novo
    # A probabilidade depende da escolha feita para gt e da temperatura Tt.
    prob = min(exp((f_objt(mu=te1,x=x)-
                      f_objt(mu=teta[(i-1),],x=x))/temp(i)),1)
    # va uniforme para testar se aceitamos o novo valor de teta ou nao
    u = runif(1)
    # matriz taxa de aceitacao de novos valor de teta ao longo da iteracao
    accept[i] = (u <= prob)
    #forma de calcular uma bernoulli com p=prob
    # se u <= prob, entao armazenamos o novo valor de teta, senao mantemos o teta inicial
    teta[i,] = (u <= prob)*te1 + (u > prob)*teta[(i-1),]
    # armazena o valor da funcao objetivo com parametro teta escolhido
    f.eval[i] =  f_objt(mu=teta[i,],x=x)
    mu=teta[i,]
  }

  # taxa de aceitacao
  # Se d for muito grande, a taxa de aceitação será muito baixa.
  tx_accept <- mean(accept[-1]) #taxa de aceitacao. regra de bolso entre 0.2 e 0.6
  #encontra o máximo
  pos <-  which.max(f.eval)
  # valor do parametro teta otimizado
  teta.ot  <-  teta[pos,]

  # lista de saidas
  # teta.ot: parametros teta otimizado
  # tx_accept: taxa de aceitacao de novos valores de teta
  # teta: matriz de valores percorridos do parametro teta
  # f.eval: vetor de valores  percorridos da funcao objetivo
  output <- list("par.opt" = teta.ot,
                 "tx_accept" = tx_accept,
                 "matriz.teta" = teta,
                 "f_objt.teta"=f.eval)
  return(output)
}

#' algoritmo EM para misturas de normais com variancia fixa E variável
#'
#'Executa o algoritmo EM para otimização numerica de uma função de máxima-verossimilhança
#'para misturas de normais com variancia fixa E variável e obtém a probabilidade
#'de uma observação pertencer a mistura 1.
#'
#' @param parini initial mean parameters 'mu'
#' @param sigmas0 initial variance matrix parameters 'sigmas0'
#' @param prop proportion of each distribution mix
#' @param amostra sample data
#' @param cc0 convergence criterion
#' @param rep0 number of replications
#' @param sigmavar variable/unknown variance matrix. default=FALSE
#'
#' @return
#' lista de saidas
#'
#' p: parametros mu's otimizado
#'
#' rep: numero de repeticoes/iteracoes realizadas pelo algoritmo
#'
#' valor_loglike: valor final da funcao loglikelihood
#'
#' sigmas0: parametros sigma otimizado
#'
#' pi1: parametro pi final (probabilidade a observação pertencer a mistura 1)
#'
#' @export
#'
#' @examples
#'
#'sample = matrix(rnorm(40*3),nrow=40,ncol=3)
#'mu_ini <- c(3,.5,-2,0,2,-0.5)
#'em_multivar_normal_mix(parini = mu_ini,amostra = sample)
#'
#'mu_ini <- c(3,.5,-2,0,2,-0.5)
#'mu <- list(mu1=c(0,0,0),mu2=c(3,3,3))
#'rvarcov1 <- matrix(rWishart(1,3,diag(3)),ncol=3,nrow=3)
#'rvarcov2 <- matrix(rWishart(1,3,diag(3)),ncol=3,nrow=3)
#'sigmas0 <- list(rvarcov1,rvarcov2)
#'set.seed(1234)
#'sample <- rmultivar_normal_mix(n=40,mu=mu,sigma = sigmas0)
#'em_multivar_normal_mix(parini = mu_ini,amostra = sample,sigmavar=TRUE)

em_multivar_normal_mix <- function(parini=rep(0,6),sigmas0=list(diag(3),diag(3)),
                                   prop=c(.5,.5),amostra,
                                   cc0=0.0001,rep0=10000,sigmavar=FALSE){

  # argumentos
  # parini: parametros de mu's iniciais
  # sigmas0: parametros da matriz de covariancia e variancia iniciais
  # prop: proporcao de cada mistura de distribuicoes
  # amostra: dados amostrais
  # cc0: criterio de convergencia
  # rep0: criterio limite de numero de repeticoes
  # sigmavar: parametro para setor se matriz de variancia e covariancia é fixa ou variavel. Padrão = Fixa

  ### Inicia EM
  # matriz dos parametros
  p=matrix(ncol=1000,nrow=length(parini))

  # valores iniciais de mu
  p0=t(parini)
  # armazena os valores iniciais de mu na matriz de parametros mu
  p[,1]=p0
  #valor inicial do criterio de convergencia
  cc=1
  #variavel numero de repeticoes
  rep=1
  #vetor de parametros encontrados no passo M
  p1=vector()

  while(cc>cc0 & rep<rep0){#loop executado enquanto nao atingir os criterios: repeticao ou convergencia
    #etapa E
    # definicao dos 'pi's
    # para os casos de matriz de covariancias variaveis ou fixas
    if(sigmavar == TRUE){
      pi1=prop[1]*mixtools::dmvnorm(amostra,mu = p0[1:3],
                        sigma = sigmas0[[1]])/(prop[1]*mixtools::dmvnorm(amostra,mu = p0[1:3],
                                                             sigma = sigmas0[[1]]) +
                                                 prop[2]*mixtools::dmvnorm(amostra,mu = p0[4:6],
                                                               sigma = sigmas0[[2]]))
    }else{
      pi1=prop[1]*mixtools::dmvnorm(amostra,mu = p0[1:3])/(prop[1]*mixtools::dmvnorm(amostra,mu = p0[1:3]) +
                                                   prop[2]*mixtools::dmvnorm(amostra,mu = p0[4:6]))
    }


    #etapa M

    # mu's estimados a partir dos pi's da etapa E
    muhat1=p1[1:3]=colSums(pi1*amostra)/sum(pi1)
    muhat2=p1[4:6]=colSums((1-pi1)*amostra)/sum((1-pi1))

    # sigmas estimados a partir dos pi's da etapa E, caso covariancia variavel
    if(sigmavar == TRUE){
      sigmahat1=(1/nrow(amostra))*sum(pi1)*t(amostra-muhat1)%*%(amostra-muhat1)
      sigmahat2=(1/nrow(amostra))*sum(1-pi1)*t(amostra-muhat2)%*%(amostra-muhat2)
    }

    # contador de repeticao do algoritmo
    rep=rep+1
    # verifica atendimento do criterio de convergencia
    # caso sigma fixo: soma da diferenca quadratica entre os parametros mu's
    # caso sigma variável: diferenca entre as funcoes objetivo (loglikelihood)
    if(sigmavar == FALSE){
      cc=sum((p0-p1)^2)}else{
        #alternativa de cc = f1-f0
        sigmashat <- list(sigmahat1,sigmahat2)
        cc=loglike_multivar_normal_mix(p1,sigmashat,amostra)-
          loglike_multivar_normal_mix(p0,sigmas0,amostra)
        # atualiza os sigmas estimados que serao utilizados para calculo dos pi's  na etapa E novamente
        sigmas0=sigmashat
      }
    # atualiza os mu's estimados que serao utilizados para calculo dos pi's na etapa E novamente
    p0=p1
    # armazena os mu's estimados na matriz de parametros mu
    p[,rep]=p0
  }

  p=p[,1:rep]

  valor_loglike=loglike_multivar_normal_mix(p0,sigmas0,amostra)

  # lista de saidas
  # p: parametros mu's otimizado
  # rep: numero de repeticoes/iteracoes realizadas pelo algoritmo
  # valor_loglike: valor final da funcao loglikelihood
  # sigmas0: parametros sigma otimizado
  # pi1: parametro pi final (probabilidade a observação pertencer a mistura 1)
  output <- list("par.opt" =   p,
                 "n repeticoes" = rep,
                 "valor_loglike" = valor_loglike,
                 "sigma.opt" = sigmas0,
                 "pi.final"=pi1
                 )
  return(output)
}

#' algoritmo de Newton Raphson caso multivariado
#'
#' @param varini vector of initial x points
#' @param f1 gradient vector for function f
#' @param f2 hessian matrix for function f
#' @param cc convergence criterion
#'
#' @return
#' vector of critical x values along iteration
#'
#' @export
#'
#' @examples
#'
#' f = function(z) {
#' x = z[1]
#' y = z[2]
#' lambda = z[3]
#' (x^4+y^2+4*x*y+lambda*(x^2+y^2-1))
#' }
#'
#' f1 = function(z){
#'  x = z[1]
#'  y = z[2]
#'  lambda = z[3]
#'  rbind(4*x^3+4*y+2*x*lambda,2*y+4*x+2*y*lambda,
#'        x^2+y^2-1)
#' }
#'
#' f2 = function(z){
#'  x = z[1]
#'  y = z[2]
#'  lambda = z[3]
#'  A = matrix(0,ncol=3,nrow=3)
#'  A[1,1] = 12*x^2+2*lambda
#'  A[1,2] = A[2,1] = 4
#'  A[1,3] = A[3,1] = 2*x
#'  A[2,2] = 2 + 2*lambda
#'  A[2,3] = A[3,2] = 2*y
#'  A[3,3] = 0
#'  return(A)
#' }
#'
#' z = x0 = rbind(3,1,0.5)
#' nraphsonmv(varini = z,f1,f2)

nraphsonmv <- function(varini,f1,f2,cc=0.0001){
  # argumentos #
  #varini = vetor inicial x
  #f1 = vetor gradiente de f
  #f2 = matriz hessiana de f
  #cc = criterio de convergencia

  conta = 0 #contagem do nro de interacoes
  x0 = varini #valor inicial de x0
  historico = c(varini,f(varini)) #para armazenar valores percorridos pelo algoritmo
  cc1= cc+1 #inicializacao da variavel teste #cc = criterio de convergencia
  while(cc1>cc){
    x1 = x0 - solve(f2(x0))%*%(f1(x0)) #solve(f2(x0)) inversa da hessiana
    cc1 = sum((x1-x0)^2)
    x0 = x1
    conta=conta+1
    historico <-  cbind(historico,c(x1,f(x1))) #registra em uma matriz todos o valores percorridos pelo algoritmo
    #print(paste("iteração",conta,"x = ",x1[1],", y=",x1[2],
    #            "e lambda=",x1[3],"com f(x)=",f(x1))) #registra as iteracoes
  }
  rownames(historico) <-  c("x","y","lambda","f(x)") # nomear as variaveis
  colnames(historico) <-  paste("i",1:(conta+1),sep="") #i indice da ordem de iteração
  varcrit <- c(historico[,conta+1]) #armazena os pontos criticos e o valor de f nesses pontos
  return(varcrit)
}
