#' @title High dimensional competing risk analysis by univariate accelerated failure time model
#'        with mediation analysis
#' @param m Starting column number from where high dimensional variates to be selected.
#' @param n Ending column number till where high dimensional variates to be selected.
#' @param survdur "Column/Variable name" consisting duration of survival.
#' @param event "Column/Variable name" consisting survival event.
#' @param ths A numeric between 0 to 100.
#' @param b Number of MCMC iterations to burn.
#' @param d Number of draws.
#' @param data High dimensional data containing survival observations and high dimensional covariates.
#' @description Given the dimension of variables and survival information risks the function
#' filters significant variables, allowing the user to fit univariate AFT model. Further, it performs mediation
#' analysis among the significant variables and provides handful variables with their alpha.a values
#' which are mediator model exposure coefficients and beta.a coefficients.
#' @return Data frame containing the beta and alpha values of active variables among the significant variables.
#' @import survival
#' @import hdbm
#' @import schoolmath
#' @export
#' @examples
#' ##
#' data(hnscc)
#' unihdaftma(m=8,n=80,survdur="os",event="death",ths=0.5,b=1000,d=10,data=hnscc2)
#' ##
unihdaftma <- function(m,n,survdur,event,ths,b,d,data){
  thresh<-ths
  burn<-b
  draws<-d
  Surv<-survdur
  Event<-event
  nbatch<-length(m:n)/5
  sq<-seq(m,n,5)
  hrt <- matrix(nrow=0,ncol=3)
  etrt <- matrix(nrow=0,ncol=3)
  pv <- matrix(nrow=0,ncol=1)
  colnames(pv)<-c("Pvalue")
  varn <- c(names(data)[m:n])

  for(i in m:n){
    wv <- WeibullReg(Surv(get(Surv),get(Event)) ~ data[,i], data=data)
    hr <- wv$HR
    etr <- wv$ETR
    pval <- wv$summary$table[2,"p"]

    hrt <- rbind(hrt,hr)
    etrt <- rbind(etrt,etr)
    pv <- rbind(pv,pval)
  }
  est <- cbind(hrt,etrt,pv)
  est  <- data.frame(est)
  rownames(est)<-varn



  est <- est[order(est$Pvalue),] #to filter or not??
  estf<-est[est$Pvalue<=0.05,]
  selvar <- c(rownames(estf))


  data2<-subset(data,select=selvar)

  M <- data.matrix(data2)

  #define parameters for hdbm
  Y<-M[,1] #data[,Event] #response variable
  A<-M[,2] #exposure variable taken as first and second column from the selected variable matrix

  C <- matrix(1, nrow(data2), 1)
  beta.m  <- rep(0, ncol(data2))
  alpha.a <- rep(0, ncol(data2))

  hdbm.out <- hdbm(Y,A,M, C, C, beta.m, alpha.a,
                   burnin = burn, ndraws = draws)

  active <- which(colSums(hdbm.out$r1 * hdbm.out$r3) > thresh) ######################
  Activevariables <- colnames(M)[active] ####################
  mbeta.a<-mean(hdbm.out$beta.a) ###################
  colm.beta.m<-apply(hdbm.out$beta.m, 2, mean) #######################
  colm.alpha.m<-apply(hdbm.out$alpha.a, 2, mean) ######################

  if(length(Activevariables)==0){
    print("No active variables")
    print("Number of active variables is 0")
  }
  if(length(Activevariables)!=0){
    dumean<-matrix(nrow = 0, ncol = 2)
    for(i in active){
      du.data<-data.frame(colm.beta.m[i],colm.alpha.m[i])
      dumean<-rbind(dumean,du.data)
    }
    act.var.means <- data.frame(dumean)
    colnames(act.var.means)<-c("Colmeans.beta.m","Colmeans.alpha.m")
    act.var.means <- data.frame(Activevariables,act.var.means)

    act.results<-list('Active variabels and their beta and alpha means'= act.var.means)
    return(act.results)
  }

}

utils::globalVariables(c("WeibullReg"))
