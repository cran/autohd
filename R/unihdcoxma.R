#' @title High dimensional univariate cox proportional hazard model.
#' @param m Starting column number from where high dimensional variates to be selected.
#' @param n Ending column number till where high dimensional variates to be selected.
#' @param survdur "Column/Variable name" consisting duration of survival.
#' @param event "Column/Variable name" consisting survival event.
#' @param ths A numeric between 0 to 100.
#' @param b Number of MCMC iterations to burn.
#' @param d Number of draws for the iterations.
#' @param data High dimensional data containing survival observations and high dimensional covariates.
#' @param sig Level of significance pre-determined by the user.
#' @description Given the dimension of variables and survival information risks the function
#' filters significant variables, allowing the user to fit univariate COx PH model. Further, it performs mediation
#' analysis among the significant variables and provides handful variables with their alpha.a values
#' which are mediator model exposure coefficients and beta.a coefficients.
#' @return Data frame containing the beta and alpha values of active variables among the significant variables.
#' @import survival
#' @import hdbm
#' @import schoolmath
#' @export
#'
#' @examples
#' ##
#' data(hnscc)
#' unihdcoxma(m=8,n=105,survdur="os",event="death",sig=0.5,ths=0.02,b=1000,d=10,data=hnscc2)
#' ##
unihdcoxma <- function(m,n,survdur,event,sig,ths,b,d,data){
  siglevel<-sig
  thresh<-ths
  burn<-b
  draws<-d
  Event<-event
  Surv<-survdur
  nbatch<-length(m:n)/5
  sq<-seq(m,n,5)
  hrpres1 <- matrix(nrow=0,ncol=2)

  Variables<-c(colnames(data)[m:n])
  for(i in m:n){

    model1 <- coxph(Surv(get(Surv),get(Event)) ~ data[,i], data=data)
    sumr <- summary(model1)
    sumrcoeff1<-round(sumr$coefficients[,2],2)
    sumrcoeff2<-round(sumr$coefficients[,5],4)
    resdata1 <- data.frame(sumrcoeff1,sumrcoeff2)
    colnames(resdata1)<-c("HR","Pvalue")

    hrpres1<-rbind(hrpres1,resdata1)
  }
  hrpres <- data.frame(Variables,hrpres1)


  hrpres <-hrpres[order(hrpres$Pvalue),] #to filter or not??
  hrpres<-hrpres[hrpres$Pvalue<=siglevel,]
  selvar <- c(hrpres$Variables)


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

