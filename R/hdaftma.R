#' @title High dimensional multivariate accelerated failure time model with bayesian
#'        mediation analysis
#' @param m Starting column number from where high dimensional variates to be selected.
#' @param n Ending column number till where high dimensional variates to be selected.
#' @param survdur "Column/Variable name" consisting duration of survival.
#' @param event "Column/Variable name" consisting survival event.
#' @param ths A numeric between 0 to 100.
#' @param b Number of MCMC iterations to burn.
#' @param sig Level of significance.
#' @param d Number of draws for the iterations.
#' @param data High dimensional data containing survival observations and high dimensional covariates.
#' @description Given the dimension of variables and survival information the function filters significant variables
#' by fitting AFT model. Further, it performs mediation analysis among the signifiant
#' variables and provides handful variables with their alpha.a values which are mediator model exposure coefficients
#' and beta.a coefficients.
#' @return Data frame containing the beta and alpha values of active variables among the significant variables.
#' @import survival
#' @import hdbm
#' @import schoolmath
#' @import SurvRegCensCov
#' @export
#'
#' @examples
#' ##
#' hdaftma(m=8,n=80,survdur="os",event="death",sig=0.05,ths=0.02,b=10,d=10,data=hnscc2)
#' ##
hdaftma <- function(m,n,survdur,event,ths,sig,b,d,data){
  thresh<-ths
  burn<-b
  draws<-d
  Surv<-survdur
  Event<-event
  siglevel<-sig
  nbatch<-length(m:n)/5
  sq<-seq(m,n,5)
  hrpres <- matrix(nrow=0,ncol=3)
  hrt <- matrix(nrow=0,ncol=3)
  etrt <- matrix(nrow=0,ncol=3)
  pv <- matrix(nrow=0,ncol=1)
  colnames(pv)<-c("Pvalue")
  cnames<-c(colnames(data)[m:n])

  for(i in 1:nbatch){
    m1=sq[i]
    n1=m1+4
    wv <- WeibullReg(Surv(get(Surv),get(Event)) ~ data[,m1]+data[,m1+1]+data[,m1+2]+data[,m1+3]+data[,m1+4], data=data)
    hr <- wv$HR
    etr <- wv$ETR
    pval<-matrix(nrow=0,ncol=1)
    p1 <- wv$summary$table[2,"p"][1]
    p2 <- wv$summary$table[2:3,"p"][2]
    p3 <- wv$summary$table[2:4,"p"][3]
    p4 <- wv$summary$table[2:5,"p"][4]
    p5 <- wv$summary$table[2:6,"p"][5]
    pval <- rbind(p1,p2,p3,p4,p5)
    rownames(pval)<-NULL

    colnames(pval)<-c("Pvalue")


    hrt <- rbind(hrt,hr)
    etrt <- rbind(etrt,etr)
    pv <- rbind(pv,pval)

  }
  if(is.decimal(nbatch)==TRUE){
    if((n-sq[nbatch+1])==0){
      m2=sq[nbatch+1]

      wv <- WeibullReg(Surv(get(Surv),get(Event)) ~ data[,m2], data =data)
      hr <- wv$HR
      etr <- wv$ETR
      pval<-matrix(nrow=0,ncol=1)
      p1 <- wv$summary$table[2,"p"][1]
      pval <- rbind(p1)
      rownames(pval)<-NULL

      colnames(pval)<-c("Pvalue")

      hrt <- rbind(hrt,hr)
      etrt <- rbind(etrt,etr)
      pv <- rbind(pv,pval)
    }
  }
  if(is.decimal(nbatch)==TRUE){
    if((n-sq[nbatch+1])==1){
      m2=sq[nbatch+1]
      n2=m2+(n-sq[nbatch+1])



      wv <- WeibullReg(Surv(get(Surv),get(Event)) ~ data[,m2]+data[,m2+1], data =data)
      hr <- wv$HR
      etr <- wv$ETR
      pval<-matrix(nrow=0,ncol=1)
      p1 <- wv$summary$table[2,"p"][1]
      p2 <- wv$summary$table[2:3,"p"][2]
      pval <- rbind(p1,p2)
      rownames(pval)<-NULL

      colnames(pval)<-c("Pvalue")

      hrt <- rbind(hrt,hr)
      etrt <- rbind(etrt,etr)
      pv <- rbind(pv,pval)


    }
  }
  if(is.decimal(nbatch)==TRUE){
    if((n-sq[nbatch+1])==2){
      m2=sq[nbatch+1]
      n2=m2+(n-sq[nbatch+1])



      wv <- WeibullReg(Surv(get(Surv),get(Event)) ~ data[,m2]+data[,m2+1]+data[,m2+2], data =data)
      hr <- wv$HR
      etr <- wv$ETR
      pval<-matrix(nrow=0,ncol=1)
      p1 <- wv$summary$table[2,"p"][1]
      p2 <- wv$summary$table[2:3,"p"][2]
      p3 <- wv$summary$table[2:4,"p"][3]
      pval <- rbind(p1,p2,p3)
      rownames(pval)<-NULL

      colnames(pval)<-c("Pvalue")

      hrt <- rbind(hrt,hr)
      etrt <- rbind(etrt,etr)
      pv <- rbind(pv,pval)

    }
  }
  if(is.decimal(nbatch)==TRUE){
    if((n-sq[nbatch+1])==4){
      m2=sq[nbatch+1]
      n2=m2+(n-sq[nbatch+1])



      wv <- WeibullReg(Surv(get(Surv),get(Event)) ~ data[,m2]+data[,m2+1]+data[,m2+2]+data[,m2+3], data =data)
      hr <- wv$HR
      etr <- wv$ETR
      pval<-matrix(nrow=0,ncol=1)
      p1 <- wv$summary$table[2,"p"][1]
      p2 <- wv$summary$table[2:3,"p"][2]
      p3 <- wv$summary$table[2:4,"p"][3]
      p4 <- wv$summary$table[2:5,"p"][4]
      pval <- rbind(p1,p2,p3,p4)
      rownames(pval)<-NULL

      colnames(pval)<-c("Pvalue")

      hrt <- rbind(hrt,hr)
      etrt <- rbind(etrt,etr)
      pv <- rbind(pv,pval)
    }
  }
  if(is.decimal(nbatch)==TRUE){
    if((n-sq[nbatch+1])==4){
      m2=sq[nbatch+1]
      n2=m2+(n-sq[nbatch+1])


      wv <- WeibullReg(Surv(get(Surv),get(Event)) ~ data[,m2]+data[,m2+1]+data[,m2+2]+data[,m2+3]+data[,m2+4], data =data)
      hr <- wv$HR
      etr <- wv$ETR
      pval<-matrix(nrow=0,ncol=1)
      p1 <- wv$summary$table[2,"p"][1]
      p2 <- wv$summary$table[2:3,"p"][2]
      p3 <- wv$summary$table[2:4,"p"][3]
      p4 <- wv$summary$table[2:5,"p"][4]
      p5 <- wv$summary$table[2:6,"p"][5]
      pval <- rbind(p1,p2,p3,p4,p5)
      rownames(pval)<-NULL

      colnames(pval)<-c("Pvalue")

      hrt <- rbind(hrt,hr)
      etrt <- rbind(etrt,etr)
      pv <- rbind(pv,pval)

    }
  }

  est <- cbind(cnames,hrt,etrt,pv)
  est  <- data.frame(est)
  rownames(est)<-NULL

  est <- est[order(est$Pvalue),] #to filter or not??
  estf<-est[est$Pvalue<=siglevel,]
  selvar <- c(estf[,"cnames"])


  data2<-subset(data,select=selvar)

  M <- data.matrix(data2)

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


