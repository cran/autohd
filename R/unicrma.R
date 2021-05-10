#' @title High dimensional competing risk analysis by univariate accelerated failure time model with mediation analysis
#'            with Weibull distribution
#' @param m Starting column number from where high dimensional variates to be selected.
#' @param n Ending column number till where high dimensional variates to be selected.
#' @param survdur "Column/Variable name" consisting duration of survival.
#' @param event "Column/Variable name" consisting survival event.
#' @param t A numeric between 0 to 100.
#' @param sig Level of significance pre-determined by the user.
#' @param b Number of MCMC iterations to burn.
#' @param d Number of draws for the iterations.
#' @param data High dimensional data containing survival observations and high dimensional covariates.
#' @description Given the dimension of variables and survival information including the competing risks the function
#' filters significant variables,allowing the user to fit univariate AFT model. Further, it performs mediation
#' analysis among the significant variables and provides handful variables with their alpha.a values
#' which are mediator model exposure coefficients and beta.a coefficients.
#' @return Data frame containing the beta and alpha values of active variables among the significant variables.
#' @import survival
#' @import hdbm
#' @import schoolmath
#' @import SurvRegCensCov
#' @export
#'
#' @examples
#' ##
#' data(hnscc2)
#' unicrma(m=8,n=100,survdur="os",event="death2",sig=0.05,t=20,b=10,d=10,data=hnscc2)
#' ##
unicrma <- function(m,n,survdur,event,sig,t,b,d,data){
  thresh<-t
  burn<-b
  draws<-d
  siglevel<-sig
  OS<-survdur
  #Event<-event
  Event<-event
   data1 <- subset(data, get(Event) != 2)
  data2 <- subset(data, get(Event) != 1)
  if(Event!="statusDeath"){
    names(data2)[names(data2) == Event] <- "statusDeath"
  }
  data2$statusDeath[data2$statusDeath == 2] <- 1 #death1 = 2 in mydata2 is assigned "1"
  dataA <- data1
  dataB <- data2



  nbatch<-length(m:n)/5
  sq<-seq(m,n,5)
  hrt1 <- matrix(nrow=0,ncol=3)
  etrt1 <- matrix(nrow=0,ncol=3)
  pv1 <- matrix(nrow=0,ncol=1)
  colnames(pv1)<-c("Pvalue")
  varn <- c(names(data)[m:n])
  hrt2 <- matrix(nrow=0,ncol=3)
  etrt2 <- matrix(nrow=0,ncol=3)
  pv2 <- matrix(nrow=0,ncol=1)
  colnames(pv2)<-c("Pvalue")


  for(i in m:n){
    wv1 <- WeibullReg(Surv(get(OS),get(Event)) ~ dataA[,i], data=dataA)
    wv2 <- WeibullReg(Surv(get(OS),statusDeath) ~ dataB[,i], data=dataB)
    hr1 <- wv1$HR
    etr1 <- wv1$ETR
    pval1 <- wv1$summary$table[2,"p"]

    hrt1 <- rbind(hrt1,hr1)
    etrt1 <- rbind(etrt1,etr1)
    pv1 <- rbind(pv1,pval1)

    hr2 <- wv2$HR
    etr2 <- wv2$ETR
    pval2 <- wv2$summary$table[2,"p"]

    hrt2 <- rbind(hrt2,hr2)
    etrt2 <- rbind(etrt2,etr2)
    pv2 <- rbind(pv2,pval2)
  }
  est1 <- cbind(varn,hrt1,etrt1,pv1)
  est1  <- data.frame(est1)
  rownames(est1)<-NULL
  est2 <- cbind(varn,hrt2,etrt2,pv2)
  est2  <- data.frame(est2)
  rownames(est2)<-NULL



  est1 <- est1[order(est1$Pvalue),] #to filter or not??
  estf1<-est1[est1$Pvalue<=siglevel,]
  selvar1 <- c(est1[,"varn"])
  est2 <- est2[order(est2$Pvalue),] #to filter or not??
  estf2<-est2[est2$Pvalue<=siglevel,]
  selvar2 <- c(est2[,"varn"])

  comselvar <- merge(estf1, estf2, by ="varn")
  if(nrow(comselvar)==0){
    print("There are no significant variables common among death due to competing risk and preogression of disease.
          Select more varaibles.")
  }
  selvar <- c(comselvar[,"varn"])

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
utils::globalVariables(c("Status","WeibullReg"))
