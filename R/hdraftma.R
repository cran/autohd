#' @title High dimensional competing risk analysis using multivariate accelerated failure time model
#'         with mediation analysis.
#' @param m Starting column number from where high dimensional variates to be selected.
#' @param n Ending column number till where high dimensional variates to be selected.
#' @param survdur "Column/Variable name" consisting duration of survival.
#' @param event "Column/Variable name" consisting survival event.
#' @param ths A numeric between 0 to 100.
#' @param b Number of MCMC iterations to burn.
#' @param d Number of draws for the iterations.
#' @param data High dimensional data containing survival observations and high dimensional covariates.
#' @param sig Level of significance pre-determined by the user.
#' @description Given the dimension of variables and survival information including the competing risks the function
#' filters significant variables, allowing the user to fit multivariate AFT model. Further, it performs mediation
#' analysis among the significant variables and provides handful variables with their alpha.a values
#' which are mediator model exposure coefficients and beta.a coefficients.
#' @return Data frame containing the beta and alpha values of active variables among the significant variables.
#' @import survival
#' @import hdbm
#' @import schoolmath
#' @import SurvRegCensCov
#'
#' @export
#'
#' @examples
#' ##
#' hdraftma(m=8,n=100,survdur="os",event="death2",sig=0.1,ths=0.02,b=10,d=10,data=hnscc2)
#' ##
hdraftma <- function(m,n,survdur,event,sig,ths,b,d,data){
  thresh<-ths
  draws<-d
  burn<-b
  Event<-event
  siglevel<-sig
  Surv<-survdur
  siglevel<-sig
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
  hrpres <- matrix(nrow=0,ncol=3)
  hrt <- matrix(nrow=0,ncol=3)
  etrt <- matrix(nrow=0,ncol=3)
  pv <- matrix(nrow=0,ncol=1)
  colnames(pv)<-c("Pvalue")
  cnames<-c(colnames(data)[m:n])

  for(i in 1:nbatch){
    m1=sq[i]
    n1=m1+4



    wv <- WeibullReg(Surv(get(Surv),get(Event)) ~ dataA[,m1]+dataA[,m1+1]+dataA[,m1+2]+dataA[,m1+3]+dataA[,m1+4], data=dataA)

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

      wv <- WeibullReg(Surv(get(Surv),get(Event)) ~ dataA[,m2], data =dataA)
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



      wv <- WeibullReg(Surv(get(Surv),get(Event)) ~ dataA[,m2]+dataA[,m2+1], data =dataA)
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



      wv <- WeibullReg(Surv(get(Surv),get(Event)) ~ dataA[,m2]+dataA[,m2+1]+dataA[,m2+2], data =dataA)
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
    if((n-sq[nbatch+1])==3){
      m2=sq[nbatch+1]
      n2=m2+(n-sq[nbatch+1])



      wv <- WeibullReg(Surv(get(Surv),get(Event)) ~ dataA[,m2]+dataA[,m2+1]+dataA[,m2+2]+dataA[,m2+3], data =dataA)
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


      wv <- WeibullReg(Surv(get(Surv),get(Event)) ~ dataA[,m2]+dataA[,m2+1]+dataA[,m2+2]+dataA[,m2+3]+dataA[,m2+4], data =dataA)
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




  hrpresa <- matrix(nrow=0,ncol=3)
  hrt1 <- matrix(nrow=0,ncol=3)
  etrt1 <- matrix(nrow=0,ncol=3)
  pv1 <- matrix(nrow=0,ncol=1)
  colnames(pv1)<-c("Pvalue")
  cnames<-c(colnames(data)[m:n])

  for(i in 1:nbatch){
    m1=sq[i]
    n1=m1+4



    wv1 <- WeibullReg(Surv(get(Surv),statusDeath) ~ dataB[,m1]+dataB[,m1+1]+dataB[,m1+2]+dataB[,m1+3]+dataB[,m1+4], data=dataB)

    hr1 <- wv1$HR
    etr1 <- wv1$ETR
    pval1<-matrix(nrow=0,ncol=1)
    p11 <- wv1$summary$table[2,"p"][1]
    p21 <- wv1$summary$table[2:3,"p"][2]
    p31 <- wv1$summary$table[2:4,"p"][3]
    p41 <- wv1$summary$table[2:5,"p"][4]
    p51 <- wv1$summary$table[2:6,"p"][5]
    pval1 <- rbind(p11,p21,p31,p41,p51)
    rownames(pval1)<-NULL

    colnames(pval1)<-c("Pvalue")


    hrt1 <- rbind(hrt1,hr1)
    etrt1 <- rbind(etrt1,etr1)
    pv1<- rbind(pv1,pval1)

  }
  if(is.decimal(nbatch)==TRUE){
    if((n-sq[nbatch+1])==0){
      m2=sq[nbatch+1]

      wv1 <- WeibullReg(Surv(get(Surv),statusDeath) ~ dataB[,m2], data =dataB)
      hr1 <- wv1$HR
      etr1 <- wv1$ETR
      pval1<-matrix(nrow=0,ncol=1)
      p11 <- wv1$summary$table[2,"p"][1]
      pval1 <- rbind(p11)
      rownames(pval1)<-NULL

      colnames(pval1)<-c("Pvalue")

      hrt1 <- rbind(hrt1,hr1)
      etrt1 <- rbind(etrt1,etr1)
      pv1 <- rbind(pv1,pval1)
    }
  }
  if(is.decimal(nbatch)==TRUE){
    if((n-sq[nbatch+1])==1){
      m2=sq[nbatch+1]
      n2=m2+(n-sq[nbatch+1])



      wv1 <- WeibullReg(Surv(get(Surv),statusDeath) ~ dataB[,m2]+dataB[,m2+1], data =dataB)
      hr1 <- wv1$HR
      etr1 <- wv1$ETR
      pval1<-matrix(nrow=0,ncol=1)
      p11 <- wv1$summary$table[2,"p"][1]
      p21 <- wv1$summary$table[2:3,"p"][2]
      pval1 <- rbind(p11,p21)
      rownames(pval1)<-NULL

      colnames(pval1)<-c("Pvalue")

      hrt1 <- rbind(hrt1,hr1)
      etrt1 <- rbind(etrt1,etr1)
      pv1 <- rbind(pv1,pval1)


    }
  }
  if(is.decimal(nbatch)==TRUE){
    if((n-sq[nbatch+1])==2){
      m2=sq[nbatch+1]
      n2=m2+(n-sq[nbatch+1])



      wv1 <- WeibullReg(Surv(get(Surv),statusDeath) ~ dataB[,m2]+dataB[,m2+1]+dataB[,m2+2], data =dataB)
      hr1 <- wv1$HR
      etr1 <- wv1$ETR
      pval1<-matrix(nrow=0,ncol=1)
      p11 <- wv1$summary$table[2,"p"][1]
      p21 <- wv1$summary$table[2:3,"p"][2]
      p31 <- wv1$summary$table[2:4,"p"][3]
      pval1 <- rbind(p11,p21,p31)
      rownames(pval1)<-NULL

      colnames(pval1)<-c("Pvalue")

      hrt1 <- rbind(hrt1,hr1)
      etrt1 <- rbind(etrt1,etr1)
      pv1 <- rbind(pv1,pval1)

    }
  }
  if(is.decimal(nbatch)==TRUE){
    if((n-sq[nbatch+1])==3){
      m2=sq[nbatch+1]
      n2=m2+(n-sq[nbatch+1])



      wv1 <- WeibullReg(Surv(get(Surv),statusDeath) ~ dataB[,m2]+dataB[,m2+1]+dataB[,m2+2]+dataB[,m2+3], data =dataB)
      hr1 <- wv1$HR
      etr1 <- wv1$ETR
      pval1<-matrix(nrow=0,ncol=1)
      p11 <- wv1$summary$table[2,"p"][1]
      p21 <- wv1$summary$table[2:3,"p"][2]
      p31 <- wv1$summary$table[2:4,"p"][3]
      p41 <- wv1$summary$table[2:5,"p"][4]
      pval1 <- rbind(p11,p21,p31,p41)
      rownames(pval1)<-NULL

      colnames(pval1)<-c("Pvalue")

      hrt1 <- rbind(hrt1,hr1)
      etrt1 <- rbind(etrt1,etr1)
      pv1 <- rbind(pv1,pval1)
    }
  }
  if(is.decimal(nbatch)==TRUE){
    if((n-sq[nbatch+1])==4){
      m2=sq[nbatch+1]
      n2=m2+(n-sq[nbatch+1])


      wv1 <- WeibullReg(Surv(get(Surv),statusDeath) ~ dataB[,m2]+dataB[,m2+1]+dataB[,m2+2]+dataB[,m2+3]+dataB[,m2+4], data =dataB)
      hr1 <- wv1$HR
      etr1 <- wv1$ETR
      pval1<-matrix(nrow=0,ncol=1)
      p11 <- wv1$summary$table[2,"p"][1]
      p21 <- wv1$summary$table[2:3,"p"][2]
      p31 <- wv1$summary$table[2:4,"p"][3]
      p41 <- wv1$summary$table[2:5,"p"][4]
      p51 <- wv1$summary$table[2:6,"p"][5]
      pval1 <- rbind(p11,p21,p31,p41,p51)
      rownames(pval1)<-NULL

      colnames(pval1)<-c("Pvalue")

      hrt1 <- rbind(hrt1,hr1)
      etrt1 <- rbind(etrt1,etr1)
      pv1 <- rbind(pv1,pval1)

    }
  }
  est1 <- cbind(cnames,hrt1,etrt1,pv1)
  est1  <- data.frame(est1)
  rownames(est1)<-NULL

  est1 <- est1[order(est1$Pvalue),] #to filter or not??
  estf1<-est1[est1$Pvalue<=siglevel,]
  selvar1 <- c(estf1[,"cnames"])


  comselvar <- merge(estf, estf1, by ="cnames")
  if(nrow(comselvar)==0){
    print("There are no significant variables common among death due to competing risk and preogression of disease.
          Select more varaibles.")
  }
  selvarcm <- c(comselvar[,"cnames"])

  data2<-subset(data,select=selvarcm)



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

utils::globalVariables(c("WeibullReg","Status","siglevel","Pvalue"))
