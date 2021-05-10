#' @title High dimensional missing data imputation and performing mediation analysis with multivariate
#'        cox proportional hazard model. It works in a multivariate setup.
#' @param m Starting column number from where high dimensional variates to be selected.
#' @param n Ending column number till where high dimensional variates to be selected.
#' @param survdur "Column/Variable name" consisting duration of survival.
#' @param event "Column/Variable name" consisting survival event.
#' @param time "Column/Variable name" consisting time of repeated observations.
#' @param sig Level of significance pre-determined by the user.
#' @param ths A numeric between 0 to 100.
#' @param b Number of MCMC iterations to burn.
#' @param d Number of draws for the iterations.
#' @param data High dimensional data containing survival observations with multiple covariates.
#' @description Given the dimension of variables and survival information the function
#' performs imputations using missForest function and filters significant variables,
#' allowing the user to fit multivariate CoxPH model with 5 variables. Further, it performs mediation
#' analysis among the significant variables and provides handful variables with their alpha.a values
#' which are mediator model exposure coefficients and beta.a coefficients.
#'
#' @return Data frame containing the beta and alpha values of active variables among the significant variables.
#' @import survival
#' @import hdbm
#' @import schoolmath
#' @import missForest
#' @export
#'
#' @examples
#' ##
#' \dontrun{
#' imphdcox(m=11,n=25,survdur="OS",event="event",time="Visit",sig=0.5,ths=20,b=10,d=10,data=srdata)
#' ##
#' }
imphdcox <- function(m,n,survdur,event,time,sig,ths,b,d,data){
  siglevel<-sig
  thresh<-ths
  burn<-b
  draws<-d
  Event<-event
  Surv<-survdur
  data <- subset(data,select = c(get(Surv),get(Event),get(time),m:n))
  m=4
  n=ncol(data)
  data.imp <- missForest(data) #imputed tmc longit
  data <- data.frame(data.imp$ximp)

  data12 <- data
  t9<-c(unique(data[,time]))

  tlast9 <- t9[length(t9)]
  data <- subset(data, get(time) == tlast9)

  nbatch<-length(m:n)/5
  sq<-seq(m,n,5)
  hrpres <- matrix(nrow=0,ncol=3)

  for(i in 1:nbatch){
    m1=sq[i]
    n1=m1+4
    cnames<-c(colnames(data)[m1:n1])
    model1 <- coxph(Surv(get(Surv),get(Event)) ~ data[,m1]+data[,m1+1]+data[,m1+2]+data[,m1+3]+data[,m1+4], data=data)
    sumr <- summary(model1)
    sumrcoeff<-sumr$coefficients[,c(2,5)]
    resdata1 <- data.frame(cnames,sumrcoeff)
    colnames(resdata1)<-c("Variables","HR","Pvalue")
    rownames(resdata1)<- NULL
    hrpres<-rbind(hrpres,resdata1)
  }
  if(is.decimal(nbatch)==TRUE){
    if((n-sq[nbatch+1])==0){
      m2=sq[nbatch+1]
      cnames<-c(colnames(data)[m2])
      model1 <- coxph(Surv(get(Surv),get(Event)) ~ data[,m2], data =data)
      sumr <- summary(model1)
      sumrcoeff<-sumr$coefficients[,c(2,5)]
      resdata1 <- data.frame(cnames,sumrcoeff)
      colnames(resdata1)<-c("Variables","HR","Pvalue")
      rownames(resdata1)<- NULL
      hrpres<-rbind(hrpres,resdata1)
    }
  }
  if(is.decimal(nbatch)==TRUE){
    if((n-sq[nbatch+1])==1){
      m2=sq[nbatch+1]
      n2=m2+(n-sq[nbatch+1])
      cnames<-c(colnames(data)[m2:n2])
      model1 <- coxph(Surv(get(Surv),get(Event)) ~ data[,m2]+data[,m2+1], data =data)
      sumr <- summary(model1)
      sumrcoeff<-sumr$coefficients[,c(2,5)]
      resdata1 <- data.frame(cnames,sumrcoeff)
      colnames(resdata1)<-c("Variables","HR","Pvalue")
      rownames(resdata1)<- NULL
      hrpres<-rbind(hrpres,resdata1)

    }
  }
  if(is.decimal(nbatch)==TRUE){
    if((n-sq[nbatch+1])==2){
      m2=sq[nbatch+1]
      n2=m2+(n-sq[nbatch+1])


      cnames<-c(colnames(data)[m2:n2])
      model1 <- coxph(Surv(get(Surv),get(Event)) ~ data[,m2]+data[,m2+1]+data[,m2+2], data =data)
      sumr <- summary(model1)
      sumrcoeff<-sumr$coefficients[,c(2,5)]
      resdata1 <- data.frame(cnames,sumrcoeff)
      colnames(resdata1)<-c("Variables","HR","Pvalue")
      rownames(resdata1)<- NULL
      hrpres<-rbind(hrpres,resdata1)

    }
  }
  if(is.decimal(nbatch)==TRUE){
    if((n-sq[nbatch+1])==3){
      m2=sq[nbatch+1]
      n2=m2+(n-sq[nbatch+1])


      cnames<-c(colnames(data)[m2:n2])
      model1 <- coxph(Surv(get(Surv),get(Event)) ~ data[,m2]+data[,m2+1]+data[,m2+2]+data[,m2+3], data =data)
      sumr <- summary(model1)
      sumrcoeff<-sumr$coefficients[,c(2,5)]
      resdata1 <- data.frame(cnames,sumrcoeff)
      colnames(resdata1)<-c("Variables","HR","Pvalue")
      rownames(resdata1)<- NULL
      hrpres<-rbind(hrpres,resdata1)

    }
  }
  if(is.decimal(nbatch)==TRUE){
    if((n-sq[nbatch+1])==4){
      m2=sq[nbatch+1]
      n2=m2+(n-sq[nbatch+1])

      cnames<-c(colnames(data)[m2:n2])
      model1 <- coxph(Surv(get(Surv),get(Event)) ~ data[,m2]+data[,m2+1]+data[,m2+2]+data[,m2+3]+data[,m2+4], data =data)
      sumr <- summary(model1)
      sumrcoeff<-sumr$coefficients[,c(2,5)]
      resdata1 <- data.frame(cnames,sumrcoeff)
      colnames(resdata1)<-c("Variables","HR","Pvalue")
      rownames(resdata1)<- NULL
      hrpres<-rbind(hrpres,resdata1)

    }
  }

  hrpres <-hrpres[order(hrpres$HR),]
  hrpres1<-hrpres[hrpres$Pvalue <= siglevel,]
  selvar <- c(hrpres1$Variables)



  data2<-data12
  t1<-c(unique(data2[,time]))

  tlast <- t1[length(t1)]
  data2a <- subset(data2, get(time) == tlast)
  data2b <- subset(data2a,select=selvar)


  M <- data.matrix(data2b)

  #define parameters for hdbm
  Y<-M[,2] #data2a[,Event] #response variable
  A<-M[,2] #exposure variable taken as first and second column from the selected variable matrix
  C <- matrix(1, nrow(data2b), 1)
  beta.m  <- rep(0, ncol(data2b))
  alpha.a <- rep(0, ncol(data2b))


  hdbm.out <- hdbm(Y,A,M, C, C, beta.m, alpha.a,
                   burnin = burn, ndraws = draws)

  active <- which(colSums(hdbm.out$r1 * hdbm.out$r3) > thresh) ######################
  Activevariables <- colnames(M)[active] ####################
  mbeta.a<-mean(hdbm.out$beta.a) ###################
  colm.beta.m<-apply(hdbm.out$beta.m, 2, mean) #######################
  colm.alpha.m<-apply(hdbm.out$alpha.a, 2, mean) ######################

  if(length(active)==0){
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
    colnames(act.var.means)<-c("colmeans.beta.m","colmeans.alpha.m")
    act.var.means <- data.frame(Activevariables,act.var.means)

    act.results<-list('Active variabels and their beta and alpha means'= act.var.means)
    return(act.results)
  }

}

