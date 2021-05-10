#' @title High dimensional missing data imputation and performing the mediation analysis with univariate cox proportional modeling.
#'
#' @param m Starting column number from where high dimensional variates to be selected.
#' @param n Ending column number till where high dimensional variates to be selected.
#' @param survdur "Column/Variable name" consisting duration of survival.
#' @param sig Level of significance pre-determined by the user.
#' @param event "Column/Variable name" consisting survival event.
#' @param t "Column/Variable name" consisting time of repeated observations.
#' @param ths A numeric between 0 to 1.
#' @param b Number of MCMC iterations to burn.
#' @param d Number of draws for the iterations.
#' @param data High dimensional data containing survival observations with multiple covariates.
#' @description Given the dimension of variables and survival information the function
#' performs imputations using missForest function and filters significant variables,
#' allowing the user to fit univariate CoxPH model. Further, it performs mediation
#' analysis among the significant variables and provides handful variables with their alpha.a values
#' which are mediator model exposure coefficients and beta.a coefficients.
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
#' impunicox(m=11,n=25,survdur="OS",event="event",t="Visit",sig=.2,ths=0.02,b=10,d=10,data=srdata)
#' ##
#' }
############################example####################################
#
impunicox <- function(m,n,survdur,event,t,sig,ths,b,d,data){
  burn<-b
  draws<-d
  time<-t
  Surv<-survdur
  Event<-event
  data <- subset(data,select = c(get(Surv),get(Event),get(time),m:n))
  siglevel<-sig
  thresh<-ths
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


  #hrpres <-hrpres[order(hrpres$Pvalue),] #to filter or not??
  hrpres <-hrpres[order(hrpres$HR),] #to filter or not??
  hrpres<-hrpres[hrpres$Pvalue<=siglevel,]
  selvar <- c(hrpres$Variables)




  t1<-c(unique(data12[,time]))

  tlast <- t1[length(t1)]
  data2a <- subset(data12, get(time) == tlast)
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
  thresh=thresh*100
  active <- which(colSums(hdbm.out$r1 * hdbm.out$r3) > thresh*100) ######################
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

