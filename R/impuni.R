#' @title High dimensional missing data imputation and performing the mediation analysis with Bayesian
#'        univariate cox proportional modeling.
#'
#' @param m Starting column number from where high dimensional variates to be selected.
#' @param n Ending column number till where high dimensional variates to be selected.
#' @param Survdur "Column/Variable name" consisting duration of survival.
#' @param event "Column/Variable name" consisting survival event.
#' @param time "Column/Variable name" consisting time of repeated observations.
#' @param lcr "Leftcensoring information"
#' @param t A numeric threshold value between 0 to 1.
#' @param i Number of MCMC iteration to perform in obtaining posterior estimates of HR by CoxPH.
#' @param b Number of MCMC iterations to burn.
#' @param d Number of draws for the iterations.
#' @param data High dimensional data containing survival observations with multiple covariates.
#' @description Given the dimension of variables and survival information the function
#' performs imputations using missForest function and filters significant variables,
#' allowing the user to do univariate survival analysis with higher number of iterations. Further, it performs mediation
#' analysis among the significant variables and provides handful variables with their alpha.a values
#' which are mediator model exposure coefficients and beta.a coefficients.
#' @return Data frame containing the beta and alpha values of active variables among the significant variables.
#' @import survival
#' @import hdbm
#' @import schoolmath
#' @import missForest
#' @import ICBayes
#' @import icenReg
#' @export
#'
#' @examples
#' \dontrun{
#' ##
#' impuni(m=8,n=25,Survdur="os",event="death",lcr=,t=0.02,i=6,b=10,d=10,data=hnscc)
#' ##
#' }
############################example####################################
#
impuni <- function(m,n,Survdur,event,time,lcr=NULL,t,i,b,d,data){
  Event<-event
  if(is.null(lcr)==FALSE){
  data1 <- subset(data,select=c(get(Survdur),get(Event),get(lcr),m:n))
  m=4
  n=ncol(data1)
  }
  if(is.null(lcr)==TRUE){
    data1 <- subset(data,select=c(get(Survdur),get(Event),m:n))
    m=3
    n=ncol(data1)
  }

  burn<-b
  draws<-d
  iter<-i
  thresh<-t
  if(sum(is.na(data1))>0){
  data.imp <- missForest(data1) #imputed tmc longit
  data1 <- data.frame(data.imp$ximp)
  }

  nbatch<-length(m:n)/5
  sq<-seq(m,n,5)
   Variables<-c(colnames(data1)[m:n])
  data12<-data1

    #model1 <- coxph(Surv(get(Survdur),get(Event)) ~ data[,i], data=data)
  survintMC <- function(m,n,Leftcensor=NULL,OS,Death,iter,data){

    data1 <- data
    data2 <- data
    colnames(data2) <- NULL
    hrt <- matrix(ncol = 1)
    hrci <- matrix(ncol = 2)
    variables <- matrix(ncol = 1)
    burn=(iter/2)
    lf<-data1[,Leftcensor]
    if(length(Leftcensor)==0)
    {
      lftcen<-c(rep(0,nrow(data)))
    }
    if(length(Leftcensor)==1)
    {
      lftcen<-lf
    }
    for(i in m:n){
      breastICB <- ICBayes(model = "case2ph",
                           L = lftcen, R = data1[,OS], status = data1[,Death],
                           xcov = data2[, i], x_user = c(0, 1),
                           knots = seq(0.1, 60.1, length = 4),
                           grids = seq(0.1, 60.1, by = 1),
                           niter = iter, burnin = (iter/2)
      )
      ngrid <- length(breastICB$S0_m)
      #plot(breastICB$grids, breastICB$S_m[1:ngrid], type = "l",
          # lty = "solid", xlab = "Survival times (months)", main = names(data1)[i],
          # ylab = "Estimated survival distributions", ylim = c(0, 1))
      #lines(breastICB$grids, breastICB$S_m[(ngrid+1):(2*ngrid)],
           # lty = "dashed")
      HR <- exp(breastICB$coef)
      HR.CI <- exp(breastICB$coef_ci)

      hrt<-rbind(hrt,HR)
      hrci<-rbind(hrci,HR.CI)
      variables<-rbind(variables,names(data1)[i])
    }
    survintMCout <- data.frame(variables,hrt,hrci)
    survintMCout  <- survintMCout[-1,]
    colnames(survintMCout ) <- c("Variables","HR","LCL","UCL")
    rownames(survintMCout ) <- NULL
    return(survintMCout)
  }

    model3 <- survintMC(m,n,Leftcensor=lcr,OS=Survdur,Death=Event,iter=6,data = data1)
    mod.res <- data.frame(model3)
    mod.res.nts <- matrix(nrow=0,ncol=4)
    mod.res.sp <- matrix(nrow=0,ncol=4)
    mod.res.sn <- matrix(nrow=0,ncol=4)
    mod.res.ntsn <- matrix(nrow=0,ncol=4)
    hrpres <- matrix(nrow=0,ncol=4)
    for(i in 1:nrow(mod.res)){
      z1 <- (1-mod.res[i,"LCL"])
      z2 <- (mod.res[i,"LCL"]-1)
      if((z1>0 & z2>0)==TRUE){
        mod.res.nts <- rbind(mod.res.nts,mod.res[i,])
      }
      if((z1>0 & z2<0)==TRUE){
        mod.res.sn <- rbind(mod.res.sn,mod.res[i,])
      }
      if((z1<0 & z2>0)==TRUE){
        mod.res.sp <- rbind(mod.res.sp,mod.res[i,])
      }
      if((z1<0 & z2<0)==TRUE){
        mod.res.ntsn <- rbind(mod.res.ntsn,mod.res[i,])
      }
    }



    if((nrow(mod.res.sp)!=0) & (nrow(mod.res.sn)!=0)){
      hrpres<-rbind(mod.res.sp,mod.res.sn)
    }
    if((nrow(mod.res.sp)!=0) & (nrow(mod.res.sn)==0)){
      hrpres<-mod.res.sp
    }
    if((nrow(mod.res.sp)==0) & (nrow(mod.res.sn)!=0)){
      hrpres<-mod.res.sn
    }

    if(nrow(hrpres)==0){print("HR Estiamtes cannot be calculated")}

    selvar <- c(hrpres$Variables)

  data2b <- subset(data1,select=selvar)


  M <- data.matrix(data2b)

  #define parameters for hdbm
  Y<-M[,1] #data2a[,Event] #response variable
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
utils::globalVariables(c("survintMC","hrpres","lines"))
