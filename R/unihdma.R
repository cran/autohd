#' @title High dimensional survival analysis using Bayesian univariate cox prportional
#'        hazard with mediation analysis
#' @param m Starting column number from where high dimensional variates to be selected.
#' @param n Ending column number till where high dimensional variates to be selected.
#' @param survdur "Column/Variable name" consisting duration of survival.
#' @param event "Column/Variable name" consisting survival event.
#' @param LC "Initial time of getting in to the study.
#' @param t A numeric between 0 to 100.
#' @param i Number of MCMC iteration to perform in obtaining posterior extimates of HR by CoxPH.
#' @param b Number of MCMC iterations to burn.
#' @param d Number of draws for the iterations.
#' @param data High dimensional data containing survival observations and high dimensional covariates.
#'
#' @description Given the dimension of variables and survival information the function filters significant variables
#' allowing the user to perform survival anlaysis with high number of iterations. Further, it performs mediation analysis among the signifiant
#' variables and provides handful variables with their alpha.a values which are mediator model exposure coefficients
#' and beta.a coefficients.
#' @return Data frame containing the beta and alpha values of active variables among the significant variables.
#' @import survival
#' @import hdbm
#' @import schoolmath
#' @import ICBayes
#' @import icenReg
#' @export
#'
#' @examples
#' ##
#' unihdma(m=8,n=15,survdur="os",event="death",LC="leftcensoring",t=0.02,i=6,b=10,d=10,data=hnscc2)
#' ##
unihdma <- function(m,n,survdur,event,LC=NULL,t,i,b,d,data){
   thresh<-t
   iter<-i
   burn<-b
   draws<-d
   Surv<-survdur
   Event<-event
   nbatch<-length(m:n)/5
   sq<-seq(m,n,5)
   Variables<-c(colnames(data)[m:n])
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
   model1 <- survintMC(m,n,Leftcensor=LC,OS=Surv,Death=Event,iter=iter,data = data)
   mod.res <- data.frame(model1)
   mod.res.nts <- matrix(nrow=0,ncol=4)
   mod.res.sp <- matrix(nrow=0,ncol=4)
   mod.res.sn <- matrix(nrow=0,ncol=4)
   mod.res.ntsn <- matrix(nrow=0,ncol=4)
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

  active <- which(colSums(hdbm.out$r1 * hdbm.out$r3) > thresh*100) ######################
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

utils::globalVariables(c("survintMC","lines"))
