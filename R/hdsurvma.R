#' High dimensional survival analysis using SurvMCmulti with mediation analysis
#'
#' @param m Starting column number from where high dimensional variates to be selected.
#' @param n Ending column number till where high dimensional variates to be selected.
#' @param Surv "Column/Variable name" consisting duration of survival.
#' @param Event "Column/Variable name" consisting survival event.
#' @param ths A numeric between 0 to 100.
#' @param chn Number of MCMC chains to perform survival analysis.
#' @param i Number of MCMC iterations to perform survival analysis.
#' @param adp Number of MCMC adaptations to perform survival analysis.
#' @param b Number of MCMC iterations to burn.
#' @param d Number of draws.
#' @param data High dimensional data containing survival observations and high dimensional covariates.
#' @description Given the dimension of variables and survival information the function filters significant variables,
#' allowing the user to perform survival analysis with high number of iterations. Further, it performs mediation analysis among the signifiant
#' variables and provides handful variables with their alpha.a values which are mediator model exposure coefficients
#' and beta.a coefficients.
#'
#' @return Data frame containing the beta and alpha values of active variables among the significant variables.
#' @import survival
#' @import hdbm
#' @import schoolmath
#' @import mlr3
#' @import rjags
#' @export

#' @examples
#' \dontrun{
#' data2 <- hnscc2[1:50,]
#' hdsurvma(m=8,n=15,Surv="os",Event="death",ths=0.02,chn=4,i=10,adp=100,b=10,d=10,data=data2)
#' }
hdsurvma <- function(m,n,Surv,Event,ths,chn,i,adp,b,d,data){
  thresh<-ths
  draws<-d
  burn<-b
  iter<-i
  adapt<-adp
  chains<-chn
  nbatch<-length(m:n)/5
  sq<-seq(m,n,5)
  dic.phd <- matrix(nrow=0,ncol=1)
  vnum <- matrix(nrow=0,ncol=1)

  survjagsmulti <- function(var1=NULL,var2=NULL,var3=NULL,var4=NULL,
                          var5=NULL,Time,Event,chains,adapt,iter,data)
  {
    if(Time!="OS"){
      names(data)[names(data) == Time] <- "OS"
    }
    if(Event!="Death"){
      names(data)[names(data) == Event] <- "Death"
    }

    multi <- c(var1,var2,var3,var4,var5)


    if(length(multi)==5){
      x1=data[,multi[1]]
      x2=data[,multi[2]]
      x3=data[,multi[3]]
      x4=data[,multi[4]]
      x5=data[,multi[5]]
      x = t(rbind(x1,x2,x3,x4,x5))
      p=5
    }
    if(length(multi)==4){
      x1=data[,multi[1]]
      x2=data[,multi[2]]
      x3=data[,multi[3]]
      x4=data[,multi[4]]
      x = t(rbind(x1,x2,x3,x4))
      p=4
    }
    if(length(multi)==3){
      x1=data[,multi[1]]
      x2=data[,multi[2]]
      x3=data[,multi[3]]
      x = t(rbind(x1,x2,x3))
      p=3
    }
    if(length(multi)==2){
      x1=data[,multi[1]]
      x2=data[,multi[2]]
      x = t(rbind(x1,x2))
      p=2
    }


    data<-data[order(data$OS),]
    var1 <- colnames(data)
    nr <- nrow(data)

    data1 <- subset(data, Death == 1) #subsetting data with death status = 1
    u <- unique(data1$OS) #creating a vector with unique values of OS

    #adding a condition for censoring time vector to include the last censored patient when censoring = 0

    if ((data$Death[nrow(data)])==0){
      u1<-c(u,data$OS[nrow(data)])
    } else {
      u1 <- u
    }
    u2 <- sort(u1)
    t.len<-(length(u2)-1)


    datafi <- list(x=x,obs.t=data$OS,t=u2,T=t.len,N=nrow(data),fail=data$Death,eps=1E-10,p=p)


    model_jags <- "
  data{
    # Set up data
  for(i in 1:N) {
    for(j in 1:T) {
    Y[i,j] <- step(obs.t[i] - t[j] + eps)
    dN[i, j] <- Y[i, j] * step(t[j + 1] - obs.t[i] - eps) * fail[i]
    }
  }
  }

  # Model
  model{
  for(i in 1:N){
    betax[i,1] <- 0
    for(k in 2:(p+1)){
      betax[i,k] <- betax[i,k-1] + beta[k-1]*x[i,k-1]
    }
  }
  for(j in 1:T) {
    for(i in 1:N) {
    dN[i, j] ~ dpois(Idt[i, j]) # Likelihood
    Idt[i, j] <- Y[i, j] * exp(betax[i,p+1]) * dL0[j] # Intensity
    }
    dL0[j] ~ dgamma(mu[j], c)
    mu[j] <- dL0.star[j] * c # prior mean hazard
  }
  c <- 0.001
  r <- 0.1
  for (j in 1 : T) {
    dL0.star[j] <- r * (t[j + 1] - t[j])
  }
  for(k in 1:p){
    beta[k] ~ dnorm(0.0,0.000001)
  }
  }"

    params <- c("beta","dL0")

    inits <-  function(){list( beta = rep(0,p), dL0 = rep(0.0001,bigt))}

    jags <- jags.model(textConnection(model_jags),
                       data = datafi,
                       n.chains = chains,
                       n.adapt = adapt)


    samps <- coda.samples(jags, params, n.iter=iter)
    s1 <- summary(samps)
    quan <- s1$quantiles[c(1:p),]
    stats <- s1$statistics[c(1:p),c(1:2)]
    results <- data.frame(stats,quan)
    expresults <- exp(results)
    names(expresults)<-c("Posterior Means","SD","2.5%","25%","50%","75%","97.5%")
    Variables <- multi
    expresults <- cbind(Variables,expresults)
    rownames(expresults) <- NULL

    d= dic.samples(jags, n.iter=iter)
    meandeviance <- round(sum(d$deviance),2)
    output <- list( 'Posterior Estimates' = expresults, 'DIC' = meandeviance)
    return(output)
  }

  survjags <- function(m,n,Time,Event,chains,adapt,iter,data)
  {
    if(Time!="OS"){
      names(data)[names(data) == Time] <- "OS"
    }
    if(Event!="Death"){
      names(data)[names(data) == Event] <- "Death"
    }


    data<-data[order(data$OS),]
    var1 <- colnames(data)
    nr <- nrow(data)

    data1 <- subset(data, Death == 1) #subsetting data with death status = 1
    u <- unique(data1$OS) #creating a vector with unique values of OS

    #adding a condition for censoring time vector to include the last censored patient when censoring = 0

    if ((data$Death[nrow(data)])==0){
      u1<-c(u,data$OS[nrow(data)])
    } else {
      u1 <- u
    }
    u2 <- sort(u1)
    t.len<-(length(u2)-1)



    model_jags <- "
  data{
    # Set up data
  for(i in 1:N) {
    for(j in 1:T) {
    Y[i,j] <- step(obs.t[i] - t[j] + eps)
    dN[i, j] <- Y[i, j] * step(t[j + 1] - obs.t[i] - eps) * fail[i]
    }
  }
  }

  # Model
  model{
  for(i in 1:N){
    betax[i,1] <- 0
    for(k in 2:(p+1)){
      betax[i,k] <- betax[i,k-1] + beta[k-1]*x[i,k-1]
    }
  }
  for(j in 1:T) {
    for(i in 1:N) {
    dN[i, j] ~ dpois(Idt[i, j]) # Likelihood
    Idt[i, j] <- Y[i, j] * exp(betax[i,p+1]) * dL0[j] # Intensity
    }
    dL0[j] ~ dgamma(mu[j], c)
    mu[j] <- dL0.star[j] * c # prior mean hazard
  }
  c <- 0.001
  r <- 0.1
  for (j in 1 : T) {
    dL0.star[j] <- r * (t[j + 1] - t[j])
  }
  for(k in 1:p){
    beta[k] ~ dnorm(0.0,0.000001)
  }
  }"

    params <- c("beta","dL0")

    inits <-  function(){list( beta = rep(0,p), dL0 = rep(0.0001,bigt))}

    x2 <- rep(0,nrow(data))

    q <- matrix(nrow=0,ncol=5)
    s <- matrix(nrow=0,ncol=2)
    di <- matrix(nrow=0,ncol=1)
    for(i in m:n){
      x1 <- data[(1:nrow(data)),i]
      x = t(rbind(x1,x2))

      datafi <- list(x=x,obs.t=data$OS,t=u2,T=t.len,N=nrow(data),fail=data$Death,eps=1E-10,p=2)

      jags <- jags.model(textConnection(model_jags),
                         data = datafi,
                         n.chains = chains,
                         n.adapt = adapt)


      samps <- coda.samples(jags, params, n.iter=iter)
      s1 <- summary(samps)
      stats <- s1$statistics[1,c(1:2)]
      s <- rbind(s,stats)
      quan <- s1$quantiles[1,]
      q <- rbind(q,quan)
      d = dic.samples(jags, n.iter=iter)
      meandeviance <- round(sum(d$deviance),2)
      di <- rbind(di,meandeviance)
    }
    results <- cbind(s,q)
    expresults <- exp(results)

    Variables <- names(data)[m:n]
    expresults <- data.frame(Variables,expresults,di)
    colnames(expresults)<-c("Variable","Posterior Means","SD","2.5%","25%","50%","75%","97.5%","DIC")
    rownames(expresults) <- NULL
    return(expresults)
  }

  if(nbatch>=1){
  for(i in 1:nbatch){
    m1=sq[i]
    n1=m1+4
    cnames<-c(colnames(data)[m1:n1])

    a1<-names(data[m1])
    a2<-names(data[m1+1])
    a3<-names(data[m1+2])
    a4<-names(data[m1+3])
    a5<-names(data[m1+4])

    model1<-survjagsmulti(var1=a1,var2=a2,var3=a3,var4=a4,var5=a5,Time=Surv,Event=Event,chains=chains,adapt=adapt,iter=iter,data=data)
    #model1 <- survMCmulti(a1,a2,a3,a4,a5,duration = Surv,event=Event,chains=chains,iter=iter,data=data)
    m1m1<-model1$DIC
    dic.phd <- rbind(dic.phd,m1m1)
    vnum<-rbind(vnum,m1)
  }

  rownames(vnum)<-NULL
  mod.fil<-cbind(dic.phd,vnum)
  colnames(mod.fil)<-c("DIC","Variable.Number")
  }
  dic.phd <- matrix(nrow=0,ncol=1)
  vnum <- matrix(nrow=0,ncol=1)
  if(is.decimal(nbatch)==TRUE){
    if((n-sq[nbatch+1])==0){
      m2=sq[nbatch+1]
      cnames<-c(colnames(data)[m2])
      #model1 <- coxph(Surv(get(Surv),get(Event)) ~ data[,m2], data =data)
      b1<-names(data[m2])

      model1 <- survjags(m=m2,n=m2,Time=Surv,Event=Event,chains=chains,adapt=adapt,iter=iter,data=data)
      #model1<-survMCmulti(var1=b1,Time=Surv,Event=Event,chains=chains,adapt=adapt,iter=iter,data=data)
      #model1 <- survMCmulti(b1,duration = Surv,event=Event,chains=chains,iter=iter,data=data)
      m1m1<-model1$DIC
      vnum<-rbind(vnum,m2)
      rownames(vnum)<-NULL
      dic.phd <- rbind(dic.phd,m1m1)
      mod.fil1<-cbind(dic.phd,vnum)
      colnames(mod.fil1)<-c("DIC","Variable.Number")
      mod.fil<-rbind(mod.fil,mod.fil1)
      mod.fil<-data.frame(mod.fil)
    }
  }

  dic.phd <- matrix(nrow=0,ncol=1)
  vnum <- matrix(nrow=0,ncol=1)
  if(is.decimal(nbatch)==TRUE){
    if((n-sq[nbatch+1])==1){
      m2=sq[nbatch+1]
      n2=m2+(n-sq[nbatch+1])


      cnames<-c(colnames(data)[m2:n2])
      #model1 <- coxph(Surv(get(Surv),get(Event)) ~ data[,m2]+data[,m2+1], data =data)
      c1<-names(data[m2])
      c2<-names(data[m2+1])
      model1<-survjagsmulti(var1=c1,var2=c2,Time=Surv,Event=Event,chains=chains,adapt=adapt,iter=iter,data=data)
      #model1 <- survMCmulti(c1,c2,duration = Surv,event=Event,chains=chains,iter=iter,data=data)
      m1m1<-model1$DIC

      vnum<-rbind(vnum,m2)
      rownames(vnum)<-NULL
      dic.phd <- rbind(dic.phd,m1m1)
      mod.fil1<-cbind(dic.phd,vnum)
      colnames(mod.fil1)<-c("DIC","Variable.Number")
      mod.fil<-rbind(mod.fil,mod.fil1)
      mod.fil<-data.frame(mod.fil)
    }
  }
  dic.phd <- matrix(nrow=0,ncol=1)
  vnum <- matrix(nrow=0,ncol=1)
  if(is.decimal(nbatch)==TRUE){
    if((n-sq[nbatch+1])==2){
      m2=sq[nbatch+1]
      n2=m2+(n-sq[nbatch+1])


      cnames<-c(colnames(data)[m2:n2])
      #model1 <- coxph(Surv(get(Surv),get(Event)) ~ data[,m2]+data[,m2+1]+data[,m2+2], data =data)
      d1<-names(data[m2])
      d2<-names(data[m2+1])
      d3<-names(data[m2+2])
      model1<-survjagsmulti(var1=d1,var2=d2,var3=d3,Time=Surv,Event=Event,chains=chains,adapt=adapt,iter=iter,data=data)
      #model1 <- survMCmulti(d1,d2,d3,duration = Surv,event=Event,chains=chains,iter=iter,data=data)
      m1m1<-model1$DIC
      vnum<-rbind(vnum,m2)
      rownames(vnum)<-NULL
      dic.phd <- rbind(dic.phd,m1m1)
      mod.fil1<-cbind(dic.phd,vnum)
      colnames(mod.fil1)<-c("DIC","Variable.Number")
      mod.fil<-rbind(mod.fil,mod.fil1)
      mod.fil<-data.frame(mod.fil)

    }
  }
  dic.phd <- matrix(nrow=0,ncol=1)
  vnum <- matrix(nrow=0,ncol=1)
  if(is.decimal(nbatch)==TRUE){
    if((n-sq[nbatch+1])==3){
      m2=sq[nbatch+1]
      n2=m2+(n-sq[nbatch+1])


      cnames<-c(colnames(data)[m2:n2])
      #model1 <- coxph(Surv(get(Surv),get(Event)) ~ data[,m2]+data[,m2+1]+data[,m2+2]+data[,m2+3], data =data)
      e1<-names(data[m2])
      e2<-names(data[m2+1])
      e3<-names(data[m2+2])
      e4<-names(data[m2+3])
      model1<-survjagsmulti(var1=e1,var2=e2,var3=e3,var4=e4,Time=Surv,Event=Event,chains=chains,adapt=adapt,iter=iter,data=data)
      #model1 <- survMCmulti(e1,e2,e3,e4,duration = Surv,event=Event,chains=chains,iter=iter,data=data)
      m1m1<-model1$DIC
      vnum<-rbind(vnum,m2)
      rownames(vnum)<-NULL
      dic.phd <- rbind(dic.phd,m1m1)
      mod.fil1<-cbind(dic.phd,vnum)
      colnames(mod.fil1)<-c("DIC","Variable.Number")
      mod.fil<-rbind(mod.fil,mod.fil1)
      mod.fil<-data.frame(mod.fil)
    }
  }
  dic.phd <- matrix(nrow=0,ncol=1)
  vnum <- matrix(nrow=0,ncol=1)
  if(is.decimal(nbatch)==TRUE){
    if((n-sq[nbatch+1])==4){
      m2=sq[nbatch+1]
      n2=m2+(n-sq[nbatch+1])

      cnames<-c(colnames(data)[m2:n2])
      #model1 <- coxph(Surv(get(Surv),get(Event)) ~ data[,m2]+data[,m2+1]+data[,m2+2]+data[,m2+3]+data[,m2+4], data =data)
      f1<-names(data[m2])
      f2<-names(data[m2+1])
      f3<-names(data[m2+2])
      f4<-names(data[m2+3])
      f5<-names(data[m2+4])
      model1<-survjagsmulti(var1=f1,var2=f2,var3=f3,var4=f4,var5=f5,Time=Surv,Event=Event,chains=chains,adapt=adapt,iter=iter,data=data)
      #model1 <- survMCmulti(f1,f2,f3,f4,f5,duration = Surv,event=Event,chains=chains,iter=iter,data=data)
      m1m1<-model1$DIC
      vnum<-rbind(vnum,m2)
      rownames(vnum)<-NULL
      dic.phd <- rbind(dic.phd,m1m1)
      mod.fil1<-cbind(dic.phd,vnum)
      colnames(mod.fil1)<-c("DIC","Variable.Number")
      mod.fil<-rbind(mod.fil,mod.fil1)
      mod.fil<-data.frame(mod.fil)
    }
  }
  mod.fil<-mod.fil[order(mod.fil$DIC),]


  per <- round(nrow(mod.fil)*0.1,0)
  mod.fil<-mod.fil[c(1:per),]
  nvar <- c(mod.fil[,"Variable.Number"])
  lvar<-nvar[length(nvar)]

  vari <- matrix(nrow=0,ncol=1)
  colnames(vari)<-c("Variabels.Filtered")

  for(i in nvar){
    if((n-i)>=5){
      vari <- rbind(vari,names(data[i]))
      vari <- rbind(vari,names(data[i+1]))
      vari <- rbind(vari,names(data[i+2]))
      vari <- rbind(vari,names(data[i+3]))
      vari <- rbind(vari,names(data[i+4]))
      }
    if((((n-i)<5) & (n-i)==0)==TRUE){
          vari <- rbind(vari,names(data[i]))
          }
    if((((n-i)<5) & (n-i)==1)==TRUE){
        vari <- rbind(vari,names(data[i]))
        vari <- rbind(vari,names(data[i+1]))
        }
    if((((n-i)<5) & (n-i)==2)==TRUE){
        vari <- rbind(vari,names(data[i]))
        vari <- rbind(vari,names(data[i+1]))
        vari <- rbind(vari,names(data[i+2]))
        }
    if((((n-i)<5) & (n-i)==3)==TRUE){
        vari <- rbind(vari,names(data[i]))
        vari <- rbind(vari,names(data[i+1]))
        vari <- rbind(vari,names(data[i+2]))
        vari <- rbind(vari,names(data[i+3]))
        }
    if((((n-i)<5) & (n-i)==4)==TRUE){
        vari <- rbind(vari,names(data[i]))
        vari <- rbind(vari,names(data[i+1]))
        vari <- rbind(vari,names(data[i+2]))
        vari <- rbind(vari,names(data[i+3]))
        vari <- rbind(vari,names(data[i+4]))}
  }


  selvar <- c(vari[,1])
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
utils::globalVariables(c("Death","N","step","obs.t","eps","fail","x1","x2","x3","x4","bigt","x5","dL0","pow","p","coda.samples","dic.samples","jags.model"))
