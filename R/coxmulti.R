#' @title Multivariate cox proportional hazard data analysis
#'
#' @param C1 Covar1
#' @param C2 Covar2
#' @param C3 Covar3
#' @param C4 Covar4
#' @param C5 Covar5
#' @param survdur "Column/Variable name" consisting duration of survival.
#' @param event "Column/Variable name" consisting survival event.
#' @param data High dimensional data containing survival observations and covariates.
#' @description Given the dimensions of the variables and survival informations. The function
#' performs multivariate Cox PH by taking 5 variables at a time.
#'
#' @return Data set containing the survival estimates and Pvalue.
#' @import survival
#' @import schoolmath
#' @import missForest
#' @export
#'
#' @examples
#' ##
#' coxmulti(C1="PGC",C2="C7",C3="HPN",C4="DDC",C5=NULL,survdur="os",event="death",data=hnscc2)
coxmulti <- function(C1=NULL,C2=NULL,C3=NULL,C4=NULL,C5=NULL,survdur,event,data){
  covar1<-C1
  covar2<-C2
  covar3<-C3
  covar4<-C4
  covar5<-C5
  Surv<- survdur
  Event<-event
  data1 <- data
   multi <- c(covar1,covar2,covar3,covar4,covar5)

  data2a <- subset(data, select=c(get(Surv),get(Event)))
  data2b <- subset(data, select=c(multi))
  data2 <- data.frame(data2a,data2b)
  if(sum(is.na(data2))==0){
    data<-data2
  }

  if(sum(is.na(data2))>=1){
    data.imp <- missForest(data2) #imputed
    data <- data.frame(data.imp$ximp)
  }

   if(length(multi)==5){
    model1 <- coxph(Surv(get(Surv),get(Event)) ~ data[,multi[1]]+data[,multi[2]]+data[,multi[3]]+data[,multi[4]]+data[,multi[5]], data=data)
    sumr <- summary(model1)
    sumrcoeff<-sumr$coefficients[,c(2,5)]
    hr <- sumr$conf.int[,-2]

    p1 <- round(sumr$coefficients[,5][1],2)
    p2 <- round(sumr$coefficients[,5][2],2)
    p3 <- round(sumr$coefficients[,5][3],2)
    p4 <- round(sumr$coefficients[,5][4],2)
    p5 <- round(sumr$coefficients[,5][5],2)
    pv <- rbind(p1,p2,p3,p4,p5)
    cn <- multi
    fd <- cbind(cn,hr,pv)
    fadata <- data.frame(fd)
    rownames(fadata)<-NULL
    colnames(fadata)<-c("Variables","HR","LCL","UCL","Pvalue")
    }
    if(length(multi)==4){
      model1 <- coxph(Surv(get(Surv),get(Event)) ~ data[,multi[1]]+data[,multi[2]]+data[,multi[3]]+data[,multi[4]], data=data)
      sumr <- summary(model1)
      sumrcoeff<-sumr$coefficients[,c(2,5)]
      hr <- sumr$conf.int[,-2]

      p1 <- round(sumr$coefficients[,5][1],2)
      p2 <- round(sumr$coefficients[,5][2],2)
      p3 <- round(sumr$coefficients[,5][3],2)
      p4 <- round(sumr$coefficients[,5][4],2)
      pv <- rbind(p1,p2,p3,p4)
      cn <- multi
      fd <- cbind(cn,hr,pv)
      fadata <- data.frame(fd)
      rownames(fadata)<-NULL
      colnames(fadata)<-c("Variables","HR","LCL","UCL","Pvalue")
      }
    if(length(multi)==3){
      model1 <- coxph(Surv(get(Surv),get(Event)) ~ data[,multi[1]]+data[,multi[2]]+data[,multi[3]], data=data)
      sumr <- summary(model1)
      sumrcoeff<-sumr$coefficients[,c(2,5)]
      hr <- sumr$conf.int[,-2]

      p1 <- round(sumr$coefficients[,5][1],2)
      p2 <- round(sumr$coefficients[,5][2],2)
      p3 <- round(sumr$coefficients[,5][3],2)
      pv <- rbind(p1,p2,p3)
      cn <- multi
      fd <- cbind(cn,hr,pv)
      fadata <- data.frame(fd)
      rownames(fadata)<-NULL
      colnames(fadata)<-c("Variables","HR","LCL","UCL","Pvalue")
      }
    if(length(multi)==2){
      model1 <- coxph(Surv(get(Surv),get(Event)) ~ data[,multi[1]]+data[,multi[2]], data=data)
      sumr <- summary(model1)
      sumrcoeff<-sumr$coefficients[,c(2,5)]
      hr <- sumr$conf.int[,-2]

      p1 <- round(sumr$coefficients[,5][1],2)
      p2 <- round(sumr$coefficients[,5][2],2)
      pv <- rbind(p1,p2)
      cn <- multi
      fd <- cbind(cn,hr,pv)
      fadata <- data.frame(fd)
      rownames(fadata)<-NULL
      colnames(fadata)<-c("Variables","HR","LCL","UCL","Pvalue")
      }
    if(length(multi)==1){
      model1 <- coxph(Surv(get(Surv),get(Event)) ~ data[,multi[1]], data=data)
      sumr <- summary(model1)
      sumrcoeff<-sumr$coefficients[,c(2,5)]
      hr <- sumr$conf.int[,-2]

      p1 <- round(sumr$coefficients[,5][1],2)
      pv <- rbind(p1,p2)
      cn <- multi
      fd <- cbind(cn,hr,pv)
      fadata <- data.frame(fd)
      rownames(fadata)<-NULL
      colnames(fadata)<-c("Variables","HR","LCL","UCL","Pvalue")
      }
  return(fadata)
  }

utils::globalVariables(c("C1","C2","C3","C4","C5"))
