#' @title High dimensional univariate cox proportional hazard analysis.
#'
#' @param m Starting column number form where study variables of high dimensional data will get selected.
#' @param n Ending column number till where study variables of high dimensional data will get selected.
#' @param sig Level of significance pre-determined by the user
#' @param survdur Column name of survival duration event, a string value. i.e. "os"
#' @param event Column name of survival event, a string value. i.e "death"
#' @param data High dimensional data containing the survival, progression and genomic observations.
#' @description Given the dimension of variables and survival information the function performs univariate Cox PH.
#' @return Data set containing the list of selected variables with HR, LCL,UCL and Pvalues through survival analysis.
#' @import survival
#' @import utils
#' @import missForest
#' @importFrom  Rdpack reprompt
#' @export
#' @examples
#' ##
#' data(hnscc2)
#' hidimcox(m=8,n=50,survdur="os",event="death",sig=0.05,data=hnscc2)
#' ##
hidimcox <- function(m,n,survdur,event,sig,data){
  siglevel<-sig
  OS<-survdur
  Event<-event
  data1 <- data
  data2 <- subset(data, select=c(get(OS),get(Event),m:n))
  if(sum(is.na(data2))==0){
    data<-data2
  }

  if(sum(is.na(data2))>=1){
    data.imp <- missForest(data2) #imputed
    data <- data.frame(data.imp$ximp)
  }

  m=3
  n=ncol(data)

  if((nrow(data)/n-m+1)>=10){
  da<-matrix(nrow = 0, ncol = 3) #creating dummy matrix
  dap<-matrix(nrow = 0, ncol = 1)

  for(i in m:n){
    coxmod <- coxph(Surv(get(OS), get(Event)) ~ data[,i], data = data)
    dsum <- summary(coxmod)
    da <- rbind(da,dsum$conf.int[,-2])
    dap <- rbind(dap,round(dsum$coefficients[,5],4))
  }
  }
  dat <- data.frame(da,dap)
  nm <- colnames(data[m:n])
  dan <- cbind(nm,dat)
  colnames(dan)<- c("Variables","HR","LCL","UCL","pvalue")
  rownames(dan)<-NULL
  estimatesresults <- dan
  sig.p.genes <- subset(dan, pvalue < siglevel)
  sig.variables <- c(sig.p.genes$Variables)

  report1 <- list()
  report1[[1]] <- "Dataset containing HR estimates, C.I and Pvalue of variables in column m to n"
  report1[[2]] <- dan

  report2 <- list()
  report2[[1]] <- "Dataset containging HR estimates, C.I and Pvalue of variables with Pvalue less than equal to siglevel"
  report2[[2]] <- sig.variables

  report3 <- list()
  report3[[1]] <- "List of significant genes"
  report3[[2]] <- sig.variables

  output<-list(report1,report2)
  return(output)
}
utils::globalVariables(c("Status","pvalue"))
