#' Chapman Robson catch curve survival estimator
#' 
#' Computes the Chapman-Robson survival estimator and it's standard error from
#' age structure data
#' 
#' 
#' @param ages vector of individual ages, coded such that the first
#' fully-recruited age class is coded as zero.
#' @return list with elemets S_hat and SE(S_hat)
#' @author Michael Scroggie
#' @examples
#' 
#' data(geocrinia)
#' chapman(geocrinia$age)
#' 
#' randata<-rgeom(1000, 1-0.3)
#' chapman(randata)
#' 
chapman<-function(ages){
  X<-sum(ages)
  n<-length(ages)
  surv<-X/(X+n-1)
  var_surv<-surv*(surv-((X-1)/(n+X-2)))
  se_surv<-sqrt(var_surv)
  out<-list(surv, se_surv)
  names(out)<-c("S_hat", "SE(S_hat)")
  return(out)
}
