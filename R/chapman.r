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