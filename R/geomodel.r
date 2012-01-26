#' Linear models for catch curve data, using the geometric probability
#' distribution
#' 
#' fits linear models to catch curve data with covariates using maximum
#' likelihood
#' 
#' 
#' @param formula formula relating the age data to associated covariates
#' @param data dataframe containing the variables in the formula
#' @return A list with elements "formula", "coefficients", "vcov", "logLik",
#' "df", "AICc"
#' @author Michael Scroggie
#' @examples
#' 
#' data(geocrinia)
#' geomodel(age~sex, data=geocrinia)
#' geomodel(age~1, data=geocrinia)
#' 
geomodel<-function(formula=NULL, data=NULL){
start.params<-rep(0, times=ncol(model.matrix(as.formula(formula), data=data)) )
modframe<-model.frame(as.formula(formula), data=data)
mod.mat<- model.matrix(as.formula(formula), data=data)
response<- modframe[,1]
 func<-function(x, p) {p*(1-p)^x}
geomloglik<-function(par, form=formula) {
  linpred<- mod.mat%*%par
  pred<-exp(linpred)/(1+exp(linpred))
  loglik<-sum(dgeom(response, 1-pred, log=T))*-1
  loglik
}
opt.params<-optim(par=start.params, fn=geomloglik, hessian=T, method="CG", 
                  control=list( trace=F))
Estimate<-opt.params$par
names(Estimate)<-dimnames(mod.mat)[[2]]
Vcov<-solve(opt.params$hessian)
Std.err<-sqrt(diag(Vcov))
coeffs<-cbind(Estimate, Std.err)
logLik<-opt.params$value
df<-length(opt.params$par)
nn<-nrow(data)
AICc<- (2*logLik)+(2*df)+(2*df*(df+1)/(nn-df-1))
out<-list(formula, coeffs, Vcov, logLik, df, AICc)
names(out)<-c("formula", "coefficients", "vcov", "logLik", "df", "AICc")
class(out)<-c("geommodel")
return(out)
}
