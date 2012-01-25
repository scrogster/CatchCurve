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
Std.err<-sqrt(diag(solve(opt.params$hessian)))
Vcov<-solve(opt.params$hessian)
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
