upsilon <- function(x, mediator, dv, conf.level = 0.95, bootstrap = TRUE, bootstrap.package = "lavaan", bootstrap.type = "ordinary", B = 1000, boot.data.out=FALSE, ...){
data <- data.frame(x=x,mediator=mediator,dv=dv)
med.model <- paste0(paste0(colnames(data)[2],'~',colnames(data)[1]),' \n ',paste0(colnames(data)[3],'~',colnames(data)[1],'+',colnames(data)[2]))
med.res <- lavaan::sem(med.model,data)
upsilon <- lavaan::coef(med.res)[1]**2*lavaan::coef(med.res)[3]**2*lavaan::inspectSampleCov(med.model,data)$cov['x','x']/lavaan::inspectSampleCov(med.model,data)$cov['dv','dv']
adj.upsilon <- (lavaan::coef(med.res)[1]**2-lavaan::vcov(med.res)[1,1])*(lavaan::coef(med.res)[3]**2-lavaan::vcov(med.res)[3,3])*lavaan::inspectSampleCov(med.model,data)$cov['x','x']/lavaan::inspectSampleCov(med.model,data)$cov['dv','dv']

if(bootstrap == TRUE){
if(bootstrap.package == 'lavaan'){
lavaan.med.boot.fun <- function(out){
data.boot <- lavaan::lavInspect(out,what = "data")
colnames(data.boot) <- lavaan::lavNames(out)
upsilon.boot <- lavaan::coef(out)[1]**2*lavaan::coef(out)[3]**2*lavaan::inspectSampleCov(med.model,data.boot)$cov['x','x']/lavaan::inspectSampleCov(med.model,data.boot)$cov['dv','dv']
adj.upsilon.boot <- (lavaan::coef(out)[1]**2-lavaan::vcov(out)[1,1])*(lavaan::coef(out)[3]**2-lavaan::vcov(out)[3,3])*lavaan::inspectSampleCov(med.model,data.boot)$cov['x','x']/lavaan::inspectSampleCov(med.model,data.boot)$cov['dv','dv']
as.vector(c(upsilon.boot,adj.upsilon.boot))
}
cat('Bootstrapping may take several minutes \n \n')
boot.med.res <- lavaan::bootstrapLavaan(med.res,R=B,type=bootstrap.type,FUN=lavaan.med.boot.fun)

upsES.out <- data.frame(est=c(upsilon,adj.upsilon),
lcl=c(ups.LCL=quantile(boot.med.res[,1],probs=(1-conf.level)/2),
adj.ups.UCL=quantile(boot.med.res[,2],probs=(1-conf.level)/2)),
ucl=c(ups.LCL=quantile(boot.med.res[,1],probs=1-(1-conf.level)/2),
adj.ups.UCL=quantile(boot.med.res[,2],probs=1-(1-conf.level)/2)),
row.names=c('Upsilon','Adj Upsilon'))
colnames(upsES.out) <- c('Estimate',paste0(conf.level*100,'% ',bootstrap.type,' LCL'),paste0(conf.level*100,'% ',bootstrap.type,' UCL'))
return(upsES.out)
}
if(bootstrap.package == 'boot'){
boot.med.boot.fun <- function(out,i){
data.boot <- out[i,]
colnames(data.boot) <- colnames(out)
out.boot <- lavaan::sem(med.model,data.boot)
upsilon.boot <- lavaan::coef(out.boot)[1]**2*lavaan::coef(out.boot)[3]**2*lavaan::inspectSampleCov(med.model,data.boot)$cov['x','x']/lavaan::inspectSampleCov(med.model,data.boot)$cov['dv','dv']
adj.upsilon.boot <- (lavaan::coef(out.boot)[1]**2-lavaan::vcov(out.boot)[1,1])*(lavaan::coef(out.boot)[3]**2-lavaan::vcov(out.boot)[3,3])*lavaan::inspectSampleCov(med.model,data.boot)$cov['x','x']/lavaan::inspectSampleCov(med.model,data.boot)$cov['dv','dv']
as.vector(c(upsilon.boot,adj.upsilon.boot))
}
cat('Bootstrapping may take several minutes \n \n')
boot.med.res <- boot::boot(data,boot.med.boot.fun,R=B)
upsES.out <- data.frame(est=c(upsilon,adj.upsilon),
lcl=c(ups.LCL=boot::boot.ci(boot.med.res,conf=conf.level,type=bootstrap.type,index=1)[[4]][4],
adj.ups.UCL=boot::boot.ci(boot.med.res,conf=conf.level,type=bootstrap.type,index=2)[[4]][4]),
ucl=c(ups.LCL=boot::boot.ci(boot.med.res,conf=conf.level,type=bootstrap.type,index=1)[[4]][5],
adj.ups.UCL=boot::boot.ci(boot.med.res,conf=conf.level,type=bootstrap.type,index=2)[[4]][5]),
row.names=c('Upsilon','Adj Upsilon'))
colnames(upsES.out) <- c('Estimate',paste0(conf.level*100,'% ',bootstrap.type,' LCL'),paste0(conf.level*100,'% ',bootstrap.type,' UCL'))
if(boot.data.out==TRUE){
return(list(upsilon=upsES.out,bootstrap.data=boot.med.res))
} else return(upsES.out)
}
} else {
upsES.out <- data.frame(Est=c(upsilon,adj.upsilon),row.names=c('Upsilon','Adj Upsilon'))
colnames(upsES.out) <- 'Estimate'
return(upsES.out)
}
}
