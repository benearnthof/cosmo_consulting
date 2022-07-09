## ----echo=FALSE,eval=TRUE,message=FALSE----
library(spatstat)
library(mvtnorm)
library(lattice)
library(mgcv)
library(pixmap)
library(numDeriv)
library(fields)
library(INLA)
set.seed(123)
if (file.exists("myinit.R")) source("myinit.R")
inla.setOption(num.threads="1:1")
inla.setOption(smtp="taucs")
opts_chunk$set(message=FALSE, warning=FALSE, tidy=FALSE,
               fig.path='figures/scale-model/')

## ----echo=TRUE,eval=FALSE----------------
#  formula = y ~  f(idx, model = "rw2", scale.model = TRUE, hyper=..)

## ----echo=TRUE,eval=FALSE----------------
#  f(idx, model="..", group=g, control.group = list(model="rw2", scale.model=TRUE))

## ----echo=F,eval=TRUE--------------------
ss.rw1.sdref = function(values,tau) { 
    n = length(values)
    stopifnot(!missing(n)) 
    stopifnot(n > 3) 
    y = NA
    
    formula = y ~ -1 + f(values, 
            model="rw1", 
            constr =TRUE, 
            values = values, 
            diagonal = 1e-10, 
            hyper = list(prec = list(
                                 initial = log(tau), 
                                 fixed=TRUE))) 
    r =inla(formula,data=data.frame(y,values),family = "gaussian",
            control.family = list(fixed =
                    TRUE), 
            control.compute=list(return.marginals=F),verbose=F)
    sd=r$summary.random$values$sd 
    return (sd) 
    }

ss.rw2.sdref= function(values,tau) { 
    n=length(values)
    stopifnot(!missing(n)) 
    stopifnot(n > 3)
    y = NA #idx =1:n
    
    A1 = matrix(1, n, 1) 
    A2 = matrix(0, n, 1) 
    for(i in 1:n) { 
        A2[i, ]= i 
    }
    
    A = rbind(c(A1), c(A2)) 
    e = rep(0, 2) 
    extraconstr = list(A = A, e= e)
    
    formula = y ~ -1 + f( values, model="rw2", constr = FALSE,
            values=values, diagonal = 1e-10, extraconstr = extraconstr, hyper
            = list(prec = list(initial = log(tau), fixed=TRUE)))
    
    r = inla(formula, data = data.frame(y, values), family =
            "gaussian", control.family = list(fixed = TRUE),
            control.compute=list(return.marginals=F),verbose=F)
    
    sd=r$summary.random$values$sd 
    return (sd) 
}

n=100 
sd.rw1=ss.rw1.sdref(1:n,1) 
sd.rw2=ss.rw2.sdref(1:n,1)
##sd.rw1.scaled=ss.rw1.sdref(1:n,1) 
##sd.rw2.scaled=ss.rw2.sdref(1:n,1)

## ----fig.keep='high',echo=FALSE,out.width='0.45\\linewidth'----
plot(sd.rw1,xlab="Node i", ylab="Marginal standard deviation") 

## ----fig.keep='high',echo=FALSE,out.width='0.45\\linewidth'----
plot(sd.rw2,xlab="Node i",ylab="Marginal standard deviation") 

## ----echo=TRUE,eval=F--------------------
#  U = sqrt(b/qgamma(alpha,a,1))

## ----echo=FALSE,eval=TRUE----------------
func.u = function(a,b,alpha,sigma.ref) { 
    upper.limit = sqrt(b)*sigma.ref/sqrt(qgamma(alpha,a,1)) 
    return(upper.limit) 
}
a=1 
b=5*10^{-5} 
alpha=0.001 
upper.limit=func.u(a,b,alpha,1) 

## ----echo=F,eval=TRUE--------------------
func<-function(x,t)
{
    f.x=x/t+sin(2*pi* x/t)-0.5
    return(f.x)
}

f.est <- function(y,x,a,b,option)
{
    formula = y~1+f(x,model="rw2",hyper=list(prec=list(prior="loggamma",param=c(a,b))),scale.model=option)
    result = inla(formula,family="gaussian",data=data.frame(y,x))
    f.est=result$summary.random$x$mean
    return(f.est)
}

## ----echo=F,eval=TRUE--------------------
### Same values of the function, same observations, but  defined on different intervals ###
set.seed(89236)
x2=0:n
x1=x2/100
x3=x2*10

sigma=0.5
f.x=func(x2,length(x2))
y=f.x+rnorm(length(f.x),0,sigma)

a=1
b=0.00005

mat1=cbind(f.x,f.est(y,x1,a,b,T),f.est(y,x1,a,b,F))
mat2=cbind(f.x,f.est(y,x2,a,b,T),f.est(y,x2,a,b,F))
mat3=cbind(f.x,f.est(y,x3,a,b,T),f.est(y,x3,a,b,F))

# Generalized marginal variances  
v1=exp(mean(log(ss.rw2.sdref(x1,1)^2)))
v2=exp(mean(log(ss.rw2.sdref(x2,1)^2)))
v3=exp(mean(log(ss.rw2.sdref(x3,1)^2)))

u1=func.u(a,b,alpha,sqrt(v1))
u2=func.u(a,b,alpha,sqrt(v2))
u3=func.u(a,b,alpha,sqrt(v3))

## ----echo=TRUE,eval=F--------------------
#  formula = y ~  f(x, model = "rw2", scale.model = TRUE, hyper=...)

## ----echo=TRUE,eval=F--------------------
#  result = inla(formula, family = "gaussian", data = data.frame(y,x))

## ----echo=TRUE,eval=F--------------------
#  result$summary.random$x$mean

## ----fig.keep='high',echo=FALSE,eval=TRUE,out.width='0.28\\linewidth'----
                                        # Estimate of f using unscaled and scaled model
r=range(c(mat1,mat2,mat3))
matplot(x1,mat1,type="lll",lty=c(1,2,4),col=c(1,4,2),xlab="x",ylab="f(x)",ylim=r,lwd=3) 
                                        #legend("topright",c("true f","scale.model=F","scale.model=T"),lty=c(1,2,4),col=c(1,2,4))
                                        #points(x,y)

## ----fig.keep='high',echo=FALSE,eval=TRUE,out.width='0.28\\linewidth'----
matplot(x2,mat2,type="lll",lty=c(1,2,4),col=c(1,4,2),xlab="x",ylab="f(x)",ylim=r,lwd=3)

## ----fig.keep='high',echo=FALSE,eval=TRUE,out.width='0.28\\linewidth'----
matplot(x3,mat3,type="lll",lty=c(1,2,4),col=c(1,4,2),xlab="x",ylab="f(x)",ylim=r,lwd=3)

## ----echo=TRUE,eval=F--------------------
#  hyper.prec = list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))

## ----echo=TRUE,eval=F--------------------
#  data(Munich)
#  g = system.file("demodata/munich.graph", package="INLA")
#  
#  ## Note that here we what to have an estimator of the effect of year
#  ## also the for years where we have no observation, therefore we give a
#  ## vector with all possible values assumed by the covariate year, that
#  ## is seq(1918,2001)
#  
#  formula = rent ~
#     f(location, model = "besag", graph = g, initial = 1, hyper = hyper.prec) +
#     f(year, model = "rw2", values = seq(1918,2001), hyper = hyper.prec) +
#     f(floor.size, model = "rw2", hyper = hyper.prec) +
#     Gute.Wohnlage + Beste.Wohnlage + Keine.Wwv + Keine.Zh +
#     Kein.Badkach  + Besond.Bad + Gehobene.Kueche +
#     zim1 + zim2 + zim3 + zim4 + zim5 + zim6 -1

## ----echo=TRUE,eval=F--------------------
#  mod  =  inla(formula, data = Munich, control.fixed=list(prec=1))

## ----echo=FALSE,eval=TRUE----------------
data(Munich)
g = system.file("demodata/munich.graph", package="INLA")

func.munich <- function(a,b,option)
{
    hyper.prec=list(prec=list(param=c(a,b)))
    formula= rent ~ f(location,model="besag",graph=g,initial=1, hyper=hyper.prec, scale.model=option) +
        f(year,model="rw2",values=seq(1918,2001), hyper=hyper.prec, scale.model=option) +
            f(floor.size,model="rw2", hyper=hyper.prec, scale.model=option) +
                Gute.Wohnlage + Beste.Wohnlage + Keine.Wwv + Keine.Zh +
                    Kein.Badkach  + Besond.Bad + Gehobene.Kueche +
                        zim1 + zim2 + zim3 + zim4 + zim5 + zim6 -1
    
    result = inla(formula, data=Munich, control.fixed=list(prec=1))
}

## ----echo=FALSE,eval=TRUE----------------
a=1
b=0.01
u.munich=func.u(a,b,alpha,1)
mod=func.munich(a,b,FALSE)
mod.scaled=func.munich(a,b,TRUE)

## ----fig.keep='high',echo=FALSE,eval=TRUE,out.width='0.45\\linewidth'----
mat=cbind(mod$summary.random$year$mean,mod.scaled$summary.random$year$mean)
matplot(mod$summary.random$year$ID,mat,type="ll",lty=c(2,1),col=c(2,1),xlab="Year of construction",ylab="",lwd=2.5)

## ----fig.keep='high',echo=FALSE,eval=TRUE,out.width='0.45\\linewidth'----
mat=cbind(mod$summary.random$floor.size$mean,mod.scaled$summary.random$floor.size$mean)
matplot(mod$summary.random$floor.size$ID,mat,type="ll",lty=c(2,1),col=c(2,1),xlab="Floor size",ylab="",lwd=2.5)

## ----echo=FALSE,eval=TRUE----------------
alpha=0.001
#u.vec=c(0.01,1,10,100)
u.vec=c(0.001,0.1,10,30)
a.vec=c(1,1,1,1)
b.vec=u.vec^2*qgamma(alpha,a.vec,1)

option=TRUE
c1=func.munich(a.vec[1],b.vec[1],option)
c2=func.munich(a.vec[2],b.vec[2],option)
c3=func.munich(a.vec[3],b.vec[3],option)
c4=func.munich(a.vec[4],b.vec[4],option)

## ----fig.keep='high',echo=FALSE,eval=TRUE,out.width='0.45\\linewidth'----
mat=cbind(c1$summary.random$year$mean,c2$summary.random$year$mean,c3$summary.random$year$mean,c4$summary.random$year$mean)
matplot(mod$summary.random$year$ID,mat,type="llll",lty=c(1,2,3,4),col=c(1,2,3,4),xlab="Year of construction",ylab="",lwd=2.5)
legend("topleft",paste("U= ",u.vec),lty=1:4,col=1:4,lwd=2.5)

## ----fig.keep='high',echo=FALSE,eval=TRUE,out.width='0.45\\linewidth'----
mat=cbind(c1$summary.random$floor.size$mean,c2$summary.random$floor.size$mean,c3$summary.random$floor.size$mean,c4$summary.random$floor.size$mean)
matplot(mod$summary.random$floor.size$ID,mat,type="llll",lty=c(1,2,3,4),col=c(1,2,3,4),xlab="Floor size",ylab="",lwd=2.5)
legend("topright",paste("U= ",u.vec),lty=1:4,col=1:4,lwd=2.5)

## ----echo=TRUE,eval=FALSE----------------
#  inla.setOption(scale.model.default = TRUE)

