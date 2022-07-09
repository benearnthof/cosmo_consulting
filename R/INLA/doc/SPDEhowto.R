## ----setting,include=FALSE----------------------------------------------------
library(knitr)
knit_hooks$set(pars = function(before, options, envir) {
    if (before) graphics::par(options$pars)
})
opts_chunk$set(message=FALSE, warning=FALSE, tidy=FALSE,
               fig.path='figures/SPDEhowto/')
library(lattice);   library(gridExtra);  library(INLA);  library(plyr) 
set.seed(123)
if (file.exists("myinit.R")) source("myinit.R")
inla.setOption(num.threads="1:1")
inla.setOption(smtp="taucs")
lcall <- inla.getOption('inla.call')

## ----locations,tidy=FALSE-----------------------------------------------------
n <- 200; coo <- matrix(runif(2*n), n) 
s2u <- .5; k <- 10; r <- 2/k ## RF params.
R <- s2u*exp(-k*as.matrix(stats::dist(coo))) 

## ----rmvnorm,tidy=FALSE-------------------------------------------------------
u <- drop(rnorm(n)%*%chol(R)) 

## ----noise,tidy=FALSE---------------------------------------------------------
x <- runif(n);  beta <- 1:2;  s2e <- 0.2
lin.pred <- beta[1] + beta[2]*x + u 
y <- lin.pred + rnorm(n, 0, sqrt(s2e)) 

## ----fmesh, tidy=FALSE, pars=list(mar=c(0,0,0.7,0)), out.width='0.45\\textwidth', fig.align='center'----
mesh <- inla.mesh.2d(coo, cutoff=r/10, 
  max.edge=c(r/4, r/2), offset=c(r/2, r)) 
plot(mesh, asp=1);  points(coo, col='red') 

## ----projector, tidy=FALSE----------------------------------------------------
A <- inla.spde.make.A(mesh=mesh, loc=coo)

## ----spde,tidy=FALSE----------------------------------------------------------
spde <- inla.spde2.pcmatern(
 mesh=mesh, alpha=1.5, 
 prior.range=c(0.2, 0.5),#P(range<0.2)=0.5
 prior.sigma=c(1, 0.5)) ## P(sigma>1)=0.5 

## ----stack-estimation,tidy=FALSE----------------------------------------------
stk.e <- inla.stack(tag='est', ## tag id
  data=list(y=y),  ## response
  A=list(1, A), ## two projection matrices
  effects=list(## two elements: 
    data.frame(b0=1, x=x), ## covariate
    idx.u=1:spde$n.spde)) ## RF index 

## ----fitting------------------------------------------------------------------
pcprec <- list(hyper=list(theta=list(
  prior='pc.prec',param=c(1,.1))))
mf <- y ~ 0 + b0 + x + f(idx.u, model=spde) 
res <- inla(mf, control.family=pcprec,
 data=inla.stack.data(stk.e), ## data 
 control.compute=list(return.marginals.predictor=TRUE),
 control.predictor=list(compute=TRUE, 
  A=inla.stack.A(stk.e)))# full projector

## ----fixed-summary------------------------------------------------------------
round(res$summary.fixed, 4) 

## ----sdlik--------------------------------------------------------------------
pmd.s2e <- inla.tmarginal(function(x) sqrt(1/x), ## inverse and square root
  res$marginals.hyperpar$'Precision for the Gaussian observations')

## ----rfpars, pars=list(mfrow=c(1,3), mar=c(3,3,0.3,0.3), mgp=c(2,0.5,0)), fig.width=7.5, fig.height=2.5, out.width='0.99\\textwidth', out.height='0.33\\textwidth', fig.align='center', fig.pos='h', fig.keep='last'----
plot(pmd.s2e, type='l', ylab='Density',  xlab=expression(sigma[e]))
abline(v=sqrt(s2e), col=2) ## add the 'true' value
plot(res$marginals.hy[[3]], type='l', xlab=expression(sigma[u]), yla='Density')
abline(v=sqrt(s2u), col=2) ## add the 'true' value
plot(res$marginals.hy[[2]], type='l', xlab='range nominal', ylab='Density')
abline(v=r, col=2) ## add the 'true' value

## ----predicts, pars=list(mfrow=c(1,1), mar=c(1.5,3,0.1,0.1), mgp=c(2,0.5,0)), fig.width=10, fig.height=2.5, out.width='0.99\\textwidth', fig.align='center', fig.pos='h', fig.keep='last'----
idx.obs <- inla.stack.index(stk.e, tag='est')$data ## index in the stack 
order.eta <- order(res$summary.fitted.values$mean[idx.obs]) 
plot(y[order.eta], pch=19, ylab='y')
segments(1:n, res$summary.fitted.val$'0.025quant'[idx.obs][order.eta], 
         1:n, res$summary.fitted.val$'0.975quant'[idx.obs][order.eta]) 

## ----project-grid-------------------------------------------------------------
nxy <- c(200, 200)
gproj <- inla.mesh.projector(mesh,  xlim=0:1, ylim=0:1, dims=nxy)
g.mean <- inla.mesh.project(gproj, res$summary.random$idx.u$mean)
g.sd <- inla.mesh.project(gproj, res$summary.random$idx.u$sd)

## ----fgrid, tidy=FALSE, fig.width=9.7, fig.height=4.5, out.width='0.97\\textwidth', out.height='0.45\\textwidth', fig.pos='h'----
library(lattice);     library(gridExtra) 
trellis.par.set(regions=list(col=terrain.colors(16))) 
grid.arrange(levelplot(g.mean, scales=list(draw=F), xlab='', ylab='', main='mean'), 
             levelplot(g.sd, scal=list(draw=F), xla='', yla='', main='sd'), nrow=1)

## ----target-loc---------------------------------------------------------------
tcoo <- rbind(c(0.3,0.9), c(0.5,0.5), c(0.7,0.3))
dim(Ap <- inla.spde.make.A(mesh=mesh, loc=tcoo)) 
x0 <- c(0.5, 0.5, 0.5)

## ----prediction-stack---------------------------------------------------------
stk.pred <- inla.stack(tag='pred', A=list(Ap, 1), data=list(y=NA), ## response as NA
  effects=list(idx.u=1:spde$n.spde, data.frame(x=x0, b0=1))) ## all idx.u
stk.full <- inla.stack(stk.e, stk.pred) ## join the data and prediction scenario
p.res <- inla(mf, data=inla.stack.data(stk.full), ## supply the full data 
  control.compute=list(return.marginals.predictor=TRUE),
  control.predictor=list(compute=TRUE, A=inla.stack.A(stk.full)), ## full 
  control.mode=list(theta=res$mode$theta, restart=FALSE))## use mode already found

## ----prdind-------------------------------------------------------------------
pred.ind <- inla.stack.index(stk.full, tag='pred')$data
round(p.res$summary.fitted.val[pred.ind,], 4)

## ----ppred, tidy=FALSE, fig.width=9.9, fig.height=3.5, out.width='0.99\\textwidth', out.height='0.35\\textwidth', fig.pos='h', fig.keep='last'----
ypost <- p.res$marginals.fitted.values[pred.ind]
names(ypost) <- paste('y', seq_along(ypost), sep='_');   
xyplot(y~x | .id, ldply(ypost), panel='llines', xlab='y', ylab='Density')

## ----echo=FALSE--------------------------
options(width=43)

## ----marginals-funcs, comment=NA---------
apropos('marginal')

## ----marginals-examples------------------
## inla.mmarginal(ypost[[1]]) ## mode
inla.qmarginal(c(0.15, 0.7), 
               ypost[[1]]) ## quantiles
inla.pmarginal(inla.qmarginal(
    0.3, ypost[[1]]), ypost[[1]]) 

## ----configs, warning=FALSE, message=FALSE----
res$.args$control.compute$config <- TRUE
res <- inla.rerun(res)

## ----echo=FALSE,results=FALSE------------
inla.setOption(inla.call=lcall) 

## ----samples-----------------------------
samples <- inla.posterior.sample(n=2e3, result=res, add.names=FALSE)

## ----xnames------------------------------
xnames <- rownames(samples[[1]]$latent) ### collect the names
idx <- lapply(c('b0', 'x', 'idx.u'), function(nam) ## for each effect
    which(substr(xnames, 1, nchar(nam))==nam)) ## find the index

## ----samplesmat--------------------------
mat.samples <- sapply(samples, function(spl) 
    c(bo=spl$latent[idx[[1]]], x=spl$latent[idx[[2]]], u=spl$latent[idx[[3]]]))

## ----y3sample----------------------------
eta3.sample <- as.matrix(cbind(b0=1, x=0.5, u=Ap)%*%mat.samples)

## ----comparey3, fig.width=9.9, fig.height=3.5, out.width='0.99\\textwidth', fig.pos='h', fig.keep='last'----
par(mfrow=c(1,3), mar=c(3,3,1,1), mgp=c(2,1,0))
for (j in 1:3) {
    hist(eta3.sample[j,], freq=FALSE, xlab=paste0('y',j), main='', col=gray(0.7))
    lines(ypost[[j]], lwd=2, col='blue')
}

## ----linpredsample-----------------------
eta.g.samp <- as.matrix(cbind(b0=1, x=0.5, 
        s=gproj$proj$A)%*%mat.samples)

## ----meansdgrid,warning=FALSE,message=FALSE----
library(matrixStats)
ss <- list(mean=rowMeans(eta.g.samp), 
 sd=rowSds(eta.g.samp), 
 p3=rowMeans(eta.g.samp>3), 
 q95=rowQuantiles(eta.g.samp, probs=0.95))

## ----predmaps, eval=FALSE----------------
#  library(fields); for (s in c(1,3,4,2))
#    image.plot(list(x=gproj$x, y=gproj$y,
#        z=matrix(ss[[s]], nxy[1])), asp=1,
#      legend.mar=3, main=names(ss)[s])
#  points(coo, pch=4, cex=0.5)

## ----predmapsview,echo=FALSE, pars=list(mfrow=c(4,1), mar=c(1.5,2,1,0.5), mgp=c(1,0.5,0)), fig.width=3.75, fig.height=11.7, out.width='0.49\\textwidth', fig.pos='h', fig.align='center'----
library(fields); for (s in c(1,3,4,2)) 
  image.plot(list(x=gproj$x, y=gproj$y, 
      z=matrix(ss[[s]], nxy[1])), asp=1,
    legend.mar=3, main=names(ss)[s]) 
points(coo, pch=4, cex=0.5)

