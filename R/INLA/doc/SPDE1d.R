## ----setting,include=FALSE----------------------------------------------------
library(knitr);  library(INLA);  library(fields) 
set.seed(123)
if (file.exists("myinit.R")) source("myinit.R")
inla.setOption(num.threads="1:1")
inla.setOption(smtp="taucs")
opts_chunk$set(message=FALSE, warning=FALSE, tidy=FALSE,
fig.path="figures/SPDE1d/")
knit_hooks$set(small.mar = function(before, options, envir) {
    if (before) par(mar = c(0.1, 0.1, .1, .1))  # smaller margin on top and right
})
knit_hooks$set(mar3311 = function(before, options, envir) {
    if (before) par(mar = c(3, 3, 1, 1), mgp=c(2, 0.7, 0))  # smaller margin on top and right
})

## ----url----------------------------------------------------------------------
u0 <- paste0("http://www.yr.no/place/Norway/S%C3%B8r-Tr%C3%B8ndelag/",
             "Trondheim/Trondheim/detailed_statistics.html")
### browseURL(u0) ### to visit the web page 

## ----read,eval=FALSE----------------------------------------------------------
#  d0 <- readLines(u0) ### read it as text (Done at 26 September 2016)

## ----echo=FALSE,results="hide"------------------------------------------------
load("TRD20160926weather.RData")

## ----index--------------------------------------------------------------------
i <- grep("<tr>", d0) ### index for each table line
i <- i[i>grep("<tbody>", d0)[2]] ### select those for the second table

## ----data---------------------------------------------------------------------
if (Sys.getlocale("LC_TIME")!="C")
    Sys.setlocale("LC_TIME", "C")
dates <- as.Date(d0[i+1], format="      <th>%b %d, %Y</th>")
tmed <- as.numeric(gsub("<td>", "", gsub("..</td>", "", d0[i+4])))
(n <- length(dates)) ### it is daily over last 13 months

## ----visualize, fig.width=9, fig.height=3, out.width="0.99\\textwidth", fig.align="center", fig.pos="h"----
pd <- pretty(c(dates, max(dates+30)), n=13)
par(mfrow=c(1,1), mar=c(3,3,0.5,2), mgp=c(2,.7,0), las=2, xaxs="i")
plot(dates, tmed, type="l", lwd=2,
     axes=FALSE, xlab="day", ylab="Temperature")
abline(h=0)
abline(h=3*(-8:9), v=pd, lty=3, col=gray(.5))
box()
axis(2, 3*(-8:9)); axis(4, 3*(-8:9))
axis(1, pd, months(pd, TRUE))

## ----mesh---------------------------------------------------------------------
coo <- as.numeric(dates-min(dates)) ## have numeric temporal coordinates
mesh <- inla.mesh.1d(loc=seq(min(coo), max(coo), by=7), ## knots (7 days)
                 boundary=c("neumann", "neumann"), ## boundaries 
                 degree=2)  ### basis function degree

## ----projector, pars=list(mar=c(0,0,1,0)), fig.width=8, fig.height=2, out.width="0.99\\textwidth", fig.align="center", fig.pos="h"----
A <- inla.spde.make.A( ## projector creator
    mesh=mesh, ## provide the mesh
    loc=coo) ### locations where to project the field
dim(A) ## an "n" by "m" projector matrix
summary(rowSums(A)) ### each line sums up to one
summary(colSums(A)) ### "how many" observations per knot

## ----SPDE model---------------------------------------------------------------
spde <- inla.spde2.pcmatern( ## model for the precision 
    mesh=mesh, ## mesh supplied
    alpha=2, ## smoothness parameter
    prior.range = c(1, 0.01), ## P(range < 1) = 0.01
    prior.sigma = c(1, 0.5)) ## P(sigma > 1) = 0.5

## ----datastack----------------------------------------------------------------
stk.e <- inla.stack( ## stack creator
  data=list(y=tmed),  ## response
  effects=list(## two elements:
    data.frame(b0=rep(1, n)), ## regressor part
    i=1:spde$n.spde),  ## RF index
  A=list(## projector list of each effect
    1, ## for the covariates
    A), ## for the RF
  tag="est") ## tag

## ----fitting------------------------------------------------------------------
formula <- y ~ 0 + b0 + ## fixed part
  f(i, model=spde) ## RF term
res <- inla( ## main function in INLA package
  formula, ## model formula
  data=inla.stack.data(stk.e), ## dataset
  control.predictor=list( ## inform projector needed in SPDE models
    A = inla.stack.A(stk.e), compute=TRUE)) ## projector from the stack data

## ----fixed--------------------------------------------------------------------
round(res$summary.fixed, 4) 

## ----hypeparsummary-----------------------------------------------------------
round(res$summary.hyperpar, 4)

## ----nugget-------------------------------------------------------------------
m.prec <- res$marginals.hyperpar$"Precision for the Gaussian observations" ## the marginal
post.s2e <- inla.tmarginal(## function to compute a tranformation 
  function(x) sqrt(1/x), ## inverse transformation and square root
  m.prec) ## marginal to be applied

## ----parameters, fig.width=9, fig.height=3, out.width="0.99\\textwidth", fig.align="center", fig.pos="h"----
par(mfrow=c(1,3), mar=c(3,3,0.3,0.3), mgp=c(2,0.5,0))
plot(post.s2e, type="l", ylab="Density", 
     xlab=expression(sigma[e]^2))
plot(res$marginals.hyperpar[[2]], type="l",
     xlab="range", ylab="Density")
plot(res$marginals.hyperpar[[3]], ty="l", 
     xlab=expression(sigma[s]), yla="Density")

## ----predicted, fig.width=9, fig.height=3, out.width="0.99\\textwidth", fig.align="center", fig.pos="h"----
par(mfrow=c(1,1), mar=c(3,3,0.3,2), mgp=c(2,0.5,0), las=2, xaxs="i")
id <- inla.stack.index(stk.e, tag="est")$data
plot(dates, tmed, type="l", axes=FALSE, ylab="Temperature", lwd=2)
for (j in 3:5)
  lines(dates, res$summary.fitted.values[id, j], lty=2, lwd=2)
box(); axis(2, 3*(-8:9)); axis(4, 3*(-8:9))
axis(1, pd, months(pd, T))
abline(h=0)
abline(h=3*(-8:9), v=pd, lty=3, col=gray(.5))

## ----other--------------------------------------------------------------------
tmax <- as.numeric(gsub("<td>", "", gsub("..</td>", "", d0[i+2])))
tmin <- as.numeric(gsub("<td>", "", gsub("..</td>", "", d0[i+3])))
tnormal <- as.numeric(gsub("<td>", "", gsub("..</td>", "", d0[i+5])))
prec <- as.numeric(gsub("<td>", "", gsub("mm</td>", "", d0[i+6])))
wind <- as.numeric(gsub("<td>", "", gsub("m/s</td>", "", d0[i+10])))
wmax <- as.numeric(gsub("<td>", "", gsub("m/s</td>", "", d0[i+9])))

## ----visualize3, fig.width=9, fig.height=6, out.width="0.99\\textwidth", fig.align="center", fig.pos="h"----
par(mfrow=c(3,1), mar=c(0.1,3,0.1,2), mgp=c(2,.7,0), las=2, xaxs="i")
plot(dates, tmed, type="l", ylim=range(tmin, tmax, na.rm=TRUE), 
     axes=FALSE, xlab="", ylab="Temperature", col="green")
lines(dates, tmin, col="blue")
lines(dates, tmax, col="red")
lines(dates, tnormal)
legend(dates[which.min(tmin)], par()$usr[4], c("normal", "max.", "aver.", "min."), 
      col=1:4, lty=1, ncol=2, xjust=0.5, bty="n")
abline(h=5*(-5:6), v=pd, lty=3, col=gray(.5))
box(); axis(2, 5*(-5:6)); axis(4, 5*(-5:6))

plot(dates, prec, type="l", axes=FALSE, xlab="")
box(); axis(2); axis(4)
abline(v=pd, h=10*(1:4), lty=3, col=gray(0.5))

par(mar=c(3, 3, 0.1, 2), new=FALSE)
plot(dates, wind, type="l", axes=FALSE, xlab="",
     ylim=range(wind, wmax, na.rm=TRUE))
lines(dates, wmax, col=2)
box(); axis(2); axis(4)
abline(v=pd, h=5*(1:3), lty=3, col=gray(0.5))
axis(1, pd, months(pd, TRUE))

## ----anomalia, fig.width=9, fig.height=3, out.width="0.99\\textwidth", fig.align="center", fig.pos="h"----
par(mar=c(3, 3, 0.1, 2), mgp=c(2,0.7,0), las=2, xaxs="i")
plot(dates, tmed-tnormal, type="l", axes=FALSE, 
     xlab="", ylab="Deviation from normal temperature")
box(); axis(2); axis(4)
abline(h=5*(-2:2), v=pd, lty=2, col=gray(0.5))
axis(1, pd, months(pd, TRUE))

