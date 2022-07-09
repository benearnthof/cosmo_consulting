## ----setup, include=FALSE-----------------------------------------------------
library(INLA)
set.seed(123)
if (file.exists("myinit.R")) source("myinit.R")
inla.setOption(num.threads="1:1")
inla.setOption(smtp="taucs")
library(knitr)
library(rmarkdown)
knitr::opts_chunk$set(echo=TRUE, cache=FALSE, message=FALSE, warning=FALSE)
knitr::opts_chunk$set(fig.path="figures/rgeneric/")

## ----eval=FALSE---------------------------------------------------------------
#  model = inla.rgeneric.define(rmodel, ...)

## ----eval=FALSE---------------------------------------------------------------
#  y ~ ... + f(idx, model=model, ...)

## ----eval=FALSE---------------------------------------------------------------
#  inla.rgeneric.ar1.model = function(
#          cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
#                  "log.prior", "quit"),
#          theta = NULL)
#  {
#      # for reference and potential storage for objects to
#      # cache, this is the environment of this function
#      # which holds arguments passed as `...` in
#      # `inla.rgeneric.define()`.
#  	envir = parent.env(environment())
#  
#      graph = function(){ <to be completed> }
#      Q = function() { <to be completed> }
#      mu = function() { <to be completed> }
#      log.norm.const = function() { <to be completed> }
#      log.prior = function() { <to be completed> }
#      initial = function() { <to be completed> }
#      quit = function() { <to be completed> }
#  
#      # sometimes this is useful, as argument 'graph' and 'quit'
#      # will pass theta=numeric(0) (or NULL in R-3.6...) as
#  	# the values of theta are NOT
#  	# required for defining the graph. however, this statement
#      # will ensure that theta is always defined.
#      if (!length(theta)) theta = initial()
#  
#      val = do.call(match.arg(cmd), args = list())
#      return (val)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  model = inla.rgeneric.define(inla.rgeneric.ar1.model, n = 100)

## ----eval=FALSE---------------------------------------------------------------
#  interpret.theta = function() {
#      return(list(prec = exp(theta[1L]),
#                  rho = 2 * exp(theta[2L])/(1 + exp(theta[2L])) - 1))
#  }

## ----eval=FALSE---------------------------------------------------------------
#  graph = function() {
#      return (Q())
#  }

## ----eval=FALSE---------------------------------------------------------------
#  Q = function() {
#      p = interpret.theta()
#      Q = p$prec/(1 - p$rho^2) *
#  	    toeplitz(c(1 + p$rho^2, -p$rho, rep(0, n - 2L)))
#      Q[1, 1] = Q[n, n] = p$prec/(1 - p$rho^2)
#      return (inla.as.sparse(Q))
#  }

## ----eval=FALSE---------------------------------------------------------------
#  Q = function() {
#      p = interpret.theta()
#      i = c(1L, n, 2L:(n - 1L), 1L:(n - 1L))
#      j = c(1L, n, 2L:(n - 1L), 2L:n)
#      x = p$prec/(1 - p$rho^2) *
#  	    c(1L, 1L, rep(1 + p$rho^2, n - 2L),
#  	      rep(-p$rho, n - 1L))
#       return (sparseMatrix(i = i, j = j, x = x, giveCsparse = FALSE))
#  }

## ----eval=FALSE---------------------------------------------------------------
#  mu = function() {
#      return(numeric(0))
#  }

## ----eval=FALSE---------------------------------------------------------------
#  log.norm.const = function() {
#     p = interpret.theta()
#     prec.i  = p$prec / (1.0 - p$rho^2)
#     val = n * (- 0.5 * log(2*pi) + 0.5 * log(prec.i)) +
#           0.5 * log(1.0 - p$rho^2)
#     return (val)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  log.norm.const = function() {
#     return (numeric(0))
#  }

## ----eval=FALSE---------------------------------------------------------------
#  log.prior = function() {
#      p = interpret.theta()
#      val = dgamma(p$prec, shape = 1, rate = 1, log=TRUE) + theta[1L] +
#            dnorm(theta[2L], mean = 0, sd = 1, log=TRUE)
#      return (val)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  initial = function() {
#      return (rep(1, 2))
#  }

## ----eval=FALSE---------------------------------------------------------------
#  quit = function() {
#      return (invisible())
#  }

## -----------------------------------------------------------------------------
n = 100
rho=0.9
x = arima.sim(n, model = list(ar = rho)) * sqrt(1-rho^2)
y = x + rnorm(n, sd = 0.1)
model = inla.rgeneric.define(inla.rgeneric.ar1.model, n=n)
formula = y ~ -1 + f(idx, model=model)
r = inla(formula, data = data.frame(y, idx = 1:n))

## -----------------------------------------------------------------------------
fformula = y ~ -1 +
	f(idx, model = "ar1",
	  hyper = list(prec = list(prior = "loggamma", param = c(1,1)),
      rho = list(prior = "normal", param = c(0,1))))
rr = inla(fformula, data = data.frame(y, idx = 1:n))

## -----------------------------------------------------------------------------
plot(inla.smarginal(rr$marginals.hyperpar[[2]]), 
	 type="l", lwd=5, col="red", xlab="stdev", ylab="density")
lines(inla.tmarginal(exp, r$internal.marginals.hyperpar[[2]]), 
     col="yellow")

## -----------------------------------------------------------------------------
plot(inla.smarginal(rr$marginals.hyperpar[[3]]), 
	 type="l", lwd=5, col="red", xlab="rho", ylab="density")
lines(inla.tmarginal(function(x) 2*exp(x)/(1+exp(x))-1,
                     r$internal.marginals.hyperpar[[3]]), 
     col="yellow")

## -----------------------------------------------------------------------------
round(rbind(native = rr$cpu.used,  
            rgeneric = r$cpu.used), digits = 3)

## -----------------------------------------------------------------------------
inla.rgeneric.iid.model

## ---- eval=FALSE--------------------------------------------------------------
#  ## In this example we do linear regression using 'rgeneric'.
#  ## The regression model is y = a + b*x + noise,  and we
#  ## define 'a + b*x + tiny.noise' as a latent model.
#  ## The dimension is length(x) and number of hyperparameters
#  ## is 2 ('a' and 'b').

## -----------------------------------------------------------------------------
rgeneric.linear.regression =
    function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", 
                     "log.prior", "quit"),
             theta = NULL)
{
	envir = parent.env(environment())

    ## artifical high precision to be added to the mean-model
    prec.high = exp(15)
    
    interpret.theta = function() {
        return(list(a = theta[1L], b = theta[2L]))
    }
    
    graph = function() {
        G = Diagonal(n = length(x), x=1)
        return(G)
    } 
    
    Q = function() {
        Q = prec.high * graph()
        return(Q)
    }
    
    mu = function() {
        par = interpret.theta()
        return(par$a + par$b * x)
    }

    log.norm.const = function() {
        return(numeric(0))
    }

    log.prior = function() {
        par = interpret.theta()
        val = (dnorm(par$a, mean=0, sd=1, log=TRUE) +
               dnorm(par$b, mean=0, sd=1, log=TRUE))
        return(val)
    }

    initial = function() {
        return(rep(0, 2))
    }
   
    quit = function() {
        return(invisible())
    }

    val = do.call(match.arg(cmd), args = list())
    return(val)
}

## -----------------------------------------------------------------------------
a = 1
b = 2
n = 50
x = rnorm(n)
eta = a + b*x
s = 0.25
y = eta + rnorm(n, sd=s)

rgen = inla.rgeneric.define(model = rgeneric.linear.regression, x=x)
r = inla(y ~ -1 + f(idx, model=rgen),
         data = data.frame(y, idx = 1:n))
rr = inla(y ~ 1 + x,
          data = data.frame(y, x),
          control.fixed = list(prec.intercept = 1, prec = 1))

## -----------------------------------------------------------------------------
plot(inla.smarginal(r$marginals.hyperpar[['Theta1 for idx']]),
     type="l", lwd=5, col="red", xlab="Intercept", ylab="density")
lines(inla.smarginal(rr$marginals.fixed$'(Intercept)'), col="yellow")

## -----------------------------------------------------------------------------
plot(inla.smarginal(r$marginals.hyperpar[['Theta2 for idx']]),
     type="l", lwd=5, col="red", xlab="Slope", ylab="density")
lines(inla.smarginal(rr$marginals.fixed$x), col="yellow")

## -----------------------------------------------------------------------------
r = inla.hyperpar(r)

## -----------------------------------------------------------------------------
plot(inla.smarginal(r$marginals.hyperpar[['Theta1 for idx']]),
     type="l", lwd=5, col="red", xlab="Intercept", ylab="density")
lines(inla.smarginal(rr$marginals.fixed$'(Intercept)'), col="yellow")

## -----------------------------------------------------------------------------
plot(inla.smarginal(r$marginals.hyperpar[['Theta2 for idx']]),
     type="l", lwd=5, col="red", xlab="Slope", ylab="density")
lines(inla.smarginal(rr$marginals.fixed$x), col="yellow")

## -----------------------------------------------------------------------------
rgeneric.test = function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"),
    theta = NULL)
{
    envir = parent.env(environment())

    graph = function() {
        return(matrix(1, n, n))
    }

    Q = function() {
        R <- matrix(sin(1:n^2), n, n)
        R <- R %*% t(R)
        diag(R) <- diag(R)+1
        Q <- exp(theta[1]) * R
        return(Q)
    }

    mu = function() return (numeric(0))

    log.norm.const = function() {
        return (numeric(0))
    }

    log.prior = function() {
        return (dgamma(exp(theta[1]), shape = 1, rate = 1, log=TRUE) + theta[1])
    }

    initial = function() {
        return(4)
    }

    if (!length(theta)) theta = initial()
    val = do.call(match.arg(cmd), args = list())

    return (val)
}

## -----------------------------------------------------------------------------
n = 200
s = .1
Q <- rgeneric.test("Q", theta = 0)
library(mvtnorm)
S <- solve(as.matrix(Q))
S <- (S + t(S))/2
x <- drop(rmvnorm(1, sigma = S))
y <- x + rnorm(n, sd = s)
cont.family = list(hyper = list(prec = list(initial=log(1/s^2), fixed=TRUE)))

r1 = inla(y ~ -1 + f(idx, model="generic", Cmatrix = Q,
                     hyper = list(prec = list(prior = "loggamma", param = c(1, 1)))), 
          data = data.frame(y = y, idx = 1:n), control.family = cont.family)
ld <- 0.5 * log(det(as.matrix(Q)))
r1$mlik <- r1$mlik + ld ## see the documentation for why

model2 = inla.rgeneric.define(rgeneric.test, n=n, optimize = FALSE)
r2 = inla(y ~ -1 + f(idx, model=model2), 
          data = data.frame(y = y, idx = 1:n), control.family = cont.family)

## -----------------------------------------------------------------------------
r2$mlik - r1$mlik

## -----------------------------------------------------------------------------
rgeneric.test.opt.1 = function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"),
    theta = NULL)
{
    envir = parent.env(environment())

    if (!exists("cache.done", envir = envir)) {
        R <- matrix(sin(1:n^2), n, n)
        R <- R %*% t(R)
        diag(R) <- diag(R)+1
        R.logdet <- log(det(R))
        R <- inla.as.sparse(R)
        idx <- which(R@i <= R@j)
        R@i <- R@i[idx]
        R@j <- R@j[idx]
        R@x <- R@x[idx]
        assign("R", R, envir = envir)
        norm.const <- -n/2 * log(2*pi) + 0.5 * R.logdet
        assign("norm.const", norm.const, envir = envir)
        assign("cache.done", TRUE, envir = envir)
    }

    graph = function() {
        return (R)
    }

    Q = function() {
        return(exp(theta[1]) * R)
    }

    mu = function() return (numeric(0))

    log.norm.const = function() {
        return (norm.const + n/2 * theta[1])
    }

    log.prior = function() {
        return (dgamma(exp(theta[1]), shape = 1, rate = 1, log=TRUE) + theta[1])
    }

    initial = function() {
        return(4)
    }

    if (!length(theta)) theta = initial()
    val = do.call(match.arg(cmd), args = list())

    return (val)
}

## -----------------------------------------------------------------------------
A=matrix(1:9,3,3)
A
inla.as.sparse(A)@x

## -----------------------------------------------------------------------------
rgeneric.test.opt.2 = function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"),
    theta = NULL)
{
    envir = parent.env(environment())

    if (!exists("cache.done", envir = envir)) {
        R <- matrix(sin(1:n^2), n, n)
        R <- R %*% t(R)
        diag(R) <- diag(R)+1
        R.logdet <- log(det(R))
        R <- inla.as.sparse(R)
        idx <- which(R@i <= R@j)
        R@i <- R@i[idx]
        R@j <- R@j[idx]
        R@x <- R@x[idx]
        assign("R", R, envir = envir)
        norm.const <- -n/2 * log(2*pi) + 0.5 * R.logdet
        assign("norm.const", norm.const, envir = envir)
        assign("cache.done", TRUE, envir = envir)
    }

    graph = function() {
        return (R)
    }

    Q = function() {
        ## since R was created with 'inla.sparse.matrix' above, the indices are sorted in a
        ## spesific order. This ordering is REQUIRED for R@x to be interpreted correctly.
        return(exp(theta[1]) * R@x)
    }

    mu = function() return (numeric(0))

    log.norm.const = function() {
        return (norm.const + n/2 * theta[1])
    }

    log.prior = function() {
        return (dgamma(exp(theta[1]), shape = 1, rate = 1, log=TRUE) + theta[1])
    }

    initial = function() {
        return(4)
    }

    if (!length(theta)) theta = initial()
    val = do.call(match.arg(cmd), args = list())

    return (val)
}

## -----------------------------------------------------------------------------
model3 = inla.rgeneric.define(rgeneric.test.opt.1, n=n, optimize = FALSE)
r3 = inla(y ~ -1 + f(idx, model=model3), 
          data = data.frame(y = y, idx = 1:n), control.family = cont.family)

model4 = inla.rgeneric.define(rgeneric.test.opt.2, n=n, optimize = TRUE)
r4 = inla(y ~ -1 + f(idx, model=model4), 
          data = data.frame(y = y, idx = 1:n), control.family = cont.family)

## -----------------------------------------------------------------------------
print(r2$mlik - r1$mlik)
print(r3$mlik - r1$mlik)
print(r4$mlik - r1$mlik)

print(rbind(native = r1$cpu[2], 
            rgeneric.plain = r2$cpu[2], 
            rgeneric.cache = r3$cpu[2], 
            rgeneric.optimze = r4$cpu[2]))

## ----eval=FALSE---------------------------------------------------------------
#  gcc -Wall -fpic -g -O -c -o cgeneric-demo.o cgeneric-demo.c
#  gcc -shared -o cgeneric-demo.so cgeneric-demo.o

## ----eval=FALSE---------------------------------------------------------------
#  cmodel <- inla.cgeneric.define(model = "inla_cgeneric_iid_model",
#                                 shlib = "cgeneric-demo.so", n = n)

## ----eval=FALSE---------------------------------------------------------------
#  rc <- inla(
#      y ~ -1 + f(idx, model = cmodel),
#      data = data.frame(y, idx = 1:n),
#      control.family = list(hyper = list(prec = list(initial = 12, fixed = TRUE))))

## ----eval=FALSE---------------------------------------------------------------
#  n <- 100
#  y <- rnorm(n)
#  
#  r <- inla(
#      y ~ -1 + f(idx, model = "iid", hyper = list(prec = list(param = c(1, 1)))),
#      data = data.frame(y, idx = 1:n),
#      control.family = list(hyper = list(prec = list(initial = 12, fixed = TRUE))))
#  
#  rmodel <- inla.rgeneric.define(inla.rgeneric.iid.model, n = n)
#  rr <- inla(
#      y ~ -1 + f(idx, model = rmodel),
#      data = data.frame(y, idx = 1:n),
#      control.family = list(hyper = list(prec = list(initial = 12, fixed = TRUE))))
#  
#  cmodel <- inla.cgeneric.define(model = "inla_cgeneric_iid_model",
#                                 shlib = "cgeneric-demo.so", n = n)
#  rc <- inla(
#      y ~ -1 + f(idx, model = cmodel),
#      data = data.frame(y, idx = 1:n),
#      control.family = list(hyper = list(prec = list(initial = 12, fixed = TRUE))))
#  
#  print(cbind(r$mlik, rr$mlik-r$mlik, rc$mlik-r$mlik))
#  print(cbind(r$cpu[2], rr$cpu[2], rc$cpu[2]))

## ----eval=FALSE---------------------------------------------------------------
#  gcc -Wall -fpic -g3 -O   -c -o cgeneric-demo.o cgeneric-demo.c
#  gcc -shared -o cgeneric-demo.so cgeneric-demo.o
#  
#  > n <- 5
#  > cmodel <- inla.cgeneric.define(model = "inla_cgeneric_ar1_model",
#  +                                shlib = "cgeneric-demo.so", n = n)
#  > inla.cgeneric.q(cmodel)
#  $theta
#  [1] 1 1
#  
#  $graph
#  5 x 5 sparse Matrix of class "dgCMatrix"
#  
#  [1,] 1 1 . . .
#  [2,] 1 1 1 . .
#  [3,] . 1 1 1 .
#  [4,] . . 1 1 1
#  [5,] . . . 1 1
#  
#  $Q
#  5 x 5 sparse Matrix of class "dgCMatrix"
#  
#  [1,]  3.45640494 -1.59726402  .           .           .
#  [2,] -1.59726402  4.19452805 -1.59726402  .           .
#  [3,]  .          -1.59726402  4.19452805 -1.59726402  .
#  [4,]  .           .          -1.59726402  4.19452805 -1.59726402
#  [5,]  .           .           .          -1.59726402  3.45640494
#  
#  $mu
#  numeric(0)
#  
#  $log.prior
#  [1] -3.13722036
#  
#  $log.norm.const
#  [1] -1.61423464

## ----eval=FALSE---------------------------------------------------------------
#  typedef enum {
#  	INLA_CGENERIC_VOID = 0,
#  	INLA_CGENERIC_Q,
#  	INLA_CGENERIC_GRAPH,
#  	INLA_CGENERIC_MU,
#  	INLA_CGENERIC_INITIAL,
#  	INLA_CGENERIC_LOG_NORM_CONST,
#  	INLA_CGENERIC_LOG_PRIOR,
#  	INLA_CGENERIC_QUIT
#  } inla_cgeneric_cmd_tp;

## ----eval=FALSE---------------------------------------------------------------
#  #include "cgeneric.h"
#  #define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))
#  #define SQR(x) ((x)*(x))
#  
#  double *inla_cgeneric_iid_model(inla_cgeneric_cmd_tp cmd,
#  	double *theta, inla_cgeneric_data_tp * data)
#  {
#  	// this reimplement `inla.rgeneric.iid.model` using cgeneric
#  
#  	double *ret = NULL, prec = (theta ? exp(theta[0]) : NAN),
#  		lprec = (theta ? theta[0] : NAN);
#  
#  	assert(!strcasecmp(data->ints[0]->name, "n")); // this will always be the case
#  	int N = data->ints[0]->ints[0];		       // this will always be the case
#  	assert(N > 0);
#  
#  	switch (cmd) {
#  	case INLA_CGENERIC_VOID:
#  	{
#  		assert(!(cmd == INLA_CGENERIC_VOID));
#  		break;
#  	}
#  
#  	case INLA_CGENERIC_GRAPH:
#  	{
#  		ret = Calloc(2 + 2 * N, double);
#  		ret[0] = N;				       /* dimension */
#  		ret[1] = N;				       /* number of (i <= j) */
#  		for (int i = 0; i < N; i++) {
#  			ret[2 + i] = i;			       /* i */
#  			ret[2 + N + i] = i;		       /* j */
#  		}
#  		break;
#  	}
#  
#  	case INLA_CGENERIC_Q:
#  	{
#  		if (1) {
#  			// optimized format
#  			ret = Calloc(2 + N, double);
#  			ret[0] = -1;			       /* code for optimized output */
#  			ret[1] = N;			       /* number of (i <= j) */
#  			for (int i = 0; i < N; i++) {
#  				ret[2 + i] = prec;
#  			}
#  		} else {
#  			// plain format, but the optimized format above is better to use
#  			ret = Calloc(2 + 3 * N, double);
#  			ret[0] = N;
#  			ret[1] = N;
#  			for (int i = 0; i < N; i++) {
#  				ret[2 + i] = i;		       /* i */
#  				ret[2 + N + i] = i;	       /* j */
#  				ret[2 + 2 * N + i] = prec;     /* Q_ij */
#  			}
#  		}
#  		break;
#  	}
#  
#  	case INLA_CGENERIC_MU:
#  	{
#  		ret = Calloc(1, double);
#  		ret[0] = 0;
#  		break;
#  	}
#  
#  	case INLA_CGENERIC_INITIAL:
#  	{
#  		ret = Calloc(2, double);
#  		ret[0] = 1;
#  		ret[1] = 4.0;
#  		break;
#  	}
#  
#  	case INLA_CGENERIC_LOG_NORM_CONST:
#  	{
#  		ret = Calloc(1, double);
#  		ret[0] = N * (-0.9189385332 + 0.5 * lprec);
#  		break;
#  	}
#  
#  	case INLA_CGENERIC_LOG_PRIOR:
#  	{
#  		// prec ~ gamma(1,1)
#  		ret = Calloc(1, double);
#  		ret[0] = -prec + lprec;
#  		break;
#  	}
#  
#  	case INLA_CGENERIC_QUIT:
#  	default:
#  		break;
#  	}
#  
#  	return (ret);
#  }

## ----eval=FALSE---------------------------------------------------------------
#  #include "cgeneric.h"
#  #define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))
#  #define SQR(x) ((x)*(x))
#  
#  double *inla_cgeneric_ar1_model(inla_cgeneric_cmd_tp cmd, double *theta,
#  	inla_cgeneric_data_tp * data)
#  {
#  	// this reimplement `inla.rgeneric.ar1.model` using cgeneric
#  
#  	double *ret = NULL, prec, lprec, rho, rho_intern;
#  
#  	if (theta) {
#  		lprec = theta[0];
#  		prec = exp(lprec);
#  		rho_intern = theta[1];
#  		rho = 2.0 * exp(rho_intern) / (1.0 + exp(rho_intern)) - 1.0;
#  	} else {
#  		prec = lprec = rho = rho_intern = NAN;
#  	}
#  
#  	assert(!strcasecmp(data->ints[0]->name, "n")); // this will always be the case
#  	int N = data->ints[0]->ints[0];		       // this will always be the case
#  	assert(N > 0);
#  	
#  	switch (cmd) {
#  	case INLA_CGENERIC_VOID:
#  	{
#  		assert(!(cmd == INLA_CGENERIC_VOID));
#  		break;
#  	}
#  
#  	case INLA_CGENERIC_GRAPH:
#  	{
#  		int m = N + N - 1, offset, i, k;
#  		ret = Calloc(2 + 2 * m, double);
#  
#  		offset = 2;
#  		ret[0] = N;				       /* dimension */
#  		ret[1] = m;				       /* number of (i <= j) */
#  		for (k = i = 0; i < N; i++) {
#  			ret[offset + k] = i;		       /* i */
#  			ret[offset + m + k++] = i;	       /* j */
#  			if (i < N - 1) {
#  				ret[offset + k] = i;	       /* i */
#  				ret[offset + m + k++] = i + 1; /* j */
#  			}
#  		}
#  		break;
#  	}
#  
#  	case INLA_CGENERIC_Q:
#  	{
#  		double param = prec / (1.0 - SQR(rho));
#  		int m = N + N - 1;
#  		int offset, i, k;
#  		ret = Calloc(2 + m, double);
#  
#  		// use optimized format.
#  		// The order of Q_ij s are then predetermined as the upper triangular of Q:
#  		//
#  		// for(i=0; i < n; i++)
#  		//     for(j=i; j<n; j++)
#  		//        ...
#  		//
#  		// but for only those (i,j)s that is defined in _GRAPH, of course
#  
#  		offset = 2;
#  		ret[0] = -1;
#  		ret[1] = m;
#  		for (i = k = 0; i < N; i++) {
#  			ret[offset + k++] = param * (i == 0 || i == N - 1 ? 1.0 : (1.0 + SQR(rho)));
#  			if (i < N - 1) {
#  				ret[offset + k++] = -param * rho;
#  			}
#  		}
#  		break;
#  	}
#  
#  	case INLA_CGENERIC_MU:
#  	{
#  		ret = Calloc(1, double);
#  		ret[0] = 0;
#  		break;
#  	}
#  
#  	case INLA_CGENERIC_INITIAL:
#  	{
#  		ret = Calloc(3, double);
#  		ret[0] = 2;
#  		ret[1] = 1.0;
#  		ret[2] = 1.0;
#  		break;
#  	}
#  
#  	case INLA_CGENERIC_LOG_NORM_CONST:
#  	{
#  		double prec_innovation = prec / (1.0 - SQR(rho));
#  		ret = Calloc(1, double);
#  		ret[0] = N * (-0.5 * log(2.0 * M_PI) +
#  			0.5 * log(prec_innovation)) + 0.5 * log(1.0 - SQR(rho));
#  		break;
#  	}
#  
#  	case INLA_CGENERIC_LOG_PRIOR:
#  	{
#  		ret = Calloc(1, double);
#  		ret[0] = -prec + lprec - 0.5 * log(2.0 * M_PI) - 0.5 * SQR(rho_intern);
#  		break;
#  	}
#  
#  	case INLA_CGENERIC_QUIT:
#  	default:
#  		break;
#  	}
#  
#  	return (ret);
#  }

## ----eval=FALSE---------------------------------------------------------------
#  typedef struct {
#  	int n_ints;
#  	inla_cgeneric_vec_tp **ints;
#  
#  	int n_doubles;
#  	inla_cgeneric_vec_tp **doubles;
#  
#  	int n_chars;
#  	inla_cgeneric_vec_tp **chars;
#  
#  	int n_mat;
#  	inla_cgeneric_mat_tp **mats;
#  
#  	int n_smat;
#  	inla_cgeneric_smat_tp **smats;
#  } inla_cgeneric_data_tp;

## ----eval=FALSE---------------------------------------------------------------
#  	assert(!strcasecmp(data->ints[0]->name, "n")); // this will always be the case
#  	int N = data->ints[0]->ints[0];		       // this will always be the case

## ----eval=FALSE---------------------------------------------------------------
#  typedef struct
#  {
#  	char *name;
#  	int len;
#  	int *ints;
#  	double *doubles;
#  	char *chars;
#  }
#  	inla_cgeneric_vec_tp;

## ----eval=FALSE---------------------------------------------------------------
#  /*
#   *      matrix storage, stored column by column, like
#   *      > matrix(1:6,2,3)
#   *      [,1] [,2] [,3]
#   *      [1,]    1    3    5
#   *      [2,]    2    4    6
#   *      > c(matrix(1:6,2,3))
#   *      [1] 1 2 3 4 5 6
#   */
#  typedef struct {
#  	char *name;
#  	int nrow;
#  	int ncol;
#  	double *x;
#  } inla_cgeneric_mat_tp;

## ----eval=FALSE---------------------------------------------------------------
#  /*
#   * sparse matrix format, stored used 0-based indices, like
#   *
#   *      > A <- inla.as.sparse(matrix(c(1,2,3,0,0,6),2,3))
#   *      > A
#   *      2 x 3 sparse Matrix of class "dgTMatrix"
#   *      [1,] 1 3 .
#   *      [2,] 2 . 6
#   *      > cbind(i=A@i, j=A@j, x=A@x)
#   *            i j x
#   *      [1,]  0 0 1
#   *      [2,]  1 0 2
#   *      [3,]  0 1 3
#   *      [4,]  1 2 6
#   */
#  typedef struct {
#  	char *name;
#  	int nrow;
#  	int ncol;
#  	int n;						       /* number of triplets (i,j,x) */
#  	int *i;
#  	int *j;
#  	double *x;
#  } inla_cgeneric_smat_tp;

