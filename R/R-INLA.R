source("global_library.R")


map1 <- readFITS(file = "data/map_rs_2.fits")

map1cols <-map1$col
temps <- map1cols[]
I <- unlist(temps)

df_temp <- as.data.frame(I)

sky <- as.CMBDataFrame(df_temp, nside = 4, ordering = "nested")
#sky <- ordering(sky, new.ordering = "ring")
# sky_cartesian <- CMBDataFrame(sky, coords  = "cartesian")
# sky_spherical <- CMBDataFrame(sky, coords = "spherical")
#CMB_temp_nested <- ordering(CMB_temp, new.ordering = "nested")
sky[31,]
#colnames(CMB_temp_nested) <- "I"
#colnames(sky_cartesian$unlist_temp) <- "I"

#plot(CMB_temp_nested) ...
plot(sky)

# plot windows slices
window(sky_spherical) <- CMBWindow(theta = c(pi/2,pi/2,pi/2.1, pi/2.1), phi = c(0,pi/100,pi/100,0))
plot(sky_spherical)


window(sky_cartesian) <- CMBWindow(x = 0, y= 3/5, z = 4/5 ,r = 0.8, set.minus = TRUE)
plot(sky_cartesian)

## Plot Cov
Cov <- covCMB(sky, max.dist = 3, num.bins = 10)
plot(Cov)
maxDist(sky)

# sub_sky <- sampleCMB(sky, sample.size = 100000)
# Cov <- covCMB(sub_sky, max.dist = 0.04, num.bins = 50)
# plot(Cov)
# maxDist(sub_sky)

df1 <- sampleCMB(sky, sample.size = 100000)
cov <- covCMB(df1, max.dist  = 0.03, num.bins = 20)
cov$v
plot(cov)


varcmb <- variogramCMB(df1, max.dist = 0.1, num.bins = 30)
ols <- variofitCMB(varcmb, fix.nug=FALSE, wei="equal", cov.model= "matern")
plot(varcmb)
lines(ols, lty=2)

# visualize Neighbourhood:
ns <- 512
rand <- rnorm(12 * ns ^ 2)
cmbdf <- CMBDataFrame(nside = ns, I = rand, ordering = "nested")
w1 <- window(cmbdf, in.pixels = c(1,9))
w2 <- window(cmbdf, in.pixels = c(2,10))
w3 <- window(cmbdf, in.pixels = c(4,12))
w1_sky <- window(sky, in.pixels = c(1,9))
w2_sky <- window(sky, in.pixels = c(2,10))
w3_sky <- window(sky, in.pixels = c(4,12))
plot(w1, col = "blue", back.col = "white", xlab = '', ylab = '', zlab = '')
plot(w2, col = "green", add = TRUE)
plot(w3, col = "orange", add = TRUE)
plot(w1_sky, col = "blue", back.col = "white", xlab = '', ylab = '', zlab = '')
plot(w2_sky, col = "green", add = TRUE)
plot(w3_sky, col = "orange", add = TRUE)
displayPixelBoundaries(nside = 2,ordering = "nested",incl.labels = 1:48,col ="red")


demoNeighbours <- function(p,j) {
  neighbours(p, j)
  displayPixels(boundary.j = j, j = j, plot.j = j+3,
                spix = neighbours(p, j),
                boundary.col = "gray",
                boundary.lwd = 1,
                incl.labels = neighbours(p, j),
                col = "blue",
                size = 3)
  rcosmo::displayPixelBoundaries(nside = 1, col = "blue", lwd = 3)
}

demoNeighbours(1,9)

# Definiere Struktur der Nachbarschaftsmatrix W für Auflösung 'res'
res <- 2
n_side <- 2^(res)
pixel <- 12*n_side^2
neighbourhood <- list()

for(i in seq_len(12*(n_side)^2)){
  temp_neigh <- neighbours(i,res)
  neighbourhood[[i]] <- head(temp_neigh, -1)
}

#neighbours(1,9)

vector <- rep(0,12*(n_side)^2)
temp_W <- replicate(12*(n_side)^2, vector, FALSE)

for(i in seq_len(12*(n_side)^2)){
  temp_W[[i]][neighbourhood[[i]]] <- 1
}

W <- as.matrix(as.data.table(temp_W)) # TODO: Die Präzisionsmatrix kann auch platzsparender  gespeichert werden wie in 'Bayesian Inference with INLA' beschrieben. Notwendig für höhere Auflösungen! 

# Um die Struktur visualiseren zu können wird Q ausgerechnet (geht aber nicht in die Berechnung # im inla Modell ein)
Q <- (diag((W%*%rep(1,pixel))[,1])-W) - Diagonal(pixel)
g = inla.read.graph(Q)
Q_inla = inla.graph2matrix(g)


image(Q_inla)


# Definiere Struktur der Nachbarschaftsmatrix W für ausgewähltes Fenster für Auflösung 'res'
res <- 4
n_side <- 2^(res)
pixel <- 12*n_side^2
neighbourhood <- list()


for(i in seq_len(nrow(sky_spherical))){
  temp_neigh <- neighbours(i,res)
  neighbourhood[[i]] <- head(temp_neigh, -1)
}

vector <- rep(0,nrow(sky_spherical))
temp_W <- replicate(nrow(sky_spherical), vector, FALSE)

for(i in seq_len(nrow(sky_spherical))){
  temp_W[[i]][neighbourhood[[i]]] <- 1
}

W <- as.matrix(as.data.table(temp_W)) # TODO: Die Präzisionsmatrix kann auch platzsparender  gespeichert werden wie in 'Bayesian Inference with INLA' beschrieben. Notwendig für höhere Auflösungen! 

# Um die Struktur visualiseren zu können wird Q ausgerechnet (geht aber nicht in die Berechnung # im inla Modell ein)
Q <- (diag((W%*%rep(1,pixel))[,1])-W) - Diagonal(pixel)
g = inla.read.graph(Q)
Q_inla = inla.graph2matrix(g)


image(Q_inla)



# Umordnung zur Bandmatrix
# g.file = inla.graph2matrix(g)
# inla.spy(g.file,  reordering = inla.qreordering(g,reordering = "band"))


# Folgende Simulation dient als Dummy, da Einlesen der Fits Datei noch nicht funktioniert!
# Simulation erfogt für Auflösung 1. Nachbarschaftsstruktur wird nicht! simuliert, ist also nur # eine sehr grobe Simulation!

# nLags = 48 # number of lags (size of region)
# # fake, uncorrelated observations
# set.seed(1)
# X = rlnorm(nLags)/10000
# ###############################################
# # fake sigma... correlated decreases distance.
# sigma = diag(nLags)
# corr = 0.05
# sigma <- corr ^ abs(row(sigma)-col(sigma))
# ###############################################
# # Y is autocorrelated...
# sky_sim <- t(X %*% chol(sigma))
# set.seed(1)
# sky_sim <- sky_sim + abs(rnorm(48,0,2)/10000)



# scale W matrix for being well defined: 
e.values <- eigen(W)$values
rho.min <- min(e.values)
rho.max <- max(e.values)
W <- W / rho.max

# this follows the implementation of a car modell in chapter 11.3 in Bayesian Infrence with 
# INLA.
interpret.theta <- function() {
  return(
    list(prec = exp(theta[1L]),
         rho = 1 / (1 + exp(-theta[2L])))
  )
}
graph <- function(){
  require(Matrix)
  
  return(Diagonal(nrow(W), x = 1) + W)
}


Q <- function() {
  require(Matrix)
  
  param <- interpret.theta()
  
  return(param$prec * (Diagonal(nrow(W), x = 1) - param$rho * W) ) # TODO: Erweiterung mit                                                                         #       theta_3
}


mu = function()
{
  return(numeric(0))
}


log.norm.const <- function() {
  param <- interpret.theta()
  n <- nrow(W)
  
  Q <- param$prec * (Diagonal(nrow(W), x = 1) - param$rho * W)
  
  res <- n * (-0.5 * log(2 * pi)) +
    0.5 * Matrix::determinant(Q, logarithm = TRUE)
  
  return(res)
}


log.prior <- function() {
  param = interpret.theta()
  
  res <- dgamma(param$prec, 1, 5e-05, log = TRUE) + log(param$prec) +
    log(1) + log(param$rho) + log(1 - param$rho)
  
  return(res)
}


initial <- function() {
  return(rep(0, 0))
}

quit <- function() {
  return(invisible())
}


'inla.rgeneric.CAR.model' <- function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
            "log.prior", "quit"),
    theta = NULL) {
  
  #Internal function
  interpret.theta <- function() {
    return(
      list(prec = exp(theta[1L]),
           rho = 1 / (1 + exp(-theta[2L])))
    )
  }
  
  graph <- function(){
    require(Matrix)
    
    return(Diagonal(nrow(W), x = 1) + W)
  }
  
  Q <- function() {
    require(Matrix)
    
    param <- interpret.theta()
    
    return(param$prec * (Diagonal(nrow(W), x = 1) - param$rho * W) )
  }
  
  mu <- function()
  {
    return(numeric(0))
  }
  
  log.norm.const <- function() {
    return(numeric(0))
    
  }
  
  log.prior <- function() {
    param = interpret.theta()
    
    res <- dgamma(param$prec, 1, 5e-05, log = TRUE) + log(param$prec) +
      log(1) + log(param$rho) + log(1 - param$rho) 
    
    return(res)
  }
  
  initial <- function() {
    return(c(0, 0))
  }
  
  quit <- function() {
    return(invisible())
  }
  
  res <- do.call(match.arg(cmd), args = list())
  return(res)
}

CAR.model <- inla.rgeneric.define(inla.rgeneric.CAR.model, W = W)

# add index for latent effect
dt_sky <-as.data.table(cbind(seq_len(nrow(sky)), sky))
colnames(dt_sky) <- c("id_x", "obs")
df_sky <- as.data.frame(dt_sky)
df_sky$obs <- df_sky$obs + 1 # wie soll mit den negativen Werten umgegangen werden?

# set car formular 
f.car <- obs ~ 1 + f(id_x, model = CAR.model)

# run inla
m.car <- inla(f.car, data = df_sky, family = "lognormal")

ggplot(data = df_sky, aes(x = id_x, y = obs)) + geom_line()

summary(m.car)


marg.prec <- inla.tmarginal(exp, m.car$marginals.hyperpar[[1]])
marg.rho <- inla.tmarginal(function(x) { 1/(1 + exp(-x))},
                           m.car$marginals.hyperpar[[2]])

df_marg.prec <- as.data.frame(marg.prec)
df_marg.rho <- as.data.frame(marg.rho)


ggplot(df_marg.prec) + 
  geom_line(aes(x = x, y = y)) +
  ylab (expression(paste(pi, "(", tau, " | ", bold(y), ")")))

ggplot(df_marg.rho) + 
  geom_line(aes(x = x, y = y)) +
  ylab (expression(paste(pi, "(", rho, " | ", bold(y), ")")))


ggplot(as.data.frame(m.car$marginals.hyperpar[[1]])) + 
  geom_line(aes(x = x, y = y)) +
  ylab (expression(paste(pi, "(", tau, " | ", bold(y), ")")))

ggplot(as.data.frame(m.car$marginals.hyperpar[[2]])) + 
  geom_line(aes(x = x, y = y)) +
  ylab (expression(paste(pi, "(", rho, " | ", bold(y), ")")))


# compare different fitting strategies
m.strategy <- lapply(c("gaussian", "simplified.laplace", "laplace"), 
                     function(st) {
                       return(lapply(c("ccd", "grid", "eb"), function(int.st) {
                         inla(f.car, data = df_sky,family = "lognormal",
                              control.inla = list(strategy = st, int.strategy = int.st),
                              control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE))
                       }))
                     })

## get table with metrics for diferent fitting strategies

metrics_table <- data.frame(matrix(data = 0, ncol = 4))
colnames(metrics_table) <- c("model", "dic", "waic", "cpo")
names <- c("gaussian", "simplified.laplace", "laplace")
op <- c("ccd", "grid", "eb")
num <- 0
for(i in seq_len(length(c("gaussian", "simplified.laplace", "laplace")))){
  for(j in seq_len(length(c("ccd", "grid", "eb")))){
    num <- num + 1
    dic_temp <- m.strategy[[i]][[j]]$dic$dic
    waic_temp <- m.strategy[[i]][[j]]$waic$waic
    cpo_temp <- -sum(log(m.strategy[[i]][[j]]$cpo$cpo))
    metrics_table[num,1] <-paste0(names[i]," with ", op[j])
    metrics_table[num,2] <- dic_temp
    metrics_table[num,3] <- waic_temp
    metrics_table[num,4] <- cpo_temp
  }
}

# set constants for table
n = 0
c_col = c("#1e3048", "#274060", "#2f5375", "#4073a0", "#5088b9")
c_col_light_blue = c("#edf2fb", "#e2eafc", "#d7e3fc", "#ccdbfd", "#c1d3fe")
c_container_width = px(800)
c_table_width = px(650)
c_rn = 30
c_save = TRUE
c_format = "png"

# get table with metrics
metrics_table %>% gt()


#Darstellung der posteriori marginals der unterschiedlichen Strategien in Modell Repräsentation
  

num <- 0
post_marg.rho_list <- list()
post_marg.prec_list <- list()
colnames(metrics_table) <- c("model", "dic", "waic", "cpo")
names <- c("gaussian", "simplified.laplace", "laplace")
op <- c("ccd", "grid", "eb")
# compare marginal posteriors of hyperparameters
for(i in seq_len(length(c("gaussian", "simplified.laplace", "laplace")))){
  for(j in seq_len(length(c("ccd", "grid", "eb")))){
    num <- num + 1
    # marg.prec_temp <- inla.tmarginal(exp, m.strategy[[i]][[j]]$marginals.hyperpar[[1]])
    # marg.rho_temp <- inla.tmarginal(function(x) { 1/(1 + exp(-x))},
    #                                 m.strategy[[i]][[j]]$marginals.hyperpar[[2]])
    
    
    marg.prec_temp <- m.strategy[[i]][[j]]$marginals.hyperpar[[1]]
    marg.rho_temp <- m.strategy[[i]][[j]]$marginals.hyperpar[[2]]
    
    df_marg.rho_temp <- as.data.frame(marg.rho_temp)
    df_marg.prec_temp <- as.data.frame(marg.prec_temp)
    
    df_marg.rho_temp <- cbind(df_marg.rho_temp, fit.st = names[i], int.st = op[j])
    df_marg.prec_temp <- cbind(df_marg.prec_temp, fit.st = names[i], int.st = op[j])
    
    post_marg.rho_list[[num]] <- df_marg.rho_temp
    post_marg.prec_list[[num]] <- df_marg.prec_temp
  }
}

df_post_marg.rho <- rbindlist(post_marg.rho_list)
df_post_marg.prec <- rbindlist(post_marg.prec_list)


ggplot(data = df_post_marg.rho, aes(x = x, y = y)) + geom_line(aes(linetype=fit.st, color=int.st)) + ylab("marginal post") +
  xlab("rho")

ggplot(data = df_post_marg.prec, aes(x = x, y = y)) + geom_line(aes(linetype=fit.st, color=int.st)) + ylab("marginal post") +
  xlab("prec")
