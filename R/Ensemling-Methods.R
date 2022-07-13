source("global_library.R")
library("gstat")

# read sky in Resolution 9:

map1 <- readFITS(file = "data/map_rs_4.fits")

map1cols <-map1$col
temps <- map1cols[]
I <- unlist(temps)

df_temp <- as.data.frame(I)

sky_rs9 <- as.CMBDataFrame(df_temp, nside = 16, ordering = "nested")


# get window slice from full sky
# Define window as one pixel at a certain resolution of the sky
# for example: 1 window is one pixel at resolution 3
# that means the window contains 4086 pixels if we slice it out from
# the full sky of resolution 9
w1_r5 <- window(sky_rs9, in.pixels = c(1), in.pixels.res = 2)
plot(w1_r5)

## Inference for single sky patches

# I. map interpreted as areal data

## Plot Cov
varcmb <- variogramCMB(w1_r5, max.dist = 0.1, num.bins = 10)
ols <- variofitCMB(varcmb, fix.nug=FALSE, wei="equal", cov.model= "matern")
plot(varcmb)
lines(ols, lty=2)

# Definiere Struktur der Nachbarschaftsmatrix W für Auflösung 9 im Patch

res <- 4  # Auflösung full-sky
# n_side <- 2^(res)
# pixel <- 12*n_side^2
n_pixel <- nrow(w1_r5)
neighbourhood <- list()

for(i in 1:n_pixel){
  temp_neigh <- neighbours(i,res)
  neighbourhood[[i]] <- head(temp_neigh, -1)
}

vector <- rep(0,n_pixel)
temp_W <- replicate(n_pixel, vector, FALSE)

for(i in 1:n_pixel){
  temp_W[[i]][neighbourhood[[i]]] <- 1
}

W <- as.matrix(as.data.table(temp_W))
W <- W[1:n_pixel,]

Q <- (diag((W%*%rep(1,n_pixel))[,1])-W) - Diagonal(n_pixel)
g = inla.read.graph(Q)
Q_inla = inla.graph2matrix(g)


image(Q_inla)


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
dt_sky <-as.data.table(cbind(seq_len(nrow(w1_r5)), w1_r5))
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




# II. geostatistical approach (SPDE Approach)
# get coordinates
sky_rs9_cart <- CMBDataFrame(sky_rs9, coords  = "cartesian")
sky_rs9_sph <- CMBDataFrame(sky_rs9, coords = "spherical")

w1_r5_cart <- window(sky_rs9_cart, in.pixels = c(1), in.pixels.res = 5)
plot(w1_r5_cart)

df_w1_r5 <- as.data.frame(w1_r5)
df_sky_rs9 <- as.data.frame(sky_rs9_cart)


# df_w1_r5$x2 <- df_w1_r5$z/df_w1_r5$x
# df_w1_r5$y2 <- df_w1_r5$y/df_w1_r5$x
# df_w1_r5$x2 <- df_w1_r5$x2/max(df_w1_r5$x2)
# df_w1_r5$y2 <- df_w1_r5$y2/max(df_w1_r5$y2)
# df_w1_r5_red <-df_w1_r5[,-c(1:3)]

coordinates(df_w1_r5) <- ~x+y+z

proj4string(df_sky_rs9) <- CRS(as.character(NA))
gridded(df_sky_rs9) = TRUE

loc_w1_r5 <- as(df_w1_r5, "SpatialPoints")

mat_w1_r5 <- as.matrix(df_w1_r5)
mat_sky_r9 <-  as.matrix(df_sky_rs9)

mesh_w1_r5 <- inla.mesh.2d(loc = df_sky_rs9)

mesh_w1_r5 <- inla.mesh.2d(loc = df_w1_r5, max.edge = 0.01761868)

max(diff(mat_w1_r5[,2]))

plot(mesh_w1_r5,1/100)
lines(pts, col = 3, with = 2)


## Inference for whole Ensemble



