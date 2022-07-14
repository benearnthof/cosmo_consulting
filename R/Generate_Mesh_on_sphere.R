source("global_library.R")
library("gstat")
library("fields")

# load map
map2 <- readFITS(file = "data/map_rs_3.fits")

# transform map into CMBDataFrame
map2cols <-map2$col
temps <- map2cols[]
I <- unlist(temps)

df_temp <- as.data.frame(I)

sky_rs2 <- as.CMBDataFrame(df_temp, nside = 8, ordering = "nested")


# following approach on Website chapter 7.3.3
# get map with cartesian coordinates
sky_rs2_cart <- CMBDataFrame(sky_rs2, coords  = "cartesian")

loc.cartesian = inla.mesh.map(loc.longlat, projection = "longlat")

# transform CMBDataFrame into DataFrame
df_sky_rs2 <- as.data.frame(sky_rs2_cart)
nrow(sky_rs2_cart)
# transform DataFrame into Matrix
mat_sky_r2 <-  as.matrix(df_sky_rs2)

# get maximal difference
max(diff(mat_sky_r2[,2]))

# generate spatial object
coordinates(df_sky_rs2) <- ~x+y+z

loc.cartesian <- as(df_sky_rs2, "SpatialPoints")

# bnd <- inla.nonconvex.hull(loc.cartesian, convex = 0.1)
# mesh <- inla.mesh.2d(boundary = bnd, cutoff = 0.05, max.edge = c(0.1))
# plot(mesh, rgl = TRUE)


# create inla.mesh
mesh_sky_r2 <- inla.mesh.2d(loc = loc.cartesian, max.edge = c(0.8843076)) #TODO: welchen cutoff, offset, max.edge verwenden?

plot(mesh_sky_r2)

#Mapping between meshes and continuous space:
#computes the sparse weight matrices needed to map between the internal 
#representation of weights for basis functions and the values of the resulting 
#functions and fields:

A = inla.spde.make.A(mesh_sky_r2, loc = loc.cartesian)


# available SPDE Models:
inla.spde.models()

spde = inla.spde2.matern(mesh_sky_r2, alpha = 2)


# Here, sigma0 is the feld standard deviation and range0 is the spatial range for theta = 0,
# and B.tau and B.kappa are matrices storing the parameter basis functions
sigma0 <- 1
size <- min(c(diff(range(mesh_sky_r2$loc[, 1])), diff(range(mesh_sky_r2$loc[, 2]))))
range0 <- size/5
kappa0 <- sqrt(8)/range0
tau0 <- 1/(sqrt(4 * pi) * kappa0 * sigma0)
spde <- inla.spde2.matern(mesh_sky_r2, B.tau = cbind(log(tau0), -1, +1),
                            B.kappa = cbind(log(kappa0), 0, -1), theta.prior.mean = c(0, 0),
                            theta.prior.prec = c(0.1, 1))

# Helper function to obtain precision of the constructed model, with standard deviation 
# a factor 3 larger than the prior median value and range equal to the prior median:
Q <- inla.spde.precision(spde, theta=c(log(3), 0))
image(Q)

# The following code then generates two samples from the model,
x <- inla.qsample(n = 2, Q)
x <- inla.qsample(n = 2, Q, constr = spde$f$extraconstr)


# Obtaining covariances is a much more costly operation, but the function 
COV_N <- inla.qinv(Q) 
image(COV_N)
# can quickly calculate all covariances between neighbours (in the Markov sense), including the
# marginal variances. Finally, 
inla.qsolve(Q,b) 
# uses the internal R-INLA methods for solving a linear system involving Q.



# the empirically derived range expression allows
#for construction of a model with known range and variance (= 1) for (theta1; theta2) = (0; 0), 
#TODO: Definiere ragnge0 heuristisch aus Daten-analyse mit rcosmo, Wie kann das noch erweitert werden
# Wie kan CAR(1)  (für d = 2, alpha = 1 on regular lattice) Prozess auf Nachbarn höherer Ordnung erweitert werden?
#FIND: Experimental helper functions for constructing
# parameterisations and priors are included in the package.

plot(mesh_sky_r2, rgl = TRUE)
lines(mesh_sky_r2$segm$bnd, mesh$loc, add = FALSE)


# For plotting fields defined on meshes, one option is to use the rgl option, which supports
# specifying colour values for the nodes, producing an interpolated field plot, and optionally
# draw triangle edges and vertices in the same plot:
plot(mesh_sky_r2, rgl = TRUE, col = x[,1], 
     color.palette = function(n) cividis(n, alpha = 1,begin =  0, end = 1),
     draw.edges = FALSE, draw.segments = TRUE, draw.vertices = FALSE)


# How to model that there is actually korrelation everywhere on the sphere?:
# Page 12 "Bayesian Spatial Modelling with R INLA"
# Models with range larger than the domain size are usually indistinguishable from intrinsic
# random fields, which can be modelled by fixing kappa to zero (or rather some small positive
# value) with B.tau = cbind(log(tau0), 1) and B.kappa = cbind(log(small), 0). Note
# that the sum-to-zero constraints often used for lattice based intrinsic Markov models is inappropriate
# due to the irregular mesh structure, and a weighted sum-to-zero constraint is
# needed to reproduce such models. The option constr = TRUE to the inla.spde.matern()
# call can be used to apply an integrate-to-zero constraint instead, which uses the triangle geometry
# to calculate the needed weights.

# Funktioniert noch nicht!
sigma0 <- 1
size <- min(c(diff(range(mesh_sky_r2$loc[, 1])), diff(range(mesh_sky_r2$loc[, 2]))))
kappa0_small <- sqrt(8)/size
tau0 <- 1/(sqrt(4 * pi) * kappa0_small * sigma0)
spde <- inla.spde2.matern(mesh_sky_r2, B.tau = cbind(log(tau0), 1),
                          B.kappa = cbind(log(kappa0_small), 0), theta.prior.mean = c(0, 0),
                          theta.prior.prec = c(0.1, 1), constr = TRUE)


proj2a = inla.mesh.projector(mesh_sky_r2, projection = "longlat", dims = c(361, 181))
proj2b = inla.mesh.projector(mesh_sky_r2, projection = "mollweide", dims = c(361, 181))


## Bayesian Inference: 
A = inla.spde.make.A(mesh_sky_r2,
                     loc = loc.cartesian,
                     index = rep(1:768, times = 1),
                     repl = rep(1, each = 768) )


mesh.index = inla.spde.make.index(name = "field",
                                  n.spde = spde$n.spde,
                                  n.repl = 1)
y <- c(sky_rs2_cart[,4])

st.est = inla.stack(data = list(y = y),
                    A = list(A),
                    effects = list(c(mesh.index, list(intercept = 1))),
                    tag = "est")


summary(st.est)
