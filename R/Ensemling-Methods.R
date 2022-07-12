source("global_library.R")
library("gstat")

# read sky in Resolution 9:

map1 <- readFITS(file = "data/map-f1z1.fits")

map1cols <-map1$col
temps <- map1cols[]
I <- unlist(temps)

df_temp <- as.data.frame(I)

sky_rs9 <- as.CMBDataFrame(df_temp, nside = 512, ordering = "nested")


# get coordinates
sky_rs9_cart <- CMBDataFrame(sky_rs9, coords  = "cartesian")
sky_rs9_sph <- CMBDataFrame(sky_rs9, coords = "spherical")

# get window slice from full sky
# Define window as one pixel at a certain resolution of the sky
# for example: 1 window is one pixel at resolution 3
# that means the window contains 4086 pixels if we slice it out from
# the full sky of resolution 9
w1_r5 <- window(sky_rs9_cart, in.pixels = c(1), in.pixels.res = 5)
plot(w1_r5)

## Inference for single sky patch

# geostatistical approach
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



