source("global_library.R")
library("gstat")

# load map
map2 <- readFITS(file = "data/map_rs_3.fits")

# transform map into CMBDataFrame
map2cols <-map2$col
temps <- map2cols[]
I <- unlist(temps)

df_temp <- as.data.frame(I)

sky_rs2 <- as.CMBDataFrame(df_temp, nside = 8, ordering = "nested")

# get map with cartesian coordinates
sky_rs2_cart <- CMBDataFrame(sky_rs2, coords  = "cartesian")

# transform CMBDataFrame into DataFrame
df_sky_rs2 <- as.data.frame(sky_rs2_cart)

# transform DataFrame into Matrix
mat_sky_r2 <-  as.matrix(df_sky_rs2)

# get maximal difference
max(diff(mat_sky_r2[,1]))

# generate spatial object
coordinates(df_sky_rs2) <- ~x+y+z

# create inla.mesh
mesh_sky_r2 <- inla.mesh.2d(loc = df_sky_rs2, max.edge = 0.8843076)

plot(mesh_sky_r2)
