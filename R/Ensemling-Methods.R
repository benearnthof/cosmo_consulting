source("global_library.R")

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
window(sky_rs9_sph) <- CMBWindow(theta = c(pi/2,pi/2,pi/2.1, pi/2.1), phi = c(0,pi/100,pi/100,0))
plot(sky_rs9_sph)
window(sky_rs9, in.pixels)