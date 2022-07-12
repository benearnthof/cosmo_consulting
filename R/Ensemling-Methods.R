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
# Define window as one pixel at a certain resolution of the sky
# for example: 1 window is one pixel at resolution 3
# that means the window contains 4086 pixels if we slice it out from
# the full sky of resolution 9
w1 <- window(sky_rs9, in.pixels = c(1), in.pixels.res = 3)
plot(w1)

pixel_list <- seq(1:100)
lapply(pixel_list, function(x) neighbours(x,9))
