# install.packages(c("sf", "rnaturalearth", "rnaturalearthdata"))

library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

esp_prov <- ne_states(country = "Spain", returnclass = "sf")

provincias_aragon <- esp_prov[esp_prov$name %in% c("Navarra", "Huesca", "Zaragoza", "Teruel"), ]

aragon <- st_union(provincias_aragon)

aragon <- st_sf(geometry = aragon)

#aragon_utm <- st_transform(aragon, 25830)

cellsize <- .1  # ~ 10 km

grid_aragon <- st_make_grid(
  aragon, #_utm
  cellsize = cellsize,
  square = TRUE
)

grid_aragon <- st_sf(geometry = grid_aragon)

grid_aragon_clip <- st_intersection(grid_aragon, aragon) # _utm

plot(st_geometry(grid_aragon_clip), border = "grey")

centroides <- st_centroid(grid_aragon_clip)
centroides$X <- st_coordinates(centroides)[,1]
centroides$Y <- st_coordinates(centroides)[,2]

plot_grid_with_significance <- function(
    grid_sf, values, lower, upper,
    title = "Map with significance",
    value_range = NULL,
    nonsig_method = c("hatch", "shade")
) {
  
  library(sf)
  library(ggplot2)
  library(dplyr)
  
  nonsig_method <- match.arg(nonsig_method)
  
  if (!inherits(grid_sf, "sf"))
    stop("grid_sf must be an sf object")
  if (any(lengths(list(values, lower, upper)) != nrow(grid_sf)))
    stop("values, lower, and upper must match number of grid cells")
  
  nonsig <- (lower <= 0 & upper >= 0)
  
  grid_sf$estimate <- values
  grid_sf$nonsig  <- nonsig
  
  p <- ggplot(grid_sf) +
    geom_sf(aes(fill = estimate), color = NA)
  
  if (is.null(value_range)) {
    p <- p + scale_fill_viridis_c(option = "plasma")
  } else {
    p <- p + scale_fill_viridis_c(option = "plasma",
                                  limits = value_range,
                                  oob = scales::squish)
  }
  
  if (nonsig_method == "hatch") {
    p <- p +
      geom_sf(data = subset(grid_sf, nonsig),
              fill = NA, color = "white", size = 0.2)
  }
  
  if (nonsig_method == "shade") {
    p <- p +
      geom_sf(data = subset(grid_sf, nonsig),
              fill = "grey70", alpha = 0.5, color = NA)
  }
  
  p +
    labs(fill = "Estimate", title = title) +
    theme_minimal()
}

### fit models for quantiles 05, 50, 95


### obtain kriged values
krigep05.1
krigep50.1
krigep95.1

### plots
plot_grid_with_significance(
  grid_aragon_clip, colMeans(krigep05.1) / sd(df0$year),
  lower = apply(krigep05.1, 2, quantile, prob = c(0.05)) / sd(df0$year),
  upper = apply(krigep05.1, 2, quantile, prob = c(0.95)) / sd(df0$year),
  value_range = c(-0.2,1.2)) + geom_point(aes(x=lon, y=lat), data = coords[,c("lon", "lat")])
plot_grid_with_significance(
  grid_aragon_clip, colMeans(krigep50.1) / sd(df0$year),
  lower = apply(krigep50.1, 2, quantile, prob = c(0.05)) / sd(df0$year),
  upper = apply(krigep50.1, 2, quantile, prob = c(0.95)) / sd(df0$year),
  value_range = c(-0.2,1.2)) + geom_point(aes(x=lon, y=lat), data = coords[,c("lon", "lat")])
plot_grid_with_significance(
  grid_aragon_clip, colMeans(krigep95.1) / sd(df0$year),
  lower = apply(krigep95.1, 2, quantile, prob = c(0.05)) / sd(df0$year),
  upper = apply(krigep95.1, 2, quantile, prob = c(0.95)) / sd(df0$year),
  value_range = c(-0.2,1.2)) + geom_point(aes(x=lon, y=lat), data = coords[,c("lon", "lat")])
