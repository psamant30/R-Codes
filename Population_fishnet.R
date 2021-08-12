create_population_fishnet <- function(pop_data, cover_area, fishnet_size, fishnet_units='mi'){
  # check that units are correct
  stopifnot(fishnet_units=='mi')
  
  # proceed with grid
  print(paste('Creating fishnet grid of dimensions', fishnet_size, fishnet_units, 'by', fishnet_size, fishnet_units, sep=' '))
  grid_width_m = mi_to_m(fishnet_size)
  print(paste('Fishnet width is', grid_width_m, 'meters', sep=' '))
  net <- st_make_grid(x=cover_area, cellsize = grid_width_m, what='polygons', square=TRUE) 
  
  # convert the fishnet raster into a sf polygon dataset and add a "net_id" column
  print('Rasterizing fishnet.')
  net_agg <- st_as_sf(net) %>% tibble::rowid_to_column(., "net_id")
  print(paste('Fishnet contains', as.character(length(unique(net_agg$net_id))), 'cells', sep=' '))
  
  # list of net cells IDs that intersect with the cover area (if cover area is irregular)
  print('Intersecting fishnet with population data.')
  net_intersect <- st_intersects(cover_area, net_agg) 
  net_cover_area <- net_agg[unique(unlist(net_intersect)),]
  net_cover_area$net_area <- st_area(net_cover_area$x)
  net_blocks_intersect <- st_intersection(pop_data, net_cover_area)
  
  # group by fishnet cell and calc block stats.
  print('Calculating block and fishnet intersection areas.')
  net_blocks_intersect <- net_blocks_intersect %>%
    # calculate the area of each tract/fishet intersection
    mutate(intersect_area = st_area(net_blocks_intersect) ) %>%
    # group by the fishnet ID
    group_by(net_id) %>%
    # calculate the fraction of each block group covered by a fishnet cell
    # then assume the population is proportional to the overlap area
    mutate(cnt = n(), pcnt_of_block = intersect_area/block_area, intersect_pop = estimate * pcnt_of_block) %>%
    arrange(net_id)
  
  print('Summarizing fishnet population sizes.')
  fishnet_pop <- net_blocks_intersect %>%
    group_by(net_id) %>%
    summarise(net_pop = sum(intersect_pop)) 
  
  # drop units (required for plotting)
  fishnet_pop$unitless_net_pop <- units::drop_units(fishnet_pop$net_pop)
  
  return(fishnet_pop)
}
