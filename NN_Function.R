cumberland_1m <- st_read(dsn = 'Data/cumberland_one_mile_fishnet.shp', crs)

crs(cumberland_1m)
st_cumberland <- sf::st_transform(
  cumberland_1m 
  , crs = nbr_crs
)








NN_point_features <- function(var_list, fishnet, k){
  NN_results <- foreach(i = length(var_list),
                        print(i),
                        .export=c('nn_function'),
                        .packages=c('raster', 'sf', 'dplyr', "FNN", "tibble", "tidyr")) %dopar% { 
                          feature <- var_list
                          print(feature)
                          fishnet_centroid_XY <- st_coordinates(st_centroid(fishnet))
                          dat <- var_list[[i]] 
                          print(dat)
                          if(nrow(dat) >= k){
                            net_NN <- nn_function(fishnet_centroid_XY,
                                                  st_coordinates(dat)[,1:2], k) %>%
                              mutate(feature_name = paste0("NN_",dat$PrimarySIC6_Description),
                                     net_id = fishnet$net_id) %>%
                              left_join(., fishnet, by = "net_id") %>%
                              #rename("value" = value.x) %>%
                              #dplyr::select(-value.y) %>%
                              st_as_sf()
                          #} else {
                            #net_NN <- data.frame(value = rep(NA, nrow(fishnet))) %>%
                              #mutate(feature_name =  paste0("NN_",feature),
                                     #net_id = fishnet$net_id) %>%
                              #left_join(., fishnet, by = "net_id") 
                              
                            #st_as_sf()                       
                          #}
                        }
  #names(NN_results) <- paste0("NN_",var_list)
  return(NN_results)
                        }
}

NN_results <- NN_point_features(var_list[[1]][["PrimarySIC6_Description"]], net_Cumberland, k_nearest_neighbors)
#print(var_list[1]["PrimarySIC6_Description"])
