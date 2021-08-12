## ----setup, include=FALSE, warning=FALSE, messages=FALSE, echo=FALSE, cache=FALSE----------------
knitr::opts_chunk$set(include=FALSE, warning=FALSE, 
                      message=FALSE, echo=FALSE, cache=TRUE, fig.align="center")


## ----packages, message=FALSE, warning=FALSE, cache=FALSE, echo=FALSE-----------------------------
library("sf")            # Spatial data objects and methods
library("mapview")       # Interactive Map Viewing
library("ggmap")         # ggplot2 addon for base maps
library("cowplot")
library("spatstat")      # KDE and other spatial functions
library("raster")        # cell-based spatial operations
library("tidyverse")     # data manipulation framework
library("Hmisc")         # using cut2() functions for ggplot legends
library("fitdistrplus")  # Distribution fitting functions
library("lubridate")     # Power tools for handling dates
library("tidycensus")
library("lwgeom")
library("Hmisc")
library("hrbrthemes")
library("gridExtra")
library("patchwork")
library("spdep")         # KNN functions
library("foreach")
library("doParallel")
library("corrplot")
library("ranger")        # randomforest implimentation      
library("glmnet")        # for Ridge and Lasso Regression
library("knitr")         # for kable table
library("kableExtra")
library("FNN")           # KNN for CPS vs. NN plots
library("groupdata2")
library("htmltools")
library("viridis")
library("viridisLite")
library("maptools")


## ----logos, echo = FALSE, include = FALSE, cache = FALSE-----------------------------------------


## ----themes--------------------------------------------------------------------------------------
mapTheme <- function() {
  theme(
    plot.title = element_text(size = 14, family = "sans", face = "plain", hjust = 0),
    plot.subtitle=element_text(size = 11, family = "sans", hjust = 0),
    plot.caption=element_text(size = 10, family = "sans", face = "italic", hjust = 0),
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.title = element_text(size = 10, family = "sans"),
    legend.text = element_text(size = 9, family = "sans"),
    panel.border = element_blank()
  )
}

plotTheme <- function() {
  theme(
    plot.title = element_text(size = 14, family = "sans", face = "plain", hjust = 0),
    plot.subtitle=element_text(size = 11, family = "sans", hjust = 0),
    plot.caption=element_text(size = 10, family = "sans", face = "italic", hjust = 0), 
    axis.title.x = element_text(size = 10, family = "sans", face = "plain", hjust = 1, vjust = -0.5),
    axis.title.y = element_text(size = 10, family = "sans", face = "plain", hjust = 1, vjust = 1),
    axis.text = element_text(size = 9, family = "sans", face = "plain"),
    panel.background = element_blank(),
    panel.grid.minor = element_line(colour = "gray"),
    panel.grid.major = element_line(colour = "gray"),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 10, family = "sans"),
    legend.text = element_text(size = 9, family = "sans"),
    axis.line = element_blank()
  )
}


## ----options-------------------------------------------------------------------------------------
mapviewOptions(basemaps = c("Stamen.TonerLite", "OpenStreetMap.DE"))
base_dir = "Data"
fishnet_grid_dim = 1000
k_direction = 8 # 4 = rook, 8 = queen
k_nearest_neighbors = 5
# Either k (e.g. 5 or 10) or "LOOCV"
n_folds = "LOOCV"
# threshold quntile for statArea grouping
stat_area_quantile = 0.60
# Number of simulations for CPS vs. NN
simulations = 1000
# Number of neighbors for CPS vs. NN
k = 5
# random seed
set.seed(717)


## ----SOURCE--------------------------------------------------------------------------------------
source('R/FUNCTIONS_VAPAP.R', echo = TRUE, keep.source = TRUE)
source('R/FEA_CREATE_VARIABLES.R', echo = TRUE, keep.source = TRUE)

setwd('~/predict-align-prevent') # set working dir as repo directory
source('R/FUNCTIONS_VAPAP.R', echo = TRUE, keep.source = TRUE)

library("sf")         # version 0.9-8
library("mapview")    # version 2.9.0
library("ggmap")      # version 3.3.0 
library("raster")     # version 3.4-5
library("tidyverse")  # version 1.3.0

#######################################################
## CREATE AN EMPTY RASTER                            ##
#######################################################

# all New Jersey counties in a single geojson file
nbr <- read_sf("https://opendata.arcgis.com/datasets/5f45e1ece6e14ef5866974a7b57d3b95_1.geojson", quiet=FALSE) %>%
  st_transform('ESRI:102311')  %>% # NAD_1983_HARN_StatePlane_New_Jersey_FIPS_2900
  filter(COUNTY=='CUMBERLAND') # subset to just cumberland county

# turn into empty raster which will later contain heatmap data
nbr_diss <- nbr %>%
  mutate(dissolve = 1) %>%
  # get rid of slivers
  st_buffer(., dist = 0.1) %>%
  group_by(dissolve) %>%
  summarise()

nbr_rast_SP <- raster(as(nbr_diss, "Spatial"), nrows = 2000, ncol = 2000)

#######################################################
## DOWNLOAD THE BASEMAP                              ##
#######################################################

# re-download the geojson file, filter, and use to get basemap
nbr <- read_sf("Data/County_Boundaries_of_NJ.geojson", quiet=FALSE) %>%
  st_transform('ESRI:102311')  %>% # NAD_1983_HARN_StatePlane_New_Jersey_FIPS_2900
  filter(COUNTY=='CUMBERLAND')

# here I define the bounding box for collecting the basemap by Cumberland county; while the original example
# restricted the area to only the bounding box (plus a small buffer) that contained the maltreatment data
nbr_bbox <- unname(st_bbox(ll(nbr))) # note the ll() function is devined in R/FUNCTIONS_VAPAP.R
print(nbr_bbox)
# the original version of the code used get_map(), a wrapper to get_googlemap(), get_openstreetmap(), and get_stamenmap()
# I had more success using the get_stamenmap() function directly
cps_base_map   <- get_stamenmap(bbox = nbr_bbox, maptype = "toner", force=TRUE)
plot(cps_base_map)


## ----neighborhoods-------------------------------------------------------------------------------
# Richmond Neighborhoods

nbr <- read_sf("https://opendata.arcgis.com/datasets/5f45e1ece6e14ef5866974a7b57d3b95_1.geojson")  %>%
  st_sf() %>%
  st_transform('ESRI:102311')%>%
  

nbr_diss <- nbr1 %>%
  mutate(dissolve = 1) %>%
  # get rid of slivers
  st_buffer(., dist = 0.1) %>%
  group_by(dissolve) %>%
  summarise()

nbr_rast_SP <- raster(as(nbr_diss, "Spatial"), nrows = 2000, ncol = 2000)


## ----basemap-------------------------------------------------------------------------------------
nbr <- read_sf("https://opendata.arcgis.com/datasets/5f45e1ece6e14ef5866974a7b57d3b95_1.geojson") %>%
  st_transform(crs = 'ESRI:102311')

nbr_bbox <- unname(st_bbox(ll(st_buffer(var_list[[1]][["geometry"]], dist = 0.1))))
print(nbr_bbox)
cps_base_map   <- get_map(location = unname(st_bbox(nbr_bbox)),
                          source = "stamen",
                          maptype = "toner")
plot(cps_base_map)
### get CPS_Accepted values (add 1 column for dissolving)
cps_dissolve <- var_list[[1]] %>%
  mutate(value = 1) %>%
  #print(cps_dissolve)
  dplyr::select(value)
  print(cps_dissolve)
st_geometry(cps_dissolve)
st_geometry(net)
cps_dissolve= st_transform(cps_dissolve, crs='NAD83') 
## ----fishnet-------------------------------------------------------------------------------------
net <- st_make_grid(nbr1, cellsize = fishnet_grid_dim)
print(net)

nbr2 <- nbr1 %>% 
  st_make_grid(n = fishnet_grid_dim)


library(tidycensus)
library(tidyverse)
options(tigris_use_cache = TRUE)
dev.new(width=100, height=50)
plot(nbr2)
#Use census block groups - population - 
census_api_key('d68574340a0011206cf42d4331abca4bf602d73a', overwrite = TRUE, install = TRUE)
nbr1 <- get_acs(state = "NJ", county = "Cumberland", geography = "tract", 
                  variables = "S0601_C01_001", geometry = TRUE)%>%
     st_transform(crs = '+proj=gnom +lat_0=39.55 +lon_0=-75.35')
nbr2 <- st_as_sf(nbr1)
# count CPS incidents per net cell - really just to get net raster into sf polygon format

net_agg <- aggregate(cps_dissolve, nbr2, sum) %>%
  tibble::rowid_to_column(.,"net_id")

# list of net cells IDs that intersect with Richmond
net_intersect <- st_intersects(nbr1, net_agg) 
# extract Richmonds net cells based on intersect ID
net_NJ <- net_agg[unique(unlist(net_intersect)),]
net_hood <- st_join(net_NJ, nbr1, largest = TRUE)
listw <- nb2listw(poly2nb(as(net_NJ, "Spatial"), queen = TRUE))


## ----population_data-----------------------------------------------------------------------------
census_api_key("d68574340a0011206cf42d4331abca4bf602d73a")
vars10 <- c("P001001") # total population (correct, I checked the web)
## get total 2010 census pop for blocks & calculate area
NJ_block <- get_decennial(geography = "block", variables = vars10, year = 2010,
summary_var = "P001001", state = 34, county = 011, geometry = TRUE) %>%
st_transform(crs ='ESRI:102311')

NJ_block <- get_decennial(
  geography = "block",
  state = "NJ",
  county = "Cumberland",
  year = 2010,
  variables = c("P001001"),
  geometry = TRUE,
  output = "wide"
)
# calc area
NJ_block <- st_read("Data/County_Boundaries_of_NJ.shp")    #was having issues with above code when I knit so pulled in from the PAP project file
NJ_block <- st_transform(nbr1)
NJ_block <- NJ_block %>%
  mutate(acre = as.numeric(st_area(NJ_block)*2.29568e-5),
         acre = units::set_units(acre, acre), 
         pop_acre_rate = estimate / acre) 
plot(NJ_block)

## ----t_intersection------------------------------------------------------------------------------
st_geometry(incidents_shp)
fishnet_CPS_var= st_transform(fishnet_CPS_var, crs='WGS 84 ') 
net_blocks_intersect <- st_intersection(NJ_block, net_NJ)

# group by cell and calc block stats.
# group by cell and calc block stats.
net_blocks_intersect <- net_blocks_intersect %>%
  mutate(intersect_area_acres = as.numeric(st_area(net_blocks_intersect)*2.29568e-5)) %>%
  group_by(net_id) %>%
  dplyr::mutate(cnt = n(),
         pcnt_of_block = intersect_area_acres/acre,
         intersect_pop = estimate * pcnt_of_block) %>%
  arrange(net_id)

plot(net_blocks_intersect)

library(moments)
print(skewness(net_blocks_intersect$pop_acre_rate))
plot(hist(net_blocks_intersect$pop_acre_rate))
r <- raster(net_blocks_intersect)
plot(r)

r1 <- as.data.frame(r, xy = TRUE) 
str(r)
q = ggplot() +
  geom_raster(data = r1,
              aes(x = x, y = y)) + 
  coord_quickmap()
plot(r1)
## ---- summarise_pop------------------------------------------------------------------------------
#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ   # <-  zeros or no zeros!!!!

fishnet_pop <- net_blocks_intersect %>% # xcc
  group_by(net_id) %>%
  dplyr::summarise(net_pop = sum(intersect_pop)) %>%
  filter(net_pop > 0)   # <-  zeros or no zeros!!!!


plot(fishnet_pop)
print(skewness(fishnet_pop$net_pop))
plot(hist(fishnet_pop$net_pop))

x <- fishnet_pop$net_pop
h<-hist(x, breaks=10, col="red", xlab="Miles Per Gallon",
        main="Histogram with Normal Curve")
xfit<-seq(min(x),max(x),length=40)
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="blue", lwd=2)

hist(starbucks)


dev.off()
options(device = "RStudioGD")
######### MAKE NET AND RATE FOR ALL CPS VARS
print(names(var_list[[1]][["Intake.Service"]]))
CPS_vars <- grep("CPS-",var_list[[1]][["Intake.Service"]], value = TRUE)
CPS_agg <- NULL
for(i in seq_along(CPS_vars)){
  var_name <- paste0("net_",CPS_vars[i])
  cat(var_name,"\n")
  
  CPS_dat <- var_list[[1]] %>%
    mutate(value = 1) %>%
    #dplyr::select(value)
  CPS_dat= st_transform(CPS_dat, crs='NAD83')
  #CPS_dat= st_transform(CPS_dat, crs='ESRI:102311')
  fishnet_CPS_var <- aggregate(CPS_dat,fishnet_pop, sum) %>%
    st_drop_geometry() %>%
    mutate(Feature = var_name) %>%
    dplyr::select(Feature,value)
  
  CPS_agg <- rbind(CPS_, fishnet_CPS_var)
}

plot(fp)
fp =st_make_grid(fishnet_pop)
fp <- lapply(fishnet_pop$geometry, Polygon)
plot(fp)
print(NROW(CPS_dat))
is.list(fishnet_pop$geometry)

CPS_vars <- unique(var_list[[1]][["Intake.Service"]])
CPS_agg <- NULL
for(i in seq_along(CPS_vars)){
  var_name <- paste0("net_",CPS_vars[i])
  cat(var_name,"\n")
  #print(var_list[[1]][["Intake.Service"]][["CPS-Family"]])
  CPS_dat <- filter(var_list[[1]], !grepl(CPS_vars[i],"Intake.Service"))%>%
    mutate(value = 1) %>%
    dplyr::select(value)
  CPS_dat= st_transform(CPS_dat, crs='')
  print(CPS_dat)
  fishnet_CPS_var <- st_join(incidents_shp, fishnet_pop, join = st_within) %>%
    #st_drop_geometry() %>%
    mutate(Feature = var_name) %>%
    dplyr::select(Feature,estimate)
  fishnet_CPS_var=st_as_sf(fishnet_CPS_var)
  fishnet_CPS_var= st_transform(fishnet_CPS_var, crs='WGS 84')
  CPS_agg <- rbind(incidents_shp, fishnet_CPS_var)
}
plot(CPS_agg)
CPS_dat <- CPS_dat %>%
  mutate(id = rep(seq(1:nrow(fishnet_pop)),length(CPS_vars))) %>%
  spread(Feature, value) %>%
  dplyr::select(-id) %>%
  mutate(geometry = fishnet_pop$geometry) %>%
  st_as_sf()
#### Spatial join of fishnet_pop and fishnet_cps to then calculate rate for all CPS features
fishnet_pop_cps <- st_join(fishnet_pop, CPS_dat, join = st_equals) %>%
  mutate_at(vars(paste0("net_",CPS_vars)), funs(rate = ./(net_pop/100)))  %>% # cps per 100 person
  rename_at(vars( contains( "_rate")), funs(paste("rate", gsub("net_|_rate", "", .), sep = "_"))) %>% 
  replace(is.na(.), 0) # replace NA with zero

fishnet_coords <- fishnet_pop_cps %>%
  st_centroid() %>%
  st_coordinates() %>%
  as.matrix()
library(rgdal)
library(rgeos)


hist(net_blocks_intersect$pop_acre_rate)
plot(cps$geometry)
plot(net_blocks_intersect$geometry)
cpsA= st_transform(cps, crs='ESRI:102311')
naa= st_transform(net_blocks_intersect, crs='ESRI:102311')
int <- as_tibble(st_intersection(cpsA$geometry, naa$geometry))


starbucks  <- as.ppp(cpsA)
plot(starbucks)


dev.new(width=1000, height=100, unit="in")
Q <- quadratcount(starbucks, nx= 200, ny=100)
plot(starbucks, pch=20, cols="grey70", main=NULL)  # Plot points
plot(Q)

print(Q)
Q.d <- intensity(Q)

# Plot the density
plot(intensity(Q, image=TRUE), main=NULL, las=1)  # Plot density raster
plot(starbucks, pch=20, cex=0.6, col=rgb(0,0,0,.5), add=TRUE)  # Add points


s  <-read_sf("Data/County_Boundaries_of_NJ.shp", quiet=FALSE)
w  <- as.owin(s)
w.km <- rescale(w, 1000)

marks(starbucks)  <- NULL
ma= as.owin(cps_base_map)
Window(starbucks) <- cps_base_map

starbucks.km <- rescale(starbucks, 1000, "km")

Q   <- quadratcount(starbucks.km, nx= 6, ny=3)
Q.d <- intensity(Q)
plot(intensity(Q, image=TRUE), main=NULL, las=1)  # Plot density raster
plot(starbucks.km, pch=20, cex=0.6, col=rgb(0,0,0,.5), add=TRUE)  # Add point
K3 <- density(starbucks.km, kernel = "disc", sigma=50) # Using a 50km bandwidth
plot(K3, main=NULL, las=1)
contour(K3, add=TRUE)
## ----CPS_COUNT_BY_FISHNET_plot-------------------------------------------------------------------
fishnet_pop_cps_cut <- fishnet_pop_cps %>%
  mutate(net_CPS_Accepted = ifelse(is.na(net_CPS_Accepted), 0, net_CPS_Accepted)) %>% 
  make_cuts(., "net_CPS_Accepted", cuts = "breaks", n_breaks = 10)

CPS_COUNT_BY_FISHNET_PLOT <- ggmap(cps_base_map) +
  geom_sf(data = ll(fishnet_pop_cps_cut), aes(fill = cut_val), inherit.aes = FALSE, color = NA, alpha = 0.8) +
  labs(title = "CPS count per\nfishnet cell") +
  scale_fill_viridis_d(na.value = NA, option = "D", direction = 1, name = "CPS Count") +
  mapTheme() +
  theme(plot.title = element_text(size = 14, family = "sans", face = "plain", hjust = 0),
        plot.subtitle=element_text(size = 11, family = "sans", hjust = 0),
        plot.caption=element_text(size = 10, family = "sans", face = "italic", hjust = 0),
        axis.line = element_blank(),
        legend.title = element_text(size = 10, family = "sans"),
        legend.text = element_text(size = 9, family = "sans"))


## ----CPS_RATE_BY_FISHNET_plot--------------------------------------------------------------------
fishnet_pop_cps_rate_cut <- fishnet_pop_cps %>%
  mutate(rate_CPS_Accepted = ifelse(is.na(rate_CPS_Accepted), 0, rate_CPS_Accepted)) %>% 
  make_cuts(., "rate_CPS_Accepted", cuts = "breaks", n_breaks = 10)

CPS_RATE_BY_FISHNET_PLOT <- ggmap(cps_base_map) +
  geom_sf(data = ll(fishnet_pop_cps_rate_cut), aes(fill = cut_val), inherit.aes = FALSE, color = NA, alpha = 0.8) +
  labs(title = "Child Protective Service rate\nper 100 people") +
  scale_fill_viridis_d(na.value = NA, option = "D", direction = 1, name = "CPS Rate\nper 100") +
  mapTheme() +
  theme(plot.title = element_text(size = 14, family = "sans", face = "plain", hjust = 0),
        plot.subtitle=element_text(size = 11, family = "sans", hjust = 0),
        plot.caption=element_text(size = 10, family = "sans", face = "italic", hjust = 0),
        axis.line = element_blank(),
        legend.title = element_text(size = 10, family = "sans"),
        legend.text = element_text(size = 9, family = "sans"))


## ----CPS_Count_Month_Year_table------------------------------------------------------------------
CPS_Counts_Year_table  <- table(lubridate::year(var_list[[3]]$Intake.RcvdDate))
CPS_Counts_Month_table <- table(lubridate::month(var_list[[3]]$Intake.RcvdDate))


## ----CPS_HIST_BY_DATE_plot-----------------------------------------------------------------------
CPS_by_year <- lubridate::year(var_list[[3]]$Intake.RcvdDate) %>%
  data.frame(year = .)
CPS_HIST_BY_DATE <- ggplot(CPS_by_year, aes(x = year)) +
  geom_histogram() +
  plotTheme()
plot(CPS_HIST_BY_DATE)

## ----CPS_POINT_BY_MONTH_plot---------------------------------------------------------------------

months <- c(1,2,3,4,5,6,7,8,9,10,11,12)
cps <- var_list[[2]] %>%
  mutate(year  = lubridate::year(Intake.RcvdDate),
         month = lubridate::month(Intake.RcvdDate),
         month = months[month],
         month = fct_relevel(month, months))
print(cps)


library(magrittr) #for the pipe
cps<- cps%>%
  dplyr::mutate(X = sf::st_coordinates(.)[,1],
                Y = sf::st_coordinates(.)[,2])

data= data.frame((ll(cps)), year = cps$year)
print(data)

chi_bb <- c(
  left = -75.5,
  bottom = 39.15,
  right = -74.7427,
  top = 39.6
)

cps_base_map <- get_stamenmap(
  bbox = chi_bb,
  zoom = 12
)
map1 <- ggmap(cps_base_map)

projection(cps_base_map) <- CRS("+WSG84")
plot(map1)

Statmap <- ggmap(cps_base_map) + geom_point(size=0.7)+
  geom_polygon(data=cps,aes(x = X, y = Y), 
                 color = "white")
plot(Statmap)

map = map1+
  ggplot(data=cps,
         aes(x=X,y=Y)) +
  geom_point(size=2)

plot(map)

m = map1+
  geom_point(data = cps, aes(x = X, y = Y), size=1, alpha=0.5) + scale_y_continuous(limits=c(39, 39.70))+scale_x_continuous(limits = c(-75.5, -76))+
  theme(legend.position = "bottom")

m = map1 + geom_point(data=cps, aes(x=X, y=Y),size=1)
plot(m)

my_sf1 <- st_as_sf(aa, coords = c('X', 'Y'))
my_sf <- st_set_crs(my_sf1, crs = 'ESRI:102311')
ggplot(my_sf1) + 
  geom_sf(aes())

aa <-cps[ , c("X", "Y")]    

mapview(my_sf1)
CPS_POINT_BY_MONTH_plot <- map1+ 
  ggplot(cps , aes(Y, X, color = as.factor(year)), size=1.5, alpha = 0.8)+
  geom_polygon(data = data , aes(x=Y, y=X))
  
  
CPS_POINT_BY_MONTH_plot <- ggmap(cps_base_map) + 
  ggplot(data=cps , aes(X, Y, color = as.factor(year)), size=1.5, alpha = 0.8)+
  geom_polygon() +
  scale_color_viridis_d(name = "Year") +
  labs(title = "CPS Accepted in NJ by Year",
       caption = "source: **************") +
  facet_wrap(~year) +
  mapTheme() +
  theme(
    legend.key = element_rect(fill = "white"),
    strip.text = element_text(face = "plain", size = 11),
    legend.position = c(0.85, 0.25) # or "none
  )
plot(CPS_POINT_BY_MONTH_plot)
## ----CPS_KDE_BY_YEAR_plot------------------------------------------------------------------------
variable = "year"
values <- unique(cps[[variable]])
year_dat <- list()
brks <- 9
window_cps <- get_window(cps, buff_dist = 10000)
for(i in seq_along(values)){
  dat <- filter(cps, !!as.name(variable) == values[i])
  points.ppp <- as.ppp(st_coordinates(ll(dat)),window_cps)
  densityRaster <- raster(density(points.ppp, scalekernel=TRUE, sigma = 0.005))
  dens_data <- gplot_data(densityRaster, maxpixels = 2500) %>%
    mutate(!!as.name(variable) := values[i])
  year_dat[[i]] <- dens_data
}
year_dat <- do.call(rbind, year_dat)

CPS_KDE_BY_YEAR_plot <- ggmap(cps_base_map) +
  geom_tile(data = year_dat, 
            aes(x,y,fill = as.factor(ntile(value,brks)), 
                group = !!as.name(variable)), alpha=0.8) +
  scale_fill_viridis_d(name = variable) +
  labs(title = "CPS accepted in NJ by year",
       caption = "Figure 5.2") +
  facet_wrap(vars(!!as.name(variable))) +
  mapTheme() +
  theme(
    legend.key = element_rect(fill = "white"),
    strip.text = element_text(face = "plain", size = 11, hjust = 0),
    strip.background = element_rect(fill = "white"),
    legend.position = "none"
  )

plot(CPS_KDE_BY_YEAR_plot)
## ----CPS_TREND_BY_MONTH_YEAR_plot----------------------------------------------------------------
CPS_by_year_month <- st_drop_geometry(var_list[[1]]) %>%
  mutate(month = lubridate::month(Intake.RcvdDate),
         year  = lubridate::year(Intake.RcvdDate))%>%
  dplyr::select(month, year) %>%
  group_by(month, year) %>%
  dplyr::mutate(m_count = n()) %>%
  distinct() %>%
  ungroup()

CPS_TREND_BY_MONTH_YEAR_plot <- ggplot(CPS_by_year_month, aes(x = year, y = m_count)) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3)) +
  labs(y="Incidents per month") +
  plotTheme()
plot(CPS_TREND_BY_MONTH_YEAR_plot)

## ----CPS_LINE_AGG_BY_MOTNH_plot------------------------------------------------------------------
CPS_agg_by_month <- st_drop_geometry(var_list[[1]]) %>%
  mutate(month = lubridate::month(Intake.RcvdDate),
         year  = lubridate::year(Intake.RcvdDate))%>%
  group_by(month) %>%
  dplyr::summarise(count = n())

CPS_LINE_AGG_BY_MOTNH_plot <- ggplot(CPS_agg_by_month, aes(x = month, y = count)) +
  scale_x_continuous(breaks = seq(1,12), labels = seq(1,12)) +
  geom_line() +
  plotTheme()

plot(CPS_LINE_AGG_BY_MOTNH_plot)
## ----CPS_LINE_NORMALIZED_plot--------------------------------------------------------------------
CPS_normalized_by_month <- st_drop_geometry(var_list[[1]]) %>%
  mutate(month = lubridate::month(Intake.RcvdDate),
         year  = lubridate::year(Intake.RcvdDate)) %>%
  group_by(year, month) %>%
  dplyr::summarise(m_total = n()) %>%
  arrange(month, year) %>%
  dplyr::select(month, year, m_total) %>%
  ungroup() %>%
  group_by(month) %>%
  mutate(m_mean = mean(m_total),
         m_sd   = sd(m_total),
         m_z    = (m_total - m_mean) / m_sd)

CPS_LINE_NORMALIZED_plot <- ggplot(CPS_normalized_by_month, aes(x = as.factor(month), 
                                                                y = m_z, group = year, 
                                                                color = as.factor(year))) +
  geom_line() +
  geom_hline(yintercept = 0, color = "gray20", linetype = "dashed") +
  scale_color_viridis_d(name = "year") +
  labs(x = "month") +
  scale_y_continuous(limits = c(-2,2)) +
  plotTheme()

plot(CPS_LINE_NORMALIZED_plot)
## ----CPS_CALENDAR_plot---------------------------------------------------------------------------
CPS_agg_cal <- st_drop_geometry(var_list[[1]]) %>%
  mutate() %>%
  mutate(day = factor(weekdays(Intake.RcvdDate,T),
                      levels = rev(c("Mon", "Tue", "Wed", "Thu","Fri", "Sat", "Sun"))),
         week = week(Intake.RcvdDate),
         month = month(Intake.RcvdDate),
         year  = year(Intake.RcvdDate)) %>%
  dplyr::select(day, week, month, year) %>%
  group_by(day, week, month, year) %>%
  summarise(day_cnt = n()) %>%
  complete(day, week, month, year) 

CPS_CALENDAR_plot <- ggplot(CPS_agg_cal, aes(x = week, y = day, fill = day_cnt)) +
  viridis::scale_fill_viridis(name="Incidents",
                              option = 'C',
                              direction = 1,
                              na.value = "gray90") +
  geom_tile(color = 'white', size = 0.1) +
  facet_wrap('year', ncol = 1) +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = seq(1, 52, length = 12),
    labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  theme_ipsum_rc()


## ----CPS_COMPARE_FISHNET_GRID_SIZE_3x2_plot------------------------------------------------------
grid_seq <- c(500,1000,1500)
p_loc_l  <- vector(mode = "list", length = length(grid_seq))
p_hist_l <- vector(mode = "list", length = length(grid_seq))
for(i in seq_along(grid_seq)){
  cat(grid_seq[i], "\n")
  net_i <- st_make_grid(nbr1, cellsize = grid_seq[i])
  net_agg_i <- aggregate(cps_dissolve, net_i, sum) %>% 
    mutate(value = ifelse(is.na(value),0,value))
  
  net_intersect_i <- st_intersects(nbr1, net_agg_i) 
  # extract Richmonds net cells based on intersect ID
  net_Richmond_i <- net_agg_i[unique(unlist(net_intersect_i)),]
  
  net_Richmond_i$class <- Hmisc::cut2(net_Richmond_i$value, g = 9)
  p_loc <- ggmap(cps_base_map) +
    geom_sf(data = ll(net_Richmond_i), aes(fill = class), 
            color = NA, inherit.aes = FALSE, size = 0.5, alpha = 0.8) +
    scale_fill_viridis_d(na.value=NA,
                         name   = paste0("Values","\n[quantiles]"),
                         breaks = levels(net_agg_i$class),
                         labels = levels(net_agg_i$class)) +
    mapTheme()
  
  p_loc_l[[i]] <- p_loc
  
  p_hist <- ggplot(net_Richmond_i, aes(x=value)) +
    geom_histogram(bins = 30) +
    scale_x_continuous(limits = c(-1,100)) +
    scale_y_continuous(limits = c(0,15)) +
    labs(title = paste0("Cell Dimensions =\n",grid_seq[i]," ft sq")) +
    plotTheme()
  
  p_hist_l[[i]] <- p_hist
}

CPS_COMPARE_FISHNET_GRID_SIZE_3x2_plot <- grid.arrange(p_hist_l[[1]], p_hist_l[[2]], p_hist_l[[3]], p_loc_l[[1]], p_loc_l[[2]], p_loc_l[[3]], ncol = 3)
plot(CPS_COMPARE_FISHNET_GRID_SIZE_3x2_plot)

## ----fitdistr_gof--------------------------------------------------------------------------------
number <- as.numeric(na.omit(fishnet_pop_cps$net_CPS_Accepted))
fitp <- fitdist(number,"pois", discrete = TRUE)
fitnb <- fitdist(number,"nbinom", discrete = TRUE)
cdfcomp(list(fitp,fitnb)) # plot
gof <- gofstat(list(fitp,fitnb))


## ----AIC_LINE_FITDISTR_plot----------------------------------------------------------------------
net_cell_dims <- seq(500,5000,50)
aic_results <- matrix(nrow=length(net_cell_dims), ncol = 3)
colnames(aic_results) <- c("cell_dim","pois","nbinom")
for(i in seq_along(net_cell_dims)){
  net <- st_make_grid(nbr,cellsize=net_cell_dims[i])
  
  cps_cnt <- aggregate(cps_dissolve, net, sum)
  
  number <- as.numeric(na.omit(cps_cnt$value))
  fitp <- fitdist(number,"pois", discrete = TRUE)
  fitnb <- fitdist(number,"nbinom", discrete = TRUE)
  gof <- gofstat(list(fitp,fitnb))
  aic_results[i,1] <- net_cell_dims[i]
  aic_results[i,2] <- as.numeric(gof$bic[1])
  aic_results[i,3] <- as.numeric(gof$bic[2])
}

AIC_LINE_FITDISTR_plot <- data.frame(aic_results) %>%
  gather(dist, aic, -cell_dim) %>%
  rename("Distribution" = dist) %>% 
  mutate(Distribution = case_when(
    Distribution == "nbinom" ~ "Negative Binomial",
    Distribution == "pois"   ~ "Poisson"
  )) %>% 
  ggplot(., aes(x = cell_dim, y = aic, group = Distribution, color = Distribution)) +
  geom_line() +
  labs(y = "AIC - goodness of fit",
       x = "Fishnet Cell Dimension (feet)") +
  plotTheme()

plot(AIC_LINE_FITDISTR_plot)
#m ,./# ----risk_protective_var_list--------------------------------------------------------------------
protective_names <- c("CHILD CARE SERVICE",
                      "GROCERS-RETAIL",
                      "CHURCHES",
                      "SCHOOLS-PRE-SCHOOL/KINDERGARTEN-ACADEMIC",
                      "PHARMACIES",
                      "DENTISTS",
                      "CLINICS",
                      "YOUTH ORGANIZATIONS & CENTERS",
                      "COUNSELORS",
                      "CLERGY",
                      "DRUG ABUSE & ADDICTION INFO & TREATMENT",
                      "CRISIS INTERVENTION SERVICE",
                      "GROUP HOMES",
                      "MINISTRIES-OUT REACH",
                      "COUNSELING SERVICES",
                      "RELIGIOUS ORGANIZATIONS",
                      "ALCOHOLISM INFORMATION & TREATMENT CTRS",
                      "CHURCH ORGANIZATIONS",
                      "HOMELESS SHELTERS",
                      "COUNSELORS-LICENSED PROFESSIONAL")
risk_names <- c("BEAUTY SALONS",
                "BARBERS",
                "HOTELS & MOTELS",
                "MANICURING",
                "CONVENIENCE STORES",
                "BARS",
                "BONDS-BAIL",
                "SERVICE STATIONS-GASOLINE & OIL",
                "LAUNDRIES-SELF SERVICE",
                "TATTOOING",
                "OILS-FUEL (WHLS)",
                "PAWNBROKERS",
                "AUTOMOBILE REPAIRING & SERVICE",
                "BEER & ALE-RETAIL",
                "CAR WASHING & POLISHING",
                "LIQUORS-RETAIL",
                "TRUCK-REPAIRING & SERVICE",
                "TRAILERS-REPAIRING & SERVICE",
                "AUTOMOBILE DETAIL & CLEAN-UP SERVICE",
                "WINES-RETAIL",
                "RECREATIONAL VEHICLES-REPAIRING & SVC")
risk_var_list1 <- var_list[[1]][grep(paste(risk_names,collapse="|"), names(var_list[[1]]), value = TRUE)]
risk_var_list <-subset(var_list[[1]], PrimarySIC6_Description %in% risk_names)
protective_var_list <- var_list[[1]][grep(paste(protective_names,collapse="|"), names(var_list[[1]][["PrimarySIC6_Description"]]), value = TRUE)]
protective_var_list<-subset(var_list[[1]], PrimarySIC6_Description %in% protective_names)
## ----risk_points_KDE_compute---------------------------------------------------------------------
boundary.w <- as.owin(as_Spatial(st_cumberland))

class(st_cumberland)

county_flat <- st_transform(st_cumberland, crs = nbr_crs)
county_owin <- as.owin(as_Spatial(county_flat))

st_as_sf(.,coords=c("long","lat"),crs=nbr_crs)

risk_plot_dat <- list()
brks <- 9
window_cps <- crs(cps, buff_dist = 10000)


e <- extent(39.15, 39.6, -75.5, -74.74270)

## create a window for spatstat

WIN <- as.owin(e)
X.ppp <- as(X, "ppp") 

SP.win <- as(e, "SpatialPolygons")
W <- as(SP.win, "owin")

US <- sf::st_read ("Data/County_Boundaries_of_NJ.shp")
SP.win1 <- as(US, "Spatial")
USw <- as(SP.win1,"owin")


library(maptools)
for(i in seq_along(risk_var_list$PrimarySIC6_Description)){
  var_dat <- risk_var_list$geometry[i]
  print(var_dat)
  test = st_coordinates(ll(var_dat))
  points.ppp <- as.ppp(test,W, fatal=TRUE)
  densityRaster <- raster(density(points.ppp, scalekernel=TRUE, sigma = 0.005))
  dens_data <- gplot_data(densityRaster, maxpixels = 2500) %>%
    mutate(variable = names(risk_var_list)[i])
  risk_plot_dat[[i]] <- dens_data
}

plot(risk_plot_dat)
risk_plot_dat <- do.call(rbind, risk_plot_dat)
plot(risk_plot_dat)
# one-liner to extract all 'geometry' cols from list and rbind
risk_compile <- sf::st_as_sf(data.table::rbindlist(lapply(risk_var_list, '[', "geometry")))
risk.points.ppp <- as.ppp(st_coordinates(ll(risk_var_list['geometry'])),boundary.w)
risk_densityRaster <- raster(density(risk.points.ppp, scalekernel=TRUE, sigma = 0.005))
risk_aggregate_plot_data <- gplot_data(risk_densityRaster, maxpixels = 2500) %>%
  mutate(variable = "Risk")%>%
  st_as_sf(risk_densityRaster, crs=nbr_crs)
spplot(r, sp.layout=bounds)
plot(density(risk.points.ppp, 0.01))
den = density(incidents_grp_net$net_id, weights = incidents_grp_net$incident_count/sum(incidents_grp_net$incident_count))
#incidents_grp_net<- incidents_grp_net>%>
plot(den)
risk_den_raster <- raster(den)
## ----RISK_KDE_FACET_PLOT-------------------------------------------------------------------------
write.csv(incidents_grp_net, "Data/Incidents_grp_net.csv")

risk_plot_dat1 <- risk_plot_dat%>%st_as_sf(.,coords=c("x","y"),crs=nbr_crs) 

library("leaflet")
library("data.table")
library("sp")
library("rgdal")
# library("maptools")
library("KernSmooth")
library("raster")
palRaster <- colorNumeric("Spectral", domain = risk_densityRaster@data@values)
plot(palRaster())
## Leaflet map with raster
leaflet() %>% addTiles() %>% 
  addRasterImage(risk_densityRaster, 
                 colors = palRaster, 
                 opacity = .8) %>%
  addLegend(pal = palRaster, 
            values = risk_densityRaster@data@values, 
            title = "Kernel Density of Points")


risk_plot_dat1  <- risk_plot_dat  %>% 
  st_as_sf(.,coords=c("x","y"),crs=nbr_crs) 
  st_transform(nbr_crs) %>% # <- may be unnecessary
  as('Spatial') %>% 
  maptools::elide(.,shift = c(39.4, -75.0)) %>% 
  as('sf') %>% 
  st_set_crs(nbr_crs) %>% 
  st_transform(nbr_crs) 

  ggplot() +
    geom_raster(data = risk_plot_dat1 , aes(x = x, y = y, fill = values)) + 
    ggtitle("Classified Elevation Map - NEON Harvard Forest Field Site")

RISK_KDE_FACET_PLOT <- ggmap(cps_base_map) +
geom_tile(data = risk_plot_dat, 
          aes(y,x,fill = as.factor(ntile(value,brks)), 
          ), alpha=0.4) +
  #facet_wrap(~variable)
  labs(title = "Spatial density of risk factors",
       caption = "Figure 5.4") +
  mapTheme() +
  theme(
    legend.key = element_rect(fill = "white"),
    strip.text = element_text(face = "plain", size = 11, hjust = 0),
    legend.position = "none",
    strip.background = element_rect(fill = "white")
  )
plot(RISK_KDE_FACET_PLOT)

risk_aggregate_plot_data1 <- st_transform(risk_aggregate_plot_data1, crs = st_crs(nbr_crs))
crs(risk_aggregate_plot_data1)
RISK_KDE_PLOT <- ggmap(cps_base_map) +
  #coord_sf(crs = st_crs("WGS84")) +
  ggplot(data = risk_aggregate_plot_data, 
            aes(y,x, fill = as.factor(ntile(value,brks)), 
                ), alpha=0.4) +
  scale_fill_viridis_d() +
  #+facet_wrap(~,scales="free_y") 
  #facet_wrap(~variable) +
  mapTheme() +
  theme(
    legend.key = element_rect(fill = "white"),
    strip.text = element_text(face = "plain", size = 11),
    legend.position = "none"
  )


plot(RISK_KDE_PLOT)
crs(risk_aggregate_plot_data)

p <- ggmap(cps_base_map) +ggplot(risk_aggregate_plot_data, aes(y, x, fill= value, text=variable)) + 
  geom_tile() +
  theme_ipsum()
plot(p)

p <- ggmap(cps_base_map) +ggplot(risk_aggregate_plot_data, aes(x=y, y=x, fill =value)) +
    geom_raster() +
    scale_fill_viridis()

plot(p)

kdeout <-cps %>% 
  with( 
    MASS::kde2d(Y, X, n = 91859,
                lims = c(
                  scales::expand_range(range(Y), .20),
                  scales::expand_range(range(X), .20)
                )
    )
  )

kde_df <- kdeout %>% 
  .[c("x", "y")] %>% 
  cross_df() %>% 
  rename("lon" = "x", "lat" = "y") %>% 
  mutate(density = as.vector(kdeout$z)) 

ggplot() +
  geom_contour_filled(aes(lon, lat, z = density), kde_df)

ggmap(, extent = "device") +
  geom_contour_filled(aes(lon, lat, z = density), kde_df, alpha = .65)

risk_aggregate_plot_data1 <- risk_aggregate_plot_data%>%st_as_sf(.,coords=c("x","y"),crs=nbr_crs)
## ----protective_points_KDE_compute---------------------------------------------------------------
protective_plot_dat <- list()
window_cps <- get_window(cps, buff_dist = 10000)
for(i in seq_along(protective_var_list)){
  var_dat <- protective_var_list[[i]]
  #print(var_dat)
  points.ppp <- as.ppp(st_coordinates(ll(var_dat)),window_cps)
  densityRaster <- raster(density(points.ppp, scalekernel=TRUE, sigma = 0.005))
  dens_data <- gplot_data(densityRaster, maxpixels = 2500) %>%
    mutate(variable = names(protective_var_list)[i])
  print(dens_data)
  protective_plot_dat[[i]] <- dens_data
}
protective_plot_dat <- do.call(rbind, protective_plot_dat)
print(protective_plot_dat)
# one-liner to extract all 'geometry' cols from list and rbind
protective_compile <- sf::st_as_sf(data.table::rbindlist(lapply(protective_var_list, '[', "geometry")))
protective.points.ppp <- as.ppp(st_coordinates(ll(protective_compile)),window_cps)
protective_densityRaster <- raster(density(protective.points.ppp, scalekernel=TRUE, sigma = 0.005))
protective_aggregate_plot_data <- gplot_data(protective_densityRaster, maxpixels = 2500) %>%
  mutate(variable = "Protective")


## ----PROTECTIVE_KDE_FACET_PLOT-------------------------------------------------------------------
PROTECTIVE_KDE_FACET_PLOT <- ggmap(cps_base_map) +
  # geom_point(data = data.frame(st_coordinates(ll(cps)),
  #                              month = cps[[variable]]),
  #            aes(x=X, y=Y), size = 1, color = "gray30", alpha = 0.75) +
  geom_tile(data = protective_plot_dat, 
            aes(x,y,fill = as.factor(ntile(value,brks)), 
                group = variable), alpha=0.8) +
  scale_fill_viridis_d(name = variable) +
  facet_wrap(~variable) +
  labs(title = "Spatial density of protective factors",
       caption = "Figure 5.3") +
  mapTheme() +
  theme(
    legend.key = element_rect(fill = "white"),
    strip.text = element_text(face = "plain", size = 11, hjust = 0),
    legend.position = "none",
    strip.background = element_rect(fill = "white")
  )

PROTECTIVE_KDE_PLOT <- ggmap(cps_base_map) +
  geom_tile(data = protective_aggregate_plot_data, 
            aes(x,y,fill = as.factor(ntile(value,brks)), 
                group = variable), alpha=0.6) +
  scale_fill_viridis_d(name = variable) +
  #facet_wrap(~variable) +
  mapTheme() +
  theme(
    legend.key = element_rect(fill = "white"),
    strip.text = element_text(face = "plain", size = 11),
    legend.position = "none"
  )


## ----spatial_weights_lattice---------------------------------------------------------------------
fishnet_knn <- knn2nb(knearneigh(fishnet_coords, k_direction))
fishnet_Weights <- nb2listw(fishnet_knn, style="W")
localMorans  <- as.data.frame(localmoran(fishnet_pop_cps$net_CPS_Accepted, fishnet_Weights))
globalMorans <- moran.mc(fishnet_pop_cps$net_CPS_Accepted, fishnet_Weights, nsim=999)


## ----GLOBAL_MORANS_PERMUTATION_plot--------------------------------------------------------------
GLOBAL_MORANS_PERMUTATION_plot <- ggplot(data.frame(res = globalMorans$res)[1:999,,0], aes(res)) + 
  geom_histogram(binwidth = 0.01) +
  geom_vline(aes(xintercept = globalMorans$statistic), colour = "red",size=1) +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(title="Observed and permuted Moran's I", 
       x = "Simulated Moran's I Value") +
  plotTheme()


## ----Morans_i_p_join-----------------------------------------------------------------------------
fishnet_pop_cps_morans <- fishnet_pop_cps
fishnet_pop_cps_morans$Ii <- localMorans$Ii
fishnet_pop_cps_morans$pvalue <- localMorans$`Pr(z > 0)`
fishnet_pop_cps_morans <- cbind(fishnet_coords, fishnet_pop_cps_morans)


## ----MORANS_I_P_plot-----------------------------------------------------------------------------
fishnet_pop_cps_morans_cut <- make_cuts(fishnet_pop_cps_morans, "net_CPS_Accepted",
                                        cuts = "breaks", n_breaks = 10)

plot_cps <- ggmap(cps_base_map) +
  geom_sf(data = ll(fishnet_pop_cps_morans_cut), aes(fill = cut_val),
          color = NA, inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_viridis_d(na.value=NA, name = "Maltreatment\nEvents") +
  labs(title = "Panel 1",
       subtitle = "CPS count by fishnet") +
  mapTheme() +
  theme(plot.title = element_text(size = 14, family = "sans", face = "plain", hjust = 0),
        plot.subtitle=element_text(size = 11, family = "sans", hjust = 0),
        plot.caption=element_text(size = 10, family = "sans", face = "italic", hjust = 0),
        axis.line = element_blank(),
        legend.title = element_text(size = 10, family = "sans"),
        legend.text = element_text(size = 9, family = "sans"))

Ii_cut <- fishnet_pop_cps_morans %>%
  mutate(Ii_cut_val = as.character(Hmisc::cut2(.$Ii, 
                                               cuts = as.numeric(quantile(round(fishnet_pop_cps_morans$Ii,2), 
                                                                          na.rm=T, p = seq(0,1,0.25))))))
plot_Ii <- ggmap(cps_base_map) +
  geom_sf(data = ll(Ii_cut), aes(fill = Ii_cut_val),
          color = NA, inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_viridis_d(na.value=NA, name = "I value", option = "D") +
  labs(title = "Panel 2",
       subtitle = "Local Moran's I value") +
  mapTheme() +
  theme(plot.title = element_text(size = 14, family = "sans", face = "plain", hjust = 0),
        plot.subtitle=element_text(size = 11, family = "sans", hjust = 0),
        plot.caption=element_text(size = 10, family = "sans", face = "italic", hjust = 0),
        axis.line = element_blank(),
        legend.title = element_text(size = 10, family = "sans"),
        legend.text = element_text(size = 9, family = "sans"))

p_cut <- fishnet_pop_cps_morans %>%
  mutate(pval_cut = ifelse(pvalue > 0.05, "Not\nSignificant", "Significant"))

plot_p <- ggmap(cps_base_map) +
  geom_sf(data = ll(p_cut), aes(fill = pval_cut),
          color = NA, inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_viridis_d(na.value=NA, name = "p-value", option = "D") +
  labs(title = "Panel 3",
       subtitle = "Stastically significant\nmaltreatment clusters",
       caption = "Figure 5.5") +
  mapTheme() +
  theme(plot.title = element_text(size = 14, family = "sans", face = "plain", hjust = 0),
        plot.subtitle=element_text(size = 11, family = "sans", hjust = 0),
        plot.caption=element_text(size = 10, family = "sans", face = "italic", hjust = 0),
        axis.line = element_blank(),
        legend.title = element_text(size = 10, family = "sans"),
        legend.text = element_text(size = 9, family = "sans"))

MORANS_I_P_plot <- cowplot::plot_grid(plot_cps, plot_Ii, plot_p, ncol =1, align = "hv", axis = "lrbt")
#cowplot::plot_grid(plot_cps, plot_Ii, plot_p,rel_widths = c(0.9,0.9,0.9),ncol = 1, align = "v")


## ----start_parallel_backend----------------------------------------------------------------------
cl <- makePSOCKcluster(24)
registerDoParallel(cl)


## ----aggregate_count_features--------------------------------------------------------------------
agg_results <- Aggregate_points_Features(var_list[[]], net_Cumberland)


## ----euclidean_distance_features-----------------------------------------------------------------
ED_results <- Euclidean_point_features(var_list, cumberland_rast_SP, nbr_diss,net_Cumberland)


## ----nearest_neighbor_features-------------------------------------------------------------------
NN_results <- NN_point_features(var_list, net_Cumberland, k_nearest_neighbors)


## ----stop_parallel_backend-----------------------------------------------------------------------
stopCluster(cl)


## ----sf1_features_download-----------------------------------------------------------------------
#v15 <- load_variables(2010, "sf1", cache = TRUE)
#View(v15)
vars_sf1 <- c("P001001",  # Total Population,
              "H001001", # housing units
              "H003003",  # vacant housing,
              "H004002",  # owned and mortgaged,
              "H004004",  # renter occupied
              "H003002",  # occupied housing units
              "H013001",  # Houshold size
              "P036001",  # Population in families
              "P036002",  # Population under 18yo
              "P012003",  # Male under 5yo
              "P012004",  # Male under 5-9yo
              "P012005",  # Male under 10-14yo
              "P012006",  # Male under 15-17yo
              "P012027",  # Female under 5yo
              "P012028",  # Female under 5-9yo
              "P012029",  # Female under 10-14yo
              "P012030"   # Female under 15-17yo
)
vars_sf1_desc <- c("Total Pop",
                   "Housing units",
                   "Housing, vacant",
                   "Housing, owned",
                   "Housing, rented",
                   "Housing, occupied units",
                   "Houshold size",
                   "Families Pop",
                   "Pop, under 18y",
                   "Male under 5y",
                   "Male under 5-9y",
                   "Male under 10-14y",
                   "Male under 15-17y",
                   "Female under 5y",
                   "Female under 5-9y",
                   "Female under 10-14y",
                   "Female under 15-17y")
vars_names <- data.frame(variable = vars_sf1, var_name = vars_sf1_desc, stringsAsFactors = FALSE)
census_api_key('d68574340a0011206cf42d4331abca4bf602d73a', overwrite = TRUE, install = TRUE)
library(censusapi)
Sys.setenv(CENSUS_KEY="d68574340a0011206cf42d4331abca4bf602d73a")
# Reload .Renviron
readRenviron("~/.Renviron")
# Check to see that the expected key is output in your R console

key = Sys.getenv("d68574340a0011206cf42d4331abca4bf602d73a")
NJ_block_sf1 <-get_decennial(key = key, geography = "block", variables = vars_sf1, year = 2010,
      summary_var = "P001001", state = 'NJ', county = 'Cumberland County', geometry = TRUE) %>%
st_transform(crs = 'WGS84')

NJ_block_sf1<- get_acs(
  geography='block group',
  variables=vars_sf1,
  state='NJ',
  county='Cumberland County',
  year=2010,
  summary_var = "P001001",
  survey='sf1',
  show_call=TRUE,
  geometry=TRUE
) %>% st_transform(crs=nbr_crs)

richmond_block_sf1 <- st_read("richmond_block_sf1.shp") %>% 
  rename(variable = variabl,
         summary_var = smmry_v)

NJ_block_sf1 <- st_transform(nbr1) %>% 
  rename(variable = variabl,
         summary_var = smmry_v)

NJ_block <- st_read("Data/County_Boundaries_of_NJ.shp")    #was having issues with above code when I knit so pulled in from the PAP project file
NJ_block <- st_transform(nbr1)

sf1 <- NJ_block_sf1 %>%
  spread(., variable, value)


## ----sf1_features_create-------------------------------------------------------------------------
sf1_block <- sf1 %>%
  mutate(acre = as.numeric(st_area(sf1)*2.29568e-5))
sf1_block= st_transform(sf1_block, crs='WGS 84 ') 
ggmap(cps_base_map) +
  geom_sf(
    data = sf1_block, 
    aes(fill = summary_value), 
    inherit.aes = FALSE, color = NA, alpha = 0.8
  )
crs(sf1_block)
crs(net_Cumberland)
net_blocks_intersect1 <- st_intersection(sf1_block, net_Cumberland) 

# group by cell and calc block stats.
net_blocks_intersect2 <- net_blocks_intersect1 %>%
  mutate(intersect_area_acres = as.numeric(st_area(net_blocks_intersect1)*2.29568e-5)) %>%
  group_by(net_id) %>%
  mutate(cnt = n(),
         pcnt_of_block = intersect_area_acres/acre) %>%
  # intersect_pop = value * pcnt_of_block) %>%
  arrange(net_id) %>%
  mutate_at(vars(matches("^P|^H")), tibbble::lst(.* pcnt_of_block))
ggmap(cps_base_map) +
  geom_sf(
    data = net_blocks_intersect2, 
    aes(fill = pcnt_of_block ), 
    inherit.aes = FALSE, color = NA, alpha = 0.8
  )
ggmap(cps_base_map) +
  geom_sf(
    data = net_blocks_intersect2, 
    aes(fill = cnt), 
    inherit.aes = FALSE, color = NA, alpha = 0.8
  )
### summarise intersect pops to each net cell and create pop rates for some
fishnet_sf1 <- net_blocks_intersect2 %>% # xcc
  group_by(net_id) %>%
  summarise_at(vars(matches("^P|^H")), funs(sum)) %>%
  dplyr::select(-pcnt_of_block) %>%
  rename_at(vars(vars_sf1), function(x) vars_sf1_desc) %>%
  mutate(`Pop, under 5y`     = rowSums(st_drop_geometry(.[grep("5y$", names(.))])),
         `Pop, under 5-9y`   = rowSums(st_drop_geometry(.[grep("5-9y$", names(.))])),
         `Pop, under 10-14y` = rowSums(st_drop_geometry(.[grep("10-14y$", names(.))])),
         `Pop, under 15-17y` = rowSums(st_drop_geometry(.[grep("15-17y$", names(.))])),
         `Pop, under 10y`    = `Pop, under 5y` + `Pop, under 5-9y`,
         `Pop, 10-17y`       = `Pop, under 10-14y` + `Pop, under 15-17y`) %>%
  mutate_at(vars(matches("^Male|^Female|^Pop")), 
            funs(rate = divide_by(.,`Total Pop`/100)))

## cast data frame to list of variables
sf1_results <- fishnet_sf1 %>%
  gather(variable, value, -net_id, -geometry) %>%
  mutate(feature_name = paste0("SF1_",variable)) %>%
  group_by(variable) %>%
  nest() %>%
  pull(data)
names(sf1_results) <- paste0("SF1_",setdiff(colnames(fishnet_sf1), c("net_id","geometry")))

ggmap(cps_base_map) +
  geom_sf(
    data = sf1_results[[1]],
    aes(fill = geometry),
    inherit.aes = FALSE, color = NA, alpha = 0.8
  )
## ----fishnet_pop_crs_net-------------------------------------------------------------------------
fishnet_pop_cps_net <- incidents_grp_net %>%
  dplyr::select(net_id, net_pop, incident_count, den) %>%
  rename(cps_rate = den,
         cps_net  = incident_count)
#Here JUST USE INCIDENTS_NET

## ----nearest_neighbor_feature_combine------------------------------------------------------------
sanity6 <- ggmap(cps_base_map) +
  geom_sf(
    data = net_NN,
    aes(fill = feature_name),
    inherit.aes = FALSE, color = NA, alpha = 0.8
  )


features <- data.frame(net_id = net_NN$net_id, stringsAsFactors = FALSE)
for(i in  seq_along(net_NN)){
  feat_i <- net_NN%>%
    print()
    st_drop_geometry() %>%
    dplyr::select(net_id, feature_name, value) %>%
    spread(feature_name, value)
  features <- left_join(features, feat_i, by = "net_id")
}
# join features to our target of cps_rate
NN_features <- feat_i %>%
  left_join(., st_drop_geometry(fishnet_pop_cps_net), by = "net_id") 


## ----euclidean_distance_feature_combine----------------------------------------------------------
features <- data.frame(net_id = ED_results[[1]][[1]]$net_id, stringsAsFactors = FALSE)
for(i in  seq_along(ED_results[[1]])){
  feat_i <- ED_results[[1]][[i]] %>%
    st_drop_geometry() %>%
    dplyr::select(net_id, feature_name, value = mean_dist ) %>% ### mean_dist  !!!
    spread(feature_name, value)
  features <- left_join(features, feat_i, by = "net_id")
}
# join features to our target of cps_rate
ED_features <- features %>%
  left_join(., st_drop_geometry(fishnet_pop_cps_net), by = "net_id")


## ----agg_feature_combine-------------------------------------------------------------------------
features <- data.frame(net_id = agg_results[[1]]$net_id, stringsAsFactors = FALSE)
for(i in  seq_along(ED_results[[1]])){
  feat_i <- agg_results[[i]] %>%
    st_drop_geometry() %>%
    dplyr::select(net_id, feature_name, value) %>%
    spread(feature_name, value)
  features <- left_join(features, feat_i, by = "net_id")
}
# join features to our target of cps_rate
agg_features <- features %>%
  left_join(., st_drop_geometry(incidents_grp_net), by = "net_id")


## ----sf1_feature_combine-------------------------------------------------------------------------
features <- data.frame(net_id = sf1_results[[1]]$net_id, stringsAsFactors = FALSE)
for(i in  seq_along(sf1_results)){
  feat_i <- sf1_results[[i]] %>%
    st_drop_geometry() %>%
    dplyr::select(net_id, feature_name, value) %>%
    spread(feature_name, value)
  features <- left_join(features, feat_i, by = "net_id")
}
# join features to our target of cps_rate
sf1_features <- features %>%
  left_join(., st_drop_geometry(fishnet_pop_cps_net), by = "net_id")

ggmap(cps_base_map) +
  geom_sf(
    data = sf1_features,
    aes(fill = inc_per_100),
    inherit.aes = FALSE, color = NA, alpha = 0.8
  )

## ----corr_feature_remove_NA----------------------------------------------------------------------
cor_NN_features <- NN_features %>%
  mutate_all(list(replace(., is.na(.), 0))) %>%
  dplyr::select(-net_id)

cor_agg_features <- agg_features %>%
  mutate_all(funs(replace(., is.na(.), 0))) %>%
  dplyr::select(-net_id)

cor_ED_features <- ED_features %>%
  mutate(cps_rate = ifelse(is.na(cps_rate),0,cps_rate),
         net_pop = ifelse(is.na(net_pop),0,net_pop)) %>%
  na.omit() %>%
  dplyr::select(-net_id)

cor_sf1_features <- sf1_features %>%
  mutate(cps_rate = ifelse(is.na(cps_rate),0,cps_rate),
         net_pop = ifelse(is.na(net_pop),0,net_pop)) %>%
  na.omit() %>%
  dplyr::select(-net_id)


## ----combine_ALL_FEATURES------------------------------------------------------------------------
ALL_FEATURES.cor = cor(ALL_FEATURES)


ALL_FEATURES <- full_join(NN_features, agg_features, by = "net_id") %>%
  full_join(.,ED_features, by = "net_id") %>%
  full_join(.,sf1_features, by = "net_id")
all.equal(ALL_FEATURES$cps_rate.x, ALL_FEATURES$cps_rate.y, 
          ALL_FEATURES$cps_rate.x.x, ALL_FEATURES$cps_rate.y.y)

ALL_FEATURES <- full_join(NN_features, sf1_features, by = "net_id", keep=TRUE) %>%
  #full_join(.,NN_features, by = "net_id") %>%
  #full_join(.,sf1_features, by = "net_id")
all.equal(ALL_FEATURES$cps_rate.x, ALL_FEATURES$cps_rate.y
          )

NN_CPS_Accepted <- ALL_FEATURES$NN_CPS_Accepted

ALL_FEATURES <- ALL_FEATURES %>%
  #select.list(-cps_rate.y ,
               # -cps_net.y,
                #-net_pop.y) %>%
  #dplyr::select(-contains("_CPS_")) %>%
  dplyr::rename(cps_net  = cps_net.x,
                cps_rate = cps_rate.x,
                net_pop  = net_pop.x) %>%
  mutate_all(funs(replace(., is.na(.), 0)))  %>%
  dplyr::rename_all(funs(make.names(.)))
## add NN_CPS_Accepted back in to ALL_FEATURES
ALL_FEATURES$NN_CPS_Accepted <- NN_CPS_Accepted


## ----corr_line_feature_data----------------------------------------------------------------------
cps_cor_ALL <- cor(ALL_FEATURES)
All_cors <- cps_cor_ALL[,"cps_net"]

p.mat_ALL <- cor.mtest(ALL_FEATURES)$p
p.mat_ALL <- p.mat_ALL[,which(colnames(cps_cor_ALL)=="cps_net")]

cor_ALL_plot <- data.frame(feature = names(All_cors), 
                           cor = as.numeric(All_cors),
                           p_value   = p.mat_ALL) %>%
  filter(!(feature %in% c("cps_rate","cps_net","net_pop","net_cps","net_id"))) %>%
  filter(!(feature %in% grep("CPS", names(All_cors),value=T))) %>%
  arrange(desc(cor)) %>% 
  mutate(p_value = ifelse(p_value >= 0.1, "Not Significant", "Significant"))

cor_ALL_plot$feature <- factor(cor_ALL_plot$feature,
                               levels=cor_ALL_plot[order(cor_ALL_plot$cor,
                                                         decreasing=F),]$feature)


## ----CORR_LINE_POSITIVE_FEATURE_plot-------------------------------------------------------------
CORR_LINE_POSITIVE_FEATURE_plot <- ggplot(dplyr::filter(cor_ALL_plot,cor >= 0), 
                                          aes(x = feature, y = cor, color = factor(p_value))) +
  geom_segment(aes(x = feature, y = 0, xend = feature, yend = cor), color = "grey50") +
  geom_point() +
  coord_flip() +
  scale_color_discrete(name = "p-value") +
  theme_bw()
plot(CORR_LINE_POSITIVE_FEATURE_plot)

## ----CORR_LINE_NEGATIVE_FEATURE_plot-------------------------------------------------------------
CORR_LINE_NEGATIVE_FEATURE_plot <- ggplot(dplyr::filter(cor_ALL_plot,cor <= 0), 
                                          aes(x = feature, y = cor, color = factor(p_value))) +
  geom_segment(aes(x = feature, y = 0, xend = feature, yend = cor), color = "grey50") +
  geom_point() +
  coord_flip() +
  scale_color_discrete(name = "p-value") +
  theme_bw()

plot(CORR_LINE_NEGATIVE_FEATURE_plot)
## ----corr_features_strong------------------------------------------------------------------------
features_cor <- cor_ALL_plot %>%
  mutate(feature = as.character(feature)) %>%
  arrange(desc(cor)) %>%
  pull(feature)
top_n <- head(features_cor,10)
bottom_n <- tail(features_cor,10)

features_strong_cor <- ALL_FEATURES %>%
  dplyr::select(top_n, bottom_n, cps_net, cps_rate, net_pop, net_id.x) %>%
  base::identity()


## ----features_protective_all---------------------------------------------------------------------
features_protective_all <- ALL_FEATURES %>%
  dplyr::select(contains("CHILD CARE SERVICE"),
                contains("CHURCHES"),
                contains("SCHOOLS-PRE-SCHOOL/KINDERGARTEN-ACADEMI"),
                contains("PHARMACIES"),
                contains("DENTISTS"),
                contains("CLINICS"),
                contains("YOUTH ORGANIZATIONS & CENTERS"),
                contains("COUNSELORS"),
                contains("CLERGY"),
                contains("SDRUG ABUSE & ADDICTION INFO & TREATMENT"),
                contains("CRISIS INTERVENTION SERVIC"),
                contains("GROUP HOMES"),
                contains("MINISTRIES-OUT REACH"),
                contains("COUNSELING SERVICES"),
                contains("ALCOHOLISM INFORMATION & TREATMENT CTRS"),
                contains("CHURCH ORGANIZATIONS"),
                contains("HOMELESS SHELTERS"),
                contains("COUNSELORS-LICENSED PROFESSIONAL"),
                NN_CPS_Accepted,
                cps_net, cps_rate, net_pop, net_id.x)


## ----features_protective_strong------------------------------------------------------------------
features_strong_protective_names <- cor_ALL_plot %>% 
  filter(feature %in% names(features_protective_all)) %>%
  mutate(prefix = str_extract(feature, "^[^_]+(?=_)"),
         suffix = str_extract(feature, "(?<=_)[^_].*"),
         feature = as.character(feature)) %>%
  group_by(suffix) %>%
  slice(which.max(abs(cor)))

features_protective_strong <- features_protective_all %>%
  dplyr::select(features_strong_protective_names$feature,
                NN_CPS_Accepted,
                cps_net, cps_rate, net_pop, net_id.x) %>%
  base::identity()


## ----features_risk_all---------------------------------------------------------------------------
features_risk_all <- ALL_FEATURES %>%
  dplyr::select(contains("BEAUTY SALONS"),
                contains("BARBERS"),
                contains("HOTELS & MOTELS"),
                contains("RUNAWAY"),
                contains("MANICURING"),
                contains("CONVENIENCE STORES"),
                contains("BARS"),
                contains("BONDS-BAILs"),
                contains("SERVICE STATIONS-GASOLINE & OIL"),
                contains("LAUNDRIES-SELF SERVICE"),
                contains("TATTOOING"),
                contains("OILS-FUEL (WHLS)"),
                contains("PAWNBROKERS"),
                contains("AUTOMOBILE REPAIRING & SERVICE"),
                contains("BEER & ALE-RETAIL"),
                contains("CAR WASHING & POLISHING"),
                contains("LIQUORS-RETAIL"),
                contains("TRUCK-REPAIRING & SERVICE"),
                contains("RAILERS-REPAIRING & SERVICE"),
                contains("AUTOMOBILE DETAIL & CLEAN-UP SERVICE"),
                contains("WINES-RETAIL"),
                contains("ARECREATIONAL VEHICLES-REPAIRING & SVC"),
                NN_CPS_Accepted,
                cps_net, cps_rate, net_pop, net_id.x)


## ----features_risk_strong------------------------------------------------------------------------
features_risk_strong_names <- cor_ALL_plot %>%
  filter(feature %in% names(features_risk_all)) %>%
  mutate(prefix = str_extract(feature, "^[^_]+(?=_)"),
         suffix = str_extract(feature, "(?<=_)[^_].*"),
         feature = as.character(feature)) %>%
  group_by(suffix) %>%
  slice(which.max(abs(cor)))

features_risk_strong <- features_risk_all %>%
  dplyr::select(features_risk_strong_names$feature,
                NN_CPS_Accepted,
                cps_net, cps_rate, net_pop, net_id.x) %>%
  base::identity()


## ----features_census_select----------------------------------------------------------------------
features_census_select <- ALL_FEATURES %>%
  dplyr::select(SF1_Pop..under.18y,
                SF1_Families.Pop,
                SF1_Housing..rented,
                SF1_Houshold.size,
                SF1_Housing..vacant,
                cps_net, cps_rate, net_pop, net_id.x)


## ----CORR_RISK_FEATURES_plot, cache = FALSE------------------------------------------------------
features_risk_strong_plot <- features_risk_strong %>%
  dplyr::select(-net_id.x)
CORR_RISK_FEATURES_plot <- feature_corrplot(features_risk_strong_plot, "Correlation of Risk Features")


## ----CORR_PROTECTIVE_FEATURES_plot, cache = FALSE------------------------------------------------
features_protective_strong_plot <- features_protective_strong %>%
  dplyr::select(-net_id.x)
CORR_PROTECTIVE_FEATURES_plot <- feature_corrplot(features_protective_strong_plot, "Correlation of Protective Features")


## ----protective_feature_CPS_vs_NN----------------------------------------------------------------
BUSI_protective_name <- var_list[["BusinessProject"]] %>%
  filter(Classification == "PROTECTIVE") %>%
  mutate(BUSTYP = paste0("BUSI_",BUSTYP)) %>%
  pull(BUSTYP) %>%
  unique()

protective_class <- c("CHILD CARE SERVICE",
                      "GROCERS-RETAIL",
                      "CHURCHES",
                      "SCHOOLS-PRE-SCHOOL/KINDERGARTEN-ACADEMIC",
                      "PHARMACIES",
                      "DENTISTS",
                      "CLINICS",
                      "YOUTH ORGANIZATIONS & CENTERS",
                      "COUNSELORS",
                      "CLERGY",
                      "DRUG ABUSE & ADDICTION INFO & TREATMENT",
                      "CRISIS INTERVENTION SERVICE",
                      "GROUP HOMES",
                      "MINISTRIES-OUT REACH",
                      "COUNSELING SERVICES",
                      "RELIGIOUS ORGANIZATIONS",
                      "ALCOHOLISM INFORMATION & TREATMENT CTRS",
                      "CHURCH ORGANIZATIONS",
                      "HOMELESS SHELTERS",
                      "COUNSELORS-LICENSED PROFESSIONAL")

BUSI_protective <- var_list[BUSI_protective_name[BUSI_protective_name %in% names(var_list)]]
resource_protective <- var_list[protective_class[protective_class %in% names(var_list)]]

protective_list <- do.call(c, list(resource_protective, BUSI_protective))


cl <- makePSOCKcluster(detectCores()-1)
registerDoParallel(cl)


rnd_protective_p_results_total <- NULL
rnd_protective_vec_total <- NULL
protective_mean_NN <- NULL
for(j in seq_along(protective_list)){
  protective_desc <- protective_list[[j]]
  cat(names(protective_list)[j], "\n")
  # need below b/c PoliceStation feature only has 5 featuresO
  k_nearest_neighbors_i <- ifelse(nrow(protective_desc) <= k_nearest_neighbors,
                                  nrow(protective_desc)-1, k_nearest_neighbors)
  protective_cps_NN <- nn_function(st_coordinates(protective_desc),
                                   st_coordinates(cps_dissolve), k_nearest_neighbors_i)
  protective_cps_NN <- mean(protective_cps_NN$value)
  
  
  rnd_protective_results <- foreach(i = seq_len(simulations),
                                    .packages = c('sf', 'dplyr', 'FNN', 'tibble', 'tidyr'),
                                    .combine  = c) %dopar% {
                                      # cat(i,"\n")
                                      if(nrow(protective_desc) > k_nearest_neighbors_i){
                                        ## b/c low count (e.g. police station) can return zero samples
                                        protective_cps_NN_rnd <- NULL
                                        while(length(protective_cps_NN_rnd) == 0){
                                          protective_cps_NN_rnd <- sf::st_sample(nbr_diss, nrow(protective_desc))
                                        }
                                        protective_cps_NN_rnd <- protective_cps_NN_rnd %>% 
                                          st_coordinates(.) %>%
                                          nn_function(., st_coordinates(cps_dissolve), k_nearest_neighbors_i)
                                      } else {
                                        protective_cps_NN_rnd <- NA
                                      }
                                    }
  rnd_protective_p_results <- data.frame(p = map_dbl(rnd_protective_results,
                                                     function(x) 1-ecdf(x)(protective_cps_NN)),
                                         Feature = names(protective_list)[j])
  rnd_protective_p_results_total <- rbind(rnd_protective_p_results_total, rnd_protective_p_results)
  
  rnd_protective_vec <- data.frame(dist = as.numeric(map_dbl(rnd_protective_results, mean)),
                                   Feature = names(protective_list)[j])
  rnd_protective_vec_total <- rbind(rnd_protective_vec_total, rnd_protective_vec)
  
  protective_mean_NN <- rbind(protective_mean_NN, data.frame(mean = protective_cps_NN,
                                                             Feature = names(protective_list)[j]))
}

stopCluster(cl)



## ----PROTECTIVE_CPS_VS_NN_plot-------------------------------------------------------------------
rnd_protective_vec_total$Feature <- factor(rnd_protective_vec_total$Feature ,
                                           levels = as.character(arrange(protective_mean_NN, mean)$Feature))
PROTECTIVE_CPS_VS_NN_plot <- ggplot(data = rnd_protective_vec_total, aes(x = dist, group = Feature, fill = Feature)) +
  geom_histogram(bins = 50) +
  geom_vline(data = protective_mean_NN, aes(xintercept = mean), size = 2) +
  #scale_x_continuous(breaks=seq(400,4500,200), labels = seq(400,4500,200)) +
  facet_wrap(~Feature, ncol = 1, scales = "free_y") +
  labs(x = paste0("Mean NN Distance (k = ",k_nearest_neighbors,")"),
       title = "Protective factors closeness test",
       caption = "Figure 5.8") +
  scale_fill_viridis_d() +
  plotTheme() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 11, family = "sans", face = "plain", hjust = 0),
    strip.background = element_rect(fill = "white")
  )


## ----PROTECTIVE_CPS_VS_NN_PVALUE_plot------------------------------------------------------------
rnd_protective_p_results_total$Feature <- factor(rnd_protective_p_results_total$Feature ,
                                                 levels = as.character(arrange(protective_mean_NN, mean)$Feature))
PROTECTIVE_CPS_VS_NN_PVALUE_plot <- ggplot(rnd_protective_p_results_total, aes(x = p, group = Feature, fill = Feature)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = 0.5) +
  scale_x_continuous(breaks=seq(0,1,0.1), labels = seq(0,1,0.1), limits = c(0,1)) +
  facet_wrap(~Feature, ncol = 1, scales = "free_y") +
  labs(x = "Probability of True Mean NN Distance Being Less Than Simulated Distance Value") +
  scale_fill_viridis_d() +
  theme_bw() +
  theme(
    legend.position = "none"
  )



## ----risk_feature_CPS_vs_NN----------------------------------------------------------------------
BUSI_risk_name <- var_list[["BusinessProject"]] %>%
  filter(Classification == "RISK") %>% 
  mutate(BUSTYP = paste0("BUSI_",BUSTYP)) %>% 
  pull(BUSTYP) %>% 
  unique()

CrimeData_risk_name <- var_list[["CrimeData"]] %>%
  filter(OFFENSE %in% c("BEAUTY SALONS",
                        "BARBERS",
                        "HOTELS & MOTELS",
                        "MANICURING",
                        "CONVENIENCE STORES",
                        "BARS",
                        "BONDS-BAIL",
                        "SERVICE STATIONS-GASOLINE & OIL",
                        "LAUNDRIES-SELF SERVICE",
                        "TATTOOING",
                        "OILS-FUEL (WHLS)",
                        "PAWNBROKERS",
                        "AUTOMOBILE REPAIRING & SERVICE",
                        "BEER & ALE-RETAIL",
                        "CAR WASHING & POLISHING",
                        "LIQUORS-RETAIL",
                        "TRUCK-REPAIRING & SERVICE",
                        "TRAILERS-REPAIRING & SERVICE",
                        "AUTOMOBILE DETAIL & CLEAN-UP SERVICE",
                        "WINES-RETAIL",
                        "RECREATIONAL VEHICLES-REPAIRING & SVC")) %>% 
  mutate(OFFENSE = paste0("CRIME_",OFFENSE)) %>% 
  pull(OFFENSE) %>% 
  unique()

VIO_risk_name <- var_list[["Violations_III_ks"]] %>%
  filter(CodeDsrp %in% c('General Violations','Unsafe Structure',
                         'Unfit Structure'))  %>% 
  mutate(CodeNbr = paste0("VIO_", CodeNbr)) %>% 
  count(CodeNbr) %>% 
  filter(n > 50) %>% 
  pull(CodeNbr) %>% 
  unique()

BUSI_risk <- var_list[BUSI_risk_name[BUSI_risk_name %in% names(var_list)]]
CRIME_risk <- var_list[CrimeData_risk_name[CrimeData_risk_name %in% names(var_list)]]
BusStops_risk <- var_list["BusStops"]
VIO_risk <- var_list[VIO_risk_name[VIO_risk_name %in% names(var_list)]]

risk_list <- do.call(c, list(BUSI_risk, CRIME_risk, BusStops_risk, VIO_risk))

cl <- makePSOCKcluster(detectCores()-1)
registerDoParallel(cl)

rnd_risk_p_results_total <- NULL
rnd_risk_vec_total <- NULL
risk_mean_NN <- NULL
for(j in seq_along(risk_list)){
  risk_desc <- risk_list[[j]]
  cat(names(risk_list)[j], "\n")
  risk_cps_NN <- nn_function(st_coordinates(risk_desc),
                             st_coordinates(cps_dissolve), k_nearest_neighbors)
  risk_cps_NN <- mean(risk_cps_NN$value)
  
  
  rnd_risk_results <- foreach(i = seq_len(simulations),
                              .packages = c('sf', 'dplyr', 'FNN', 'tibble', 'tidyr'),
                              .combine  = c) %dopar% { 
                                cat(i,"\n")
                                if(nrow(risk_desc) >= k_nearest_neighbors){
                                  risk_cps_NN_rnd <- sf::st_sample(nbr_diss, nrow(risk_desc)) %>%
                                    st_coordinates(.) %>%
                                    nn_function(., st_coordinates(cps_dissolve), k_nearest_neighbors)
                                } else {
                                  risk_cps_NN_rnd <- NA                    
                                }
                              }
  rnd_risk_p_results <- data.frame(p = map_dbl(rnd_risk_results, 
                                               function(x) 1-ecdf(x)(risk_cps_NN)), 
                                   Feature = names(risk_list)[j])
  rnd_risk_p_results_total <- rbind(rnd_risk_p_results_total, rnd_risk_p_results)
  
  rnd_risk_vec <- data.frame(dist = as.numeric(map_dbl(rnd_risk_results, mean)),
                             Feature = names(risk_list)[j])
  rnd_risk_vec_total <- rbind(rnd_risk_vec_total, rnd_risk_vec)
  
  risk_mean_NN <- rbind(risk_mean_NN, data.frame(mean = risk_cps_NN, 
                                                 Feature = names(risk_list)[j]))
}

stopCluster(cl)


## ----RISK_CPS_VS_NN_plot-------------------------------------------------------------------------
rnd_risk_vec_total$Feature <- factor(rnd_risk_vec_total$Feature , 
                                     levels = as.character(arrange(risk_mean_NN, mean)$Feature))
RISK_CPS_VS_NN_plot <- ggplot(data = rnd_risk_vec_total, aes(x = dist, group = Feature, fill = Feature)) +
  geom_histogram(bins = 50) +
  geom_vline(data = risk_mean_NN, aes(xintercept = mean), size = 2) +
  scale_x_continuous(breaks=seq(400,4500,200), labels = seq(400,4500,200)) +
  facet_wrap(~Feature, ncol = 1, scales = "free_y") +
  labs(x = paste0("Mean NN Distance (k = ",k_nearest_neighbors,")")) +
  scale_fill_viridis_d(option = "A") +
  theme_bw() +
  theme(
    legend.position = "none"
  )


## ----RISK_CPS_VS_NN_PVALUE_plot------------------------------------------------------------------
rnd_risk_p_results_total$Feature <- factor(rnd_risk_p_results_total$Feature , 
                                           levels = as.character(arrange(risk_mean_NN, mean)$Feature))
RISK_CPS_VS_NN_PVALUE_plot <- ggplot(rnd_risk_p_results_total, aes(x = p, group = Feature, fill = Feature)) +
  geom_histogram(bins = 50) + 
  geom_vline(xintercept = 0.5) +
  scale_x_continuous(breaks=seq(0,1,0.1), labels = seq(0,1,0.1), limits = c(0,1)) +
  facet_wrap(~Feature, ncol = 1, scales = "free_y") +
  labs(x = "Probability of True Mean NN Distance Being Less Than Simulated Distance Value") +
  scale_fill_viridis_d(option = "A") +
  theme_bw() +
  theme(
    legend.position = "none"
  )


## ----violation_feature_CPS_vs_NN-----------------------------------------------------------------
vio_vars <- var_list[grep("VIO_",names(var_list), value = TRUE)]

#################################
cl <- makePSOCKcluster(32)
registerDoParallel(cl)
#################################
rnd_vio_p_results_total <- NULL
rnd_vio_vec_total <- NULL
vio_mean_NN <- NULL
for(j in seq_along(vio_vars)){
  vio_desc <- vio_vars[[j]]
  cat(names(vio_vars)[j], "\n")
  vio_cps_NN <- nn_function(st_coordinates(vio_desc),
                            st_coordinates(cps_dissolve), k)
  vio_cps_NN <- mean(vio_cps_NN$value)
  
  
  rnd_vio_results <- foreach(i = seq_len(simulations),
                             .packages = c('sf', 'dplyr', 'FNN', 'tibble', 'tidyr'),
                             .combine  = c) %dopar% { 
                               cat(i,"\n")
                               if(nrow(vio_desc) >= k){
                                 vio_cps_NN_rnd <- sf::st_sample(nbr_diss, nrow(vio_desc)) %>%
                                   st_coordinates(.) %>%
                                   nn_function(., st_coordinates(cps_dissolve), k_nearest_neighbors)
                               } else {
                                 vio_cps_NN_rnd <- NA                    
                               }
                             }
  rnd_vio_p_results <- data.frame(p = map_dbl(rnd_vio_results, 
                                              function(x) 1-ecdf(x)(vio_cps_NN)), 
                                  Feature = names(vio_vars)[j])
  rnd_vio_p_results_total <- rbind(rnd_vio_p_results_total, rnd_vio_p_results)
  
  rnd_vio_vec <- data.frame(dist = as.numeric(map_dbl(rnd_vio_results, mean)),
                            Feature = names(vio_vars)[j])
  rnd_vio_vec_total <- rbind(rnd_vio_vec_total, rnd_vio_vec)
  
  vio_mean_NN <- rbind(vio_mean_NN, data.frame(mean = vio_cps_NN, 
                                               Feature = names(vio_vars)[j]))
}

#################################
stopCluster(cl)
#################################


## ----VIOLATION_CPS_VS_NN_plot--------------------------------------------------------------------
rnd_vio_vec_total$Feature <- factor(rnd_vio_vec_total$Feature , 
                                    levels = as.character(arrange(vio_mean_NN, mean)$Feature))
VIOLATION_CPS_VS_NN_plot <- ggplot(data = rnd_vio_vec_total, aes(x = dist, group = Feature, fill = Feature)) +
  geom_histogram(bins = 50) +
  geom_vline(data = vio_mean_NN, aes(xintercept = mean), size = 2) +
  scale_x_continuous(breaks=seq(400,4500,200), labels = seq(400,4500,200)) +
  facet_wrap(~Feature, ncol = 1, scales = "free_y") +
  labs(x = paste0("Mean NN Distance (k = ",k_nearest_neighbors,")"),
       title = "Violation-related risk factors closeness test",
       caption = "Figure 5.7") +
  scale_fill_viridis_d() +
  plotTheme() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 11, family = "sans", face = "plain", hjust = 0),
    strip.background = element_rect(fill = "white")
  )


## ----VIOLATION_CPS_VS_NN_PVALUE_plot-------------------------------------------------------------
rnd_vio_p_results_total$Feature <- factor(rnd_vio_p_results_total$Feature , 
                                          levels = as.character(arrange(vio_mean_NN, mean)$Feature))
VIOLATION_CPS_VS_NN_PVALUE_plot <- ggplot(rnd_vio_p_results_total, aes(x = p, group = Feature, fill = Feature)) +
  geom_histogram(bins = 50) + 
  geom_vline(xintercept = 0.5) +
  scale_x_continuous(breaks=seq(0,1,0.1), labels = seq(0,1,0.1), limits = c(0,1)) +
  facet_wrap(~Feature, ncol = 1, scales = "free_y") +
  labs(x = "Probability of True Mean NN Distance Being Less Than Simulated Distance Value") +
  scale_fill_viridis_d(option = "A") +
  theme_bw() +
  theme(
    legend.position = "none"
  )


## ----crime_feature_CPS_vs_NN---------------------------------------------------------------------
crime_vars <- var_list[grep("CRIME_",names(var_list), value = TRUE)]

#################################
cl <- makePSOCKcluster(32)
registerDoParallel(cl)
#################################
rnd_crime_p_results_total <- NULL
rnd_crime_vec_total <- NULL
crime_mean_NN <- NULL
for(j in seq_along(crime_vars)){
  crime_desc <- crime_vars[[j]]
  cat(names(crime_vars)[j], "\n")
  crime_cps_NN <- nn_function(st_coordinates(crime_desc),
                              st_coordinates(cps_dissolve), k)
  crime_cps_NN <- mean(crime_cps_NN$value)
  
  
  rnd_crime_results <- foreach(i = seq_len(simulations),
                               .packages = c('sf', 'dplyr', 'FNN', 'tibble', 'tidyr'),
                               .combine  = c) %dopar% { 
                                 cat(i,"\n")
                                 if(nrow(crime_desc) >= k){
                                   crime_cps_NN_rnd <- sf::st_sample(nbr_diss, nrow(crime_desc)) %>%
                                     st_coordinates(.) %>%
                                     nn_function(., st_coordinates(cps_dissolve), k_nearest_neighbors)
                                 } else {
                                   crime_cps_NN_rnd <- NA                    
                                 }
                               }
  rnd_crime_p_results <- data.frame(p = map_dbl(rnd_crime_results, 
                                                function(x) 1-ecdf(x)(crime_cps_NN)), 
                                    Feature = names(crime_vars)[j])
  rnd_crime_p_results_total <- rbind(rnd_crime_p_results_total, rnd_crime_p_results)
  
  rnd_crime_vec <- data.frame(dist = as.numeric(map_dbl(rnd_crime_results, mean)),
                              Feature = names(crime_vars)[j])
  rnd_crime_vec_total <- rbind(rnd_crime_vec_total, rnd_crime_vec)
  
  crime_mean_NN <- rbind(crime_mean_NN, data.frame(mean = crime_cps_NN, 
                                                   Feature = names(crime_vars)[j]))
}

#################################
stopCluster(cl)
#################################


## ----CRIME_CPS_VS_NN_plot------------------------------------------------------------------------
rnd_crime_vec_total$Feature <- factor(rnd_crime_vec_total$Feature , 
                                      levels = as.character(arrange(crime_mean_NN, mean)$Feature))
CRIME_CPS_VS_NN_plot <- ggplot(data = rnd_crime_vec_total, aes(x = dist, group = Feature, fill = Feature)) +
  geom_histogram(bins = 50) +
  geom_vline(data = crime_mean_NN, aes(xintercept = mean), size = 2) +
  scale_x_continuous(breaks=seq(400,4500,200), labels = seq(400,4500,200)) +
  facet_wrap(~Feature, ncol = 1, scales = "free_y") +
  labs(x = paste0("Mean NN Distance (k = ",k_nearest_neighbors,")"),
       title = "Crime-related risk factors closeness test",
       caption = "Figure 5.6") +
  scale_fill_viridis_d() +
  plotTheme() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 11, family = "sans", face = "plain", hjust = 0),
    strip.background = element_rect(fill = "white")
  )


## ----CRIME_CPS_VS_NN_PVALUE_plot-----------------------------------------------------------------
rnd_crime_p_results_total$Feature <- factor(rnd_crime_p_results_total$Feature , 
                                            levels = as.character(arrange(crime_mean_NN, mean)$Feature))
CRIME_CPS_VS_NN_PVALUE_plot <- ggplot(rnd_crime_p_results_total, aes(x = p, group = Feature, fill = Feature)) +
  geom_histogram(bins = 50) + 
  geom_vline(xintercept = 0.5) +
  scale_x_continuous(breaks=seq(0,1,0.1), labels = seq(0,1,0.1), limits = c(0,1)) +
  facet_wrap(~Feature, ncol = 1, scales = "free_y") +
  labs(x = "Probability of True Mean NN Distance Being Less Than Simulated Distance Value") +
  scale_fill_viridis_d(option = "A") +
  theme_bw() +
  theme(
    legend.position = "none"
  )


## ----busifeature_CPS_vs_NN-----------------------------------------------------------------------
busi_vars <- var_list[grep("BUSI_",names(var_list), value = TRUE)]

#################################
cl <- makePSOCKcluster(32)
registerDoParallel(cl)
#################################
rnd_busi_p_results_total <- NULL
rnd_busi_vec_total <- NULL
busi_mean_NN <- NULL
for(j in seq_along(busi_vars)){
  busi_desc <- busi_vars[[j]]
  cat(names(busi_vars)[j], "\n")
  busi_cps_NN <- nn_function(st_coordinates(busi_desc),
                             st_coordinates(cps_dissolve), k)
  busi_cps_NN <- mean(busi_cps_NN$value)
  
  
  rnd_busi_results <- foreach(i = seq_len(simulations),
                              .packages = c('sf', 'dplyr', 'FNN', 'tibble', 'tidyr'),
                              .combine  = c) %dopar% { 
                                cat(i,"\n")
                                if(nrow(busi_desc) >= k){
                                  busi_cps_NN_rnd <- sf::st_sample(nbr_diss, nrow(busi_desc)) %>%
                                    st_coordinates(.) %>%
                                    nn_function(., st_coordinates(cps_dissolve), k_nearest_neighbors)
                                } else {
                                  busi_cps_NN_rnd <- NA                    
                                }
                              }
  rnd_busi_p_results <- data.frame(p = map_dbl(rnd_busi_results, 
                                               function(x) 1-ecdf(x)(busi_cps_NN)), 
                                   Feature = names(busi_vars)[j])
  rnd_busi_p_results_total <- rbind(rnd_busi_p_results_total, rnd_busi_p_results)
  
  rnd_busi_vec <- data.frame(dist = as.numeric(map_dbl(rnd_busi_results, mean)),
                             Feature = names(busi_vars)[j])
  rnd_busi_vec_total <- rbind(rnd_busi_vec_total, rnd_busi_vec)
  
  busi_mean_NN <- rbind(busi_mean_NN, data.frame(mean = busi_cps_NN, 
                                                 Feature = names(busi_vars)[j]))
}

#################################
stopCluster(cl)
#################################


## ----BUSINESS_CPS_VS_NN_plot---------------------------------------------------------------------
rnd_busi_vec_total$Feature <- factor(rnd_busi_vec_total$Feature , 
                                     levels = as.character(arrange(busi_mean_NN, mean)$Feature))
BUSINESS_CPS_VS_NN_plot <- ggplot(data = rnd_busi_vec_total, aes(x = dist, group = Feature, fill = Feature)) +
  geom_histogram(bins = 50) +
  geom_vline(data = busi_mean_NN, aes(xintercept = mean), size = 2) +
  #scale_x_continuous(breaks=seq(400,4500,200), labels = seq(400,4500,200)) +
  facet_wrap(~Feature, ncol = 1, scales = "free_y") +
  labs(x = paste0("Mean NN Distance (k = ",k_nearest_neighbors,")"),
       title = "Business-related features closeness test",
       caption = "Figure 5.8") +
  scale_fill_viridis_d() +
  plotTheme() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 11, family = "sans", face = "plain", hjust = 0),
    strip.background = element_rect(fill = "white")
  )


## ----BUSINESS_CPS_VS_NN_PVALUE_plot--------------------------------------------------------------
rnd_busi_p_results_total$Feature <- factor(rnd_busi_p_results_total$Feature , 
                                           levels = as.character(arrange(busi_mean_NN, mean)$Feature))
BUSINESS_CPS_VS_NN_PVALUE_plot <- ggplot(rnd_busi_p_results_total, aes(x = p, group = Feature, fill = Feature)) +
  geom_histogram(bins = 50) + 
  geom_vline(xintercept = 0.5) +
  scale_x_continuous(breaks=seq(0,1,0.1), labels = seq(0,1,0.1), limits = c(0,1)) +
  facet_wrap(~Feature, ncol = 1, scales = "free_y") +
  labs(x = "Probability of True Mean NN Distance Being Less Than Simulated Distance Value") +
  scale_fill_viridis_d(option = "A") +
  theme_bw() +
  theme(
    legend.position = "none"
  )


## ----features_prep-------------------------------------------------------------------------------
target_var <- "cps_net"
features_protective_strong2 <- dplyr::select(features_protective_strong, -cps_rate, -net_pop)
features_risk_strong2 <- dplyr::select(features_risk_strong, -cps_rate, -net_pop)
features_census_select2     <- dplyr::select(features_census_select, -cps_rate, -net_pop)


## ----model_data_prep-----------------------------------------------------------------------------
og_dat <- full_join(features_risk_strong, features_census_select, by = "net_id.x") %>%
  full_join(., features_protective_strong, by = "net_id.x") %>% 
  dplyr::select(-net_pop.y, -cps_net.y, -cps_rate.y,
                -net_pop.x, -cps_net.x, -cps_rate.x,
                -NN_CPS_Accepted.y) %>% 
  rename("NN_CPS_Accepted" = NN_CPS_Accepted.x)

dat    <- og_dat %>% dplyr::select(-cps_rate, -net_pop, -net_id.x) %>%
  mutate_at(vars(-cps_net), scale_this) %>%
  identity() # line ender (does nothing)

net_hood <- st_join(net_Richmond, nbr, largest = TRUE)
all.equal(net_hood$net_id, og_dat$net_id)
og_dat$.block_id <- net_hood$name


## ----neighborhood_fixed_effects------------------------------------------------------------------
hood_matrix <- model.matrix(cps_net~.block_id,og_dat)
hood_model <- lm(sqrt(og_dat$cps_net) ~ hood_matrix)
dat$hood_fixed <- predict(hood_model, type = "response")^2
og_dat$hood_fixed <- predict(hood_model, type = "response")^2


## ----create_cv_fold_tibble-----------------------------------------------------------------------
all_hoods <- length(unique(net_hood$name))
n_folds = ifelse(n_folds == "LOOCV", all_hoods, n_folds)
folds_index <- groupdata2::fold(og_dat, k = n_folds, id_col = '.block_id')$.folds

cv_tbl <- tibble(folds = seq_len(n_folds),
                 train = NA, train_y = NA, train_index = NA, train_net_id = NA,
                 test  = NA, test_y  = NA, test_index  = NA, test_net_id  = NA)
for(k in seq_len(n_folds)){
  fold_i  <- which(folds_index == k)
  cv_tbl[k,]$train         <- list(dat[-fold_i,])
  cv_tbl[k,]$test          <- list(dat[ fold_i,])
  cv_tbl[k,]$train_y       <- list(og_dat[-fold_i,target_var])
  cv_tbl[k,]$test_y        <- list(og_dat[ fold_i,target_var])
  cv_tbl[k,]$train_index   <- list(setdiff(seq(1:nrow(dat)),fold_i))
  cv_tbl[k,]$test_index    <- list(fold_i)
  cv_tbl[k,]$train_net_id  <- list(og_dat[-fold_i,"net_id"])
  cv_tbl[k,]$test_net_id   <- list(og_dat[ fold_i,"net_id"])
}


## ----NEIGHBORHOOD_FOLDS_plot---------------------------------------------------------------------
cv_sf <- left_join(og_dat, net_Richmond, by = "net_id") %>%
  st_as_sf() %>%
  dplyr::select(.block_id)
NEIGHBORHOOD_FOLDS_plot <- plot(cv_sf)


## ----Poisson_regression--------------------------------------------------------------------------
po_cv_tbl <- cv_tbl %>%
  mutate(fit   = map(train, glm_fit, 
                     formula =  paste("cps_net ~ ."), 
                     family = "poisson"),
         pred  = map2(fit, test, lm_predict, sqrt = FALSE),
         mdl_nam = "GLM - Poisson") %>% 
  score_model()
cat("Test Set MAE:",mean(po_cv_tbl$MAE),"\n")
cat("Test Set logdev:",mean(po_cv_tbl$logdev, na.rm=TRUE),"\n")


## ----POISSON_REGRESSION_FIT_plot-----------------------------------------------------------------
POISSON_REGRESSION_FIT_plot <- plot_fold_pred(po_cv_tbl$pred, po_cv_tbl$test_y, type = "fit")


## ----Random_Forest_regression--------------------------------------------------------------------
rf_cv_tbl <- cv_tbl %>%
  mutate(fit   = map(train, rf_fit, formula = "cps_net ~ .", mtry_add = 2, importance = "impurity"),
         pred  = map2(fit, test, lm_predict),
         mdl_nam = "Random Forest") %>% 
  score_model()
cat("Test Set MAE:",mean(rf_cv_tbl$MAE),"\n")
cat("Test Set logdev:",mean(rf_cv_tbl$logdev, na.rm=TRUE),"\n")


## ----Radom_Forest_variable_importance_PLOT-------------------------------------------------------
varimp_dat <- data.frame(importance = rf_cv_tbl$fit[[1]]$variable.importance) %>% 
  rownames_to_column("variable")

RF_VARIMP_PLOT <- ggplot(varimp_dat, aes(x=reorder(variable,importance), y=importance, fill=importance))+ 
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  labs(y = "Variable Importance",
       x = " ", 
       title = "Feature importance",
       subtitle = "Random Forest sub-model",
       caption = "Figure 6.4") +
  guides(fill=F)+
  scale_fill_viridis_c() +
  plotTheme()



## ----RANDOM_FOREST_FIT_plot----------------------------------------------------------------------
RANDOM_FOREST_FIT_plot <- plot_fold_pred(rf_cv_tbl$pred, rf_cv_tbl$test_y, type = "fit")


## ----spatial_error_regression--------------------------------------------------------------------
spat_durbin <- errorsarlm(sqrt(cps_net) ~ ., data = dat, listw, etype ="emixed")
spat_durbin_tbl <- tibble(
  fit   = list(spat_durbin),
  pred  = map(fit, sar_pred),
  test_y= list(dat$cps_net),
  test_net_id = list(og_dat$net_id),
  mdl_nam = "Spatial Durbin - sqrt") %>% 
  score_model()
cat("Test Set MAE:",mean(spat_durbin_tbl$MAE),"\n")
cat("Test Set logdev:",mean(spat_durbin_tbl$logdev, na.rm=TRUE),"\n")


## ----SPATIAL_ERROR_FIT_plot----------------------------------------------------------------------
SPATIAL_ERROR_FIT_plot <- plot_fold_pred(spat_durbin_tbl$pred, dat$cps_net, type = "fit")


## ----gather_OOF_predictions----------------------------------------------------------------------
po_pred_dat <- po_cv_tbl %>%
  unnest(pred) %>%
  mutate(test_y = po_cv_tbl %>% unnest(test_y) %>% pull(test_y),
         test_net_id = po_cv_tbl %>% unnest(test_net_id) %>% pull(test_net_id))

po_pred_geoplot <- model_pred_geoplot(po_pred_dat$pred,
                                      po_pred_dat$test_y,
                                      po_pred_dat$test_net_id,
                                      net_Richmond, cps_base_map, "po")

rf_pred_dat <- rf_cv_tbl %>%
  unnest(pred) %>%
  mutate(test_y = rf_cv_tbl %>% unnest(test_y) %>% pull(test_y),
         test_net_id = rf_cv_tbl %>% unnest(test_net_id) %>% pull(test_net_id))

rf_pred_geoplot <- model_pred_geoplot(rf_pred_dat$pred,
                                      rf_pred_dat$test_y,
                                      rf_pred_dat$test_net_id,
                                      net_Richmond, cps_base_map,
                                      "Random Forest")

sarlm_pred_dat <- spat_durbin_tbl %>%
  unnest(pred) %>%
  mutate(test_y = spat_durbin_tbl %>% unnest(test_y) %>% pull(test_y),
         test_net_id = spat_durbin_tbl %>% unnest(test_net_id) %>% pull(test_net_id))

sarlm_pred_geoplot <- model_pred_geoplot(sarlm_pred_dat$pred,
                                         sarlm_pred_dat$test_y,
                                         sarlm_pred_dat$test_net_id,
                                         net_Richmond, cps_base_map,
                                         "SARLM")


## ----join_model_predictions----------------------------------------------------------------------
cps_preds <- og_dat %>% 
  dplyr::select(net_id, cps_net) %>% 
  left_join(., dplyr::select(po_pred_dat,
                             net_id = test_net_id,
                             pred_lm = pred), by = "net_id") %>%
  left_join(., dplyr::select(rf_pred_dat, 
                             net_id = test_net_id,
                             pred_rf = pred), by = "net_id") %>% 
  left_join(., dplyr::select(sarlm_pred_dat, 
                             net_id = test_net_id,
                             pred_sarlm = pred), by = "net_id") %>% 
  mutate_if(is.double, round, 2)


## ----meta_model_stacking-------------------------------------------------------------------------
if(all.equal(cps_preds$net_id, net_hood$net_id)){
  cat("Predictions and spatial data are in same order, GOOD to go!", "\n")
} else {
  cat("There is a PROBLEM with order of predictions and spatial data; Likely Errors!","\n")
}

cps_preds_cv_dat <- dplyr::select(cps_preds, -net_id)
ens_cv_tbl <- tibble(folds = seq_len(n_folds),
                     train = NA, train_y = NA, train_index = NA, train_net_id = NA,
                     test  = NA, test_y  = NA, test_index  = NA, test_net_id  = NA)
for(k in seq_len(n_folds)){
  fold_i  <- which(folds_index == k)
  ens_cv_tbl[k,]$train         <- list(cps_preds_cv_dat[-fold_i,])
  ens_cv_tbl[k,]$test          <- list(cps_preds_cv_dat[ fold_i,])
  ens_cv_tbl[k,]$train_y       <- list(cps_preds_cv_dat[-fold_i,target_var])
  ens_cv_tbl[k,]$test_y        <- list(cps_preds_cv_dat[ fold_i,target_var])
  ens_cv_tbl[k,]$train_index   <- list(setdiff(seq(1:nrow(cps_preds_cv_dat)),fold_i))
  ens_cv_tbl[k,]$test_index    <- list(fold_i)
  ens_cv_tbl[k,]$train_net_id  <- list(cps_preds[-fold_i,"net_id"])
  ens_cv_tbl[k,]$test_net_id   <- list(cps_preds[ fold_i,"net_id"])
}

ens_cv_tbl <- ens_cv_tbl %>%
  mutate(fit   = map(train, rf_fit, formula = "cps_net ~ pred_rf + pred_sarlm"),
         pred  = map2(fit, test, lm_predict),
         # pred  = map(pred, round),
         mdl_nam = "Meta-Model") %>% 
  score_model()

cat("Test Set MAE:",mean(ens_cv_tbl$MAE),"\n")
cat("Test Set logdev:",mean(ens_cv_tbl$logdev),"\n")


## ----META_MODEL_FIT_plot-------------------------------------------------------------------------
META_MODEL_FIT_plot <- plot_fold_pred(ens_cv_tbl$pred, ens_cv_tbl$test_y, type = "fit") +
  labs(x = "Observed Maltreatment Counts",
       y = "Predicted Maltreatment Counts",
       title = "Predicted vs. observed maltreatment counts",
       caption = "Figure 1.7") +
  plotTheme() +
  theme(panel.border = element_blank())


## ----join_meta_model_predictions-----------------------------------------------------------------
ens_pred_dat <- ens_cv_tbl %>% 
  unnest(pred) %>% 
  mutate(test_y = ens_cv_tbl %>% unnest(test_y) %>% pull(test_y),
         test_net_id = ens_cv_tbl %>% unnest(test_net_id) %>% pull(test_net_id)) 

ens_pred_geoplot <- model_pred_geoplot(ens_pred_dat$pred, 
                                       ens_pred_dat$test_y, 
                                       ens_pred_dat$test_net_id,
                                       net_Richmond, cps_base_map, 
                                       "Meta-Model")
cps_preds2 <- cps_preds %>% 
  left_join(., dplyr::select(ens_pred_dat, 
                             net_id = test_net_id,
                             pred_ens = pred) %>% 
              mutate(pred_ens = round(pred_ens,2)), by = "net_id") 


## ----PREDICTION_MAP_plots------------------------------------------------------------------------
POISSON_MODEL_PREDICTION_MAP_plot <- cowplot::plot_grid(po_pred_geoplot[[2]] + 
                                                          labs(title = "Poisson Regression",
                                                               subtitle = "Predicted Maltreatment Count") +
                                                          mapTheme() + 
                                                          theme(panel.border = element_blank()), 
                                                        po_pred_geoplot[[1]] + 
                                                          labs(subtitle = "MAE") +
                                                          scale_fill_viridis_d(name = "MAE") +
                                                          mapTheme() +
                                                          theme(panel.border = element_blank()), 
                                                        align = "h")
RF_MODEL_PREDICTION_MAP_plot <- cowplot::plot_grid(rf_pred_geoplot[[2]] + 
                                                     labs(title = "Random Forest",
                                                          subtitle = "Predicted Maltreatment Count") +
                                                     mapTheme() + 
                                                     theme(panel.border = element_blank()), 
                                                   rf_pred_geoplot[[1]] + 
                                                     labs(subtitle = "MAE") +
                                                     scale_fill_viridis_d(name = "MAE") +
                                                     mapTheme() + 
                                                     theme(panel.border = element_blank()), 
                                                   align = "h")
SARLM_MODEL_PREDICTION_MAP_plot <- cowplot::plot_grid(sarlm_pred_geoplot[[2]] + 
                                                        labs(title = "Spatial Durbin Model",
                                                             subtitle = "Predicted Maltreatment Count") +
                                                        mapTheme() + 
                                                        theme(panel.border = element_blank()), 
                                                      sarlm_pred_geoplot[[1]] + 
                                                        labs(subtitle = "MAE") +
                                                        scale_fill_viridis_d(name = "MAE") +
                                                        mapTheme() + 
                                                        theme(panel.border = element_blank()), 
                                                      align = "h")
META_MODEL_PREDICTION_MAP_plot <- cowplot::plot_grid(ens_pred_geoplot[[2]] + 
                                                       labs(title = "Meta-Model",
                                                            subtitle = "Predicted Maltreatment Count",
                                                            caption = "Figure 6.2") +
                                                       mapTheme() + 
                                                       theme(panel.border = element_blank()), 
                                                     ens_pred_geoplot[[1]] + 
                                                       labs(subtitle = "MAE") +
                                                       scale_fill_viridis_d(name = "MAE") +
                                                       mapTheme() + 
                                                       theme(panel.border = element_blank()), 
                                                     align = "h")


## ----model_error_by_decile-----------------------------------------------------------------------
models <- bind_rows(rf_cv_tbl, spat_durbin_tbl, ens_cv_tbl, po_cv_tbl)

CV_preds_long <- models %>%
  group_by(mdl_nam) %>%
  unnest(pred, test_y) 

## map over all quantiles to get error metrics
quantile_errors <- CV_preds_long %>%
  nest(-mdl_nam) %>%
  mutate(q      = list(seq(0,1,0.01)),
         pred   = map(data, "pred"),
         test_y = map(data, "test_y")) %>%
  dplyr::select(-data) %>%
  unnest(q, .preserve = c(pred, test_y)) %>%
  filter(q != 0) %>% 
  mutate(q_dat  = pmap(list(pred, test_y, q), quantile_error),
         q_pred = map(q_dat, "pred"),
         q_obs  = map(q_dat, "obs"),
         q_RMSE = map2_dbl(q_pred, q_obs, rmse),
         q_MAE  = map2_dbl(q_pred, q_obs, mae),
         q_logdev  = map2_dbl(q_pred, q_obs, logdev_p),
         y_max  = quantile(seq(0,max(dat$cps_net)), q),
         q_cnt  = nrow(og_dat) - map_int(q_dat, nrow))

q_error_plotdat <- quantile_errors %>%
  dplyr::select(mdl_nam, q, q_RMSE, q_MAE, q_logdev)
q_cnt_plotdat <- quantile_errors %>% 
  dplyr::select(mdl_nam, q, y_max, q_cnt) %>% 
  filter(q != 0) %>%
  mutate(q_pcnt = (q_cnt / nrow(og_dat)))
q_error_mean <- q_error_plotdat %>%
  group_by(mdl_nam) %>%
  summarise(mean_RMSE = mean(q_RMSE, na.rm = TRUE),
            mean_MAE  = mean(q_MAE, na.rm = TRUE),
            mean_logdev  = mean(q_logdev, na.rm = TRUE)) %>%
  arrange(desc(mean_logdev))
print(q_error_mean)


## ----ERROR_DECILE_plots--------------------------------------------------------------------------
LOGDEV_MODEL_ERROR_BY_DECILE_plot <- ggplot(data = q_error_plotdat, aes(x=q, y=q_logdev, group = mdl_nam, color = factor(mdl_nam))) +
  geom_line(size = 1) +
  scale_color_viridis_d(name = "Model") +
  scale_y_continuous(limits=c(0,1)) +
  labs(y = "Logarithmic Score",
       caption = "Figure 6.2 - Goodness of fit by decile") +
  plotTheme() +
  theme(legend.position = "none")

MAE_MODEL_ERROR_BY_DECILE_plot <- ggplot(data = q_error_plotdat, aes(x=q, y=q_MAE, group = mdl_nam, color = factor(mdl_nam))) +
  geom_line(size = 1) +
  scale_color_viridis_d(name = "Model") +
  labs(y = "MAE") +
  plotTheme()  +
  theme(legend.position = "none")

COUNT_BY_DECILE_plot <- ggplot(data = q_cnt_plotdat, 
                               aes(x=q, y=q_cnt, group = mdl_nam, color = factor(mdl_nam))) +
  geom_line(size = 1) +
  scale_x_continuous(breaks=seq(0,1,0.1), labels = seq(0,1,0.1)) +
  scale_color_viridis_d(name = "Model") +
  labs(y = "Number of Predictions in Each Decile",
       x = "Decile") +
  plotTheme()  +
  theme(legend.position = "none")


legend <- get_legend(LOGDEV_MODEL_ERROR_BY_DECILE_plot + plotTheme() + theme(legend.position = "right"))


## ----Model_Error_Results_table-------------------------------------------------------------------
model_results <- models %>%
  dplyr::select("Model Name" = mdl_nam, R2, RMSE, MAE, logdev) %>%
  group_by(`Model Name`) %>%
  arrange(`Model Name`) %>%
  summarise(R2_mean      = mean(R2, na.rm=TRUE),
            R2_sd        = sd(R2, na.rm=TRUE),
            MAE_mean     = mean(MAE, na.rm=TRUE),
            MAE_sd       = sd(MAE, na.rm=TRUE),
            RMSE_mean    = mean(RMSE, na.rm=TRUE),
            RMSE_sd      = sd(RMSE, na.rm=TRUE),
            logdev_mean  = mean(logdev, na.rm=TRUE),
            logdev_sd    = sd(logdev, na.rm=TRUE)) 
Model_Error_Results_table <- model_results %>%
  kable(., format = "html", digits = 3) %>%
  kable_styling()

meta_log_mean <- model_results[which(model_results$`Model Name` == "Meta-Model"),"logdev_mean",drop=TRUE]
meta_log_sd <- model_results[which(model_results$`Model Name` == "Meta-Model"),"logdev_sd",drop=TRUE]
meta_log_error <- qnorm(0.975)*meta_log_sd/sqrt(nrow(ens_cv_tbl))
meta_log_error_lower <- round(meta_log_mean - meta_log_error,3)
meta_log_error_upper <- round(meta_log_mean + meta_log_error,3)

meta_MAE_mean <- model_results[which(model_results$`Model Name` == "Meta-Model"),"MAE_mean",drop=TRUE]
meta_MAE_sd <- model_results[which(model_results$`Model Name` == "Meta-Model"),"MAE_sd",drop=TRUE]
meta_MAE_error <- qnorm(0.975)*meta_MAE_sd/sqrt(nrow(ens_cv_tbl))
meta_MAE_error_lower <- round(meta_MAE_mean - meta_MAE_error,3)
meta_MAE_error_upper <- round(meta_MAE_mean + meta_MAE_error,3)


## ----aggregate_model_errors_to_neighborhood------------------------------------------------------
error_geoplot <-  net_Richmond %>%
  left_join(., ens_pred_dat, by = c("net_id" = "test_net_id"),
            feature_name = paste0("Meta-Model", "dev")) %>%
  score_model() %>%
  mutate(dev_p_inv = 1 - logdev) %>% 
  make_cuts(., "logdev", cuts = "breaks", n_breaks = 5)

# error metrics to points
error_points <- st_centroid(error_geoplot) %>%
  dplyr::select(logdev, MAE, test_y)

# aggreate mean errors to neighborhoods
neighborhood_metric_logdev <- error_points %>%
  aggregate(., nbr, mean) %>%
  dplyr::select(logdev) %>% 
  make_cuts(., "logdev")

neighborhood_metric_MAE<- error_points %>%
  aggregate(., nbr, mean) %>%
  dplyr::select(MAE) %>% 
  mutate(MAE = round(MAE,2)) %>% 
  make_cuts(., "MAE")


## ----MODEL_ERROR_BY_NEIGHBORHOOD_plots-----------------------------------------------------------
LOGDEV_BY_NEIGHBORHOOD_plot <- make_fishnet_dist_plot(neighborhood_metric_logdev, cps_base_map, legend = "right", 
                                                      direction = 1, var_name = "Deviance", 
                                                      title = "Out-of-Fold error by neighborhood") + 
  labs(caption = "Figure 6.3",
       subtitle = "Logarithmic score") +
  mapTheme()

MAE_BY_NEIGHBORHOOD_plot <- make_fishnet_dist_plot(neighborhood_metric_MAE, cps_base_map, legend = "right", 
                                                   direction = 1, var_name = "MAE") +
  labs(subtitle = "MAE") +
  mapTheme()


## ----statistical_area_download-------------------------------------------------------------------
#Below we calculate poverty and nonWhite rates by neighborhood by converting tracts to centroids and spatial joining with 
#neighbothoods statistical areas. 

#get statarea
nbr_statAreas <- read_sf("https://data.richmondgov.com/resource/8kyq-v9j2.geojson") %>%
  st_transform(crs = 102747) %>% 
  mutate(stat_area_id = id)


## ----ACS_data_download---------------------------------------------------------------------------
#download poverty and population data
tract10 <- get_acs(geography = "tract", variables = c("B02001_001","B02001_002E","B17001_002"), 
                   year = 2010, state=51, county=760, geometry=T)


## ----ACS_Rates-----------------------------------------------------------------------------------
tract10 <- tract10 %>%
  dplyr::select(variable,estimate) %>%
  as.data.frame() %>%
  spread(variable,estimate) %>%
  rename(TotalPop=B02001_001,
         NumberWhites=B02001_002,
         TotalPoverty=B17001_002) %>%
  mutate(percentNonWhite = ifelse(TotalPop > 0, ((TotalPop - NumberWhites) / TotalPop),0),
         percentPoverty  = ifelse(TotalPop > 0, TotalPoverty / TotalPop, 0),
         tract_id        = row_number()) %>%
  st_sf() %>%
  st_transform(102747) 
tract10$tract_area <- st_area(tract10)


## ----census_statistical_area_spatial_intersection------------------------------------------------
#do the spatial join, create poverty and non whites rates by district. Create a dummy for rates >= stat_area_quantile percentile
# create intersection of tract10 and statareas
nbr_statAreas.intersect <- st_intersection(tract10, nbr_statAreas)
# get % tract in statares and mulitply by pop totals from each tract
# result is the total tract pops distributed to the statarea by % of tract in statare
nbr_statAreas.spJoin <- nbr_statAreas.intersect %>% 
  mutate(intersect_area = st_area(nbr_statAreas.intersect)) %>% 
  # get % of tract and multiply totals by percent area of tract in statarea
  group_by(tract_id) %>% 
  mutate(intersect_pcnt_of_tract = as.numeric(intersect_area) / as.numeric(tract_area),
         intersect_TotalPop = round(TotalPop * intersect_pcnt_of_tract, 1),
         intersect_NumberWhites = round(NumberWhites * intersect_pcnt_of_tract, 1),
         intersect_TotalPoverty = round(TotalPoverty * intersect_pcnt_of_tract, 1)) %>%
  ungroup() %>% 
  # sum the fraction of pop totals up to statarea
  group_by(stat_area_id) %>%
  summarise(statarea_TotalPop = sum(intersect_TotalPop),
            statarea_NumberWhites = sum(intersect_NumberWhites),
            statarea_TotalPoverty = sum(intersect_TotalPoverty)) %>% 
  # make quantites of interest
  mutate(percentNonWhite = ifelse(statarea_TotalPop > 0, 
                                  ((statarea_TotalPop - statarea_NumberWhites) / statarea_TotalPop),0),
         percentPoverty = ifelse(statarea_TotalPop > 0, 
                                 statarea_TotalPoverty / statarea_TotalPop, 0))

# classify by quantile and make dummy variable
nbr_statAreas.spJoin <- nbr_statAreas.spJoin %>% 
  mutate(poverty.percentile = ifelse(percentPoverty >=
                                       quantile(nbr_statAreas.spJoin$percentPoverty, 
                                                p = stat_area_quantile, na.rm=T),"1",0),
         nonWhite.percentile = ifelse(percentNonWhite >=
                                        quantile(nbr_statAreas.spJoin$percentNonWhite, 
                                                 p = stat_area_quantile, na.rm=T),1,0))


## ----STAT_AREA_CATEGORY_plot---------------------------------------------------------------------
STAT_AREA_CATEGORY_plot <- nbr_statAreas.spJoin %>%
  dplyr::select(poverty.percentile,nonWhite.percentile) %>%
  gather(var,val,-geometry) %>%
  ggplot() +
  geom_sf(aes(fill=factor(val))) +
  facet_wrap(~var) +
  theme_bw()


## ----aggregate_model_area_to_statistical_area----------------------------------------------------
# aggreate mean errors to statareas
stat_area_metric_logdev <- error_points %>%
  aggregate(., nbr_statAreas.spJoin, mean) %>%
  dplyr::select(logdev) %>% 
  mutate(logdev = round(logdev, 3)) %>% 
  make_cuts(., "logdev")
stat_area_metric_MAE<- error_points %>%
  aggregate(., nbr_statAreas.spJoin, mean) %>%
  dplyr::select(MAE) %>% 
  mutate(MAE = round(MAE, 3)) %>% 
  make_cuts(., "MAE")

# aggregate sum of CPS incidents to statarea
stat_area_cps <- error_points %>%
  aggregate(., nbr_statAreas.spJoin, sum) %>%
  dplyr::select(test_y)
stat_area_errors <- stat_area_metric_logdev %>% 
  st_join(., stat_area_metric_MAE, join = st_equals) %>% 
  st_join(., stat_area_cps, join = st_equals) %>% 
  st_join(., nbr_statAreas.spJoin, join = st_equals)

# group by poverty and get median of statarea aggregate errors
poverty_aggregate <- stat_area_errors %>% 
  group_by(poverty.percentile) %>% 
  summarise(med_dev = round(median(logdev),3),
            med_MAE = round(median(MAE),3),
            med_CPS = sum(test_y)) %>% 
  st_drop_geometry() %>% 
  dplyr::select(poverty.percentile, med_dev, med_MAE, med_CPS)

# group by nonwhite and get median of statarea aggregate errors
nonwhite_aggregate <- stat_area_errors %>% 
  group_by(nonWhite.percentile) %>% 
  summarise(med_dev = round(median(logdev),3),
            med_MAE = round(median(MAE),3),
            med_CPS = sum(test_y)) %>% 
  st_drop_geometry() %>% 
  dplyr::select(nonWhite.percentile, med_dev, med_MAE, med_CPS)

print(poverty_aggregate)
print(nonwhite_aggregate)


## ----STATISTICAL_AREA_MODEL_ERROR_plot-----------------------------------------------------------
logdev_stat_area_plot <- make_fishnet_dist_plot(stat_area_metric_logdev, cps_base_map, 
                                                legend = "right", 
                                                direction = -1, var_name = "Deviance",
                                                title = "Out-of-Fold-Error by Statistical Area")

MAE_stat_area_plot <- make_fishnet_dist_plot(stat_area_metric_MAE, cps_base_map, 
                                             legend = "right", 
                                             direction = 1, var_name = "MAE",
                                             title = "Out-of-Fold-Error by Statistical Area")
STATISTICAL_AREA_MODEL_ERROR_plot <- cowplot::plot_grid(logdev_stat_area_plot , MAE_stat_area_plot, 
                                                        align = "h", labels = "Out-of-Fold Error by Statistical Area")


## ----prediction_to_bin_class---------------------------------------------------------------------
error_geoplot$pred_bin_class <- bin_class(error_geoplot, "pred")

p.summ <- error_geoplot %>%
  group_by(pred_bin_class) %>%
  dplyr::summarize(obs.total = sum(test_y),
                   obs.cnt = n()) %>% 
  rename(sens_group = pred_bin_class) %>%
  filter(!is.na(sens_group)) %>%
  identity()


## ----compute_KDE---------------------------------------------------------------------------------
cps_ppp <- as.ppp(st_coordinates(cps_dissolve), W = st_bbox(net_Richmond))
cps_KDE <- spatstat::density.ppp(cps_ppp)

cps_KDE_tbl <- as.data.frame(cps_KDE) %>%
  st_as_sf(coords = c("x", "y"), crs = 102747) %>%
  aggregate(., net_Richmond, mean) %>%
  mutate(net_id = net_Richmond$net_id)


## ----KDE_to_bin_class----------------------------------------------------------------------------
if(all.equal(error_geoplot$net_id, cps_KDE_tbl$net_id)){
  cat("Good to go!")
} else {
  cat("Join will be an error, Net_id index does not match")
}

error_geoplot$kde_bin_class  <- bin_class(cps_KDE_tbl, "value")

kde.summ <- error_geoplot %>%
  group_by(kde_bin_class) %>%
  dplyr::summarize(kde.total = sum(test_y),
                   kde.cnt = n()) %>% 
  rename(sens_group = kde_bin_class) %>%
  filter(!is.na(sens_group)) %>%
  identity()


## ----REALTIVE_SENSITIVITY_plot-------------------------------------------------------------------
REALTIVE_SENSITIVITY_KDE <- ggmap(cps_base_map) +
  geom_sf(data = ll(kde.summ), aes(fill = factor(sens_group)), 
          color = NA, alpha = 0.85, inherit.aes = FALSE) +
  geom_sf(data = ll(cps_dissolve), inherit.aes = FALSE, size = 1) +
  scale_fill_viridis_d(na.value = NA, option = "D", direction = 1,
                       name = "Risk\nCategory") +
  labs(title = "Risk categories from KDE",
       caption = "Figure 6.4") +
  mapTheme()

REALTIVE_SENSITIVITY_PREDICTIONS <- ggmap(cps_base_map) +
  geom_sf(data = ll(p.summ), aes(fill = factor(sens_group)), 
          color = NA, alpha = 0.85, inherit.aes = FALSE) +
  geom_sf(data = ll(cps_dissolve), inherit.aes = FALSE, size = 1) +
  scale_fill_viridis_d(na.value = NA, option = "D", direction = 1,
                       name = "Risk\nCategory") +
  labs(title = "Risk categories from meta-model",
       caption = "Figure 1.5") +
  mapTheme()


## ----REALTIVE_RISK_BARPLOT_COMPARE_plot----------------------------------------------------------
countComparisons <- merge(st_drop_geometry(p.summ), st_drop_geometry(kde.summ)) %>%
  mutate_if(is.double, round, 3) %>% 
  mutate(Category = rev(c("90% - 100%", "70% - 89%", "50% - 69%", 
                          "30% - 49%", "1% - 29%"))) %>%
  dplyr::mutate(kernelPct = round(kde.total / sum(kde.total),4),
                fittedPct = round(obs.total / sum(obs.total), 4))

countComparisonsLong <- countComparisons %>% 
  gather(Variable, Value, kernelPct:fittedPct)

REALTIVE_RISK_BARPLOT_COMPARE_plot <- ggplot(data=countComparisonsLong, aes(Category,Value)) +
  geom_bar(aes(fill = Variable), position = "dodge", stat="identity", color = NA) +
  scale_fill_viridis_d(name = " ",
                       labels=c("Meta-model", "Kernel Density")) +
  labs(x= "Predicted Risk Levels",
       y="Percent of Test Set Cases",
       title= "Goodness of fit: Spatial risk model vs. Kernel Density",
       caption = "Figure 1.6") +
  plotTheme() +
  theme(axis.line = element_blank())


## ----alignPhase----------------------------------------------------------------------------------
##In this section change "sensitivty class" to "risk category"

#the final prediction map
predMap <- 
  error_geoplot #%>%
#st_transform(102747)
#removals
removals <- 
  read.csv("C:/projects/PAP_Virginia/data/z_alignPhase/removals.csv") %>%
  st_as_sf(coords = c("X", "Y"), crs = 102747, agr = "constant")
#service visits
visits <- 
  read.csv("C:/projects/PAP_Virginia/data/z_alignPhase/healthFamilyServices.csv") %>%
  filter(!is.na(X)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") %>%
  st_transform(102747) %>%
  #some are outside of the the study area, so select by those that itersect
  .[st_union(predMap),]
#SCAN data - 'stop child abuse now' network centers. geocode. There are only 27 in the study area
scan <- read.csv("C:/projects/PAP_Virginia/data/z_alignPhase/scanPreventionResources.csv") %>%
  mutate(Street.Address = as.character(Street.Address))  %>%
  mutate(Street.Address = paste(Street.Address, "Richmond, Va."))
scan <-
  scan %>%
  bind_cols(scan,geocode(scan$Street.Address, source="dsk")) %>%
  filter(!is.na(lon)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, agr = "constant") %>%
  st_transform(102747) %>%
  .[st_union(predMap),]

#7.1 - population per risk category
popID <- ALL_FEATURES %>% 
  dplyr::select(net_id, net_pop) %>% 
  as.data.frame()

binID <- error_geoplot %>% 
  dplyr::select(net_id, pred_bin_class) %>% 
  as.data.frame()

popBin <- left_join(popID, binID, by = "net_id") %>% 
  group_by(pred_bin_class) %>% 
  summarise(sumPop = sum(net_pop)) %>%
  mutate(pctPop = sumPop/sum(sumPop),
         Category = rev(c("90% - 100%", "70% - 89%", "50% - 69%", 
                          "30% - 49%", "1% - 29%")))

popPerRiskPlot <- ggplot(data=popBin, aes(Category, sumPop)) +
  geom_bar(position = "dodge", stat="identity") +
  labs(x= "Predicted Risk Levels",
       y= "Number of People",
       title = "Population per risk category",
       caption = "Figure 2.1") +
  plotTheme()

#7.2 - poverty rate correlation with predicted maltreatment count
#take the net id and geometry for the fishnet 
geom_fishnet <- error_geoplot %>% 
  dplyr::select(net_id, geometry)

#get population and poverty info at tract level
poverty_tract <- tract10 %>% 
  dplyr::select(TotalPoverty, tract_id, geometry) %>% 
  mutate(tract_acre = as.numeric(st_area(.)*2.29568e-5),
         pov_acre_rate = TotalPoverty/tract_acre)

population_tract <- tract10 %>% 
  dplyr::select(TotalPop, tract_id, geometry) %>% 
  mutate(tract_acre = as.numeric(st_area(.)*2.29568e-5),
         pop_acre_rate = TotalPop/tract_acre)

#intersect tract poverty and fishnet
pov_tracts_intersect <- st_intersection(poverty_tract, geom_fishnet)

pov_tracts_intersect <- pov_tracts_intersect %>%
  mutate(int_area_acres = as.numeric(st_area(pov_tracts_intersect)*2.29568e-5)) %>%
  group_by(net_id) %>%
  mutate(cnt = n(),
         pcnt_of_block = int_area_acres/tract_acre,
         int_pov =TotalPoverty * pcnt_of_block) %>%
  arrange(net_id)

fishnet_poverty <- pov_tracts_intersect %>% # xcc
  group_by(net_id) %>%
  summarise(net_pov = sum(int_pov)) %>% 
  as.data.frame()

#intersect tract population and fishnet
pop_tracts_intersect <- st_intersection(population_tract, geom_fishnet)

pop_tracts_intersect <- pop_tracts_intersect %>%
  mutate(int_area_acres = as.numeric(st_area(pop_tracts_intersect)*2.29568e-5)) %>%
  group_by(net_id) %>%
  mutate(cnt = n(),
         pcnt_of_block = int_area_acres/tract_acre,
         int_pop =TotalPop * pcnt_of_block) %>%
  arrange(net_id)

fishnet_population <- pop_tracts_intersect %>% # xcc
  group_by(net_id) %>%
  summarise(net_pop = sum(int_pop)) %>% 
  as.data.frame()

pov_pop_fishnet <- left_join(fishnet_poverty, fishnet_population, by = "net_id") %>% 
  mutate(povRate = net_pov/net_pop) %>% 
  dplyr::select(-geometry.x) %>% 
  rename(geometry = geometry.y) %>% 
  st_sf() %>% 
  st_transform(102747)

#map it:
povertyRateMap <- ggmap(cps_base_map) +
  geom_sf(data=ll(pov_pop_fishnet), aes(fill=factor(ntile(povRate, 5))), inherit.aes = FALSE, alpha = 0.8, color = NA) +
  scale_fill_viridis_d(labels = as.character(round(quantile(pov_pop_fishnet$povRate,
                                                            c(.1,.2,.4,.6,.8),na.rm=T), 4)),
                       name="Poverty\nRate") +
  labs(title = "Weighted poverty rate",
       caption = "Figure 2.2") +
  mapTheme()

pov_pop_fishnet_pred <- left_join(pov_pop_fishnet, error_geoplot %>% 
                                    dplyr::select(net_id, pred) %>% 
                                    #mutate(pred = round(pred)) %>% 
                                    as.data.frame(), by = "net_id") %>% 
  filter(pred > 0)

povRatePredPlot <- ggplot(pov_pop_fishnet_pred, aes(x=povRate, y=pred)) + 
  geom_point() +
  labs(x = "Poverty Rate",
       y = "Predicted\nMaltreatment Counts",
       title = "Relationship between predicted maltreatment counts\nand poverty rate",
       caption = "Figure 2.3") +
  plotTheme()

#7.3 - map of risk categoeries with removals overlayed
removalsMap <-
  ggmap(cps_base_map) +
  geom_sf(data=ll(predMap), aes(fill=factor(pred_bin_class)), inherit.aes = FALSE, alpha = 0.8, color = NA) +
  geom_sf(data=ll(removals),aes(colour="Removals"),colour="red", inherit.aes = FALSE, size = 1) +
  scale_fill_viridis_d(name = "Risk\nCategory") +
  labs(title="Predicted risk levels and removals",
       subtitle="Removals in red",
       caption = "Figure 2.4") +
  mapTheme()

#7.4 - count of removals by risk categoery
removalsPlot <-
  removals %>%
  mutate(counter=1) %>%
  aggregate(predMap,FUN=length) %>%
  dplyr::select(counter) %>%
  mutate(counter = ifelse(is.na(counter),0,counter)) %>%
  bind_cols(predMap) %>%
  mutate(Category = case_when(pred_bin_class == 1 ~ "1% - 29%",
                              pred_bin_class == 2 ~ "30% - 49%%",
                              pred_bin_class == 3 ~ "50% - 69%",
                              pred_bin_class == 4 ~ "70% - 89%",
                              pred_bin_class == 5 ~ "90% - 100%")) %>%
  group_by(Category) %>%
  dplyr::summarize(percentCount = sum(counter)/nrow(removals)) %>%
  ggplot(aes(Category,percentCount)) +
  geom_bar(position = "dodge", stat="identity") +
  labs(x= "Predicted Risk Levels",
       y="Percent of Removals",
       title= "Percent of removals by risk category",
       caption = "Figure 2.5") +
  plotTheme()

#7.4 - map of protective land uses by type
#add the data
protectiveAlign <- 
  read.csv("C:/projects/PAP_Virginia/data/z_alignPhase/protective_align.csv") %>%
  st_as_sf(coords = c("x", "y"), crs = 102747, agr = "constant") %>%
  .[st_union(predMap),]

#map 
protectiveUsesByType <- ggmap(cps_base_map) +
  geom_sf(data=ll(st_union(predMap)), inherit.aes = FALSE, alpha = 0.6, color = NA) +
  geom_sf(data=ll(protectiveAlign), aes(colour=use), inherit.aes = FALSE) +
  scale_colour_viridis_d(name="Protective\nuses",
                         label = c("Church", "Community Center", "Fire Station", "Homeless Shelter",
                                   "Library", "Oasis", "Police Station", "School")) +
  labs(title="Protective land uses",
       caption = "Figure 2.7") +
  mapTheme()

#7.5 mea predictive count by quarter
protectiveAlign.buffers <-
  st_centroid(predMap) %>%
  aggregate(st_buffer(protectiveAlign,1320),FUN=mean) 

#map average predicted event by quarter mile buffer
quarterMileBuffer_protective <-
  ggmap(cps_base_map) +
  geom_sf(data=ll(st_union(predMap)), inherit.aes = FALSE, alpha = 0.8, color = NA) +
  geom_sf(data=ll(protectiveAlign.buffers), aes(fill=pred), inherit.aes = FALSE) +
  scale_fill_viridis_c(name = "Mean\npredicted\ncount") +
  labs(title = "Mean predicted count by quarter mile buffer",
       subtitle = "Protective uses",
       caption = "Figure 2.8") +
  guides(fill = guide_colourbar(reverse = TRUE)) +
  mapTheme()

#table of top protective places
protectiveTable <- 
  protectiveAlign.buffers %>%
  bind_cols(protectiveAlign) %>%
  group_by(use) %>%
  top_n(n=3,wt=pred) %>%
  as.data.frame() %>%
  mutate(Mean_Predicted_Count = round(pred)) %>%
  dplyr::select(use,name,address,Mean_Predicted_Count) %>%
  arrange(use,-Mean_Predicted_Count) %>%
  kable() %>% 
  kable_styling()

#7.6 - visits as a function of risk levels
homeVisits <-
  ggmap(cps_base_map) +
  geom_sf(data=ll(predMap), aes(fill=factor(pred_bin_class)), inherit.aes = FALSE, alpha = 0.8, color = NA) +
  geom_sf(data=ll(visits), inherit.aes = FALSE, colour="red") +
  scale_fill_viridis_d(name = "Risk\nCategory") +
  labs(title="Predicted risk levels and home visits",
       subtitle="Home visits in red",
       caption = "Figure 2.12") +
  mapTheme()

#7.7 visits by risk category
vistsCategory <-
  visits %>%
  mutate(counter=1) %>%
  aggregate(predMap,FUN=length) %>%
  dplyr::select(counter) %>%
  mutate(counter = ifelse(is.na(counter),0,counter)) %>%
  bind_cols(predMap) %>%
  mutate(Category = case_when(pred_bin_class == 1 ~ "1% - 29%",
                              pred_bin_class == 2 ~ "30% - 49%%",
                              pred_bin_class == 3 ~ "50% - 69%",
                              pred_bin_class == 4 ~ "70% - 89%",
                              pred_bin_class == 5 ~ "90% - 100%")) %>%
  group_by(Category) %>%
  dplyr::summarize(percentCount = sum(counter)/nrow(visits)) %>%
  ggplot(aes(Category,percentCount)) +
  geom_bar(position = "dodge", stat="identity") +
  labs(x= "Predicted Risk Levels",
       y="Percent of Visits",
       title= "Percent of visits by risk category") +
  plotTheme()

#7.8 scan centers
scanCenters <-
  ggmap(cps_base_map) +
  geom_sf(data=ll(predMap), aes(fill=factor(pred_bin_class)), inherit.aes = FALSE, alpha = 0.8, color = NA) +
  geom_sf(data=ll(scan),inherit.aes = FALSE, colour="red") +
  scale_fill_viridis_d(name = "Risk\nCategory") +
  labs(title="Predicted risk levels and SCAN Centers",
       subtitle="SCAN Centers in red",
       caption = "Figure 2.13") +
  mapTheme()

#7.9 scan buffers and risk
scan.buffers <-
  st_centroid(predMap) %>%
  aggregate(st_buffer(scan,1320),FUN=mean) 

#map average predicted event by quarter mile buffer
scanMap <-
  ggmap(cps_base_map) +
  geom_sf(data=ll(st_union(predMap)), inherit.aes = FALSE, alpha = 0.8, color = NA) +
  geom_sf(data=ll(scan.buffers), aes(fill=pred), inherit.aes = FALSE) +
  scale_fill_viridis_c(name = "Mean\npredicted\ncount") +
  labs(title = "Mean predicted count by quarter mile buffer",
       subtitle = "SCAN Centers",
       caption = "Figure 2.14") +
  guides(fill = guide_colourbar(reverse = TRUE)) +
  mapTheme()

#find top 3 Scan centers
scan.buffers %>%
  bind_cols(scan) %>%
  top_n(n=3,wt=pred) %>%
  as.data.frame() %>%
  mutate(Mean_Predicted_Count = round(pred)) %>%
  dplyr::select(SCAN.s.Trauma.Informed.Community.Network..TICN.,Street.Address,Mean_Predicted_Count) %>%
  arrange(-Mean_Predicted_Count) %>%
  kable()

#Gap
#get neighborhoods
nbr_statAreas <- read_sf("https://data.richmondgov.com/resource/8kyq-v9j2.geojson") %>%
  st_sf() %>%
  st_transform(102747) %>%
  #calculate area
  mutate(area.sq_miles = as.numeric(st_area(.) * 3.58701e-8))

#create a dummy variable field for where risk category == 5
predMap$is5th <- ifelse(predMap$pred_bin_class == 5,1,0)
#count 5th risk quintile per neighborhood
nbr_statAreas.demand <-
  st_centroid(predMap) %>%
  dplyr::select(is5th,geometry) %>%
  aggregate(nbr_statAreas,FUN=sum) %>%
  bind_cols(
    st_centroid(predMap) %>%
      mutate(counter=1) %>%
      dplyr::select(counter,geometry) %>%  
      aggregate(nbr_statAreas,FUN=sum)) %>%
  mutate(relativeRisk = ntile((is5th / counter), 100))

#calculate number of protective centers within 
nbr_statAreas.supply <-
  protectiveAlign %>%
  mutate(counter= 1) %>%
  dplyr::select(counter,geometry) %>%
  aggregate(nbr_statAreas,FUN=sum) %>%
  bind_cols(nbr_statAreas) %>%
  dplyr::select(counter,geometry,area.sq_miles) %>%
  mutate(relativeProtective = ntile((counter/area.sq_miles),100))

#put demand and supply together and look at difference
nbr_statAreas.gap <-
  bind_cols(nbr_statAreas.demand, nbr_statAreas.supply) %>%
  mutate(gap = relativeProtective - relativeRisk)

#map average predicted event by neighborhoods
bl <- colorRampPalette(c("green2","green3","green4"))(200)                      
re <- colorRampPalette(c("darkred", "red3","red1"))(200)

gapMap <-
  ggmap(cps_base_map) +
  geom_sf(data=ll(nbr_statAreas.gap), aes(fill=gap), inherit.aes = FALSE,alpha = 0.8) +
  scale_fill_gradientn(colours=c(re,"green2", bl),
                       limits = c(-100,100),
                       breaks = c(-100,100),
                       labels = c("More risk\nthan protection","More protection\nthan risk")) +
  labs(title = "Gap analysis:\nComparing relative risk to protective resources",
       subtitle = "By Neighborhood Statistical Area",
       caption = "Figure 1.8") +
  mapTheme() +
  theme(legend.position = 'bottom',
        legend.title=element_blank(),
        legend.text=element_text(size=10, face = "bold"),
        legend.justification = "center") +
  guides(fill= guide_colorbar(barwidth=15,barheight=2))


## ----childRiskAreas------------------------------------------------------------------------------
child_tract10 <- get_acs(geography = "tract", variables = c("B01001_003E", "B01001_004E", "B01001_005E", "B01001_006E",
                                                            "B01001_027E", "B01001_028E", "B01001_029E", "B01001_030E"), 
                         year = 2010, state=51, county=760, geometry=T) %>%
  dplyr::select(variable,estimate) %>%
  as.data.frame() %>%
  spread(variable,estimate) %>%
  rename(MU5=B01001_003,
         M5_9=B01001_004,
         M10_14=B01001_005,
         M15_17 = B01001_006,
         FU5=B01001_027,
         F5_9=B01001_028,
         F10_14=B01001_029,
         F15_17 = B01001_030) %>%
  mutate(sum_U18 = (MU5 + M5_9 + M10_14 + M15_17 + FU5 + F5_9 + F10_14 + F15_17),   #39600
         tract_id = row_number()) %>%
  st_sf() %>%
  st_transform(102747) 

child_tract10 <- child_tract10 %>% 
  mutate(area = st_area(child_tract10))

children_tract <- child_tract10 %>% 
  dplyr::select(sum_U18, tract_id, geometry) %>% 
  mutate(tract_acre = as.numeric(st_area(.)*2.29568e-5),
         pop_acre_rate = sum_U18/tract_acre)


#intersect with fishnet
child_tracts_intersect <- st_intersection(children_tract, geom_fishnet)

child_tracts_intersect <- child_tracts_intersect %>% 
  mutate(int_area_acres = as.numeric(st_area(child_tracts_intersect)*2.29568e-5)) %>%
  group_by(net_id) %>%
  mutate(cnt = n(),
         pcnt_of_block = int_area_acres/tract_acre,
         int_child = sum_U18 * pcnt_of_block) %>%
  arrange(net_id)

fishnet_children <- child_tracts_intersect  %>%
  group_by(net_id) %>%
  summarise(net_child = sum(int_child)) %>%    #39599.9
  as.data.frame()


childBin <- left_join(binID, fishnet_children, by = "net_id") %>% 
  filter(net_child != "NA") %>% 
  group_by(pred_bin_class) %>% 
  summarise(count = n(),
            sumChild = sum(net_child)) %>%     #39599.9
  mutate(pctChild = sumChild/sum(sumChild),
         Category = rev(c("90% - 100%", "70% - 89%", "50% - 69%", 
                          "30% - 49%", "1% - 29%")))


## ----CPSCOUNT_by_nbr-----------------------------------------------------------------------------
nbr_CPSCount <- st_join(fishnet_pop_cps, nbr, largest = TRUE) %>% 
  dplyr::select(name, geometry, net_CPS_Accepted, net_id) %>% 
  group_by(name) %>% 
  summarise(count = sum(net_CPS_Accepted)) %>% 
  filter(name != "NA")

CPSCOUNT_by_nbr_plot <- ggmap(cps_base_map) +
  geom_sf(data = ll(nbr_CPSCount), aes(fill = factor(ntile(count, 5))), color = NA, alpha = 0.8, inherit.aes = FALSE) +
  scale_fill_viridis_d(labels = as.character(quantile(nbr_CPSCount$count,
                                                      c(.1,.2,.4,.6,.8),na.rm=T)),
                       direction = 1,
                       name="Count\nQuantile\nBreaks") +
  labs(title = "Number of Child Protective Service events\nby neighborhood",
       caption = "Figure 1.1") +
  mapTheme()


## ----GLOBAL_MORANS_TEST_MAE----------------------------------------------------------------------
#matrix of coordinates
nbrCoords <- neighborhood_metric_MAE %>% 
  st_centroid() %>% 
  mutate(X = st_coordinates(.)[,1],
         Y = st_coordinates(.)[,2]) %>% 
  as.data.frame() %>% 
  dplyr::select(X, Y) %>% 
  as.matrix()


#make weights for test - neighborhood level
neighborNbrs <- knn2nb(knearneigh(nbrCoords, 5))
spatialWeights <- nb2listw(neighborNbrs, style="W")

moranTest_MAE <- moran.mc(neighborhood_metric_MAE$MAE, spatialWeights, nsim = 999)


## ----childFatality-------------------------------------------------------------------------------
deaths <- 
  read.csv("C:/projects/PAP_Virginia/data/z_alignPhase/deaths.csv") %>%
  st_as_sf(coords = c("X", "Y"), crs = 102747, agr = "constant") %>% 
  filter(ChildFatality == "Y")

fatalitiesMap <- ggmap(cps_base_map) +
  geom_sf(data=ll(predMap), aes(fill=factor(pred_bin_class)), inherit.aes = FALSE, alpha = 0.8, color = NA) +
  geom_sf(data=ll(deaths),aes(colour="ChildFatality"),colour="red", inherit.aes = FALSE, size = 1.5) +
  scale_fill_viridis_d(name = "Risk\nCategory") +
  labs(title="Predicted risk levels and child fatalities",
       subtitle="Fatalities in red",
       caption = "Figure 2.6") +
  mapTheme()

fatalitiesPlot <- deaths %>%
  mutate(counter=1) %>%
  aggregate(predMap,FUN=length) %>%
  dplyr::select(counter) %>%
  mutate(counter = ifelse(is.na(counter),0,counter)) %>%
  bind_cols(predMap) %>%
  mutate(Category = case_when(pred_bin_class == 1 ~ "1% - 29%",
                              pred_bin_class == 2 ~ "30% - 49%%",
                              pred_bin_class == 3 ~ "50% - 69%",
                              pred_bin_class == 4 ~ "70% - 89%",
                              pred_bin_class == 5 ~ "90% - 100%")) %>%
  group_by(Category) %>%
  dplyr::summarize(percentCount = sum(counter)/nrow(deaths)) %>%
  ggplot(aes(Category,percentCount)) +
  geom_bar(position = "dodge", stat="identity") +
  labs(x= "Predicted Risk Levels",
       y="Percent of Child Fatalities",
       title= "Percent of child fatalities by risk category") +
  plotTheme()


## ----churchesBuffers-----------------------------------------------------------------------------
churches <- read.csv("C:/projects/PAP_Virginia/data/z_alignPhase/protective_align.csv") %>%
  st_as_sf(coords = c("x", "y"), crs = 102747, agr = "constant") %>% 
  filter(use == "Church") %>% 
  .[st_union(predMap),]

church.buffers <-
  st_centroid(predMap) %>%
  aggregate(st_buffer(churches,1320),FUN=mean)

#map average predicted event by quarter mile buffer
churchMap <- ggmap(cps_base_map) +
  geom_sf(data=ll(st_union(predMap)), inherit.aes = FALSE, alpha = 0.8, color = NA) +
  geom_sf(data=ll(church.buffers), aes(fill=pred), inherit.aes = FALSE) +
  scale_fill_viridis_c(name = "Mean\npredicted\ncount",
                       labels = round) +
  labs(title = "Mean predicted count by quarter mile buffer",
       subtitle = "Churches and other faith organizations",
       caption = "Figure 2.9") +
  guides(fill = guide_colourbar(reverse = TRUE)) +
  mapTheme()

#find top 3 Scan centers
top_church <- church.buffers %>%
  bind_cols(churches) %>%
  top_n(n=10,wt=pred) %>%
  as.data.frame() %>%
  mutate(Mean_Predicted_Count = round(pred)) %>%
  dplyr::select(name,address,Mean_Predicted_Count) %>%
  arrange(-Mean_Predicted_Count) %>%
  kable() %>% 
  kable_styling()


## ----childcareBuffers----------------------------------------------------------------------------
childcare <- st_read("C:/projects/PAP_Virginia/data/z_alignPhase/Richmond_child_providers.shp") %>%
  dplyr::select(Loc_name, Match_addr, Type_Of_Ca, License_Ty) %>% 
  st_transform(102747) %>% 
  .[st_union(predMap),]

resourcehomes <- read.csv("C:/projects/PAP_Virginia/data/z_alignPhase/resourceHomes.csv") %>%
  filter(X != "#NUM!",
         Y != "#NUM!") %>% 
  mutate(X = as.double(as.character(X)),
         Y = as.double(as.character(Y)),
         Type = ifelse(Type%in% c("Emergency Foster Care Home", "FFC/Kinship/Relative", 
                                  "Kinship/Relative Non-Paid", "LCPA Homes", "LDSS Home"), 
                       "Individual", as.character(Type)),
         ResourceName = ifelse(Type == "Individual", "Individual", as.character(ResourceName)),
         Address = ifelse(Type == "Individual", "Redacted (HIPPA)", as.character(Address))) %>%
  st_as_sf(coords = c("X", "Y"), crs = 102747, agr = "constant") %>% 
  .[st_union(predMap),]

childcare.buffers <-
  st_centroid(predMap) %>%
  aggregate(st_buffer(childcare,1320),FUN=mean)

resourcehomes.buffers <-
  st_centroid(predMap) %>%
  aggregate(st_buffer(resourcehomes,1320),FUN=mean) 

#map average predicted event by quarter mile buffer
childcareMap <- ggmap(cps_base_map) +
  geom_sf(data=ll(st_union(predMap)), inherit.aes = FALSE, alpha = 0.8, color = NA) +
  geom_sf(data=ll(childcare.buffers), aes(fill=pred), inherit.aes = FALSE) +
  scale_fill_viridis_c(name = "Mean\npredicted\ncount") +
  labs(title = "Mean predicted count by quarter mile buffer",
       subtitle = "Licensed child care providers",
       caption = "Figure 2.10") +
  guides(fill = guide_colourbar(reverse = TRUE)) +
  mapTheme()

resourcehomeMap <- ggmap(cps_base_map) +
  geom_sf(data=ll(st_union(predMap)), inherit.aes = FALSE, alpha = 0.8, color = NA) +
  geom_sf(data=ll(resourcehomes.buffers), aes(fill=pred), inherit.aes = FALSE) +
  scale_fill_viridis_c(name = "Mean\npredicted\ncount",
                       labels = round) +
  labs(title = "Mean predicted count by quarter mile buffer",
       subtitle = "Locations certified to accept foster care children",
       caption = "Figure 2.11") +
  guides(fill = guide_colourbar(reverse = TRUE)) +
  mapTheme()

#find top 3 
top_childcare <- childcare.buffers %>% 
  bind_cols(childcare) %>%
  group_by(License_Ty) %>%
  top_n(n=3,wt=pred) %>%
  as.data.frame() %>%
  mutate(Mean_Predicted_Count = round(pred)) %>%
  dplyr::select(License_Ty, Match_addr,Mean_Predicted_Count) %>%
  rename(license_type = License_Ty,
         address = Match_addr) %>% 
  arrange(license_type,-Mean_Predicted_Count) %>%
  kable() %>% 
  kable_styling()

top_resourcehomes <- resourcehomes.buffers %>% 
  bind_cols(resourcehomes) %>%
  group_by(Type) %>%
  top_n(n=3,wt=pred) %>%
  as.data.frame() %>%
  mutate(Mean_Predicted_Count = round(pred)) %>%
  dplyr::select(ResourceName,Type,Address,Mean_Predicted_Count) %>%
  rename(type = Type,
         name = ResourceName) %>% 
  arrange(type,-Mean_Predicted_Count) %>%
  kable() %>% 
  kable_styling()


## ----cache = TRUE,echo=FALSE, warning=FALSE, include = TRUE, fig.width=11.5, fig.height=5--------
cowplot::plot_grid(CPSCOUNT_by_nbr_plot, CPS_RATE_BY_FISHNET_PLOT, align = "hv", axis = "lrbt")


## ----message=FALSE, warning=FALSE, cache=TRUE, echo=FALSE, include=TRUE, fig.width = 9-----------
#pull out rec centers
rec_centers <- var_list$CommunityCenters
rec_centers_count_map <-
  ggmap(cps_base_map) +
  geom_sf(data = ll(rec_centers), inherit.aes = FALSE, color = "red", size = 2, alphae = 0.9) +
  geom_sf(data = ll(nbr_diss), inherit.aes = FALSE, color = "black", fill = NA, size = 1) +
  labs(title = "Community center locations") +
  #scale_fill_viridis_d(na.value = NA, option = "D", direction = 1, name = "Count") +
  mapTheme()


#grab and map count
rec_centers_agg <- ALL_FEATURES %>% 
  dplyr::select(agg_CommunityCenters, net_id) %>% 
  left_join(., net_Richmond, by = "net_id") %>% 
  mutate(value = ifelse(is.na(agg_CommunityCenters), 0, agg_CommunityCenters)) %>% 
  st_sf()
rec_centers_agg <- make_cuts(rec_centers_agg, "value", cuts = "breaks", n_breaks = 5)
rec_centers_agg_map <-
  ggmap(cps_base_map) +
  geom_sf(data = ll(rec_centers_agg), aes(fill = cut_val), inherit.aes = FALSE, color = NA, alpha = 0.8) +
  # geom_sf(data = ll(rec_centers), inherit.aes = FALSE, color = "red", size = 4) +  
  #geom_sf(data = ll(nbr_diss), inherit.aes = FALSE, color = "black", fill = NA, size = 2) +
  labs(title = "Community center count\nby fishnet") +
  scale_fill_viridis_d(na.value = NA, option = "D", direction = 1, name = "Count") +
  mapTheme()

#grab and counteuclidean distance
rec_centers_ed <- ALL_FEATURES %>% 
  dplyr::select(ed_CommunityCenters, net_id) %>% 
  left_join(., net_Richmond, by = "net_id") %>% 
  mutate(value = ifelse(is.na(ed_CommunityCenters), 0, ed_CommunityCenters)) %>% 
  st_sf()
rec_centers_ed <- make_cuts(rec_centers_ed, "value",  "value", cuts = "breaks", n_breaks = 5)
rec_centers_ed_map <-
  ggmap(cps_base_map) +
  geom_sf(data = ll(rec_centers_ed), aes(fill = cut_val), inherit.aes = FALSE, color = NA, alpha = 0.8) +
  geom_sf(data = ll(rec_centers), inherit.aes = FALSE, color = "red", size = 1, alpha = 0.9) +
  #geom_sf(data = ll(nbr_diss), inherit.aes = FALSE, color = "black", fill = NA, size = 2) +
  labs(title = "Community center\neuclidean distance",
       caption = "Figure 1.2") +
  scale_fill_viridis_d(na.value = NA, option = "D", direction = 1, name = "NN Distance") +
  mapTheme()

#grab and map avg nn 
rec_centers_NN <- ALL_FEATURES %>% 
  dplyr::select(NN_CommunityCenters, net_id) %>% 
  left_join(., net_Richmond, by = "net_id") %>% 
  mutate(value = ifelse(is.na(NN_CommunityCenters), 0, NN_CommunityCenters)) %>% 
  st_sf()
rec_centers_NN <- make_cuts(rec_centers_NN, "value",  "value", cuts = "breaks", n_breaks = 5)
rec_centers_NN_map <-
  ggmap(cps_base_map) +
  geom_sf(data = ll(rec_centers_NN), aes(fill = cut_val), inherit.aes = FALSE, color = NA, alpha = 0.8) +
  geom_sf(data = ll(rec_centers), inherit.aes = FALSE, color = "red", size = 1, alpha = 0.9) +
  #geom_sf(data = ll(nbr_diss), inherit.aes = FALSE, color = "black", fill = NA, size = 2) +
  labs(title = "Community center average\nnearest neighbor distance") +
  scale_fill_viridis_d(na.value = NA, option = "D", direction = 1, name = "NN Distance") +
  mapTheme()

cowplot::plot_grid(rec_centers_count_map,rec_centers_agg_map,rec_centers_ed_map, rec_centers_NN_map, align="hv", axis = "rlbt", ncol=2)


## ----cache = TRUE,echo=FALSE,warning=FALSE,include = TRUE,fig.align="center"---------------------
ggmap(cps_base_map) +
  geom_sf(data = ll(p_cut), aes(fill = pval_cut),
          color = NA, inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_viridis_d(na.value=NA, name = "p-value", option = "D") +
  labs(title = "Stastically significant maltreatment clusters",
       caption = "Figure 1.3") +
  mapTheme()


## ----message=FALSE, warning=FALSE, cache=TRUE, echo=FALSE, include=TRUE, fig.height=18-----------
CRIME_CPS_VS_NN_plot


## ----cache = TRUE,echo=FALSE,warning=FALSE, include = TRUE, fig.align = "center"-----------------
REALTIVE_SENSITIVITY_PREDICTIONS


## ----cache = TRUE,echo=FALSE,warning=FALSE, include = TRUE, fig.width=8, fig.align="center"------
REALTIVE_RISK_BARPLOT_COMPARE_plot


## ----cache = TRUE,echo=FALSE,warning=FALSE, include = TRUE, fig.align="center"-------------------
META_MODEL_FIT_plot


## ----cache = TRUE, echo=FALSE,warning=FALSE, include = TRUE, fig.height = 7, fig.width = 6.5, fig.align = "center"----
gapMap


## ----cache = TRUE, echo=FALSE,warning=FALSE, include = TRUE, fig.align="center"------------------
ggmap(cps_base_map) +
  geom_sf(data=ll(st_union(predMap)), inherit.aes = FALSE, alpha = 0.8, color = NA) +
  geom_sf(data=ll(scan.buffers), aes(fill=pred), inherit.aes = FALSE) +
  scale_fill_viridis_c(name = "Mean predicted\ncount") +
  labs(title = "Mean predicted count by quarter mile buffer",
       subtitle = "SCAN Centers",
       caption = "Figure 1.9") +
  guides(fill = guide_colourbar(reverse = TRUE)) +
  mapTheme()


## ----include = TRUE, fig.width=8, fig.align="center"---------------------------------------------
popPerRiskPlot


## ----include = TRUE, fig.align="center"----------------------------------------------------------
    povertyRateMap


## ----include = TRUE, fig.width = 8, fig.align="center"-------------------------------------------
povRatePredPlot


## ----include = TRUE, fig.align="center"----------------------------------------------------------
removalsMap


## ----include = TRUE, fig.width = 8, fig.align="center"-------------------------------------------
removalsPlot


## ----include = TRUE, fig.width=14, fig.height=7--------------------------------------------------
cowplot::plot_grid(fatalitiesMap, fatalitiesPlot, align = "hv", axis = "lrbt")


## ----include = TRUE, fig.width = 10, fig.align="center"------------------------------------------
protectiveUsesByType


## ----include = TRUE, fig.width = 10, fig.align="center"------------------------------------------
quarterMileBuffer_protective


## ----include = TRUE------------------------------------------------------------------------------
protectiveTable


## ----include = TRUE, fig.width=10.5--------------------------------------------------------------
churchMap


## ----include = TRUE------------------------------------------------------------------------------
top_church


## ----include = TRUE, fig.width = 10, fig.align="center"------------------------------------------
childcareMap


## ----include = TRUE------------------------------------------------------------------------------
top_childcare


## ----include = TRUE, fig.width = 10, fig.align="center"------------------------------------------
resourcehomeMap


## ----include = TRUE------------------------------------------------------------------------------
top_resourcehomes


## ----include = TRUE, fig.width=14, fig.height=7--------------------------------------------------
cowplot::plot_grid(homeVisits, vistsCategory, align = "hv", axis = "lrbt")


## ----include = TRUE, fig.align="center"----------------------------------------------------------
scanCenters


## ----include = TRUE,fig.align="center", fig.width = 10-------------------------------------------
scanMap 


## ----include = TRUE------------------------------------------------------------------------------
scan.buffers %>%
  bind_cols(scan) %>%
  top_n(n=3,wt=pred) %>%
  as.data.frame() %>%
  mutate(Mean_Predicted_Count = round(pred)) %>%
  dplyr::select(SCAN.s.Trauma.Informed.Community.Network..TICN.,Street.Address, Mean_Predicted_Count) %>%
  arrange(-Mean_Predicted_Count) %>%
  rename("Scan Center" = SCAN.s.Trauma.Informed.Community.Network..TICN.,
         "Street Address" = Street.Address,
         "Mean Predicted Count" = Mean_Predicted_Count) %>% 
  kable() %>% 
  kable_styling()


## ----include = TRUE, fig.width=7, fig.height=6.5, fig.align="center"-----------------------------
ggmap(cps_base_map) +
  geom_sf(data=ll(nbr_statAreas.gap), aes(fill=gap), inherit.aes = FALSE,alpha = 0.8) +
  scale_fill_gradientn(colours=c(re,"green2", bl),
                       limits = c(-100,100),
                       breaks = c(-100,100),
                       labels = c("More risk\nthan protection","More protection\nthan risk")) +
  labs(title = "Gap analysis:\nComparing relative risk to protective resources",
       subtitle = "By Neighborhood Statistical Area",
       caption = "Figure 2.15") +
  mapTheme() +
  theme(legend.position = 'bottom',
        legend.title=element_blank(),
        legend.text=element_text(size=10, face = "bold"),
        legend.justification = "center") +
  guides(fill= guide_colorbar(barwidth=15,barheight=2))


## ----eval=TRUE, echo=FALSE, include = TRUE, cache=FALSE------------------------------------------
datasetsTable <- rbind(
  data.frame(Dataset = "Crime",
             Description = "Locations of reported violent and nonviolent incidents that represent risk in the city"),
  data.frame(Dataset = "Points of Interest",
             Description = "Places in Richmond that characterize environmental factors and are categorized as either risk (i.e. Pawnbrokers, Payday loan locations, ABC stores) or protective (i.e. Libraries, Community Centers, Grocery stores) features"),
  data.frame(Dataset = "Businesses",
             Description = "Businesses that either characterize nuisance land uses or protective land uses"),
  data.frame(Dataset = "Code Violations",
             Description = "Indicates the safety and overall quality of structures"),
  data.frame(Dataset = "CPS Accepted Events",
             Description = "Locations of child matreatment incidents")
)

datasetsTable  %>% 
  kable(., format = "html") %>% 
  kable_styling()


## ----cache = TRUE,echo=FALSE, warning=FALSE, include = TRUE, fig.height=10-----------------------
fishnet1 <- ggmap(cps_base_map) +
  geom_sf(data=ll(net_hood), aes(), inherit.aes = FALSE, alpha = 0.8) +
  labs(title="Fishnet - Richmond, VA") +
  mapTheme() +
  theme(axis.line = element_blank(),
        plot.title = element_text(hjust = 0, face = "plain", size = 14, family = "sans"))

CPS_RATE_BY_FISHNET_PLOT2 <- ggmap(cps_base_map) +
  geom_sf(data = ll(fishnet_pop_cps_rate_cut), aes(fill = cut_val), inherit.aes = FALSE, color = NA, alpha = 0.8) +
  labs(title = "CPS rate per 100 people",
       caption = "Figure 4.2") +
  scale_fill_viridis_d(na.value = NA, option = "D", direction = 1, name = "CPS Rate\nper 100") +
  mapTheme() +
  theme(plot.title = element_text(size = 14, family = "sans", face = "plain", hjust = 0),
    plot.subtitle=element_text(size = 11, family = "sans", hjust = 0),
    plot.caption=element_text(size = 10, family = "sans", face = "italic", hjust = 0),
    axis.line = element_blank(),
    legend.title = element_text(size = 10, family = "sans"),
    legend.text = element_text(size = 9, family = "sans"))

cowplot::plot_grid(fishnet1, CPS_COUNT_BY_FISHNET_PLOT, CPS_RATE_BY_FISHNET_PLOT2, ncol = 1, align = "hv", axis = "lrbt")


## ----message=FALSE, warning=FALSE, cache=TRUE, echo=FALSE, include=TRUE, fig.width = 8-----------
#pull out rec centers
rec_centers <- var_list$CommunityCenters
rec_centers_count_map <-
ggmap(cps_base_map) +
  geom_sf(data = ll(rec_centers), inherit.aes = FALSE, color = "red", size = 2, alphae = 0.9) +
  geom_sf(data = ll(nbr_diss), inherit.aes = FALSE, color = "black", fill = NA, size = 1) +
  labs(title = "Community center locations") +
  #scale_fill_viridis_d(na.value = NA, option = "D", direction = 1, name = "Count") +
  mapTheme()


#grab and map count
rec_centers_agg <- ALL_FEATURES %>% 
  dplyr::select(agg_CommunityCenters, net_id) %>% 
  left_join(., net_Richmond, by = "net_id") %>% 
  mutate(value = ifelse(is.na(agg_CommunityCenters), 0, agg_CommunityCenters)) %>% 
  st_sf()
rec_centers_agg <- make_cuts(rec_centers_agg, "value", cuts = "breaks", n_breaks = 5)
rec_centers_agg_map <-
ggmap(cps_base_map) +
  geom_sf(data = ll(rec_centers_agg), aes(fill = cut_val), inherit.aes = FALSE, color = NA, alpha = 0.8) +
  # geom_sf(data = ll(rec_centers), inherit.aes = FALSE, color = "red", size = 4) +  
  #geom_sf(data = ll(nbr_diss), inherit.aes = FALSE, color = "black", fill = NA, size = 2) +
  labs(title = "Community center count\nby fishnet") +
  scale_fill_viridis_d(na.value = NA, option = "D", direction = 1, name = "Count") +
  mapTheme()

#grab and counteuclidean distance
rec_centers_ed <- ALL_FEATURES %>% 
  dplyr::select(ed_CommunityCenters, net_id) %>% 
  left_join(., net_Richmond, by = "net_id") %>% 
  mutate(value = ifelse(is.na(ed_CommunityCenters), 0, ed_CommunityCenters)) %>% 
  st_sf()
rec_centers_ed <- make_cuts(rec_centers_ed, "value",  "value", cuts = "breaks", n_breaks = 5)
rec_centers_ed_map <-
ggmap(cps_base_map) +
  geom_sf(data = ll(rec_centers_ed), aes(fill = cut_val), inherit.aes = FALSE, color = NA, alpha = 0.8) +
  geom_sf(data = ll(rec_centers), inherit.aes = FALSE, color = "red", size = 1, alpha = 0.9) +
  #geom_sf(data = ll(nbr_diss), inherit.aes = FALSE, color = "black", fill = NA, size = 2) +
  labs(title = "Community center\neuclidean distance",
       caption = "Figure 4.3") +
  scale_fill_viridis_d(na.value = NA, option = "D", direction = 1, name = "NN Distance") +
  mapTheme()

#grab and map avg nn 
rec_centers_NN <- ALL_FEATURES %>% 
  dplyr::select(NN_CommunityCenters, net_id) %>% 
  left_join(., net_Richmond, by = "net_id") %>% 
  mutate(value = ifelse(is.na(NN_CommunityCenters), 0, NN_CommunityCenters)) %>% 
  st_sf()
rec_centers_NN <- make_cuts(rec_centers_NN, "value",  "value", cuts = "breaks", n_breaks = 5)
rec_centers_NN_map <-
ggmap(cps_base_map) +
  geom_sf(data = ll(rec_centers_NN), aes(fill = cut_val), inherit.aes = FALSE, color = NA, alpha = 0.8) +
  geom_sf(data = ll(rec_centers), inherit.aes = FALSE, color = "red", size = 1, alpha = 0.9) +
  #geom_sf(data = ll(nbr_diss), inherit.aes = FALSE, color = "black", fill = NA, size = 2) +
  labs(title = "Community center average\nnearest neighbor distance") +
  scale_fill_viridis_d(na.value = NA, option = "D", direction = 1, name = "NN Distance") +
  mapTheme()

cowplot::plot_grid(rec_centers_count_map,rec_centers_agg_map,rec_centers_ed_map, rec_centers_NN_map, align="hv", axis = "rlbt", ncol=2)


## ----cache = TRUE,echo=FALSE,warning=FALSE, message = FALSE, include = TRUE----------------------
ggmap(cps_base_map) +
  geom_sf(data=ll(net_hood), aes(fill=name), inherit.aes = FALSE, color = NA, alpha = 0.8) +
  scale_fill_viridis_d() +
  guides(fill=FALSE) +
  labs(title="LOGOCV Neighborhoods",
       caption = "Figure 4.4") +
  mapTheme()


## ----countEventsFishnet2, warning=FALSE, message = FALSE, cache = TRUE, echo=FALSE, include = TRUE, fig.align = "center", fig.width=8----
lubridate::year(var_list[["CPS_Accepted"]]$RDate) %>%
  data.frame(year = .) %>%
  ggplot(., aes(x = year)) +
  geom_histogram() +
  labs(title="Count of maltreatment events by year",
       caption = "Figure 5.1") +
  plotTheme()


## ----warning=FALSE, message = FALSE, cache = TRUE, echo=FALSE, include = TRUE, fig.align = "center", fig.width = 10----
CPS_KDE_BY_YEAR_plot


## ----warning=FALSE, message = FALSE, cache = TRUE, echo=FALSE, include = TRUE, fig.align = "center", fig.width=12, fig.height=12----
PROTECTIVE_KDE_FACET_PLOT


## ----warning=FALSE, message = FALSE, cache = TRUE, echo=FALSE, include = TRUE, fig.align = "center", fig.width=12, fig.height=12----
RISK_KDE_FACET_PLOT


## ----moranIExp, cache = TRUE, warning=FALSE,echo=FALSE, include = TRUE,include = TRUE, fig.height=10----
MORANS_I_P_plot


## ----cache = TRUE, warning=FALSE,echo=FALSE, include = TRUE, fig.height=18-----------------------
CRIME_CPS_VS_NN_plot


## ----cache = TRUE, warning=FALSE,echo=FALSE, include = TRUE, fig.height=14-----------------------
VIOLATION_CPS_VS_NN_plot


## ----cache = TRUE, warning=FALSE,echo=FALSE, include = TRUE, fig.height = 17---------------------
BUSINESS_CPS_VS_NN_plot


## ----corrplotStrongRisk,  cache = TRUE, warning=FALSE,echo=FALSE,fig.height=10,fig.width=10.5, include = TRUE----
#features_risk_strong_plot <- features_risk_strong %>%
# dplyr::select(-net_id)
risk_cp <- feature_corrplot(features_risk_strong_plot, "Correlation of risk features")
#CORR_RISK_FEATURES_plot


## ----corrplotStrongProtective, cache = TRUE, warning=FALSE, echo=FALSE, include = TRUE, fig.height=10,fig.width=10.5----
#features_protective_strong_plot <- features_protective_strong %>%
#dplyr::select(-net_id)
protective_cp <- feature_corrplot(features_protective_strong_plot, "Correlation of protective features")
#CORR_PROTECTIVE_FEATURES_plot


## ----include = TRUE------------------------------------------------------------------------------
Model_Error_Results_table


## ----include = TRUE, fig.align = "center"--------------------------------------------------------
plot_fold_pred(ens_cv_tbl$pred, ens_cv_tbl$test_y, type = "fit") +
  labs(x = "Observed Maltreatment Counts",
       y = "Predicted Maltreatment Counts",
       title = "Predicted vs. observed maltreatment counts",
       caption = "Figure 6.1") +
  plotTheme() +
  theme(panel.border = element_blank())


## ----include = TRUE, fig.height=12---------------------------------------------------------------
cowplot::plot_grid(POISSON_MODEL_PREDICTION_MAP_plot, RF_MODEL_PREDICTION_MAP_plot, SARLM_MODEL_PREDICTION_MAP_plot, META_MODEL_PREDICTION_MAP_plot, ncol = 1, align = "hv", axis = "lrbt")


## ----include = TRUE, fig.width=12, fig.height=5--------------------------------------------------
cowplot::plot_grid(LOGDEV_BY_NEIGHBORHOOD_plot,MAE_BY_NEIGHBORHOOD_plot,align = "h")


## ----include = TRUE, fig.height=7.5, fig.width=10------------------------------------------------
RF_VARIMP_PLOT


## ----include = TRUE------------------------------------------------------------------------------
data.frame(rbind(cbind(Poverty="Low",Log_score=.658,Count_events=668),
                 cbind(Poverty="High",Log_score=.635,Count_events=1145))) %>%
  kable() %>%
  kable_styling()

data.frame(rbind(cbind(Non_White="Low",Log_score=.757,Count_events=524),
                 cbind(Non_White="High",Log_score=.471,Count_events=1289))) %>%
  kable() %>%
  kable_styling()


## ----include = TRUE, fig.width=10----------------------------------------------------------------
REALTIVE_SENSITIVITY_PREDICTIONS2 <- ggmap(cps_base_map) +
  geom_sf(data = ll(p.summ), aes(fill = factor(sens_group)), 
          color = NA, alpha = 0.8, inherit.aes = FALSE) +
  geom_sf(data = ll(cps_dissolve), inherit.aes = FALSE, size = 1) +
  scale_fill_viridis_d(na.value = NA, option = "D", direction = 1,
                       name = "Risk\nCategory") +
  labs(title = "Risk categories from Meta-Model") +
  mapTheme()

cowplot::plot_grid(REALTIVE_SENSITIVITY_KDE, REALTIVE_SENSITIVITY_PREDICTIONS2, align = "h")


## ----include = TRUE, fig.width=8-----------------------------------------------------------------
ggplot(data=countComparisonsLong, aes(Category,Value)) +
  geom_bar(aes(fill = Variable), position = "dodge", stat="identity", color = NA) +
  scale_fill_viridis_d(name = " ",
                       labels=c("Risk & Protective\nFeatures", "Kernel Density")) +
  labs(x= "Predicted Risk Levels",
       y="Percent of Test Set Cases",
       title= "Goodness of fit: Spatial risk model vs. Kernel Density",
       caption = "Figure 6.5") +
  plotTheme()


## ----APPENDIX_1a, eval=FALSE, echo=TRUE, include = TRUE, cache=FALSE-----------------------------
## # Read data filenames
## files <-list.files(file.path(base_dir,"/data"), pattern = "*\\.xls$|*\\.csv$")
## var_list <- vector(mode = "list")
## var_names <- NULL
## # Loop over files and load
## for(i in seq_along(files)){
##   filename <- str_sub(files[i], start = 1, end = -5)
##   sf_i <- tryCatch({
##     if(tools::file_ext(files[i]) == "xls"){
##       dat <- readxl::read_xls(file.path(base_dir,"data",files[i]))
##     } else if(tools::file_ext(files[i]) == "csv"){
##       dat <- read.csv(file.path(base_dir,"data",files[i]))
##     }
##     dat %>%
##       filter(!is.na(X) | !is.na(Y)) %>%
##       # Make files spatial `sf` object and assign geographic projection
##       st_as_sf(., coords = c("X", "Y"), crs = 102747)
##   }, error = function(e){
##     cat(filename, "error = ",e$message,"\n")
##     return(e)
##   }
##   )
##   # add dataframes to list
##   if(!inherits(sf_i, "error")){
##     var_list[[length(var_list)+1]] <- sf_i
##     var_names[length(var_list)] <- filename
##   }
## }
## names(var_list) <- var_names


## ----APPENDIX_1b, eval=FALSE, echo=TRUE, include = TRUE, cache=FALSE-----------------------------
## # Get Richmond neighborhoods
## nbr <- read_sf("https://data.richmondgov.com/resource/7juf-nwis.geojson")  %>%
##   st_sf() %>%
##   st_transform(102747)
## # Create neighborhood outline
## nbr_diss <- nbr %>%
##   mutate(dissolve = 1) %>%
##   # get rid of slivers
##   st_buffer(., dist = 0.1) %>%
##   group_by(dissolve) %>%
##   summarise()
## 
## # Create neighborhood raster
## nbr_rast_SP <- raster(as(nbr_diss, "Spatial"), nrows = 2000, ncol = 2000)
## 
## # Query Population data
## vars10 <- c("P0010001") # total population
## # Get total 2010 census population for blocks & calculate area
## richmond_block <- get_decennial(geography = "block", variables = vars10, year = 2010,
##                                 summary_var = "P0010001", state = 51, county = 760, geometry = TRUE) %>%
##   st_transform(crs = 102747)
## 
## # Calculate area
## richmond_block <- richmond_block %>%
##   mutate(acre = as.numeric(st_area(richmond_block)*2.29568e-5),
##          # acre = units::set_units(acre, acre),
##          pop_acre_rate = value / acre)


## ----APPENDIX_1c, eval=FALSE, echo=TRUE, include = TRUE, cache=FALSE-----------------------------
## # set global parameter for fishnet grid dimensions
## fishnet_grid_dim = 1200
## # Fishnet creation
## net <- st_make_grid(nbr, cellsize = fishnet_grid_dim)
## 
## # Count CPS incidents per net cell
## net_agg <- aggregate(cps_dissolve, net, sum) %>%
##   tibble::rowid_to_column(.,"net_id")
## 
## # List of net cells IDs that intersect with Richmond Neighborhoods
## net_intersect <- st_intersects(nbr, net_agg)
## 
## # Extract Richmonds net cells based on intersect ID
## net_Richmond <- net_agg[unique(unlist(net_intersect)),]
## 
## # Join neighborhood attributes to fishnet grid cells
## net_hood <- st_join(net_Richmond, nbr, largest = TRUE)
## 
## # Calculate spatial neighborhood matrix
## listw <- nb2listw(poly2nb(as(net_Richmond, "Spatial"), queen = TRUE))


## ----APPENDIX_1d, eval=FALSE, echo=TRUE, include = TRUE, cache=FALSE-----------------------------
## # Compute intersection of census blocks and fishnet
## net_blocks_intersect <- st_intersection(richmond_block, net_Richmond)
## 
## # ... Calculate percent of each blocks area within each fishnet cell and divide block pop by that percent
## 
## # Summerise population by fishnet cell
## fishnet_pop <- net_blocks_intersect %>% # xcc
##   group_by(net_id) %>%
##   summarise(net_pop = sum(intersect_pop)) %>%
##   filter(net_pop > 0)   # <-  zeros or no zeros!!!!
## 
## # ... make cps_agg which is aggregate of each variable by cell
## 
## # Spatial join of fishnet_pop and fishnet_cps to then calculate rate for all CPS features
## fishnet_pop_cps <- st_join(fishnet_pop, CPS_agg, join = st_equals) %>%
##   mutate_at(vars(paste0("net_",CPS_vars)), funs(rate = ./(net_pop/100)))  %>% # cps per 100 person
##   rename_at(vars( contains( "_rate")), funs(paste("rate", gsub("net_|_rate", "", .), sep = "_"))) %>%
##   replace(is.na(.), 0) # replace NA with zero
## 


## ----APPENDIX_2a, eval=FALSE, echo=TRUE, include = TRUE, cache=FALSE-----------------------------
## # Create list of Protective variables
## protective_class <- c("CommunityCenters","FireStations",
##                       "HomelessShelters","Libraries","Parks","PointsOfInterest",
##                       "PoliceStations","PublicSchools","ResourceOASIS","SNAP_WIC",
##                       "VotingStations")
## # Add column for variable set name and assign to new list
## protective_vars <- list()
## for(i in seq_along(protective_class)){
##   dat <- var_list[[protective_class[i]]] %>%
##     mutate(feature_name = protective_class[i],
##            class = "protective") %>%
##     dplyr::select(feature_name, class)
##   protective_vars[[i]] <- dat
## }
## 
## # Pull protective features from business features
## Businesses_protective <- var_list[["BusinessProject"]] %>%
##   filter(Classification == "PROTECTIVE") %>%
##   mutate(feature_name = "BusinessProject",
##          class = "protective") %>%
##   dplyr::select(feature_name, class)
## # Add to protective feature list
## protective_vars[[length(protective_vars)+1]] <- Businesses_protective
## # Turn into a spatial sf object with columsn for feature name, protective vs. risk class, and point coordinates
## protective_vars <- do.call(rbind, protective_vars)
## # Add dataframe to varaibles list
## var_list[["Protective"]] <- protective_vars
## 
## ## Simialr to above, but for Risk variables
## risk_class <- c("BusStops")
## risk_vars <- list()
## for(i in seq_along(risk_class)){
##   dat <- var_list[[risk_class[i]]] %>%
##     mutate(feature_name = risk_class[i],
##            class = "risk") %>%
##     dplyr::select(feature_name, class)
##   risk_vars[[i]] <- dat
## }
## # Business risks
## Business_risk <- var_list[["BusinessProject"]] %>%
##   filter(Classification == "RISK") %>%
##   mutate(feature_name = "BusinessProject",
##          class = "risk") %>%
##   dplyr::select(feature_name, class)
## risk_vars[[length(risk_vars)+1]] <- Business_risk
## # Crime realted risks
## CrimeData_risk <- var_list[["CrimeData"]] %>%
##   filter(OFFENSE %in% c('DRUG/NARCOTIC VIOLATION','SIMPLE ASSAULT, DOMESTIC',
##                         'Runaway','AGGRAVATED ASSAULT DOMESTIC')) %>%
##   mutate(feature_name = "CrimeData",
##          class = "risk") %>%
##   dplyr::select(feature_name, class)
## risk_vars[[length(risk_vars)+1]] <- CrimeData_risk
## # Code violation risks
## Violations_III_ks_risk <- var_list[["Violations_III_ks"]] %>%
##   filter(CodeDsrp %in% c('General Violations','Unsafe Structure',
##                          'Unfit Structure')) %>%
##   mutate(feature_name = "Violations_III_ks_risk",
##          class = "risk") %>%
##   dplyr::select(feature_name, class)
## risk_vars[[length(risk_vars)+1]] <- Violations_III_ks_risk
## # Bind into dataframe
## risk_vars <- do.call(rbind, risk_vars)
## # Add to var_list
## var_list[["Risk"]] <- risk_vars


## ----APPENDIX_2b, eval=FALSE, echo=TRUE, include = TRUE, cache=FALSE-----------------------------
## # Functions to make thee feature types:
## # Mean Nearest Neighbor distance
## NN_point_features <- function(var_list, fishnet, k){
##   NN_results <- foreach(i = seq_along(var_list),
##                         .export=c('nn_function'),
##                         .packages=c('raster', 'sf', 'dplyr', "FNN", "tibble", "tidyr")) %dopar% {
##                           feature <- names(var_list)[i]
##                           # Create centroid of fishnet cell
##                           fishnet_centroid_XY <- st_coordinates(st_centroid(fishnet))
##                           dat <- var_list[[i]]
##                           # If more features the `k` nearest neighbors...
##                           if(nrow(dat) >= k){
##                             # get mean distance and id for `k` nearest neighbors ...
##                             # from center of each fishnet cell
##                             net_NN <- nn_function(fishnet_centroid_XY,
##                                                   st_coordinates(dat)[,1:2], k) %>%
##                               mutate(feature_name = paste0("NN_",feature),
##                                      net_id = fishnet$net_id) %>%
##                               left_join(., fishnet, by = "net_id") %>%
##                               rename("value" = value.x) %>%
##                               dplyr::select(-value.y) %>%
##                               st_as_sf()
##                           } else {
##                             # If there are not enough features to make NN, then it is NA
##                             net_NN <- data.frame(value = rep(NA, nrow(fishnet))) %>%
##                               mutate(feature_name =  paste0("NN_",feature),
##                                      net_id = fishnet$net_id) %>%
##                               left_join(., fishnet, by = "net_id") %>%
##                               rename("value" = value.x) %>%
##                               dplyr::select(-value.y) %>%
##                               st_as_sf()
##                           }
##                         }
##   names(NN_results) <- paste0("NN_",names(var_list))
##   return(NN_results)
## }
## 
## # Nearest Euclidean Distance
## Euclidean_point_features <- function(var_list, dist_raster, raster_mask, fishnet){
##   ED_results <- foreach::foreach(i = seq_along(var_list),
##                                  .combine='comb', .multicombine=TRUE,
##                                  .init=list(list(), list()),
##                                  .export=c('distanceFromPoints', 'raster_to_fishnet'),
##                                  .packages=c('raster', 'sf', 'dplyr')) %dopar% {
##                                    feature <- names(var_list)[i]
##                                    # Calculate distance raster from all point features
##                                    bs_dist <- distanceFromPoints(dist_raster,
##                                                                  sf::st_coordinates(var_list[[feature]]))
##                                    # Clip raster to Richmond City
##                                    bs_clip <- raster::mask(bs_dist, mask = as(raster_mask, "Spatial"))
##                                    # Extract raster distance values within each fishnet cell and take mean
##                                    fea_mean_dist <- raster_to_fishnet(bs_clip,fishnet,paste0("ed_",feature))
##                                    list(fea_mean_dist, bs_clip)
##                                  }
##   # Add results to list
##   dist_results <- ED_results[[1]]
##   dist_rasters <- ED_results[[2]]
##   names(dist_results) <- paste0("ED_",names(var_list))
##   names(dist_rasters) <- paste0("ED_",names(var_list))
##   return(list(dist_results, dist_raster))
## }
## 
## # Aggregare Point Count
## Aggregate_points_Features <- function(var_list, fishnet){
##   agg_results <- foreach(i = seq_along(var_list),
##                          .packages=c('raster', 'sf', 'dplyr')) %dopar% {
##                            feature <- names(var_list)[i]
##                            dat <- var_list[[i]] %>%
##                              mutate(value = 1) %>%
##                              dplyr::select(value)
##                            # Count feature points within each fishnet grid cell
##                            net_agg <- aggregate(dat, fishnet, sum) %>%
##                              mutate(feature_name = paste0("agg_",feature),
##                                     net_id = fishnet$net_id)
##                          }
##   names(agg_results) <- paste0("agg_",names(var_list))
##   return(agg_results)
## }
## 
## #Make `ALL_FEATUES` dataframe by joining all Mean Nearest Neighbor, Mean Euclidean Distance, and Aggregare Count features
## ALL_FEATURES <- full_join(NN_features, agg_features, by = "net_id") %>%
##   full_join(.,ED_features, by = "net_id") %>%
##   # Join in census data features
##   full_join(.,sf1_features, by = "net_id")
## 
## # Clean up joined fields and replace `NA` with zeros
## ALL_FEATURES <- ALL_FEATURES %>%
##   dplyr::select(-cps_rate.y, -cps_rate.x.x, -cps_rate.y.y,
##                 -cps_net.y, -cps_net.x.x, -cps_net.y.y,
##                 -net_pop.y, -net_pop.x.x, -net_pop.y.y) %>%
##   dplyr::select(-contains("_CPS_")) %>%
##   dplyr::rename(cps_net  = cps_net.x,
##                 cps_rate = cps_rate.x,
##                 net_pop  = net_pop.x) %>%
##    mutate_all(funs(replace(., is.na(.), 0)))  %>%
##   dplyr::rename_all(funs(make.names(.)))
## 


## ----APPENDIX_2c, eval=FALSE, echo=TRUE, include = TRUE, cache=FALSE-----------------------------
## # Compute correlation between all pariwise features
## cps_cor_ALL <- cor(ALL_FEATURES)
## All_cors <- cps_cor_ALL[,"cps_net"]
## # Compute p-values
## p.mat_ALL <- cor.mtest(ALL_FEATURES)$p
## p.mat_ALL <- p.mat_ALL[,which(colnames(cps_cor_ALL)=="cps_net")]
## # Prepare data for plotting
## cor_ALL_plot <- data.frame(feature = names(All_cors),
##                            cor = as.numeric(All_cors),
##                            p_value   = p.mat_ALL) %>%
##   filter(!(feature %in% c("cps_rate","cps_net","net_pop","net_cps","net_id"))) %>%
##   filter(!(feature %in% grep("CPS", names(All_cors),value=T))) %>%
##   arrange(desc(cor)) %>%
##   mutate(p_value = ifelse(p_value >= 0.1, "Not Significant", "Significant"))
## cor_ALL_plot$feature <- factor(cor_ALL_plot$feature,
##                                levels=cor_ALL_plot[order(cor_ALL_plot$cor,
##                                                          decreasing=F),]$feature)
## 
## # For wach Protective variable, extract the feature type (mean NN, mean distance, aggregate) that as the highest absolute correlation
## features_strong_protective_names <- cor_ALL_plot %>%
##   filter(feature %in% names(features_protective_all)) %>%
##   mutate(prefix = str_extract(feature, "^[^_]+(?=_)"),
##          suffix = str_extract(feature, "(?<=_)[^_].*"),
##          feature = as.character(feature)) %>%
##   group_by(suffix) %>%
##   # Highest absolute correlation
##   slice(which.max(abs(cor)))
## features_protective_strong <- features_protective_all %>%
##   dplyr::select(features_strong_protective_names$feature,
##                 NN_CPS_Accepted,
##                 cps_net, cps_rate, net_pop, net_id)
## 
## # Do the same for Protective


## ----APPENDIX_2d, eval=FALSE, echo=TRUE, include = TRUE, cache=FALSE-----------------------------
## # Join Risk and Protective features to selected census feaures to form data for regression models
## og_dat <- full_join(features_risk_strong, features_census_select, by = "net_id") %>%
##   full_join(., features_protective_strong, by = "net_id") %>%
##   dplyr::select(-net_pop.y, -cps_net.y, -cps_rate.y,
##                 -net_pop.x, -cps_net.x, -cps_rate.x,
##                 -NN_CPS_Accepted.y) %>%
##   rename("NN_CPS_Accepted" = NN_CPS_Accepted.x)
## 
## # Remove unneeded columns and center & scale all variables except the dependent variable
## dat    <- og_dat %>% dplyr::select(-cps_rate, -net_pop, -net_id) %>%
##   mutate_at(vars(-cps_net), scale_this)
## # Add neighborhood name
## net_hood <- st_join(net_Richmond, nbr, largest = TRUE)
## og_dat$.block_id <- net_hood$name


## ----APPENDIX_3a, eval=FALSE, echo=TRUE, include = TRUE, cache=FALSE-----------------------------
## # bin the countof maltreatment events
## fishnet_pop_cps_cut <- fishnet_pop_cps %>%
##   mutate(net_CPS_Accepted = ifelse(is.na(net_CPS_Accepted), 0, net_CPS_Accepted)) %>%
##   make_cuts(., "net_CPS_Accepted", cuts = "breaks", n_breaks = 10)
## # plot maltreatment event counts on top of Richmond City basemap
## CPS_COUNT_BY_FISHNET_PLOT <- ggmap(cps_base_map) +
##   geom_sf(data = ll(fishnet_pop_cps_cut), aes(fill = cut_val), inherit.aes = FALSE, color = NA) +
##   labs(title = "CPS Count per Fishnet Cell") +
##   scale_fill_viridis_d(na.value = NA, option = "D", direction = 1, name = "CPS Count") +
##   theme_bw()


## ----APPENDIX_3b, eval=FALSE, echo=TRUE, include = TRUE, cache=FALSE-----------------------------
## # Normalizes values by month
## CPS_normalized_by_month <- st_drop_geometry(var_list[["CPS_Accepted"]]) %>%
##     mutate(month = lubridate::month(RDate),
##            year  = lubridate::year(RDate)) %>%
##     group_by(year, month) %>%
##     summarise(m_total = n()) %>%
##     arrange(month, year) %>%
##     dplyr::select(month, year, m_total) %>%
##     ungroup() %>%
##     group_by(month) %>%
##     # Normalize
##     mutate(m_mean = mean(m_total),
##            m_sd   = sd(m_total),
##            m_z    = (m_total - m_mean) / m_sd)
## # Plot
## CPS_LINE_NORMALIZED_plot <- ggplot(CPS_normalized_by_month, aes(x = as.factor(month),
##                                           y = m_z, group = year,
##                                           color = as.factor(year))) +
##     geom_line() +
##     geom_hline(yintercept = 0, color = "gray20", linetype = "dashed") +
##     scale_y_continuous(limits = c(-2,2)) +
##     theme_bw()


## ----APPENDIX_3c, eval=FALSE, echo=TRUE, include = TRUE, cache=FALSE-----------------------------
## # Risk KDE facet plot
## RISK_KDE_FACET_PLOT <- ggmap(cps_base_map) +
##   geom_tile(data = risk_plot_dat,
##             aes(x,y,fill = as.factor(ntile(value,brks)),
##                 group = variable), alpha=0.5) +
##   scale_fill_viridis_d(name = variable) +
##   facet_wrap(~variable) +
##   theme(
##     axis.title=element_blank(),
##     axis.text=element_blank(),
##     axis.ticks=element_blank(),
##     legend.key = element_rect(fill = "white"),
##     strip.text = element_text(face = "bold", size = 12),
##     legend.position = "none"
##   )


## ----APPENDIX_3d, eval=FALSE, echo=TRUE, include = TRUE, cache=FALSE-----------------------------
## # Calculate fishent grid nearest neighbors using the "Queen" case (all directions)
## # convert to a neighborhood graph
## fishnet_knn <- knn2nb(knearneigh(fishnet_coords, k_direction))
## # List of neighborhood weights
## fishnet_Weights <- nb2listw(fishnet_knn, style="W")
## # Simulate Global Moran's I with `999` simulations
## globalMorans <- moran.mc(fishnet_pop_cps$net_CPS_Accepted, fishnet_Weights, nsim=999)
## 
## # Plot Global Moran's I simulated distribution and observed value
## # note: the `globalMorans` I object will have 1000 rows, but only the first 999 are simulated.
## # The 1000th row is the observed value
## GLOBAL_MORANS_PERMUTATION_plot <- ggplot(data.frame(res = globalMorans$res)[1:999,,0], aes(res)) +
##   geom_histogram(binwidth = 0.01) +
##   geom_vline(aes(xintercept = globalMorans$statistic), colour = "red",size=1) +
##   scale_x_continuous(limits = c(-1, 1)) +
##   labs(title="Observed and permuted Moran's I", x = "Simulated Moran's I Value") +
##   theme_bw()


## ----APPENDIX_3e, eval=FALSE, echo=TRUE, include = TRUE, cache=FALSE-----------------------------
## # calculate Local Moran's I for each fishnet cell
## localMorans  <- as.data.frame(localmoran(fishnet_pop_cps$net_CPS_Accepted, fishnet_Weights))
## # Moran's I join
## fishnet_pop_cps_morans        <- fishnet_pop_cps
## fishnet_pop_cps_morans$Ii     <- localMorans$Ii
## fishnet_pop_cps_morans$pvalue <- localMorans$`Pr(z > 0)`
## fishnet_pop_cps_morans        <- cbind(fishnet_coords, fishnet_pop_cps_morans)
## 
## # Bin maltreatment event counts
## fishnet_pop_cps_morans_cut <- make_cuts(fishnet_pop_cps_morans, "net_CPS_Accepted",
##                                         cuts = "breaks", n_breaks = 10)
## # Plot maltreatment events (assign to a variable)
## plot_cps <- ggmap(cps_base_map) +
##   geom_sf(data = ll(fishnet_pop_cps_morans_cut), aes(fill = cut_val),
##           color = NA, inherit.aes = FALSE, alpha = 0.9) +
##   scale_fill_viridis_d(na.value=NA, name = "Maltreatment Events", option = "D" ) +
##   theme_void() +
##   theme(
##     legend.position = "right"
##   )
## 
## # Bin Local Moran's I statistic
## Ii_cut <- fishnet_pop_cps_morans %>%
##   mutate(Ii_cut_val = as.character(Hmisc::cut2(.$Ii,
##                                                cuts = as.numeric(quantile(round(fishnet_pop_cps_morans$Ii,2),
##                                                                           na.rm=T, p = seq(0,1,0.25))))))
## # plot binned Local Moran's I statistic (assign to a variable)
## plot_Ii <- ggmap(cps_base_map) +
##   geom_sf(data = ll(Ii_cut), aes(fill = Ii_cut_val),
##           color = NA, inherit.aes = FALSE, alpha = 0.9) +
##   scale_fill_viridis_d(na.value=NA, name = "Local Moran's I", option = "D")+
##   theme_void()+
##   theme(
##     legend.position = "right"
##   )
## 
## # Bin Local Moran's I p-value
## p_cut <- fishnet_pop_cps_morans %>%
##   mutate(pval_cut = ifelse(pvalue > 0.05, "Not Significant", "Significant"))
## # Plot binned p-value (assign to a variable)
## plot_p <- ggmap(cps_base_map) +
##   geom_sf(data = ll(p_cut), aes(fill = pval_cut),
##           color = NA, inherit.aes = FALSE, alpha = 0.9) +
##   scale_fill_viridis_d(na.value=NA, name = "p-value", option = "D")    +
##   theme_void()+
##   theme(
##     legend.position = "right"
##   )
## # use `cowplot` to put plots together
## MORANS_I_P_plot <- cowplot::plot_grid(plot_cps, plot_Ii, plot_p,
##                                       rel_widths = c(0.9,0.9,0.9),
##                                       ncol = 1, align = "v")


## ----APPENDIX_4a, eval=FALSE, echo=TRUE, include = TRUE, cache=FALSE-----------------------------
## # neighborhood fixed effect model
## hood_matrix <- model.matrix(cps_net~.block_id,og_dat)
## hood_model <- lm(sqrt(og_dat$cps_net) ~ hood_matrix)
## dat$hood_fixed <- predict(hood_model, type = "response")^2
## og_dat$hood_fixed <- predict(hood_model, type = "response")^2
## # Create CV neighborhood fold index
## all_hoods <- length(unique(net_hood$name))
## n_folds = ifelse(n_folds == "LOOCV", all_hoods, n_folds)
## folds_index <- groupdata2::fold(og_dat, k = n_folds, id_col = '.block_id')$.folds
## # create tibble with all CV folds and assocaited data
## cv_tbl <- tibble(folds = seq_len(n_folds),
##                  train = NA, train_y = NA, train_index = NA, train_net_id = NA,
##                  test  = NA, test_y  = NA, test_index  = NA, test_net_id  = NA)
## for(k in seq_len(n_folds)){
##   fold_i  <- which(folds_index == k)
##   cv_tbl[k,]$train         <- list(dat[-fold_i,])
##   cv_tbl[k,]$test          <- list(dat[ fold_i,])
##   cv_tbl[k,]$train_y       <- list(og_dat[-fold_i,target_var])
##   cv_tbl[k,]$test_y        <- list(og_dat[ fold_i,target_var])
##   cv_tbl[k,]$train_index   <- list(setdiff(seq(1:nrow(dat)),fold_i))
##   cv_tbl[k,]$test_index    <- list(fold_i)
##   cv_tbl[k,]$train_net_id  <- list(og_dat[-fold_i,"net_id"])
##   cv_tbl[k,]$test_net_id   <- list(og_dat[ fold_i,"net_id"])


## ----APPENDIX_4b, eval=FALSE, echo=TRUE, include = TRUE, cache=FALSE-----------------------------
## # Helper function for GLM fit
## glm_fit <- function(dat, formula, family){
##   glm_model <- glm(formula, data = dat, family = family)
##   return(glm_model)
## }
## # HElper funtion for Random Forest fit
## rf_fit <- function(dat, formula, mtry_add = 0, importance = "none"){
##   mtry <- floor(sqrt(ncol(dat)-1))+mtry_add
##   rf_model <- ranger(formula, data = dat,
##                      mtry = mtry,
##                      splitrule = "variance",
##                      importance = importance,
##                      num.trees = 500,
##                      min.node.size = 10)
##   return(rf_model)
## }
## # Helper function for model scoring
## score_model <- function(dat){
##   dat <- dat %>%
##     mutate(R2     = map2_dbl(pred, test_y, r_squared),
##            MAE    = map2_dbl(pred, test_y, mae),
##            MAAPE  = map2_dbl(pred, test_y, maape),
##            RMSE   = map2_dbl(pred, test_y, rmse),
##            logdev = map2_dbl(pred, test_y, logdev_p))
##   return(dat)
## }
## # Fit the Poisson GLM model
## po_cv_tbl <- cv_tbl %>%
##   mutate(fit   = map(train, glm_fit,
##                      formula =  paste("cps_net ~ ."),
##                      family = "poisson"),
##          pred  = map2(fit, test, lm_predict, sqrt = FALSE),
##          mdl_nam = "GLm - Poisson") %>%
##   score_model()
## # Fit the Random Forest model
## rf_cv_tbl <- cv_tbl %>%
##   mutate(fit   = map(train, rf_fit, formula = "cps_net ~ .", mtry_add = 2, importance = "impurity"),
##          pred  = map2(fit, test, lm_predict),
##          mdl_nam = "Random Forest") %>%
##   score_model()
## # Fit the Spatial Durbin model
## spat_durbin <- errorsarlm(sqrt(cps_net) ~ ., data = dat, listw, etype ="emixed")
## spat_durbin_tbl <- tibble(
##   fit   = list(spat_durbin),
##   pred  = map(fit, sar_pred),
##   test_y= list(dat$cps_net),
##   test_net_id = list(og_dat$net_id),
##   mdl_nam = "Spatial Durbin - sqrt") %>%
##   score_model()


## ----APPENDIX_4c, eval=FALSE, echo=TRUE, include = TRUE, cache=FALSE-----------------------------
## # Extract predictions from Possion model (same is done for other two models)
## po_pred_dat <- po_cv_tbl %>%
##   unnest(pred) %>%
##   mutate(test_y = po_cv_tbl %>% unnest(test_y) %>% pull(test_y),
##          test_net_id = po_cv_tbl %>% unnest(test_net_id) %>% pull(test_net_id))
## # Plot prediction map
## po_pred_geoplot <- model_pred_geoplot(po_pred_dat$pred,
##                                       po_pred_dat$test_y,
##                                       po_pred_dat$test_net_id,
##                                       net_Richmond, cps_base_map, "po")
## # Join all predictions from three models into single data frame
## cps_preds <- og_dat %>%
##   dplyr::select(net_id, cps_net) %>%
##   left_join(., dplyr::select(po_pred_dat,
##                              net_id = test_net_id,
##                              pred_lm = pred), by = "net_id") %>%
##   left_join(., dplyr::select(rf_pred_dat,
##                              net_id = test_net_id,
##                              pred_rf = pred), by = "net_id") %>%
##   left_join(., dplyr::select(sarlm_pred_dat,
##                              net_id = test_net_id,
##                              pred_sarlm = pred), by = "net_id") %>%
##   mutate_if(is.double, round, 2)


## ----APPENDIX_4d, eval=FALSE, echo=TRUE, include = TRUE, cache=FALSE-----------------------------
## # Create CV neighborhood folds for meta-model
## cps_preds_cv_dat <- dplyr::select(cps_preds, -net_id)
## ens_cv_tbl <- tibble(folds = seq_len(n_folds),
##                      train = NA, train_y = NA, train_index = NA, train_net_id = NA,
##                      test  = NA, test_y  = NA, test_index  = NA, test_net_id  = NA)
## for(k in seq_len(n_folds)){
##   fold_i  <- which(folds_index == k)
##   ens_cv_tbl[k,]$train         <- list(cps_preds_cv_dat[-fold_i,])
##   ens_cv_tbl[k,]$test          <- list(cps_preds_cv_dat[ fold_i,])
##   ens_cv_tbl[k,]$train_y       <- list(cps_preds_cv_dat[-fold_i,target_var])
##   ens_cv_tbl[k,]$test_y        <- list(cps_preds_cv_dat[ fold_i,target_var])
##   ens_cv_tbl[k,]$train_index   <- list(setdiff(seq(1:nrow(cps_preds_cv_dat)),fold_i))
##   ens_cv_tbl[k,]$test_index    <- list(fold_i)
##   ens_cv_tbl[k,]$train_net_id  <- list(cps_preds[-fold_i,"net_id"])
##   ens_cv_tbl[k,]$test_net_id   <- list(cps_preds[ fold_i,"net_id"])
## }
## # Fit meta-model with Random Forest
## ens_cv_tbl <- ens_cv_tbl %>%
##   mutate(fit   = map(train, rf_fit, formula = "cps_net ~ pred_rf + pred_sarlm"),
##          pred  = map2(fit, test, lm_predict),
##          # pred  = map(pred, round),
##          mdl_nam = "Meta-Model") %>%
##   score_model()


## ----APPENDIX_4e, eval=FALSE, echo=TRUE, include = TRUE, cache=FALSE-----------------------------
## # Helper function for quantile error
## quantile_error <- function(pred,obs,quant){
##   preds <- data.frame(pred = pred, obs = obs) %>%
##     filter(quantile(seq(0,max(obs)), quant)>obs)
##   return(preds)
## }
## # Join/bind model prediction tables
## models <- bind_rows(rf_cv_tbl, spat_durbin_tbl, ens_cv_tbl, po_cv_tbl)
## # Unnest predictions by model
## CV_preds_long <- models %>%
##   group_by(mdl_nam) %>%
##   unnest(pred, test_y)
## ## Map over all quantiles to get error metrics
## quantile_errors <- CV_preds_long %>%
##   nest(-mdl_nam) %>%
##   mutate(q      = list(seq(0,1,0.01)),
##          pred   = map(data, "pred"),
##          test_y = map(data, "test_y")) %>%
##   dplyr::select(-data) %>%
##   unnest(q, .preserve = c(pred, test_y)) %>%
##   filter(q != 0) %>%
##   mutate(q_dat  = pmap(list(pred, test_y, q), quantile_error),
##          q_pred = map(q_dat, "pred"),
##          q_obs  = map(q_dat, "obs"),
##          q_RMSE = map2_dbl(q_pred, q_obs, rmse),
##          q_MAE  = map2_dbl(q_pred, q_obs, mae),
##          q_logdev  = map2_dbl(q_pred, q_obs, logdev_p),
##          y_max  = quantile(seq(0,max(dat$cps_net)), q),
##          q_cnt  = nrow(og_dat) - map_int(q_dat, nrow))
## 
## # Map over all predictions grouped by model to calculate mean and sd for error metrics
## model_results <- models %>%
##   dplyr::select("Model Name" = mdl_nam, R2, RMSE, MAE, logdev) %>%
##   group_by(`Model Name`) %>%
##   arrange(`Model Name`) %>%
##   summarise(R2_mean      = mean(R2, na.rm=TRUE),
##             R2_sd        = sd(R2, na.rm=TRUE),
##             MAE_mean     = mean(MAE, na.rm=TRUE),
##             MAE_sd       = sd(MAE, na.rm=TRUE),
##             RMSE_mean    = mean(RMSE, na.rm=TRUE),
##             RMSE_sd      = sd(RMSE, na.rm=TRUE),
##             logdev_mean  = mean(logdev, na.rm=TRUE),
##             logdev_sd    = sd(logdev, na.rm=TRUE))
## Model_Error_Results_table <- model_results %>%
##   kable(., format = "html", digits = 3) %>%
##   kable_styling()


## ----APPENDIX_4f, eval=FALSE, echo=TRUE, include = TRUE, cache=FALSE-----------------------------
## # Download statarea
## nbr_statAreas <- read_sf("https://data.richmondgov.com/resource/8kyq-v9j2.geojson") %>%
##   st_transform(crs = 102747) %>%
##   mutate(stat_area_id = id)
## 
## # Download poverty and population data
## tract10 <- get_acs(geography = "tract", variables = c("B02001_001","B02001_002E","B17001_002"),
##                    year = 2010, state=51, county=760, geometry=T)
## # Aggregate data and make percentages
## tract10 <- tract10 %>%
##   dplyr::select(variable,estimate) %>%
##   as.data.frame() %>%
##   spread(variable,estimate) %>%
##   rename(TotalPop=B02001_001,
##          NumberWhites=B02001_002,
##          TotalPoverty=B17001_002) %>%
##   mutate(percentNonWhite = ifelse(TotalPop > 0, ((TotalPop - NumberWhites) / TotalPop),0),
##          percentPoverty  = ifelse(TotalPop > 0, TotalPoverty / TotalPop, 0),
##          tract_id        = row_number()) %>%
##   st_sf() %>%
##   st_transform(102747)
## tract10$tract_area <- st_area(tract10)
## 
## # Aggreate mean errors to statareas
## stat_area_metric_logdev <- error_points %>%
##   aggregate(., nbr_statAreas.spJoin, mean) %>%
##   dplyr::select(logdev) %>%
##   mutate(logdev = round(logdev, 3)) %>%
##   make_cuts(., "logdev")
## stat_area_metric_MAE<- error_points %>%
##   aggregate(., nbr_statAreas.spJoin, mean) %>%
##   dplyr::select(MAE) %>%
##   mutate(MAE = round(MAE, 3)) %>%
##   make_cuts(., "MAE")
## 
## # Aggregate sum of CPS incidents to statarea
## stat_area_cps <- error_points %>%
##   aggregate(., nbr_statAreas.spJoin, sum) %>%
##   dplyr::select(test_y)
## stat_area_errors <- stat_area_metric_logdev %>%
##   st_join(., stat_area_metric_MAE, join = st_equals) %>%
##   st_join(., stat_area_cps, join = st_equals) %>%
##   st_join(., nbr_statAreas.spJoin, join = st_equals)
## 
## # Group by poverty and get median of statarea aggregate errors
## poverty_aggregate <- stat_area_errors %>%
##   group_by(poverty.percentile) %>%
##   summarise(med_dev = round(median(logdev),3),
##             med_MAE = round(median(MAE),3),
##             med_CPS = sum(test_y)) %>%
##   st_drop_geometry() %>%
##   dplyr::select(poverty.percentile, med_dev, med_MAE, med_CPS)


## ----APPENDIX_4g, eval=FALSE, echo=TRUE, include = TRUE, cache=FALSE-----------------------------
## # Function to create sensitivity classes by binning and cutting quantities
## bin_class <- function(dat, bin_col = "pred",
##                       quantile_labels = 100, break_vec = c(-1, 30, 50, 70, 90, 100)){
##   if(is(dat, "sf")){
##     dat <- st_drop_geometry(dat)
##   }
##   pred_bin <- as.numeric(.bincode(dat[,bin_col]+1e-8, # wiggle factor to get above zero
##                                   breaks = quantile(dat[,bin_col],
##                                                     seq(0,1, by=(1/quantile_labels)),
##                                                     na.rm = TRUE,
##                                                     labels = seq(1,quantile_labels,1))))
##   pred_bin_class <- as.numeric(cut(pred_bin,
##                                    breaks = break_vec,
##                                    na.rm  = TRUE,
##                                    labels = seq(1,length(break_vec)-1,1)))
##   pred_bin_class <- ifelse(is.na(pred_bin_class), length(break_vec)-1, pred_bin_class)
## }
## 
## error_geoplot$pred_bin_class <- bin_class(error_geoplot, "pred")
## # Count recorded incidents by meta-model prediction sensitivity classes
## p.summ <- error_geoplot %>%
##   group_by(pred_bin_class) %>%
##   dplyr::summarize(obs.total = sum(test_y),
##                    obs.cnt = n()) %>%
##   rename(sens_group = pred_bin_class) %>%
##   filter(!is.na(sens_group)) %>%
##   identity()
## 
## # Compute KDE
## cps_ppp <- as.ppp(st_coordinates(cps_dissolve), W = st_bbox(net_Richmond))
## cps_KDE <- spatstat::density.ppp(cps_ppp)
## # COnvert KDE to fishnet grid
## cps_KDE_tbl <- as.data.frame(cps_KDE) %>%
##   st_as_sf(coords = c("x", "y"), crs = 102747) %>%
##   aggregate(., net_Richmond, mean) %>%
##   mutate(net_id = net_Richmond$net_id)
## 
## error_geoplot$kde_bin_class  <- bin_class(cps_KDE_tbl, "value")
## # Count incident counts by KDE sensitvity classes
## kde.summ <- error_geoplot %>%
##   group_by(kde_bin_class) %>%
##   dplyr::summarize(kde.total = sum(test_y),
##                    kde.cnt = n()) %>%
##   rename(sens_group = kde_bin_class) %>%
##   filter(!is.na(sens_group)) %>%
##   identity()
## 
## # Compile counts of incidents from KDE and predictions then create percentages
## countComparisons <- merge(st_drop_geometry(p.summ), st_drop_geometry(kde.summ)) %>%
##   mutate_if(is.double, round, 3) %>%
##   mutate(Category = rev(c("90% - 100%", "70% - 89%", "50% - 69%",
##                           "30% - 49%", "1% - 29%"))) %>%
##   dplyr::mutate(kernelPct = round(kde.total / sum(kde.total),4),
##                 fittedPct = round(obs.total / sum(obs.total), 4))
## 
## countComparisonsLong <- countComparisons %>%
##   gather(Variable, Value, kernelPct:fittedPct)
## # Bat plot of results
## REALTIVE_RISK_BARPLOT_COMPARE_plot <- ggplot(data=countComparisonsLong, aes(Category,Value)) +
##   geom_bar(aes(fill = Variable), position = "dodge", stat="identity", color = NA) +
##   scale_fill_viridis_d(name = " ",
##                        labels=c("Risk & Protective\nFeatures", "Kernel Density")) +
##   labs(x= "Predicted Risk Levels",
##        y="Percent of Test Set Cases",
##        title= "Goodness of Fit: Spatial Risk Model vs. Kernel Density hotspot",
##        caption = "Figure 1.6 - Kernel density comparison") +
##   plotTheme() +
##   theme(axis.line = element_blank())


## ----APPENDIX_3_OBJECTS, eval=TRUE, echo=FALSE, include = TRUE, cache=FALSE----------------------
# Objects
objects <- data.frame("Name" = setdiff(ls(), lsf.str()), stringsAsFactors = FALSE) %>% 
  mutate(Class    = purrr::map_chr(Name, ~ class(get(.x))[1]),
         Size_kB  = purrr::map_dbl(Name, ~ round(object.size(get(.x))*0.001,2)),
         Rows     = ifelse(Class %in% c("data.frame","tibble"), 
                           purrr::map(Name, ~ nrow(get(.x))), NA),
         Columns  = ifelse(Class %in% c("data.frame","tibble"), 
                           purrr::map(Name, ~ ncol(get(.x))), NA),
         Length   = ifelse(Class %in% c("numeric","integer","character"), 
                           purrr::map(Name, ~ length(get(.x))), NA)) %>% 
  arrange(Name)

objects %>%
  kable(., format = "html", digits = 3) %>%
  kable_styling()


## ----APPENDIX_3_FUNCTIONS, eval=TRUE, echo=FALSE, include = TRUE, cache=FALSE--------------------
# Objects
funs <-  data.frame("Name" = ls(), stringsAsFactors = FALSE) %>% 
  mutate(Class    = purrr::map_chr(Name, ~ class(get(.x))[1])) %>% 
  filter(Class == "function")  %>% 
  arrange(Name)

funs %>%
  kable(., format = "html", digits = 3) %>%
  kable_styling()


