#title: Figures for analysis of common nighthawk disturbance response
#author: Elly C. Knight

library(tidyverse)
library(gridExtra)
library(ggmap)
library(cowplot)
library(sf)
library(raster)
library(ggspatial)
library(ggsn)
library(magick)

my.theme <- theme_classic() +
  theme(text=element_text(size=12, family="Arial"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        plot.title=element_text(size=12, hjust = 0.5))

ggmap_rast <- function(map){
  map_bbox <- attr(map, 'bb') 
  .extent <- extent(as.numeric(map_bbox[c(2,4,1,3)]))
  my_map <- raster(.extent, nrow= nrow(map), ncol = ncol(map))
  rgb_cols <- setNames(as.data.frame(t(col2rgb(map))), c('red','green','blue'))
  red <- my_map
  values(red) <- rgb_cols[['red']]
  green <- my_map
  values(green) <- rgb_cols[['green']]
  blue <- my_map
  values(blue) <- rgb_cols[['blue']]
  stack(red,green,blue)
}

#FIGURE 1: STUDY AREA####

#Wrangle
dat <- read.csv("PDTGDataWrangled&Cleaned.csv")

sites <- dat %>% 
  dplyr::select(StationKey, DepYear, Latitude, Longitude, disturbance, time, pine, wetland, types, count) %>% 
  unique()

#Part A. Disturbance study sites----

#TO DO: NEEDS STUDY AREA BOUNDARY####

#A1. Maps----
#Get background data
register_google(key="AIzaSyCta9P4x7jGNELznpwlx07VZkkLVk3FP4M")

center <- dat %>% 
  summarize(long=mean(Longitude),
            lat=mean(Latitude))

map <- get_map(center, zoom=7, force=TRUE, maptype="satellite", color="color")

map_attributes <- attributes(map)

map_transparent <- matrix(adjustcolor(map, 
                                                alpha.f = 0.8), 
                                    nrow = nrow(map))
attributes(map_transparent) <- map_attributes

#trying to reproject imagery
map.rast <- ggmap_rast(map_transparent)
crs(map.rast) <- 4326
map.rast.10tm <- map.rast %>% 
  projectRaster(crs=3400) %>% 
  as.data.frame(xy=TRUE)
  
#Fire map
fire <- raster("/Volumes/ECK004/PDTG/Analysis/TIFs/fire_sa.tif") %>% 
  aggregate(fact=30) %>% 
  projectRaster(crs=4326) %>% 
  as.data.frame(xy=TRUE) %>% 
  dplyr::filter(!is.na(fire_sa))

sites.fire <- sites %>% 
  dplyr::filter(disturbance=="fire")

#Add study area perimeter
sa <- read_sf("/Volumes/ECK004/GIS/Projects/PDTG/PDTG_studyarea_V1.shp") %>% 
  st_make_valid() %>% 
  st_transform(crs=4326)

extent(fire)
sa.crop <- st_intersection(sa, st_set_crs(st_as_sf(as(raster::extent(-114.9958, -109.9933, 54.58757, 57.92461), "SpatialPolygons")), st_crs(sa)))

write_sf(sa.crop, "StudyAreaCropped.shp")

sa.crop <- readOGR("StudyAreaCropped.shp") %>% 
  fortify() %>% 
  mutate(ID = row_number())
plot(sa.crop)

ggplot() +
  geom_polygon(data=sa.crop, aes(long, lat), fill=NA, colour="black")

map.fire <- ggmap(map) +
    geom_raster(data=fire, aes(x=x, y=y, fill=fire_sa), na.rm=TRUE) +
    geom_polygon(data=sa.crop, aes(long, lat), fill=NA, colour="black") +
    geom_point(aes(x = Longitude, y = Latitude),
               shape = 21,
               colour="grey25",
               fill="grey75",
               data = sites.fire, 
               alpha = 1,
               size=3) +
  ggspatial::annotation_north_arrow(location = "br",
                                    style = ggspatial::north_arrow_orienteering(fill = c("grey75", "grey25"), line_col = "grey75", text_col="grey75"),
                                    height=unit(1, "cm"),
                                    width=unit(1, "cm")) +
  ggsn::scalebar(x.min = -115, x.max = -112.6, 
                 y.min = 54.6, y.max = 56, 
                 transform=TRUE, model="WGS84",
                 dist=50, dist_unit="km",
                 box.fill=c("grey75", "grey25"),
                 box.color="grey75",height=0.1,
                 st.bottom=FALSE, st.dist=0.08, st.size=3, st.color="grey75") +
  scale_fill_viridis_c(name="Year of\ndisturbance", option="A", limits=c(1935, 2015)) +  
  my.theme +
  coord_sf(datum = NA) +
  xlab("") +
  ylab("") +
  xlim(c(min(sites$Longitude)-0.1, max(sites$Longitude)+0.1)) +
  ylim(c(min(sites$Latitude)-0.1, max(sites$Latitude)+0.1)) +
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm"))
map.fire

#Harvest map
cc <- raster("/Volumes/ECK004/PDTG/Analysis/TIFs/harvest_sa.tif") %>% 
  aggregate(fact=30) %>% 
  projectRaster(crs=4326) %>% 
  as.data.frame(xy=TRUE) %>% 
  dplyr::filter(!is.na(harvest_sa))

sites.cc <- sites %>% 
  dplyr::filter(disturbance=="cc")

map.cc <- ggmap(map) +
  geom_raster(data=cc, aes(x=x, y=y, fill=harvest_sa), na.rm=TRUE) +
  geom_point(aes(x = Longitude, y = Latitude),
             shape = 21,
             colour="grey25",
             fill="grey75",
             data = sites.cc, 
             alpha = 1,
             size=3) +
  scale_fill_viridis_c(name="Year of\ndisturbance", option="A", limits=c(1935, 2015)) +   
  my.theme +
  coord_sf(datum = NA) +
  xlab("") +
  ylab("") +
  xlim(c(min(sites$Longitude)-0.1, max(sites$Longitude)+0.1)) +
  ylim(c(min(sites$Latitude)-0.1, max(sites$Latitude)+0.1)) +
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm"))

#Wellpad map
well <- raster("/Volumes/ECK004/PDTG/Analysis/TIFs/wellsmerge_sa.tif") %>% 
  aggregate(fact=30) %>% 
  projectRaster(crs=4326) %>% 
  as.data.frame(xy=TRUE) %>% 
  dplyr::filter(!is.na(wellsmerge_sa))

sites.well <- sites %>% 
  dplyr::filter(disturbance=="well")

map.well <- ggmap(map) +
  geom_raster(data=well, aes(x=x, y=y, fill=wellsmerge_sa), na.rm=TRUE) +
  geom_point(aes(x = Longitude, y = Latitude),
             shape = 21,
             colour="grey25",
             fill="grey75",
             data = sites.well, 
             alpha = 1,
             size=3) +
  scale_fill_viridis_c(name="Year of\ndisturbance", option="A", limits=c(1935, 2015)) +    
  my.theme +
  coord_sf(datum = NA) +
  xlab("") +
  ylab("") +
  xlim(c(min(sites$Longitude)-0.1, max(sites$Longitude)+0.1)) +
  ylim(c(min(sites$Latitude)-0.1, max(sites$Latitude)+0.1)) +
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm"))
#A2. Single disturbance histograms----

#Wrangle
site.1 <- sites %>% 
  dplyr::filter(types==1) %>% 
  group_by(StationKey, DepYear, Latitude, Longitude, disturbance, pine, wetland) %>% 
  summarize(time=min(time)) %>% 
  ungroup()

#Fire histogram
site.1.fire <- site.1 %>% 
  dplyr::filter(disturbance=="fire")
hist.fire <- ggplot(site.1.fire) +
  geom_histogram(aes(x=time)) +
  scale_x_continuous(limits=c(0,80), breaks=c(seq(0,80,20))) +   
  my.theme +
  theme(legend.position = "none",
        plot.margin = unit(c(0,1,0,-1), "cm")) +
  xlab("") +
  ylab("") +
  ylim(c(0,30))

#Harvest histogram
site.1.cc <- site.1 %>% 
  dplyr::filter(disturbance=="cc")
hist.cc <- ggplot(site.1.cc) +
  geom_histogram(aes(x=time)) +
  scale_x_continuous(limits=c(0,80), breaks=c(seq(0,80,20))) +   
  my.theme +
  theme(legend.position = "none",
        plot.margin = unit(c(0,1,0,-1), "cm")) +
  xlab("") +
  ylab("Number of study sites") +
  ylim(c(0,30))

#Wellpad histogram
site.1.well <- site.1 %>% 
  dplyr::filter(disturbance=="well")
hist.well <- ggplot(site.1.well) +
  geom_histogram(aes(x=time)) +
  scale_x_continuous(limits=c(0,80), breaks=c(seq(0,80,20))) +   
  my.theme +
  theme(legend.position = "none",
        plot.margin = unit(c(0,1,0,-1), "cm")) +
  xlab("Years since disturbance") +
  ylab("") +
  ylim(c(0,30))

#A3. Photos----
jpeg.fire <- image_read("/Users/ellyknight/Documents/UoA/Projects/Projects/PDTG/Analysis/postdisturbance_temporal_gradient/Figs/P6250957.jpg") %>% 
  image_fill('none') %>% 
  as.raster()
photo.fire <- ggplot() + 
  annotation_raster(jpeg.fire, 0, 1, 0, 1) +
  ylab("Fire") +
  my.theme +
  theme(axis.line.x.bottom = element_blank(),
        axis.line.y.right = element_blank(),
        axis.title.y.right = element_text(size=16),
        plot.margin = unit(c(-0.5,0,0.4,0), "cm")) +
  scale_y_continuous(position="right")
#photo.fire

jpeg.cc <- image_read("/Users/ellyknight/Documents/UoA/Projects/Projects/PDTG/Analysis/postdisturbance_temporal_gradient/Figs/P7221108.jpg") %>% 
  image_fill('none') %>% 
  as.raster()
photo.cc <- ggplot() + 
  annotation_raster(jpeg.cc, 0, 1, 0, 1) +
  ylab("Harvest") +
  my.theme +
  theme(axis.line.x.bottom = element_blank(),
        axis.line.y.right = element_blank(),
        axis.title.y.right = element_text(size=16),
        plot.margin = unit(c(-0.5,0,0.4,0), "cm")) +
  scale_y_continuous(position="right")

jpeg.well <- image_read("/Users/ellyknight/Documents/UoA/Projects/Projects/PDTG/Analysis/postdisturbance_temporal_gradient/Figs/stock-photo-88923883.jpg") %>% 
  image_fill('none') %>% 
  as.raster()
photo.well <- ggplot() + 
  annotation_raster(jpeg.well, 0, 1, 0, 1) +
  ylab("Well site") +
  my.theme +
  theme(axis.line.x.bottom = element_blank(),
        axis.line.y.right = element_blank(),
        axis.title.y.right = element_text(size=16),
        plot.margin = unit(c(-0.5,0,0.4,0), "cm")) +
  scale_y_continuous(position="right")

#A5. Make legend----
plot.legend <- ggmap(map) +
  geom_raster(data=cc, aes(x=x, y=y, fill=harvest_sa), na.rm=TRUE) +
  geom_point(aes(x = Longitude, y = Latitude),
             shape = 21,
             colour="grey25",
             fill="grey75",
             data = sites.cc, 
             alpha = 1,
             size=3) +
  scale_fill_viridis_c(name="Year of\ndisturbance", option="A", limits=c(1935, 2015)) +   
  my.theme +
  coord_sf(datum = NA) +
  xlab("") +
  ylab("") +
  xlim(c(min(sites$Longitude)-0.1, max(sites$Longitude)+0.1)) +
  ylim(c(min(sites$Latitude)-0.1, max(sites$Latitude)+0.1)) +
  theme(legend.position="right") +
  guides(fill = guide_colourbar(barwidth = 2.5, barheight = 15))
#plot.legend

legend <- get_legend(plot.legend)

#Part B. North America map####

nam <- map_data("world", region=c("Canada", 
                                  "USA", 
                                  "Mexico",
                                  "Guatemala", 
                                  "Belize", 
                                  "El Salvador",
                                  "Honduras", 
                                  "Nicaragua", 
                                  "Costa Rica",
                                  "Panama", 
                                  "Jamaica", 
                                  "Cuba", 
                                  "The Bahamas",
                                  "Haiti", 
                                  "Dominican Republic", 
                                  "Antigua and Barbuda",
                                  "Dominica", 
                                  "Saint Lucia", 
                                  "Saint Vincent and the Grenadines", 
                                  "Barbados",
                                  "Grenada",
                                  "Trinidad and Tobago")) %>% 
  dplyr::filter(!group%in%c(258:264))

nam.eq <- nam %>% 
  st_as_sf(coords=c("long", "lat"), crs=4326) %>% 
  st_transform(crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs") %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  cbind(nam)

sites.eq.center <- sites %>% 
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326) %>% 
  st_transform(crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs") %>%
  st_coordinates() %>% 
  as.data.frame() %>% 
  dplyr::summarize(X=mean(X),
            Y=mean(Y))

map.nam <- ggplot() +
  geom_polygon(data=nam.eq, aes(x=X, y=Y, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  geom_text(label="★", aes(x=X, y=Y), size=10, family = "HiraKakuPro-W3", data=sites.eq.center, colour="black") +
  xlim(c(-4000000, 3000000)) +
  coord_sf(datum = NA) +
  xlab("") +
  ylab("") +
  my.theme +
  theme(plot.margin = unit(c(0,0,0,-1), "cm"),
        legend.position="bottom")
#map.nam

#Part C: Put it together####

ggsave(plot=grid.arrange(photo.fire, photo.cc, photo.well,
                        map.fire, map.cc, map.well,
                         hist.fire, hist.cc, hist.well,
                         legend,
                         map.nam,
                         widths = c(3, 6, 4, 4),
                         heights = c(4, 2, 2, 4),
                         layout_matrix = rbind(c(11, 4, 7, 1),
                                               c(11, 5, 8, 2),
                                               c(10, 5, 8, 2),
                                               c(10, 6, 9, 3))),
       "Figs/Fig1StudyArea.jpeg", width=12, height=8, units="in", device="jpeg")

#FIGURE 2. DETECTABILITY & VEGETATION COVARIATES####

#Detectability predictions
pred.det <- read.csv("PDTGDetectionCovariatePredictions.csv") %>% 
  mutate(mod = "Detectability covariates") %>% 
  rename(fit = Det,
         upr = DetUpper,
         lwr = DetLower)
pred.det$response <- factor(pred.det$response, labels=c("Territorial habitat use", "Extraterritorial habitat use"))
pred.det$cov <- factor(pred.det$cov, levels=c("s2n2"), labels=c("StN (4.4-5.6 kHz)"))

#Vegetation predictions
pred.cov.boom <- read.csv("PDTGVegetationPredictionsBoom.csv") %>% 
  mutate(response="Territorial habitat use",
         cov="Proportion of pine forest") %>% 
  rename(x=pine,,
         fit = Occu,
         upr = OccuUpper,
         lwr = OccuLower)
pred.cov.call <- read.csv("PDTGVegetationPredictionsCall.csv") %>% 
  mutate(response="Extraterritorial habitat use",
         cov="Mean wetland probability") %>% 
  rename(x=wetland,
         fit = Occu,
         upr = OccuUpper,
         lwr = OccuLower)
pred.cov <- rbind(pred.cov.boom, pred.cov.call) %>% 
  dplyr::select(response, cov, x, fit, upr, lwr) %>% 
  mutate(mod = "Vegetation covariates")

pred <- rbind(pred.det, pred.cov)

plot.det.boom <- ggplot(subset(pred.det, response=="Territorial habitat use")) +
  geom_ribbon(aes(x=x, ymin=lwr, ymax=upr), alpha=0.4) +
  geom_line(aes(x=x, y=fit), colour="yellow3", size=1.2) +
  my.theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("StN (4.4-5.6 kHz)") +
  ylab("Probability of territorial detection")
#plot.det.boom

plot.det.call <- ggplot(subset(pred.det, response=="Extraterritorial habitat use")) +
  geom_ribbon(aes(x=x, ymin=lwr, ymax=upr), alpha=0.4) +
  geom_line(aes(x=x, y=fit), colour="yellow3", size=1.2) +
  my.theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("StN (4.4-5.6 kHz)") +
  ylab("Probability of extraterritorial detection")
#plot.det.call

plot.cov.boom <- ggplot(subset(pred.cov, response=="Territorial habitat use")) +
  geom_ribbon(aes(x=x, ymin=lwr, ymax=upr), alpha=0.4) +
  geom_line(aes(x=x, y=fit), colour="springgreen4", size=1.2) +
  my.theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("Proportion of pine forest") +
  ylab("Probability of territorial habitat use")
#plot.cov.boom

plot.cov.call <- ggplot(subset(pred.cov, response=="Extraterritorial habitat use")) +
  geom_ribbon(aes(x=x, ymin=lwr, ymax=upr), alpha=0.4) +
  geom_line(aes(x=x, y=fit), colour="steelblue3", size=1.2) +
  my.theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("Mean wetland probability") +
  ylab("Probability of extraterritorial habitat use")
#plot.cov.call

ggsave(plot=grid.arrange(plot.det.boom, plot.cov.boom, plot.det.call, plot.cov.call, ncol=2, nrow=2), 
       filename="Figs/Fig2Covariates.jpeg", width=10, height=10, units="in", device="jpeg")

#FIGURE 3. DISTURBANCE EFFECTS####

pred.1d.boom <- read.csv("PDTGSingleDisturbancePredictionsBoom.csv")
pred.1d.boom$disturbance <- factor(pred.1d.boom$disturbance, levels=c("fire", "cc", "well"), labels=c("Fire", "Harvest", "Wellpad"))

plot.1d.boom <- ggplot(pred.1d.boom) +
  geom_ribbon(aes(x=time, ymin=OccuLower, ymax=OccuUpper, group=factor(factor(pine))), colour="grey70", alpha=0.2) +
  geom_line(aes(x=time, y=Occu, colour=factor(pine)))+
  facet_wrap(~disturbance) +
  labs(x="Years since disturbance", y="Probability of territorial habitat use", colour="Proportion of pine forest")+
  xlim(c(0,80)) +
  ylim(c(0,1)) +
  my.theme +
  scale_colour_manual(values=c("darkgoldenrod1", "chartreuse4")) +
  theme(legend.position="bottom")
plot.1d.boom

ggsave(plot=plot.1d.boom, 
       filename="Figs/Fig5SingleDisturbance.jpeg", width=8, height=4, units="in", device="jpeg")

#SUMMARY STATS####

dat <- read.csv("PDTGDataWrangled&Cleaned.csv")
val <- read.csv("PDTGDataWrangledValidation.csv")

#Number of ARU recordings
nrow(dat)

#Number of ARU deployments
dat %>% 
  dplyr::select(ID) %>% 
  unique() %>% 
  nrow()

#Number of study sites
dat %>% 
  dplyr::select(StationKey) %>% 
  unique() %>% 
  nrow()

#Recordings per study site
dat %>% 
  dplyr::select(ID, fileID) %>% 
  unique() %>% 
  group_by(ID) %>% 
  summarize(n=n()) %>% 
  ungroup() %>% 
  summary()

dat %>% 
  dplyr::select(ID, fileID) %>% 
  unique() %>% 
  group_by(ID) %>% 
  summarize(n=n()) %>% 
  ungroup() %>% 
  summarize(sd=sd(n))

#Sites with CONI detections
site.coni <- dat %>% 
  group_by(StationKey) %>% 
  summarize(boom=sum(siteboom),
            call=sum(sitecall)) %>% 
  ungroup() %>% 
  mutate(boom = ifelse(boom > 0, 1, 0),
         call = ifelse(call > 0, 1, 0))
table(site.coni$boom, site.coni$call)

id.coni <- dat %>% 
  group_by(ID) %>% 
  summarize(boom=sum(siteboom),
            call=sum(sitecall)) %>% 
  ungroup() %>% 
  mutate(boom = ifelse(boom > 0, 1, 0),
         call = ifelse(call > 0, 1, 0))
table(id.coni$boom, id.coni$call)

rec.coni <- dat %>% 
  group_by(fileID) %>% 
  summarize(boom=sum(siteboom),
            call=sum(sitecall)) %>% 
  ungroup() %>% 
  mutate(boom = ifelse(boom > 0, 1, 0),
         call = ifelse(call > 0, 1, 0))
table(rec.coni$boom, rec.coni$call)

#Recognizer stats
nrow(val)
table(val$validation)
nrow(val)