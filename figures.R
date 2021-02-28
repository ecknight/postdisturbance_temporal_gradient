library(tidyverse)
library(gridExtra)
library(ggmap)
library(cowplot)
library(sf)
library(raster)

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

#FIGURE 1: STUDY AREA####

#Wrangle
setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/PDTG/Analysis")
dat <- read.csv("PDTGDataWrangled&Cleaned.csv")

sites <- dat %>% 
  dplyr::select(StationKey, DepYear, Latitude, Longitude, disturbance, time, pine, wetland, types, count) %>% 
  unique()

#Part A. Disturbance study sites----

#A1. Maps----
#Get background data
register_google(key="AIzaSyCta9P4x7jGNELznpwlx07VZkkLVk3FP4M")

center <- dat %>% 
  summarize(long=mean(Longitude),
            lat=mean(Latitude))

map <- get_map(center, zoom=7, force=TRUE, maptype="satellite", color="bw")

map_attributes <- attributes(map)

map_transparent <- matrix(adjustcolor(map, 
                                                alpha.f = 0.8), 
                                    nrow = nrow(map))
attributes(map_transparent) <- map_attributes

#Fire map
sites.fire <- sites %>% 
  dplyr::filter(disturbance=="fire")
map.fire <- ggmap(map_transparent) +
  geom_point(aes(x = Longitude, y = Latitude,
                 fill=time),
             shape = 21,
             colour="grey75",
             data = sites.fire, 
             alpha = 1,
             size=3) +
  scale_fill_viridis_c(name="Years since\ndisturbance", direction=-1) + 
  ggtitle("Fire") +
  my.theme +
  coord_sf(datum = NA) +
  xlab("") +
  ylab("") +
  xlim(c(min(sites$Longitude)-0.1, max(sites$Longitude)+0.1)) +
  ylim(c(min(sites$Latitude)-0.1, max(sites$Latitude)+0.1)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        legend.position = "none")

#Harvest map
sites.cc <- sites %>% 
  dplyr::filter(disturbance=="cc")
map.cc <- ggmap(map_transparent) +
  geom_point(aes(x = Longitude, y = Latitude,
                 fill=time),
             shape = 21,
             colour="grey75",
             data = sites.cc, 
             alpha = 1,
             size=3) +
  scale_fill_viridis_c(name="Years since\ndisturbance",
                       limits=c(0,80), direction=-1) +   
  ggtitle("Harvest") +
  my.theme +
  coord_sf(datum = NA) +
  xlab("") +
  ylab("") +
  xlim(c(min(sites$Longitude)-0.1, max(sites$Longitude)+0.1)) +
  ylim(c(min(sites$Latitude)-0.1, max(sites$Latitude)+0.1)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        legend.position = "none")

#Wellpad map
sites.well <- sites %>% 
  dplyr::filter(disturbance=="well")
map.well <- ggmap(map_transparent) +
  geom_point(aes(x = Longitude, y = Latitude,
                 fill=time),
             shape = 21,
             colour="grey75",
             data = sites.well, 
             alpha = 1,
             size=3) +
  scale_fill_viridis_c(name="Years since\ndisturbance",
                       limits=c(0,80), direction=-1) +   
  ggtitle("Wellpad") +
  my.theme +
  coord_sf(datum = NA) +
  xlab("") +
  ylab("") +
  xlim(c(min(sites$Longitude)-0.1, max(sites$Longitude)+0.1)) +
  ylim(c(min(sites$Latitude)-0.1, max(sites$Latitude)+0.1)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        legend.position = "none")

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
  geom_histogram(aes(x=time, fill=factor(time))) +
  scale_fill_viridis_d(name="Years since\ndisturbance",
                       breaks=c(seq(0,80,5)), direction=-1) +   
  my.theme +
  theme(legend.position = "none") +
  xlab("Years since disturbance") +
  ylab("Number of study sites") +
  ylim(c(0,30))

#Harvest histogram
site.1.cc <- site.1 %>% 
  dplyr::filter(disturbance=="cc")
hist.cc <- ggplot(site.1.cc) +
  geom_histogram(aes(x=time, fill=factor(time))) +
  scale_fill_viridis_d(name="Years since\ndisturbance",
                       breaks=c(seq(0,80,5)), direction=-1) +   
  my.theme +
  theme(legend.position = "none") +
  xlab("Years since disturbance") +
  ylab("Number of study sites") +
  ylim(c(0,30))

#Wellpad histogram
site.1.well <- site.1 %>% 
  dplyr::filter(disturbance=="well")
hist.well <- ggplot(site.1.well) +
  geom_histogram(aes(x=time, fill=factor(time))) +
  scale_fill_viridis_d(name="Years since\ndisturbance",
                       breaks=c(seq(0,80,5)), direction=-1) +   
  my.theme +
  theme(legend.position = "none") +
  xlab("Years since disturbance") +
  ylab("Number of study sites") +
  ylim(c(0,30))

#A3. Dual disturbance plots----

#Wrangle
sites.2.min <- sites %>% 
  dplyr::filter(types==2) %>% 
  group_by(StationKey, DepYear, Latitude, Longitude, disturbance, pine, wetland) %>% 
  summarize(time=min(time)) %>% 
  ungroup()

sites.2 <- sites.2.min %>% 
  mutate(ID = row_number()) %>% 
  pivot_wider(id_cols=c(StationKey, DepYear, pine, wetland),
              names_from=disturbance,
              names_prefix="time.",
              values_from=c(time))  %>% 
  rowwise(StationKey, DepYear, pine, wetland) %>% 
  mutate(time.min = min(time.cc, time.well, time.fire, na.rm=TRUE)) %>% 
  ungroup() %>% 
  left_join(sites.2.min %>% 
              dplyr::rename(time.min=time))

#CCFR plot
sites.2.frcc <- sites.2 %>% 
  dplyr::filter(is.na(time.well))

plot.frcc <- ggplot(sites.2.frcc) +
  geom_point(aes(x=time.fire, y=time.cc)) +
  my.theme + 
  xlab("Years since fire") +
  ylab("Years since harvest")


#CCFR plot
sites.2.ccwell <- sites.2 %>% 
  dplyr::filter(is.na(time.fire))

plot.ccwell <- ggplot(sites.2.ccwell) +
  geom_point(aes(x=time.cc, y=time.well)) +
  my.theme + 
  xlab("Years since harvest") +
  ylab("Years since wellpad clearing")

#WellFR plot
#CCFR plot
sites.2.wellfr <- sites.2 %>% 
  dplyr::filter(is.na(time.cc))

plot.wellfr <- ggplot(sites.2.wellfr) +
  geom_point(aes(x=time.well, y=time.fire)) +
  my.theme + 
  xlab("Years since wellpad clearing") +
  ylab("Years since fire")

#A5. Make legend----
plot.legend <- ggmap(map_transparent) +
  geom_point(aes(x = Longitude, y = Latitude,
                 fill=time),
             shape = 21,
             colour="grey75",
             data = sites.well, 
             alpha = 1,
             size=3) +
  scale_fill_viridis_c(name="Years since\ndisturbance",
                       limits=c(0,80), direction=-1) +
  my.theme +
  theme(legend.position="right") +
  guides(fill = guide_colourbar(barwidth = 2.5, barheight = 20))
#plot.legend

legend <- get_legend(plot.legend)

#Part B. Vegetation Covariates####

#B1. Maps----

sa <- read_sf("/Volumes/ECK004/GIS/Projects/PDTG/PDTG_studyarea_V1.shp")

#Pine
#pine <- raster("/Users/ellyknight/Documents/UoA/Projects/Projects/PDTG/Analysis/pine-200.tif") %>% 
#  aggregate(factor=30) %>% 
#  as.data.frame(xy=TRUE)

map.pine <- ggmap(map_transparent) +
#  geom_raster(data=pine, aes(x=x, y=y, fill=pine.200)) +
#  geom_sf(data=sa, fill=NA, colour="black") +
  geom_point(aes(x = Longitude, y = Latitude,
                 fill=pine*100),
             shape = 21,
             colour="grey75",
             data = sites, 
             alpha = 1,
             size=3) +
#  scale_fill_gradient(name="% pine forest\nwithin 200m",
#                      low="black", high="green") +
  scale_fill_viridis_c(name="% pine forest\nwithin 200m",
                       option="cividis",
                       direction=-1,
                       limits=c(0,100)) +
#  ggtitle("Fire") +
  my.theme +
  coord_sf(datum = NA) +
  xlab("") +
  ylab("") +
  xlim(c(min(sites$Longitude)-0.1, max(sites$Longitude)+0.1)) +
  ylim(c(min(sites$Latitude)-0.1, max(sites$Latitude)+0.1)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

#Wetland
#wetland <- raster("/Users/ellyknight/Documents/UoA/Projects/Projects/PDTG/Analysis/wetland-200.tif")

map.wetland <- ggmap(map_transparent) +
  geom_point(aes(x = Longitude, y = Latitude,
                 fill=wetland),
             shape = 21,
             colour="grey75",
             data = sites, 
             alpha = 1,
             size=3) +
  scale_fill_viridis_c(name="Mean wetland\nprobability\nwithin 200m",
                      option="plasma",
                      direction=-1,
                      limits=c(0,1)) +
  #  ggtitle("Fire") +
  my.theme +
  coord_sf(datum = NA) +
  xlab("") +
  ylab("") +
  xlim(c(min(sites$Longitude)-0.1, max(sites$Longitude)+0.1)) +
  ylim(c(min(sites$Latitude)-0.1, max(sites$Latitude)+0.1)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

#B2. Histograms----

#Pine histogram
hist.pine <- ggplot(sites) +
  geom_histogram(aes(x=pine*100, fill=factor(pine*100))) +
  scale_fill_viridis_d(name="% pine forest\nwithin 200m",
                       option="cividis",
                       direction=-1,
                       breaks=c(seq(0,100,5))) + 
  my.theme +
  theme(legend.position = "none") +
  xlab("% pine forest within 200m") +
  ylab("Number of study sites")

#Wetland histogram
hist.wetland <- ggplot(sites) +
  geom_histogram(aes(x=wetland, fill=factor(wetland))) +
  scale_fill_viridis_d(name="Mean wetland\nprobability\nwithin 200m",
                       option="plasma",
                       direction=-1,
                       breaks=c(seq(0,1,0.05))) +
  my.theme +
  theme(legend.position = "none") +
  xlab("Mean wetland probability within 200m") +
  ylab("Number of study sites")

#Part C. Study area map####

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
  st_transform(crs=102008) %>%
  st_coordinates() %>% 
  data.frame() %>% 
  cbind(nam)

sites.eq.center <- sites %>% 
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326) %>% 
  st_transform(crs=102008) %>%
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
  theme(plot.margin = unit(c(0,0,0,-1.5), "cm"),
        legend.position="bottom")

#Part D: Put it together####

ggsave(plot=grid.arrange(map.fire, map.cc, map.well,
                         hist.fire, hist.cc, hist.well,
                         plot.frcc, plot.ccwell, plot.wellfr,
                         legend,
                         map.nam,
                         widths = c(6, 3, 3, 3, 2),
                         heights = c(3, 3, 3),
                         layout_matrix = rbind(c(11, 1, 2, 3, 10),
                                               c(11, 4, 5, 6, 10),
                                               c(11, 7, 8, 9, 10))),
       "Figs/Fig1StudyArea.jpeg", width=15, height=8, units="in", device="jpeg")

#FIGURE 2. CONI SPECTROGRAM####
wav <- readWave("/Users/ellyknight/Documents/UoA/Presentations/Materials/Recordings/CONI Common nighthawk.wav")
wav.f <- bwfilter(wav, from=200, to=5400, output="Wave")
spectro <- ggspectro(wav.f, flim=c(0,8), wl=1024, ovlp=60) +
  stat_contour(geom="polygon", aes(fill=..level..), bins=30, show.legend=FALSE) +
  scale_fill_continuous(name="Amplitude\n(dB)\n", limits=c(-40,0),
                        na.value="transparent", low="white", high="black") + 
  xlab("Time (s)") +
  ylab("Frequency (kHz)") +
  my.theme
spectro

ggplot2::ggsave(plot=spectro, "/Users/ellyknight/Documents/UoA/Presentations/Materials/Recordings/CONI Common nighthawk - boom.tiff", width=6, height=8, units="in", device="tiff")

#FIGURE 3. DETECTABILITY COVARIATES####
pred.det <- read.csv("PDTGDetectionCovariatePredictions.csv") 
pred.det$response <- factor(pred.det$response, labels=c("Territorial habitat use", "Home range habitat use"))
pred.det$cov <- factor(pred.det$cov, levels=c("set", "doy", "s2n1", "s2n2"), labels=c("Hours since sunset", "Day of year", "StN (0.6-1.2 kHz)", "StN (4.4-5.6 kHz)"))

plot.det <- ggplot(pred.det) +
  geom_ribbon(aes(x=x, ymin=DetLower, ymax=DetUpper), alpha=0.4) +
  geom_line(aes(x=x, y=Det, colour=cov)) +
  scale_colour_viridis_d(option="plasma") +
  facet_grid(response~cov, scales="free") +
  my.theme +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("") +
  ylab("Probability of detection")

ggsave(plot=plot.det, filename="Figs/Fig3Detectability.jpeg", width=8, height=4, units="in", device="jpeg")

#FIGURE 4. VEGETATION COVARIATES####
pred.cov.boom <- read.csv("PDTGVegetationPredictionsBoom.csv")
pred.cov.call <- read.csv("PDTGVegetationPredictionsCall.csv")

plot.cov.boom <- ggplot(pred.cov.boom) +
  geom_line(aes(x=pine, y=Occu))+
  geom_ribbon(aes(x=pine, ymin=OccuLower, ymax=OccuUpper), colour="grey70", alpha=0.2) +
  my.theme +
  ylim(c(0,1)) + 
  xlab("Proportion of pine forest") +
  ylab("Probability of territorial habitat use")

plot.cov.call <- ggplot(pred.cov.call) +
  geom_line(aes(x=wetland, y=Occu))+
  geom_ribbon(aes(x=wetland, ymin=OccuLower, ymax=OccuUpper), colour="grey70", alpha=0.2) +
  my.theme +
  ylim(c(0,1)) +
  xlab("Mean wetland probability") +
  ylab("Probability of home range habitat use")

ggsave(plot=grid.arrange(plot.cov.boom, plot.cov.call, ncol=2), 
       filename="Figs/Fig4Vegetation.jpeg", width=8, height=4, units="in", device="jpeg")

#FIGURE 5. SINGLE DISTURBANCES####

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

#Number of study sites
dat %>% 
  dplyr::select(ID) %>% 
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
table(val$validation)
nrow(val)
309