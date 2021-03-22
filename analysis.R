#load("/Users/ellyknight/Documents/UoA/Projects/Projects/PDTG/Analysis/PDTGV16.Rdata")
#save.image("/Users/ellyknight/Documents/UoA/Projects/Projects/PDTG/Analysis/PDTGV16.Rdata")

options(scipen = 999999)

library(tidyverse)
library(raster)
library(lubridate)
library(maptools)
library(RPresence)
library(unmarked)
library(MuMIn)
library(lme4)
library(usdm)
library(GGally)
library(AICcmodavg)
library(GGally)
library(seewave)
library(sf)
library(lwgeom)
library(dggridR)
library(suncalc)
library(gridExtra)
library(data.table)

my.theme <- theme_classic() +
  theme(text=element_text(size=16, family="Arial"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        plot.title=element_text(size=14))


#TO DO: FIX MODEL NAMES FROM COV TO MOD####
#TO DO: FIX OCCUPANCY DETECTIONS FORMATTING####

#A. READ IN ALL DATA####
#1. Read in study site data----
site.1 <- read.csv("PDTG_RecordingSelection_SingleYearSingleDisturbance_V5.csv") %>% 
  dplyr::select(ProjectID, DepYear, StationKey, Latitude, Longitude) %>% 
  unique()

site.2<- read.csv("PDTG_RecordingSelection_SingleYearMultiDisturbance_V5.csv") %>% 
  dplyr::select(ProjectID, DepYear, StationKey, Latitude, Longitude) %>% 
  unique()

site.3 <- read.csv("PDTG_SiteSelection_MultiYear_V4.csv") %>% 
  dplyr::select(ProjectID, DepYear, StationKey, Latitude, Longitude) %>% 
  unique()

  
site <- rbind(site.1, site.2, site.3) %>% 
  unique()

#2. Read in recording lists----
rec.1 <- read.csv("PDTG_RecordingSelection_SingleYearSingleDisturbance_V5.csv") %>% 
  mutate(fileID = str_sub(filename, -100, -5),
         fileID = str_replace_all(fileID, "[$]", "_"),
         fileID = str_replace_all(fileID, "[+]", "-")) %>% 
  dplyr::select(StationKey, DepYear, fileID, yday, sundiff)

rec.2 <- read.csv("PDTG_RecordingSelection_SingleYearMultiDisturbance_V5.csv") %>% 
  mutate(fileID = str_sub(filename, -100, -5),
         fileID = str_replace_all(fileID, "[$]", "_"),
         fileID = str_replace_all(fileID, "[+]", "-")) %>% 
  dplyr::select(StationKey, DepYear, fileID, yday, sundiff)

rec.3 <- read.csv("PDTG_RecordingSelection_MultiYear_V6.csv") %>% 
  mutate(fileID = str_sub(filename, -100, -5),
         fileID = str_replace_all(fileID, "[$]", "_"),
         fileID = str_replace_all(fileID, "[+]", "-")) %>% 
  dplyr::select(StationKey, DepYear, fileID, yday, sundiff)
  
rec <- rbind(rec.1, rec.2, rec.3) %>% 
  dplyr::filter(sundiff >= -0.5,
                sundiff <= 1.5)
summary(rec)

#3. Read in & wrangle validated recognizer results----
val1 <- read.table("PDTG_y1d1_1_CONI0mv2_20_60_results_validated.txt",
                   sep="\t",
                   col.names = c("filepath", "start", "duration", "rsl", "quality", "score", "recognizer", "validation"))
val2 <- read.table("PDTG_y1d1_2_CONI0mv2_20_60_results_validated.txt",
                   sep="\t",
                   col.names = c("filepath", "start", "duration", "rsl", "quality", "score", "recognizer", "validation"))
val3 <- read.table("PDTG_y1d1_3_CONI0mv2_20_60_results_validated.txt",
                   sep="\t",
                   col.names = c("filepath", "start", "duration", "rsl", "quality", "score", "recognizer", "validation"))
val4 <- read.table("PDTG_y1d1_4_CONI0mv2_20_60_results_validated.txt",
                   sep="\t",
                   col.names = c("filepath", "start", "duration", "rsl", "quality", "score", "recognizer", "validation"))
val5 <- read.table("PDTG_y1d2_1_CONI0mv2_20_60_results_validated.txt",
                   sep="\t",
                   col.names = c("filepath", "start", "duration", "rsl", "quality", "score", "recognizer", "validation"))
val6 <- read.table("PDTG_y1d2_2_CONI0mv2_20_60_results_validated.txt",
                   sep="\t",
                   col.names = c("filepath", "start", "duration", "rsl", "quality", "score", "recognizer", "validation"))
val7 <- read.table("PDTG_y1d2_3_CONI0mv2_20_60_results_validated.txt",
                   sep="\t",
                   col.names = c("filepath", "start", "duration", "rsl", "quality", "score", "recognizer", "validation"))
val8 <- read.table("PDTG_y2b_3min_CONI0mv2_20_60_results_validated.txt",
                   sep="\t",
                   col.names = c("filepath", "start", "duration", "rsl", "quality", "score", "recognizer", "validation"))

val <- rbind(val1, val2, val3, val4, val5, val6, val7, val8) %>% 
  separate(filepath, into=c("f1", "f2", "f3", "f4", "treatment", "f6", "file"), sep="/", remove=FALSE) %>%
  dplyr::select(file, start, score, rsl, validation) %>% 
  mutate(validation = case_when(validation=="" ~ "n",
                                validation=="b" ~ "b",
                                validation=="n" ~ "n",
                                validation=="y" ~ "y",
                                validation=="y - missed boom" ~ "y",
                                validation=="yb" ~ "b",
                                validation=="yy" ~ "y",
                                validation=="nn" ~ "n")) %>% 
  mutate(start=as.character(start),
         start=ifelse(nchar(start)<10, paste0("00:", start), start)) %>% 
  separate(start, into=c("starthr", "startmin", "startsec"), sep=":", remove=FALSE) %>% 
  mutate(startmin=as.numeric(startmin)) %>% 
  dplyr::filter(as.numeric(startmin) < 3) #Shorten to 3 minute recordings

write.csv(val, "PDTGDataWrangledValidation.csv", row.names = FALSE)

#4. Summarize validation by recording----
val.rec <- val %>% 
  group_by(file, validation) %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  dplyr::filter(validation !="n") %>% 
  spread(key = validation, value = n) %>% 
  rename(boom = b, call = y) %>% 
  mutate(boom = ifelse(is.na(boom), 0, boom),
         call = ifelse(is.na(call), 0, call)) %>% 
  mutate(fileID = str_sub(file, -100, -5),
         fileID = str_replace_all(fileID, "[$]", "_"),
         fileID = str_replace_all(fileID, "[+]", "-"),
         fileID = str_replace_all(fileID, "_0-1_", "_"),
         fileID = str_replace_all(fileID, "[.]", "-"),
         fileID = str_replace_all(fileID, "-wav", ".wav")) %>% 
  dplyr::select(-file) %>% 
  right_join(rec) %>% 
  mutate(boom=ifelse(is.na(boom), 0, boom),
         call=ifelse(is.na(call), 0, call),
         boompres=ifelse(boom > 0, 1, 0),
         callpres=ifelse(call > 0, 1, 0))

site.val.rec <- val.rec %>% 
  dplyr::select(StationKey, DepYear) %>% 
  unique() %>% 
  inner_join(site) %>% 
  dplyr::filter(!is.na(Latitude))

#5. Read in and wrangle hard rain data-----
hr1 <- read.csv("HardRainResults_y1d1_1.csv")
hr2 <- read.csv("HardRainResults_y1d1_2.csv")
hr3 <- read.csv("HardRainResults_y1d1_3.csv")
hr4 <- read.csv("HardRainResults_y1d1_4.csv")
hr5 <- read.csv("HardRainResults_y1d2_1.csv")
hr6 <- read.csv("HardRainResults_y1d2_2.csv")
hr7 <- read.csv("HardRainResults_y1d2_3.csv")
hr8 <- read.csv("HardRainResults_y2b_1.csv")
hr9 <- read.csv("HardRainResults_y2b_2.csv")
hr10 <- read.csv("HardRainResults_y2b_3.csv")

hr <- rbind(hr1, hr2, hr3, hr4, hr5, hr6, hr7, hr8, hr9, hr10) %>% 
  mutate(segment = ifelse(is.na(segment), 0, as.numeric(segment))) %>% 
  dplyr::filter(segment < 6) %>% 
  group_by(file) %>% 
  summarize(psd.1 = mean(band.1.psd),
            psd.2 = mean(band.2.psd),
            s2n.1 = mean(band.1.s2n),
            s2n.2 = mean(band.2.s2n)) %>% 
  ungroup() %>% 
  rename(filename = file) %>% 
  mutate(filename = str_replace_all(filename, "[.]", "-"),
         filename = str_replace_all(filename, "-wav", ".wav")) %>% 
  mutate(fileID = str_sub(filename, -100, -9)) %>% 
  dplyr::select(-filename)

#B. EXTRACT SPATIAL COVARIATES####

#See what the sampling looks like with all the sampled points

#1. Buffer----
site.sf <- site.val.rec %>% 
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326) %>% 
  st_transform(crs=" +proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m
+no_defs") %>% 
  mutate(Latitude = site.val.rec$Latitude,
         Longitude = site.val.rec$Longitude) %>% 
  unique()

site.buff <- site.sf %>% 
  st_buffer(dist=200)

#2. Intersect resampled disturbance polygons with original polygon to get FID----
#fire.sa <- read_sf("/Volumes/ECK004/GIS/Projects/PDTG/fire_sa.shp") %>% 
#  st_make_valid() %>% 
#  st_set_crs(3400) %>% 
#  st_intersection(st_make_valid(read_sf("/Volumes/ECK004/GIS/Projects/PDTG/fire_sa_year.shp")))
#write_sf(fire.sa, "/Volumes/ECK004/GIS/Projects/PDTG/fire_sa_year_area.shp")

fire.sa <- read_sf("/Volumes/ECK004/GIS/Projects/PDTG/fire_sa_year_area.shp") %>% 
  st_transform(crs=crs(site.buff))

#cc.sa <- read_sf("/Volumes/ECK004/GIS/Projects/PDTG/harvest_sa.shp") %>% 
#  st_make_valid() %>% 
#  st_set_crs(3400) %>% 
#  st_intersection(st_make_valid(read_sf("/Volumes/ECK004/GIS/Projects/PDTG/harvest_sa_year.shp")))
#write_sf(cc.sa, "/Volumes/ECK004/GIS/Projects/PDTG/harvest_sa_year_area.shp")

cc.sa <- read_sf("/Volumes/ECK004/GIS/Projects/PDTG/harvest_sa_year_area.shp") %>% 
  st_transform(crs=crs(site.buff))

#wells.sa <- read_sf("/Volumes/ECK004/GIS/Projects/PDTG/wells_sa.shp") %>% 
#  st_make_valid() %>% 
#  st_set_crs(3400) %>% 
#  st_intersection(st_make_valid(read_sf("/Volumes/ECK004/GIS/Projects/PDTG/wells_sa_year.shp")))
#write_sf(wells.sa, "/Volumes/ECK004/GIS/Projects/PDTG/well_sa_year_area.shp")

wells.sa <- read_sf("/Volumes/ECK004/GIS/Projects/PDTG/well_sa_year_area.shp") %>% 
  st_transform(crs=crs(site.buff))

#3. Intersect sites with disturbance polygons----
site.fire <- data.frame(sitearea= as.numeric(st_area(st_intersection(fire.sa, site.buff)))) %>% 
  cbind(data.frame(st_intersection(fire.sa, site.buff))) %>% 
  rename(year = YEAR, class = BURN_CLASS, distarea=HECTARES_U) %>% 
  dplyr::select(StationKey, year, sitearea, distarea) %>% 
  mutate(distarea = distarea*10000,
         disturbance="fire")

site.cc <- data.frame(sitearea= as.numeric(st_area(st_intersection(cc.sa, site.buff)))) %>% 
  cbind(data.frame(st_intersection(cc.sa, site.buff))) %>% 
  rename(year = YEAR, distarea = Shape_Area) %>% 
  dplyr::select(StationKey, year, sitearea, distarea) %>% 
  mutate(disturbance="cc")

site.well <- data.frame(sitearea= as.numeric(st_area(st_intersection(wells.sa, site.buff)))) %>% 
  cbind(data.frame(st_intersection(wells.sa, site.buff))) %>% 
  rename(year = YEAR, distarea = Shape_Area) %>% 
  dplyr::select(StationKey, year, sitearea, distarea) %>% 
  mutate(disturbance="well")

#4. Extract pine and wetland values----
pine <- raster("pine-200.tif")
wetland <- raster("wetland-200.tif")

covs.sf <- data.frame(pine = raster::extract(x=pine, y=site.sf),
                      wetland = raster::extract(x=wetland, y=site.sf)) %>% 
  bind_cols(site.sf) %>% 
  dplyr::filter(!is.na(pine))

#5. Join with station key list----

#Sum disturbances of same year and add disturbance presence columns
site.covs <- rbind(site.fire, site.cc, site.well) %>% 
  left_join(covs.sf) %>% 
  group_by(ProjectID, StationKey, DepYear, Latitude, Longitude, disturbance, year, pine, wetland) %>% 
  summarize(distarea=sum(distarea),
            sitearea=sum(sitearea)) %>% 
  ungroup()

#6. Set up grid sampling----
grid <- dgconstruct(spacing=1, metric=TRUE)

site.grid <- site.sf %>% 
  mutate(cell = dgGEO_to_SEQNUM(grid, Longitude, Latitude)$seqnum) %>% 
  arrange(cell) %>% 
  right_join(site.covs) %>% 
  data.frame() %>% 
  dplyr::select(-geometry)  %>% 
  mutate(time=DepYear-year,
         timeclass=round(time/5)*5) %>% 
  dplyr::filter(time > 0)

#7. Count # of disturbances at each site----
site.dist <- site.grid %>% 
  group_by(ProjectID, StationKey, DepYear, Latitude, Longitude, disturbance, pine, wetland, cell) %>% 
  summarize(count=n()) %>% 
  group_by(ProjectID, StationKey, DepYear, Latitude, Longitude, pine, wetland, cell) %>% 
  summarize(types=n(),
            count=sum(count)) %>% 
  ungroup() %>% 
  right_join(site.grid)

#8. Put together with validated results----
dat.rec <- site.dist %>% 
  left_join(val.rec) %>% 
  left_join(hr) %>% 
  mutate(ID = paste0(StationKey, DepYear)) %>% 
  filter(yday < 213,
         !is.na(s2n.2),
         !is.na(call)) %>% #filter out days after Aug 1 
  group_by(StationKey, DepYear) %>% 
  mutate(n = row_number()) %>% 
  ungroup() %>% 
  dplyr::select(-boom, -call) %>% 
  dplyr::rename(boom=boompres, call=callpres)

dat <- dat.rec %>% 
  group_by(StationKey, DepYear) %>% 
  summarize(siteboom=ifelse(sum(boom)>0, 1, 0),
            sitecall=ifelse(sum(call)>0, 1, 0)) %>% 
  ungroup() %>% 
  inner_join(dat.rec) %>% 
  mutate(call=ifelse(siteboom==1, 0, call),
         ID = paste(StationKey, DepYear))

summary(dat)

#9. Final site list----
set.seed(1234)
site.dat <- dat %>% 
  dplyr::select(StationKey, DepYear) %>% 
  unique() %>% 
  sample_n(800) %>% 
  left_join(dat) %>% 
  dplyr::select(ProjectID, StationKey, DepYear, Latitude, Longitude, disturbance, time, sitearea, pine, wetland, cell, types, count) %>% 
  unique()

dat <- dat %>% 
  inner_join(site.dat)

write.csv(dat, "PDTGDataWrangled&Cleaned.csv", row.names = FALSE)

#10. Visualize availability----
ggplot() +
  geom_point(data=subset(dat, boom==0), aes(x=Longitude, y=Latitude), colour="black") +
  geom_point(data=subset(dat, boom==1), aes(x=Longitude, y=Latitude), colour="red") +
#  geom_hex(data=dat, aes(x=Longitude, y=Latitude)) +
  facet_wrap(~disturbance)

ggplot() +
  geom_point(data=subset(dat, call==0), aes(x=Longitude, y=Latitude), colour="black") +
  geom_point(data=subset(dat, call==1), aes(x=Longitude, y=Latitude), colour="red") +
  #  geom_hex(data=dat, aes(x=Longitude, y=Latitude)) +
  facet_grid( ~ types)

#C. DETECTION COVARIATES####
#1. Check for vif & covariation----
vif(dat %>% 
      dplyr::select(sundiff, yday, psd.1, psd.2, s2n.2, s2n.1) %>% 
      data.frame())

cor(dat %>% 
      dplyr::select(sundiff, yday, psd.1, psd.2, s2n.2, s2n.1) %>% 
      data.frame())
#Take out psd.2

vif(dat %>% 
      dplyr::select(sundiff, yday, psd.1, s2n.2, s2n.1) %>% 
      data.frame())

cor(dat %>% 
      dplyr::select(sundiff, yday, psd.1, s2n.2, s2n.1) %>% 
      data.frame())
#all good

ggpairs(dat %>% 
          dplyr::select(sundiff, yday, psd.2, s2n.2, s2n.1) %>% 
          data.frame())

#2. Format into unmarked objects----
hist.boom <- dat %>% 
  dplyr::select(ID, n, boom) %>% 
  spread(key = n, value = boom) %>% 
  arrange(ID) %>% 
  dplyr::select(-ID) %>% 
  data.frame()

hist.call <- dat%>% 
  dplyr::select(ID, n, call) %>% 
  spread(key = n, value = call) %>% 
  arrange(ID) %>% 
  dplyr::select(-ID) %>% 
  data.frame()

hist.set <- dat %>% 
  dplyr::select(ID, n, sundiff) %>% 
  spread(key = n, value = sundiff) %>% 
  arrange(ID) %>% 
  dplyr::select(-ID) %>% 
  data.frame() 

hist.Doy <- dat %>% 
  dplyr::select(ID, n, yday) %>% 
  spread(key=n, value=yday) %>% 
  arrange(ID) %>% 
  dplyr::select(-ID) %>% 
  data.frame()

hist.psd1 <- dat %>% 
  dplyr::select(ID, n, psd.1) %>% 
  spread(key=n, value=psd.1) %>% 
  arrange(ID) %>% 
  dplyr::select(-ID) %>% 
  data.frame()

hist.s2n1 <- dat %>% 
  dplyr::select(ID, n, s2n.1) %>% 
  spread(key=n, value=s2n.1) %>% 
  arrange(ID) %>% 
  dplyr::select(-ID) %>% 
  data.frame()

hist.s2n2 <- dat %>% 
  dplyr::select(ID, n, s2n.2) %>% 
  spread(key=n, value=s2n.2) %>% 
  arrange(ID) %>% 
  dplyr::select(-ID) %>% 
  data.frame()

obs.cov <- list(Set = hist.set, Doy = hist.Doy, psd1 = hist.psd1, s2n2 = hist.s2n2, s2n1 = hist.s2n1)

site.cov <- dat %>% 
  dplyr::select(ID, pine, wetland) %>% 
  unique() %>% 
  arrange(ID) %>% 
  dplyr::select(-ID) %>% 
  data.frame()

df.boom <- unmarkedFrameOccu(hist.boom, siteCovs=site.cov, obsCovs=obs.cov)
df.call <- unmarkedFrameOccu(hist.call, siteCovs=site.cov, obsCovs=obs.cov)

#3. Visualize----
#BOOM
ggplot(dat, aes(y=boom, x=sundiff)) +
  geom_jitter() +
  geom_smooth() #second order

ggplot(dat, aes(y=boom, x=yday)) +
  geom_jitter() +
  geom_smooth()

ggplot(dat, aes(y=boom, x=psd.1)) +
  geom_point() +
  geom_smooth() #third order?

ggplot(dat, aes(y=boom, x=s2n.1)) +
  geom_point() +
  geom_smooth() #second order

ggplot(dat, aes(y=boom, x=s2n.2)) +
  geom_point() +
  geom_smooth() #second order

#CALL
ggplot(dat, aes(y=call, x=sundiff)) +
  geom_jitter() +
  geom_smooth() #second order

ggplot(dat, aes(y=call, x=yday)) +
  geom_jitter() +
  geom_smooth()

ggplot(dat, aes(y=call, x=psd.1)) +
  geom_point() +
  geom_smooth() #third order?

ggplot(dat, aes(y=call, x=s2n.1)) +
  geom_point() +
  geom_smooth() #second order

ggplot(dat, aes(y=call, x=s2n.2)) +
  geom_point() +
  geom_smooth() #second order


#4. Model----
detboom <- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2) + psd1
                ~1,
                data=df.boom)
detboom.dredge <- dredge(detboom, trace=2)
detboom.dredge
detboom.dredge.tbl <- detboom.dredge %>% 
  data.frame() %>% 
  dplyr::select(df, AICc, delta, weight) %>% 
  mutate(AICc=round(AICc, 2),
         delta=round(delta, 2),
         weight=round(weight, 2)) %>% 
  head(6)
View(detboom.dredge.tbl)

detboom.use <- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
                    ~1,
                    data=df.boom) 
summary(detboom.use)


detcall <- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2) + psd1
                ~1,
                data=df.call)
detcall.dredge <- dredge(detcall, trace=2)
detcall.dredge
detcall.dredge.tbl <- detcall.dredge %>% 
  data.frame() %>% 
  dplyr::select(df, AICc, delta, weight) %>% 
  mutate(AICc=round(AICc, 2),
         delta=round(delta, 2),
         weight=round(weight, 2)) %>% 
  head(6)
View(detcall.dredge.tbl)

detcall.use <- occu(~ Doy + s2n2 + I(s2n2^2)
                    ~1,
                    data=df.call)
summary(detcall.use)

#Thought: makes sense that boom includes s2n1 and s2n2 because it includes both frequency bands for the call

#5. Predict----
newdat.det.boom.s2n1 <- expand.grid(s2n1 = seq(min(dat$s2n.1), max(dat$s2n.1), 0.05),
                      s2n2 = mean(dat$s2n.2),
                      Doy = mean(dat$yday),
                      Set = mean(dat$sundiff))
pred.det.boom.s2n1 <- newdat.det.boom.s2n1 %>% 
  cbind(predict(detboom.use, type="det", newdata=newdat.det.boom.s2n1)) %>% 
  dplyr::rename(Det=Predicted, DetSE=SE, DetLower=lower, DetUpper=upper) %>% 
  mutate(response="boom",
         cov="s2n1") %>% 
  rename(x=s2n1) %>% 
  dplyr::select(response, cov, x, Det, DetLower, DetUpper)

newdat.det.boom.s2n2 <- expand.grid(s2n2 = seq(min(dat$s2n.2), max(dat$s2n.2), 0.05),
                                    s2n1 = mean(dat$s2n.1),
                                    Doy = mean(dat$yday),
                                    Set = mean(dat$sundiff))
pred.det.boom.s2n2 <- newdat.det.boom.s2n2 %>% 
  cbind(predict(detboom.use, type="det", newdata=newdat.det.boom.s2n2)) %>% 
  dplyr::rename(Det=Predicted, DetSE=SE, DetLower=lower, DetUpper=upper) %>% 
  mutate(response="boom",
         cov="s2n2") %>% 
  rename(x=s2n2) %>% 
  dplyr::select(response, cov, x, Det, DetLower, DetUpper)

newdat.det.boom.doy <- expand.grid(s2n2 = mean(dat$s2n.2),
                                    s2n1 = mean(dat$s2n.1),
                                    Doy = seq(min(dat$yday), max(dat$yday), 1),
                                    Set = mean(dat$sundiff))
pred.det.boom.doy <- newdat.det.boom.doy %>% 
  cbind(predict(detboom.use, type="det", newdata=newdat.det.boom.doy)) %>% 
  dplyr::rename(Det=Predicted, DetSE=SE, DetLower=lower, DetUpper=upper) %>% 
  mutate(response="boom",
         cov="doy") %>% 
  rename(x=Doy) %>% 
  dplyr::select(response, cov, x, Det, DetLower, DetUpper)

newdat.det.boom.set <- expand.grid(s2n2 = mean(dat$s2n.2),
                                   s2n1 = mean(dat$s2n.1),
                                   Doy = mean(dat$yday),
                                   Set = seq(min(dat$sundiff), max(dat$sundiff), 0.05))
pred.det.boom.set <- newdat.det.boom.set %>% 
  cbind(predict(detboom.use, type="det", newdata=newdat.det.boom.set)) %>% 
  dplyr::rename(Det=Predicted, DetSE=SE, DetLower=lower, DetUpper=upper) %>% 
  mutate(response="boom",
         cov="set") %>% 
  rename(x=Set) %>% 
  dplyr::select(response, cov, x, Det, DetLower, DetUpper)

newdat.det.call.s2n2 <- expand.grid(s2n2 = seq(min(dat$s2n.2), max(dat$s2n.2), 0.05),
                                    s2n1 = mean(dat$s2n.1),
                                    Doy = mean(dat$yday),
                                    Set = mean(dat$sundiff))
pred.det.call.s2n2 <- newdat.det.call.s2n2 %>% 
  cbind(predict(detcall.use, type="det", newdata=newdat.det.call.s2n2)) %>% 
  dplyr::rename(Det=Predicted, DetSE=SE, DetLower=lower, DetUpper=upper) %>% 
  mutate(response="call",
         cov="s2n2") %>% 
  rename(x=s2n2) %>% 
  dplyr::select(response, cov, x, Det, DetLower, DetUpper)

newdat.det.call.doy <- expand.grid(s2n2 = mean(dat$s2n.2),
                                   s2n1 = mean(dat$s2n.1),
                                   Doy = seq(min(dat$yday), max(dat$yday), 1),
                                   Set = mean(dat$sundiff))
pred.det.call.doy <- newdat.det.call.doy %>% 
  cbind(predict(detcall.use, type="det", newdata=newdat.det.call.doy)) %>% 
  dplyr::rename(Det=Predicted, DetSE=SE, DetLower=lower, DetUpper=upper) %>% 
  mutate(response="call",
         cov="doy") %>% 
  rename(x=Doy) %>% 
  dplyr::select(response, cov, x, Det, DetLower, DetUpper)

pred.det <- rbind(pred.det.boom.s2n1, pred.det.boom.s2n2, pred.det.boom.set, pred.det.boom.doy, pred.det.call.s2n2, pred.det.call.doy)

write.csv(pred.det, "PDTGDetectionCovariatePredictions.csv", row.names = FALSE)

#6. Plot----
ggplot(pred.det) +
  geom_line(aes(x=x, y=Det, colour=cov)) +
  geom_ribbon(aes(x=x, ymin=DetLower, ymax=DetUpper), alpha=0.4) +
  scale_colour_viridis_d() +
  facet_grid(response~cov, scales="free") +
  my.theme

#7. Calculate detectability---
det.boom.det <- linearComb(detboom.use['det'], c(1, mean(dat$sundiff), mean(dat$yday), mean(dat$s2n.1), (mean(dat$s2n.1))^2, mean(dat$s2n.2), (mean(dat$s2n.2))^2)) %>% 
  backTransform() %>% 
  confint(level=0.0)
det.boom.det

det.call.det <- linearComb(detcall.use['det'], c(1, mean(dat$yday), mean(dat$s2n.2), (mean(dat$s2n.2))^2)) %>% 
  backTransform() %>% 
  confint(level=0.0)
det.call.det

#D. HABITAT COVARIATES####
#1. Check for vif & covariation----
vif(site.dat %>% 
      dplyr::select(pine, wetland) %>% 
      data.frame())

cor(site.dat %>% 
      dplyr::select(pine, wetland) %>% 
      data.frame())

#2. New data for prediciton----
new.dat <- data.frame(expand.grid(pine = seq(0, 1, 0.05),
                                  wetland = seq(0, 1, 0.05),
                                  s2n1 = mean(dat$s2n.1),
                                  s2n2 = mean(dat$s2n.2),
                                  Doy = mean(dat$yday),
                                  Set = mean(dat$sundiff)))

#3. Visualize----
#BOOM
ggplot(dat, aes(y=boom, x=pine)) +
  geom_jitter() +
  geom_smooth() #linear

ggplot(dat, aes(y=boom, x=wetland)) +
  geom_jitter() +
  geom_smooth() #2nd order

#CALL
ggplot(dat, aes(y=call, x=pine)) +
  geom_jitter() +
  geom_smooth() #linear?

ggplot(dat, aes(y=call, x=wetland)) +
  geom_jitter() +
  geom_smooth() #2nd order


#4. Spatial thinning----
boot <- 100 

cov.boom.aic <- list()
cov.call.aic <- list()
cov.boom.mod <- list()
cov.call.mod <- list()
cov.boom.pred <- list()
cov.call.pred <- list()
for(i in 1:boot){
  
  dat.pinewet <- site.dat %>% 
    dplyr::select(ProjectID, StationKey, DepYear, cell) %>% 
    group_by(cell) %>% 
    sample_n(1) %>% 
    ungroup() %>% 
    left_join(dat)
  
  #5. Format into unmarked objects----
  hist.pinewet.boom <- dat.pinewet %>% 
    dplyr::select(ID, n, boom) %>% 
    spread(key = n, value = boom) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.pinewet.call <- dat.pinewet%>% 
    dplyr::select(ID, n, call) %>% 
    spread(key = n, value = call) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.pinewet.set <- dat.pinewet %>% 
    dplyr::select(ID, n, sundiff) %>% 
    spread(key = n, value = sundiff) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame() 
  
  hist.pinewet.Doy <- dat.pinewet %>% 
    dplyr::select(ID, n, yday) %>% 
    spread(key=n, value=yday) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.pinewet.psd1 <- dat.pinewet %>% 
    dplyr::select(ID, n, psd.1) %>% 
    spread(key=n, value=psd.1) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.pinewet.s2n1 <- dat.pinewet %>% 
    dplyr::select(ID, n, s2n.1) %>% 
    spread(key=n, value=s2n.1) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.pinewet.s2n2 <- dat.pinewet %>% 
    dplyr::select(ID, n, s2n.2) %>% 
    spread(key=n, value=s2n.2) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  obs.cov.pinewet <- list(Set = hist.pinewet.set, Doy = hist.pinewet.Doy, psd1 = hist.pinewet.psd1, s2n2 = hist.pinewet.s2n2, s2n1 = hist.pinewet.s2n1)
  
  site.cov.pinewet <- dat.pinewet %>% 
    dplyr::select(ID, pine, wetland) %>% 
    unique() %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  df.boom.pinewet <- unmarkedFrameOccu(hist.pinewet.boom, siteCovs=site.cov.pinewet, obsCovs=obs.cov.pinewet)
  df.call.pinewet <- unmarkedFrameOccu(hist.pinewet.call, siteCovs=site.cov.pinewet, obsCovs=obs.cov.pinewet)
  
  #6. Model----
  cov1 <- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
               ~pine*wetland + pine*I(wetland^2),
               data=df.boom.pinewet)
  cov2<- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
              ~pine + wetland + I(wetland^2),
              data=df.boom.pinewet)
  cov3<- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
              ~pine,
              data=df.boom.pinewet)
  cov4<- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
              ~wetland + I(wetland^2),
              data=df.boom.pinewet)
  cov5<- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
              ~1,
              data=df.boom.pinewet)
  
  cov.boom <- cov3
  
  cov.boom.aic[[i]] <- aictab(list(cov1, cov2, cov3, cov4, cov5), sort=FALSE) %>% 
    data.frame() %>% 
    mutate(boot=i)
  
  #Save out top model (hand picked from preview of results, not top AIC)
  cov.boom.mod[[i]] <- summary(cov3)[['state']] %>% 
    mutate(boot=i,
           var=row.names(data.frame(cov3@estimates@estimates[["state"]]@estimates)))
  
  cov1 <- occu(~ Doy + s2n2 + I(s2n2^2)
               ~pine*wetland + pine*I(wetland^2),
               data=df.call.pinewet)
  cov2<- occu(~ Doy + s2n2 + I(s2n2^2)
              ~pine + wetland + I(wetland^2),
              data=df.call.pinewet)
  cov3<- occu(~ Doy + s2n2 + I(s2n2^2)
              ~pine,
              data=df.call.pinewet)
  cov4<- occu(~ Doy + s2n2 + I(s2n2^2)
              ~wetland +I(wetland^2),
              data=df.call.pinewet)
  cov5<- occu(~ Doy + s2n2 + I(s2n2^2)
              ~1,
              data=df.call.pinewet)
  
  cov.call <- cov4
  
  cov.call.aic[[i]] <- aictab(list(cov1, cov2, cov3, cov4, cov5), sort=FALSE) %>% 
    data.frame() %>% 
    mutate(boot=i)
  
  #Save out top model (hand picked from preview of results, not top AIC)
  cov.call.mod[[i]] <- summary(cov4)[['state']] %>% 
    mutate(boot=i,
           var=row.names(data.frame(cov4@estimates@estimates[["state"]]@estimates)))

  
  #7. Predict----
  cov.boom.pred[[i]] <- new.dat %>% 
    cbind(predict(cov.boom, type="state", newdata=new.dat)) %>% 
    dplyr::rename(Occu=Predicted, OccuE=SE, OccuLower=lower, OccuUpper=upper) %>% 
    cbind(predict(cov.boom, type="det", newdata=new.dat)) %>% 
    dplyr::rename(Det=Predicted, DetSE=SE, DetLower=lower, DetUpper=upper)
  
  cov.call.pred[[i]] <- new.dat %>% 
    cbind(predict(cov.call, type="state", newdata=new.dat)) %>% 
    dplyr::rename(Occu=Predicted, OccuE=SE, OccuLower=lower, OccuUpper=upper) %>% 
    cbind(predict(cov.call, type="det", newdata=new.dat)) %>% 
    dplyr::rename(Det=Predicted, DetSE=SE, DetLower=lower, DetUpper=upper)
  
  print(paste0("Finished bootstrap ", i))
  
}

#8. Collapse bootstrapped results----
#AIC----
cov.boom <- rbindlist(cov.boom.aic) %>% 
  group_by(Modnames) %>% 
  summarize(delta.mn=round(mean(Delta_AICc), 2),
            delta.sd=round(sd(Delta_AICc), 2),
            delta.ci=round(sd(Delta_AICc)/sqrt(100)*1.96, 2),
            wt.mn=round(mean(AICcWt), 2),
            wt.sd=round(sd(AICcWt), 2),
            wt.ci=round(sd(AICcWt)/sqrt(100)*1.96, 2),
            k=mean(K),
            aic.mn=round(mean(AICc), 2),
            aic.sd=round(sd(AICc), 2)) %>% 
  ungroup() %>% 
  mutate(aic=paste0(aic.mn, " (SD=", aic.sd, ")"),
         delta=paste0(delta.mn, " (SD=", delta.sd, ")"),
         wt=paste0(wt.mn, " (SD=", wt.sd, ")")) %>% 
  arrange(delta.mn)

cov.boom.tbl <- cov.boom %>%   
  dplyr::select(Modnames, k, aic, delta, wt)
View(cov.boom.tbl)

ggplot(cov.boom) +
  geom_bar(aes(x=Modnames, y=delta.mn), stat="identity") +
  geom_errorbar(aes(x=Modnames, ymin=delta.mn-delta.ci, ymax=delta.mn+delta.ci), width=0.2, position=position_dodge(0.9))
#PINE PINE PINE

cov.call <- rbindlist(cov.call.aic) %>% 
  group_by(Modnames) %>% 
  summarize(delta.mn=round(mean(Delta_AICc), 2),
            delta.sd=round(sd(Delta_AICc), 2),
            delta.ci=round(sd(Delta_AICc)/sqrt(100)*1.96, 2),
            wt.mn=round(mean(AICcWt), 2),
            wt.sd=round(sd(AICcWt), 2),
            wt.ci=round(sd(AICcWt)/sqrt(100)*1.96, 2),
            k=mean(K),
            aic.mn=round(mean(AICc), 2),
            aic.sd=round(sd(AICc), 2)) %>% 
  ungroup() %>% 
  mutate(aic=paste0(aic.mn, " (SD=", aic.sd, ")"),
         delta=paste0(delta.mn, " (SD=", delta.sd, ")"),
         wt=paste0(wt.mn, " (SD=", wt.sd, ")"))

cov.call.tbl <- cov.call %>%   
  dplyr::select(Modnames, k, aic, delta, wt) %>% 
  arrange(delta)
View(cov.call.tbl)

ggplot(cov.call) +
  geom_bar(aes(x=Modnames, y=delta.mn), stat="identity") +
  geom_errorbar(aes(x=Modnames, ymin=delta.mn-delta.ci, ymax=delta.mn+delta.ci), width=0.2, position=position_dodge(0.9))
#Wetland

#Coefficients----
cov.boom.coeff <- rbindlist(cov.boom.mod) %>% 
  dplyr::filter(var=="pine") %>% 
  summarize(est.mn=mean(Estimate),
            est.sd=sd(Estimate),
            se.mn=mean(SE),
            se.sd=sd(SE)) %>% 
  mutate(est.mn.p = plogis(est.mn),
         est.sd.p = plogis(est.sd),
         se.mn.p = plogis(se.mn),
         se.sd.p = plogis(se.sd))
cov.boom.coeff

cov.call.coeff <- rbindlist(cov.call.mod) %>% 
  dplyr::filter(var!="(Intercept)") %>% 
  group_by(var) %>% 
  summarize(est.mn=mean(Estimate),
            est.sd=sd(Estimate),
            se.mn=mean(SE),
            se.sd=sd(SE)) %>% 
  ungroup() %>% 
  mutate(est.mn.p = plogis(est.mn),
         est.sd.p = plogis(est.sd),
         se.mn.p = plogis(se.mn),
         se.sd.p = plogis(se.sd))
cov.call.coeff

#Predictions----
pred.cov.boom <- rbindlist(cov.boom.pred) %>% 
  group_by(pine) %>% 
  summarize(Occu = mean(Occu),
            OccuLower = mean(OccuLower),
            OccuUpper = mean(OccuUpper)) %>% 
  ungroup()
head(pred.cov.boom)

pred.cov.call <- rbindlist(cov.call.pred) %>% 
  group_by(wetland) %>% 
  summarize(Occu = mean(Occu),
            OccuLower = mean(OccuLower),
            OccuUpper = mean(OccuUpper)) %>% 
  ungroup()
head(pred.cov.call)

write.csv(pred.cov.boom, "PDTGVegetationPredictionsBoom.csv", row.names = FALSE)
write.csv(pred.cov.call, "PDTGVegetationPredictionsCall.csv", row.names = FALSE)

#10. Plot----
plot.cov.boom <- ggplot(pred.cov.boom) +
  geom_line(aes(x=pine, y=Occu))+
  geom_ribbon(aes(x=pine, ymin=OccuLower, ymax=OccuUpper), colour="grey70", alpha=0.2) +
  my.theme

plot.cov.call <- ggplot(pred.cov.call) +
  geom_line(aes(x=wetland, y=Occu))+
  geom_ribbon(aes(x=wetland, ymin=OccuLower, ymax=OccuUpper), colour="grey70", alpha=0.2) +
  my.theme

grid.arrange(plot.cov.boom, plot.cov.call, ncol=2)

#E. HOW TO QUANTIFY DISTURBANCE?####
#1. Wrangle covariates----
site.1t2d <- site.dat %>% 
  dplyr::filter(types==1, count>1) %>% 
  group_by(StationKey, DepYear, cell, pine, wetland, disturbance) %>% 
  mutate(n=row_number()) %>% 
  ungroup() %>% 
  pivot_wider(id_cols=c(StationKey, DepYear, cell, pine, wetland, disturbance),
              names_from=n,
              values_from=c(time, sitearea)) %>% 
  rowwise(StationKey, DepYear, cell, pine, wetland, disturbance) %>% 
  mutate(timemin=min(time_1, time_2),
         timemean=mean(time_1, time_2),
         siteareamax=max(sitearea_1, sitearea_2),
         siteareasum=sum(sitearea_1, sitearea_2)) %>% 
  ungroup() %>% 
  arrange(StationKey)

table(site.1t2d$disturbance)
nrow(site.1t2d)

#2. Check for vif & covariation----
vif(site.1t2d %>% 
      dplyr::select(timemin, timemean, siteareamax, siteareasum) %>% 
      data.frame())

cor(site.1t2d %>% 
      dplyr::select(timemin, timemean, siteareamax, siteareasum) %>% 
      data.frame())
#Take out timemean, siteareasum

ggpairs(site.1t2d %>% 
          dplyr::select(pine, disturbance, timemin, siteareamax) %>% 
          data.frame())

#3. Spatial thinning----
boot <- 100

occ2t.boom.aic <- list()
occ2t.call.aic <- list()
for(i in 1:boot){
  
  dat.1t2d <- data.frame()
  
  for(h in 1:length(disturbances)){
    
    disturbance.h <- disturbances[h]
    
    site.h <- site.1t2d %>% 
      dplyr::filter(disturbance==disturbance.h)
    
    dat.h <- site.h %>% 
      group_by(cell) %>% 
      sample_n(1) %>% 
      ungroup() %>% 
      left_join(dat)
    
    dat.1t2d <- rbind(dat.1t2d, dat.h)
  }
  
  table(dat.1t2d$boom, dat.1t2d$call)
  
  #5. Format into unmarked objects----
  hist.boom.1t2d <- dat.1t2d %>% 
    dplyr::select(ID, n, boom) %>% 
    spread(key = n, value = boom) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.call.1t2d <- dat.1t2d%>% 
    dplyr::select(ID, n, call) %>% 
    spread(key = n, value = call) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.Set.1t2d <- dat.1t2d %>%
    dplyr::select(ID, n, sundiff) %>% 
    spread(key = n, value = sundiff) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame() 
  
  hist.Doy.1t2d <- dat.1t2d %>% 
    dplyr::select(ID, n, yday) %>% 
    spread(key=n, value=yday) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.s2n1.1t2d <- dat.1t2d %>% 
    dplyr::select(ID, n, s2n.1) %>% 
    spread(key=n, value=s2n.1) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.s2n2.1t2d <- dat.1t2d %>% 
    dplyr::select(ID, n, s2n.2) %>% 
    spread(key=n, value=s2n.2) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  obs.cov.1t2d <- list(Set = hist.Set.1t2d, Doy = hist.Doy.1t2d, s2n1 = hist.s2n1.1t2d, s2n2 = hist.s2n2.1t2d)
  
  site.cov.1t2d <- dat.1t2d %>% 
    dplyr::select(ID, disturbance, pine, wetland, timemin, siteareamax) %>% 
    unique() %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  df.boom <- unmarkedFrameOccu(hist.boom.1t2d, siteCovs=site.cov.1t2d, obsCovs=obs.cov.1t2d)
  df.call <- unmarkedFrameOccu(hist.call.1t2d, siteCovs=site.cov.1t2d, obsCovs=obs.cov.1t2d)
  
  #6. Model----
  rm(list=ls(pattern="^mod"))
  #BOOM
  mod1 <- try(occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
                   ~timemin*siteareamax + pine,
                   data=df.boom))
  mod2 <- try(occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
                   ~timemin + siteareamax + pine,
                   data=df.boom))
  mod3 <- try(occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
                   ~timemin + pine,
                   data=df.boom))
  mod4 <- try(occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
                   ~siteareamax + pine,
                   data=df.boom))
  
  rm(err)
  err <- try(aictab(list(mod1, mod2, mod3, mod4), sort=FALSE))
  
  if(inherits(err, "try-error")){
    next
  }
  
  occ2t.boom.aic[[i]] <- aictab(list(mod1, mod2, mod3, mod4), sort=FALSE) %>% 
    data.frame() %>% 
    mutate(boot=i,
           nboom=boom)
  
  rm(list=ls(pattern="^mod"))
  #CALL
  mod1 <- try(occu(~ Doy + s2n2 + I(s2n2^2)
                   ~timemin*siteareamax + wetland + I(wetland^2),
                   data=df.call))
  mod2 <- try(occu(~ Doy + s2n2 + I(s2n2^2)
                   ~timemin + siteareamax + wetland + I(wetland^2),
                   data=df.call))
  mod3 <- try(occu(~ Doy + s2n2 + I(s2n2^2)
                   ~timemin + wetland + I(wetland^2),
                   data=df.call))
  mod4 <- try(occu(~ Doy + s2n2 + I(s2n2^2)
                   ~siteareamax + wetland + I(wetland^2),
                   data=df.call))
  
  rm(err)
  err <- try(aictab(list(mod1, mod2, mod3, mod4)))
  
  if(inherits(err, "try-error")) {
    next
  }
  
  occ2t.call.aic[[i]] <- aictab(list(mod1, mod2, mod3, mod4), sort=FALSE) %>% 
    data.frame() %>% 
    mutate(boot=i,
           ncall=call)
  
  print(paste0("Finished bootstrap ", i))
  
}

#7. Collapse bootstrapped results----
occ2t.boom <- rbindlist(occ2t.boom.aic) %>% 
  group_by(Modnames) %>% 
  summarize(delta.mn=round(mean(Delta_AICc), 2),
            delta.sd=round(sd(Delta_AICc), 2),
            delta.ci=round(sd(Delta_AICc)/sqrt(100)*1.96, 2),
            wt.mn=round(mean(AICcWt), 2),
            wt.sd=round(sd(AICcWt), 2),
            wt.ci=round(sd(AICcWt)/sqrt(100)*1.96, 2),
            k=mean(K),
            aic.mn=round(mean(AICc), 2),
            aic.sd=round(sd(AICc), 2)) %>% 
  ungroup() %>% 
  mutate(aic=paste0(aic.mn, " (SD=", aic.sd, ")"),
         delta=paste0(delta.mn, " (SD=", delta.sd, ")"),
         wt=paste0(wt.mn, " (SD=", wt.sd, ")")) %>% 
  arrange(delta.mn)

occ2t.boom.tbl <- occ2t.boom %>%   
  dplyr::select(Modnames, k, aic, delta, wt)
View(occ2t.boom.tbl)

ggplot(occ2t.boom) +
  geom_bar(aes(x=Modnames, y=delta.mn), stat="identity") +
  geom_errorbar(aes(x=Modnames, ymin=delta.mn-delta.ci, ymax=delta.mn+delta.ci), width=0.2, position=position_dodge(0.9))
#Minimum time

occ2t.call <- rbindlist(occ2t.call.aic) %>% 
  group_by(Modnames) %>% 
  summarize(delta.mn=round(mean(Delta_AICc), 2),
            delta.sd=round(sd(Delta_AICc), 2),
            delta.ci=round(sd(Delta_AICc)/sqrt(100)*1.96, 2),
            wt.mn=round(mean(AICcWt), 2),
            wt.sd=round(sd(AICcWt), 2),
            wt.ci=round(sd(AICcWt)/sqrt(100)*1.96, 2),
            k=mean(K),
            aic.mn=round(mean(AICc), 2),
            aic.sd=round(sd(AICc), 2)) %>% 
  ungroup() %>% 
  mutate(aic=paste0(aic.mn, " (SD=", aic.sd, ")"),
         delta=paste0(delta.mn, " (SD=", delta.sd, ")"),
         wt=paste0(wt.mn, " (SD=", wt.sd, ")"))

occ2t.call.tbl <- occ2t.call %>%   
  dplyr::select(Modnames, k, aic, delta, wt) %>% 
  arrange(delta)
View(occ2t.call.tbl)
#Minimum time

#F. DOES DISTURBANCE TYPE MATTER?####
#1. Wrangle covariates----
site.1d1t <- site.dat %>% 
  dplyr::filter(types==1) %>% 
  group_by(StationKey, DepYear, cell, Latitude, Longitude, disturbance, pine, wetland) %>% 
  summarize(time=min(time)) %>% 
  ungroup()
  
table(site.1d1t$time, site.1d1t$disturbance)

#2. Check for vif & covariation----
vif(site.1d1t %>% 
      dplyr::select(time, pine, wetland) %>% 
      data.frame())

cor(site.1d1t %>% 
      dplyr::select(time, pine, wetland) %>% 
      data.frame())
#All good

ggpairs(site.1d1t %>% 
          dplyr::select(time, disturbance, pine, wetland))

#3. New data for prediction----
hist.fire <- site.1d1t %>% 
  dplyr::filter(disturbance=="fire")
fire.max <- max(hist.fire$time)
fire.min <- min(hist.fire$time)

hist.cc <- site.1d1t %>% 
  dplyr::filter(disturbance=="cc")
cc.max <- max(hist.cc$time)
cc.min <- min(hist.cc$time)

hist.wells <- site.1d1t %>% 
  dplyr::filter(disturbance=="well")
wells.max <- max(hist.wells$time)
wells.min <- min(hist.wells$time)

new.dat <- data.frame(expand.grid(disturbance=c("fire", "well", "cc"),
                                  time=seq(min(site.1d1t$time), max(site.1d1t$time, 1)),
                                  pine = c(0.1, 0.9),
                                  wetland = c(0.1, 0.9),
                                  s2n1 = mean(dat$s2n.1),
                                  s2n2 = mean(dat$s2n.2),
                                  Doy = mean(dat$yday),
                                  Set = mean(dat$sundiff)))

#4. Spatial thinning----
boot <- 100
disturbances <- unique(site.1d1t$disturbance)

occ.1d.boom.aic <- list()
occ.1d.call.aic <- list()
occ.1d.boom.mod <- list()
occ.1d.call.mod <- list()
occ.1d.boom.pred <- list()
occ.1d.call.pred <- list()
for(g in 1:boot){
  
  dat.1d1t <- data.frame()
  
  for(h in 1:length(disturbances)){
    
    disturbance.h <- disturbances[h]
    
    site.h <- site.1d1t %>% 
      dplyr::filter(disturbance==disturbance.h)
    
    dat.h <- site.h %>% 
      group_by(cell) %>% 
      sample_n(1) %>% 
      ungroup() %>% 
      left_join(dat)
    
    dat.1d1t <- rbind(dat.1d1t, dat.h)
  }
  
  #6. Format into unmarked objects----
  hist.boom.1d1t <- dat.1d1t %>% 
    dplyr::select(ID, n, boom) %>% 
    spread(key = n, value = boom) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.call.1d1t <- dat.1d1t%>% 
    dplyr::select(ID, n, call) %>% 
    spread(key = n, value = call) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.Set.1d1t <- dat.1d1t %>% 
    dplyr::select(ID, n, sundiff) %>% 
    spread(key=n, value=sundiff) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.Doy.1d1t <- dat.1d1t %>% 
    dplyr::select(ID, n, yday) %>% 
    spread(key=n, value=yday) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.s2n1.1d1t <- dat.1d1t %>% 
    dplyr::select(ID, n, s2n.1) %>% 
    spread(key=n, value=s2n.1) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.s2n2.1d1t <- dat.1d1t %>% 
    dplyr::select(ID, n, s2n.2) %>% 
    spread(key=n, value=s2n.2) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  obs.cov.1d1t <- list(Set = hist.Set.1d1t, Doy = hist.Doy.1d1t, s2n2 = hist.s2n2.1d1t, s2n1 = hist.s2n1.1d1t)
  
  site.cov.1d1t <- dat.1d1t %>% 
    dplyr::select(ID, pine, wetland, time, disturbance) %>% 
    unique() %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  df.boom <- unmarkedFrameOccu(hist.boom.1d1t, siteCovs=site.cov.1d1t, obsCovs=obs.cov.1d1t)
  df.call <- unmarkedFrameOccu(hist.call.1d1t, siteCovs=site.cov.1d1t, obsCovs=obs.cov.1d1t)
  
  #7. Model----
  #BOOM
  occ1 <- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
               ~time*pine + time*disturbance,
               data=df.boom)
  occ2 <- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
               ~time*pine + disturbance,
               data=df.boom)
  occ3 <- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
               ~time*disturbance + pine,
               data=df.boom)
  occ4 <- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
               ~time + disturbance + pine,
               data=df.boom)
  occ5 <- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
               ~time*pine,
               data=df.boom)
  occ6 <- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
               ~time + pine,
               data=df.boom)
  occ7 <- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
               ~disturbance*pine,
               data=df.boom)
  occ8 <- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
               ~disturbance + pine,
               data=df.boom)
  occ9 <- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
               ~pine,
               data=df.boom)
  
  #Save out AIC
  occ.1d.boom.aic[[g]] <- aictab(list(occ1, occ2, occ3, occ4, occ5, occ6, occ7, occ8, occ9), sort=FALSE) %>% 
    data.frame() %>% 
    mutate(boot=g)
  
  #Save out top model (hand picked from preview of results, not top AIC)
  occ.1d.boom.mod[[g]] <- summary(occ1)[['state']] %>% 
    mutate(boot=g,
           var=row.names(data.frame(occ1@estimates@estimates[["state"]]@estimates)))
  
  #CALL
  occ1 <- occu(~ Doy + s2n2 + I(s2n2^2)
               ~time*wetland + time*I(wetland^2) + time*disturbance,
               data=df.call)
  occ2 <- occu(~ Doy + s2n2 + I(s2n2^2)
               ~time*wetland + time*I(wetland^2) + disturbance,
               data=df.call)
  occ3 <- occu(~ Doy + s2n2 + I(s2n2^2)
               ~time + wetland + I(wetland^2) + time*disturbance,
               data=df.call)
  occ4 <- occu(~ Doy + s2n2 + I(s2n2^2)
               ~time + disturbance + wetland + I(wetland^2),
               data=df.call)
  occ5 <- occu(~ Doy + s2n2 + I(s2n2^2)
               ~time*wetland + time*I(wetland^2),
               data=df.call)
  occ6 <- occu(~ Doy + s2n2 + I(s2n2^2)
               ~time + wetland + I(wetland^2),
               data=df.call)
  occ7 <- occu(~ Doy + s2n2 + I(s2n2^2)
               ~disturbance*wetland + disturbance*I(wetland^2),
               data=df.call)
  occ8 <- occu(~ Doy + s2n2 + I(s2n2^2)
               ~disturbance + wetland + I(wetland^2),
               data=df.call)
  occ9 <- occu(~ Doy + s2n2 + I(s2n2^2)
               ~wetland + I(wetland^2),
               data=df.call)
  
  #Save out AIC
  occ.1d.call.aic[[g]] <- aictab(list(occ1, occ2, occ3, occ4, occ5, occ6, occ7, occ8, occ9), sort=FALSE) %>% 
    data.frame() %>% 
    mutate(boot=g)
  
  #Save out top model (hand picked from preview of results, not top AIC)
  occ.1d.call.mod[[g]] <- summary(occ8)[['state']] %>% 
    mutate(boot=g,
           var=row.names(data.frame(occ8@estimates@estimates[["state"]]@estimates)))
  
  #8. Predict----
  occ.1d.boom.pred[[g]] <- new.dat %>% 
    cbind(predict(occboom, type="state", newdata=new.dat)) %>% 
    dplyr::rename(Occu=Predicted, OccuE=SE, OccuLower=lower, OccuUpper=upper) %>% 
    cbind(predict(occboom, type="det", newdata=new.dat)) %>% 
    dplyr::rename(Det=Predicted, DetSE=SE, DetLower=lower, DetUpper=upper) %>% 
    dplyr::filter((disturbance=="fire" & time <= fire.max & time >= fire.min) |
                    (disturbance=="well" & time <= wells.max & time >= wells.min) |
                    (disturbance=="cc" & time <= cc.max & time >= cc.min))
  
  occ.1d.call.pred[[g]] <- new.dat %>% 
    cbind(predict(occcall, type="state", newdata=new.dat)) %>% 
    dplyr::rename(Occu=Predicted, OccuE=SE, OccuLower=lower, OccuUpper=upper) %>% 
    cbind(predict(occcall, type="det", newdata=new.dat)) %>% 
    dplyr::rename(Det=Predicted, DetSE=SE, DetLower=lower, DetUpper=upper) %>% 
    dplyr::filter((disturbance=="fire" & time <= fire.max & time >= fire.min) |
                    (disturbance=="well" & time <= wells.max & time >= wells.min) |
                    (disturbance=="cc" & time <= cc.max & time >= cc.min))
  
  print(paste0("COMPLETED BOOTSTRAP ", g))
  
}


#9. Collapse bootstrapped results----
#AIC----
aic.1d.boom <- rbindlist(occ.1d.boom.aic) %>% 
  group_by(Modnames) %>% 
  summarize(delta.mn=round(mean(Delta_AICc), 2),
            delta.sd=round(sd(Delta_AICc), 2),
            delta.ci=round(sd(Delta_AICc)/sqrt(100)*1.96, 2),
            wt.mn=round(mean(AICcWt), 2),
            wt.sd=round(sd(AICcWt), 2),
            wt.ci=round(sd(AICcWt)/sqrt(100)*1.96, 2),
            k=mean(K),
            aic.mn=round(mean(AICc), 2),
            aic.sd=round(sd(AICc), 2)) %>% 
  ungroup() %>% 
  mutate(aic=paste0(aic.mn, " (SD=", aic.sd, ")"),
         delta=paste0(delta.mn, " (SD=", delta.sd, ")"),
         wt=paste0(wt.mn, " (SD=", wt.sd, ")")) %>% 
  arrange(delta.mn)

aic.1d.boom.tbl <- aic.1d.boom %>%   
  dplyr::select(Modnames, k, aic, delta, wt)
View(aic.1d.boom.tbl)

ggplot(aic.1d.boom) +
  geom_bar(aes(x=Modnames, y=delta.mn), stat="identity") +
  geom_errorbar(aes(x=Modnames, ymin=delta.mn-delta.ci, ymax=delta.mn+delta.ci), width=0.2, position=position_dodge(0.9))

aic.1d.call <- rbindlist(occ.1d.call.aic) %>% 
  group_by(Modnames) %>% 
  summarize(delta.mn=round(mean(Delta_AICc), 2),
            delta.sd=round(sd(Delta_AICc), 2),
            delta.ci=round(sd(Delta_AICc)/sqrt(100)*1.96, 2),
            wt.mn=round(mean(AICcWt), 2),
            wt.sd=round(sd(AICcWt), 2),
            wt.ci=round(sd(AICcWt)/sqrt(100)*1.96, 2),
            k=mean(K),
            aic.mn=round(mean(AICc), 2),
            aic.sd=round(sd(AICc), 2)) %>% 
  ungroup() %>% 
  mutate(aic=paste0(aic.mn, " (SD=", aic.sd, ")"),
         delta=paste0(delta.mn, " (SD=", delta.sd, ")"),
         wt=paste0(wt.mn, " (SD=", wt.sd, ")"))

aic.1d.call.tbl <- aic.1d.call %>%   
  dplyr::select(Modnames, k, aic, delta, wt) %>% 
  arrange(delta)
View(aic.1d.call.tbl)

ggplot(aic.1d.call) +
  geom_bar(aes(x=Modnames, y=delta.mn), stat="identity") +
  geom_errorbar(aes(x=Modnames, ymin=delta.mn-delta.ci, ymax=delta.mn+delta.ci), width=0.2, position=position_dodge(0.9))

#Coefficients----
mod.1d.boom <- rbindlist(occ.1d.boom.mod) %>% 
  group_by(var) %>% 
  summarize(est.mn=mean(Estimate),
            est.sd=sd(Estimate),
            se.mn=mean(SE),
            se.sd=sd(SE)) %>% 
  ungroup()
mod.1d.boom

mod.1d.call <- rbindlist(occ.1d.call.mod) %>% 
  group_by(var) %>% 
  summarize(est.mn=mean(Estimate),
            est.sd=sd(Estimate),
            se.mn=mean(SE),
            se.sd=sd(SE)) %>% 
  ungroup()
mod.1d.call

#Predictions----
pred.1d.boom <- rbindlist(occ.1d.boom.pred) %>% 
  group_by(disturbance, time, pine, wetland) %>% 
  summarize(Occu = mean(Occu),
            OccuLower = mean(OccuLower),
            OccuUpper = mean(OccuUpper)) %>% 
  ungroup()
head(pred.1d.boom)

pred.1d.call <- rbindlist(occ.1d.call.pred) %>% 
  group_by(disturbance, time, pine, wetland) %>% 
  summarize(Occu = mean(Occu),
            OccuLower = mean(OccuLower),
            OccuUpper = mean(OccuUpper)) %>% 
  ungroup()
head(pred.1d.call)

write.csv(pred.1d.boom, "PDTGSingleDisturbancePredictionsBoom.csv", row.names=FALSE)
write.csv(pred.1d.call, "PDTGSingleDisturbancePredictionsCall.csv", row.names=FALSE)

#10. Plot----
plot.1d.boom <- ggplot(pred.1d.boom) +
  geom_line(aes(x=time, y=Occu, colour=factor(pine)))+
  geom_ribbon(aes(x=time, ymin=OccuLower, ymax=OccuUpper, group=factor(factor(pine))), colour="grey70", alpha=0.2) +
  facet_wrap(~disturbance) +
  labs(x="Years since disturbance", y="Territorial habitat use", colour="% pine")+
  theme(legend.justification=c(0,0), legend.position=c(0,0)) +
  theme(legend.justification=c(0,0), legend.position=c(0,0)) +
  xlim(c(0,80)) +
  ylim(c(0,1)) +
  my.theme

plot.1d.call <- ggplot(pred.1d.call) +
  geom_line(aes(x=time, y=Occu, colour=factor(wetland)))+
  geom_ribbon(aes(x=time, ymin=OccuLower, ymax=OccuUpper, group=factor(factor(wetland))), colour="grey70", alpha=0.2) +
  facet_wrap(~disturbance) +
  labs(x="Years since disturbance", y="Home range habitat use", colour="% wetland")+
  theme(legend.justification=c(0,0), legend.position=c(0,0)) +
  theme(legend.justification=c(0,0), legend.position=c(0,0)) +
  xlim(c(0,80)) +
  ylim(c(0,1)) +
  my.theme

grid.arrange(plot.1d.boom, plot.1d.call, nrow=2)



#G. IF THERE ARE TWO TYPES, DO THEY INTERACT?####
#1. Wrangle covariates----
#Take out stations where two disturbances are the same year
site.mintime.2d <- site.dat %>% 
  dplyr::filter(types==2) %>% 
  group_by(StationKey, DepYear, cell, Latitude, Longitude, disturbance, pine, wetland) %>% 
  summarize(time=min(time)) %>% 
  ungroup()

site.2d <- site.mintime.2d %>% 
  mutate(ID = row_number()) %>% 
  pivot_wider(id_cols=c(StationKey, DepYear, cell, pine, wetland),
              names_from=disturbance,
              names_prefix="time.",
              values_from=c(time))  %>% 
  rowwise(StationKey, DepYear, cell, pine, wetland) %>% 
  mutate(time.min = min(time.cc, time.well, time.fire, na.rm=TRUE)) %>% 
  ungroup() %>% 
  left_join(site.mintime.2d %>% 
              dplyr::rename(time.min=time)) %>% 
  dplyr::rename(time.dist=disturbance) %>% 
  mutate(combo=case_when(is.na(time.cc) ~ "frwell",
                         is.na(time.well) ~ "ccfr",
                         is.na(time.fire) ~ "wellcc")) %>% 
  mutate(ID=paste(StationKey, DepYear))

#2. Visualize availability----
ggplot(dat %>% 
         dplyr::filter(types==2)) + 
  geom_jitter(aes(x=pine, y=time, colour=disturbance)) +
  geom_smooth(aes(x=pine, y=time, colour=disturbance), method="lm") + 
  facet_wrap(~disturbance)


#G1. CC vs FIRE####

#3. Check for vif & covariation----
site.ccfr <- site.2d %>%
  dplyr::filter(is.na(time.well))

vif(site.ccfr %>% 
      dplyr::select(time.cc, time.fire, time.min, pine, wetland) %>% 
      data.frame())

cor(site.ccfr %>% 
      dplyr::select(time.cc, time.fire, time.min, pine, wetland) %>% 
      data.frame())

ggpairs(site.ccfr %>% 
          dplyr::select(time.cc, time.fire, time.min, pine, wetland) %>% 
          data.frame())

#pine and time.fire strongly negatively correlated
#time.min and time.cc very strongly negatively correlated


boot <- 100

occ.ccfr.boom.aic <- list()
occ.ccfr.call.aic <- list()
occ.ccfr.boom.mod <- list()
occ.ccfr.call.mod <- list()
for(i in 1:boot){
  
  #4. Spatial thinning----
  dat.ccfr <- site.2d %>%
    dplyr::filter(is.na(time.well)) %>% 
    unique() %>% 
    group_by(cell) %>% 
    sample_n(1) %>% 
    ungroup() %>% 
    left_join(dat)
  
  table(dat.ccfr$boom, dat.ccfr$call)
  
  #5. Format into unmarked objects----
  hist.boom.ccfr <- dat.ccfr %>% 
    dplyr::select(ID, n, boom) %>% 
    spread(key = n, value = boom) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.call.ccfr <- dat.ccfr%>% 
    dplyr::select(ID, n, call) %>% 
    spread(key = n, value = call) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.Set.ccfr <- dat.ccfr %>% 
    dplyr::select(ID, n, sundiff) %>% 
    spread(key=n, value=sundiff) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.Doy.ccfr <- dat.ccfr %>% 
    dplyr::select(ID, n, yday) %>% 
    spread(key=n, value=yday) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.s2n1.ccfr <- dat.ccfr %>% 
    dplyr::select(ID, n, s2n.1) %>% 
    spread(key=n, value=s2n.1) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.s2n2.ccfr <- dat.ccfr %>% 
    dplyr::select(ID, n, s2n.2) %>% 
    spread(key=n, value=s2n.2) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  obs.cov.ccfr <- list(Set = hist.Set.ccfr, Doy = hist.Doy.ccfr, s2n1 = hist.s2n1.ccfr, s2n2 = hist.s2n2.ccfr)
  
  site.cov.ccfr <- dat.ccfr %>% 
    dplyr::select(ID, pine, wetland, time.min, time.dist, time.cc, time.fire) %>% 
    unique() %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  df.boom.ccfr <- unmarkedFrameOccu(hist.boom.ccfr, siteCovs=site.cov.ccfr, obsCovs=obs.cov.ccfr)
  df.call.ccfr <- unmarkedFrameOccu(hist.call.ccfr, siteCovs=site.cov.ccfr, obsCovs=obs.cov.ccfr)
  
  #6. Model----
  #BOOM
  occ1 <- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
               ~time.cc*time.fire*pine,
               data=df.boom.ccfr)
  occ2 <- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
               ~time.cc*pine + time.fire*pine,
               data=df.boom.ccfr)
  occ3 <- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
               ~time.cc*pine,
               data=df.boom.ccfr)
  occ4 <- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
               ~time.fire*pine,
               data=df.boom.ccfr)
  occ5 <- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
               ~time.min*pine,
               data=df.boom.ccfr)
  
  #Save out AIC
  occ.ccfr.boom.aic[[i]] <- aictab(list(occ1, occ2, occ3, occ4, occ5), sort=FALSE) %>% 
    data.frame() %>% 
    mutate(boot=i)
  
  #Save out top model (hand picked from preview of results, not top AIC)
  occ.ccfr.boom.mod[[i]] <- summary(occ4)[['state']] %>% 
    mutate(boot=i,
           var=row.names(data.frame(occ4@estimates@estimates[["state"]]@estimates)))
  
  print(paste0("COMPLETED BOOTSTRAP ", i))
  
}

#7. Collapse bootstrapped results----
#AIC----
aic.ccfr.boom <- rbindlist(occ.ccfr.boom.aic) %>% 
  group_by(Modnames) %>% 
  summarize(delta.mn=round(mean(Delta_AICc), 2),
            delta.sd=round(sd(Delta_AICc), 2),
            delta.ci=round(sd(Delta_AICc)/sqrt(100)*1.96, 2),
            wt.mn=round(mean(AICcWt), 2),
            wt.sd=round(sd(AICcWt), 2),
            wt.ci=round(sd(AICcWt)/sqrt(100)*1.96, 2),
            k=mean(K),
            aic.mn=round(mean(AICc), 2),
            aic.sd=round(sd(AICc), 2)) %>% 
  ungroup() %>% 
  mutate(aic=paste0(aic.mn, " (SD=", aic.sd, ")"),
         delta=paste0(delta.mn, " (SD=", delta.sd, ")"),
         wt=paste0(wt.mn, " (SD=", wt.sd, ")")) %>% 
  arrange(delta.mn)

aic.ccfr.boom.tbl <- aic.ccfr.boom %>%   
  dplyr::select(Modnames, k, aic, delta, wt)
View(aic.ccfr.boom.tbl)

ggplot(aic.ccfr.boom) +
  geom_bar(aes(x=Modnames, y=delta.mn), stat="identity") +
  geom_errorbar(aes(x=Modnames, ymin=delta.mn-delta.ci, ymax=delta.mn+delta.ci), width=0.2, position=position_dodge(0.9))

#Coefficients----
mod.ccfr.boom <- rbindlist(occ.ccfr.boom.mod) %>% 
  group_by(var) %>% 
  summarize(est.mn=mean(Estimate),
            est.sd=sd(Estimate),
            se.mn=mean(SE),
            se.sd=sd(SE)) %>% 
  ungroup()
mod.ccfr.boom

#G2. FIRE VS WELL####

#3. Check for vif & covariation----
site.firewell <- site.2d %>%
  dplyr::filter(is.na(time.cc))

vif(site.firewell %>% 
      dplyr::select(time.fire, time.well, time.min, pine, wetland) %>% 
      data.frame())

cor(site.firewell %>% 
      dplyr::select(time.fire, time.well, time.min, pine, wetland) %>% 
      data.frame())

ggpairs(site.firewell %>% 
          dplyr::select(time.fire, time.well, time.min, pine, wetland) %>% 
          data.frame())
#time min pretty strongly correlated with time.fire and time.well

boot <- 100

occ.firewell.boom.aic <- list()
occ.firewell.call.aic <- list()
occ.firewell.boom.mod <- list()
occ.firewell.call.mod <- list()
for(i in 1:boot){
  
  #4. Spatial thinning----
  dat.firewell <- site.2d %>%
    dplyr::filter(is.na(time.cc)) %>% 
    unique() %>% 
    group_by(cell) %>% 
    sample_n(1) %>% 
    ungroup() %>% 
    left_join(dat)
  
  table(dat.firewell$boom, dat.firewell$call)
  
  #5. Format into unmarked objects----
  hist.boom.firewell <- dat.firewell %>% 
    dplyr::select(ID, n, boom) %>% 
    spread(key = n, value = boom) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.call.firewell <- dat.firewell%>% 
    dplyr::select(ID, n, call) %>% 
    spread(key = n, value = call) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.Set.firewell <- dat.firewell %>% 
    dplyr::select(ID, n, sundiff) %>% 
    spread(key=n, value=sundiff) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.Doy.firewell <- dat.firewell %>% 
    dplyr::select(ID, n, yday) %>% 
    spread(key=n, value=yday) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.s2n1.firewell <- dat.firewell %>% 
    dplyr::select(ID, n, s2n.1) %>% 
    spread(key=n, value=s2n.1) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.s2n2.firewell <- dat.firewell %>% 
    dplyr::select(ID, n, s2n.2) %>% 
    spread(key=n, value=s2n.2) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  obs.cov.firewell <- list(Set = hist.Set.firewell, Doy = hist.Doy.firewell, s2n1 = hist.s2n1.firewell, s2n2 = hist.s2n2.firewell)
  
  site.cov.firewell <- dat.firewell %>% 
    dplyr::select(ID, pine, wetland, time.min, time.dist, time.fire, time.well) %>% 
    unique() %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  df.boom.firewell <- unmarkedFrameOccu(hist.boom.firewell, siteCovs=site.cov.firewell, obsCovs=obs.cov.firewell)
  df.call.firewell <- unmarkedFrameOccu(hist.call.firewell, siteCovs=site.cov.firewell, obsCovs=obs.cov.firewell)
  
  #6. Model----
  #BOOM
  occ1 <- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
               ~time.well*time.fire*pine,
               data=df.boom.firewell)
  occ2 <- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
               ~time.well*pine + time.fire*pine,
               data=df.boom.firewell)
  occ3 <- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
               ~time.well*pine,
               data=df.boom.firewell)
  occ4 <- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
               ~time.fire*pine,
               data=df.boom.firewell)
  occ5 <- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
               ~time.min*pine,
               data=df.boom.firewell)
  
  #Save out AIC
  occ.firewell.boom.aic[[i]] <- aictab(list(occ1, occ2, occ3, occ4, occ5), sort=FALSE) %>% 
    data.frame() %>% 
    mutate(boot=i)
  
  #Save out top model (hand picked from preview of results, not top AIC)
  occ.firewell.boom.mod[[i]] <- summary(occ4)[['state']] %>% 
    mutate(boot=i,
           var=row.names(data.frame(occ4@estimates@estimates[["state"]]@estimates)))

  print(paste0("COMPLETED BOOTSTRAP ", i))
  
}

#7. Collapse bootstrapped results----
#AIC----
aic.firewell.boom <- rbindlist(occ.firewell.boom.aic) %>% 
  group_by(Modnames) %>% 
  summarize(delta.mn=round(mean(Delta_AICc), 2),
            delta.sd=round(sd(Delta_AICc), 2),
            delta.ci=round(sd(Delta_AICc)/sqrt(100)*1.96, 2),
            wt.mn=round(mean(AICcWt), 2),
            wt.sd=round(sd(AICcWt), 2),
            wt.ci=round(sd(AICcWt)/sqrt(100)*1.96, 2),
            k=mean(K),
            aic.mn=round(mean(AICc), 2),
            aic.sd=round(sd(AICc), 2)) %>% 
  ungroup() %>% 
  mutate(aic=paste0(aic.mn, " (SD=", aic.sd, ")"),
         delta=paste0(delta.mn, " (SD=", delta.sd, ")"),
         wt=paste0(wt.mn, " (SD=", wt.sd, ")")) %>% 
  arrange(delta.mn)

aic.firewell.boom.tbl <- aic.firewell.boom %>%   
  dplyr::select(Modnames, k, aic, delta, wt)
View(aic.firewell.boom.tbl)

ggplot(aic.firewell.boom) +
  geom_bar(aes(x=Modnames, y=delta.mn), stat="identity") +
  geom_errorbar(aes(x=Modnames, ymin=delta.mn-delta.ci, ymax=delta.mn+delta.ci), width=0.2, position=position_dodge(0.9))

#Coefficients----
mod.firewell.boom <- rbindlist(occ.firewell.boom.mod) %>% 
  group_by(var) %>% 
  summarize(est.mn=mean(Estimate),
            est.sd=sd(Estimate),
            se.mn=mean(SE),
            se.sd=sd(SE)) %>% 
  ungroup()
mod.firewell.boom


#G3. WELL vs CC####

#3. Check for vif & covariation----
site.wellcc <- site.2d %>%
  dplyr::filter(is.na(time.fire))

vif(site.wellcc %>% 
      dplyr::select(time.cc, time.well, time.min, pine, wetland) %>% 
      data.frame())

cor(site.wellcc %>% 
      dplyr::select(time.cc, time.well, time.min, pine, wetland) %>% 
      data.frame())

ggpairs(site.wellcc %>% 
          dplyr::select(time.cc, time.well, time.min, pine, wetland) %>% 
          data.frame())
#time min pretty strongly correlated with time.cc and time.well
#verrrrry little pine

boot <- 100

occ.wellcc.boom.aic <- list()
occ.wellcc.call.aic <- list()
occ.wellcc.boom.mod <- list()
occ.wellcc.call.mod <- list()
for(i in 1:boot){
  
  #4. Spatial thinning----
  dat.wellcc <- site.2d %>%
    dplyr::filter(is.na(time.fire)) %>% 
    unique() %>% 
    group_by(cell) %>% 
    sample_n(1) %>% 
    ungroup() %>% 
    left_join(dat)
  
  table(dat.wellcc$boom, dat.wellcc$call)
  
  #5. Format into unmarked objects----
  hist.boom.wellcc <- dat.wellcc %>% 
    dplyr::select(ID, n, boom) %>% 
    spread(key = n, value = boom) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.call.wellcc <- dat.wellcc%>% 
    dplyr::select(ID, n, call) %>% 
    spread(key = n, value = call) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.Set.wellcc <- dat.wellcc %>% 
    dplyr::select(ID, n, sundiff) %>% 
    spread(key=n, value=sundiff) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.Doy.wellcc <- dat.wellcc %>% 
    dplyr::select(ID, n, yday) %>% 
    spread(key=n, value=yday) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.s2n1.wellcc <- dat.wellcc %>% 
    dplyr::select(ID, n, s2n.1) %>% 
    spread(key=n, value=s2n.1) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.s2n2.wellcc <- dat.wellcc %>% 
    dplyr::select(ID, n, s2n.2) %>% 
    spread(key=n, value=s2n.2) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  obs.cov.wellcc <- list(Set = hist.Set.wellcc, Doy = hist.Doy.wellcc, s2n1 = hist.s2n1.wellcc, s2n2 = hist.s2n2.wellcc)
  
  site.cov.wellcc <- dat.wellcc %>% 
    dplyr::select(ID, pine, wetland, time.min, time.dist, time.cc, time.well) %>% 
    unique() %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  df.boom.wellcc <- unmarkedFrameOccu(hist.boom.wellcc, siteCovs=site.cov.wellcc, obsCovs=obs.cov.wellcc)
  df.call.wellcc <- unmarkedFrameOccu(hist.call.wellcc, siteCovs=site.cov.wellcc, obsCovs=obs.cov.wellcc)
  
  #6. Model----
  #BOOM
  mod1 <- try(occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
               ~time.cc*time.well*pine,
               data=df.boom.wellcc))
  mod2 <- try(occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
               ~time.cc*pine + time.well*pine,
               data=df.boom.wellcc))
  mod3 <- try(occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
               ~time.cc*pine,
               data=df.boom.wellcc))
  mod4 <- try(occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
               ~time.well*pine,
               data=df.boom.wellcc))
  mod5 <- try(occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
               ~time.min*pine,
               data=df.boom.wellcc))
  
  rm(err)
  err <- try(aictab(list(mod1, mod2, mod3, mod4, mod5), sort=FALSE))
  
  if(inherits(err, "try-error")){
    next
  }
  
  #Save out AIC
  occ.wellcc.boom.aic[[i]] <- aictab(list(mod1, mod2, mod3, mod4, mod5), sort=FALSE) %>% 
    data.frame() %>% 
    mutate(boot=i)
  
  #Save out top model (hand picked from preview of results, not top AIC)
  occ.wellcc.boom.mod[[i]] <- summary(mod4)[['state']] %>% 
    mutate(boot=i,
           var=row.names(data.frame(mod4@estimates@estimates[["state"]]@estimates)))
  
  rm(list=ls(pattern="^mod"))
  
  print(paste0("COMPLETED BOOTSTRAP ", i))
  
}

#7. Collapse bootstrapped results----
#AIC----
aic.wellcc.boom <- rbindlist(occ.wellcc.boom.aic) %>% 
  group_by(Modnames) %>% 
  summarize(delta.mn=round(mean(Delta_AICc), 2),
            delta.sd=round(sd(Delta_AICc), 2),
            delta.ci=round(sd(Delta_AICc)/sqrt(100)*1.96, 2),
            wt.mn=round(mean(AICcWt), 2),
            wt.sd=round(sd(AICcWt), 2),
            wt.ci=round(sd(AICcWt)/sqrt(100)*1.96, 2),
            k=mean(K),
            aic.mn=round(mean(AICc), 2),
            aic.sd=round(sd(AICc), 2)) %>% 
  ungroup() %>% 
  mutate(aic=paste0(aic.mn, " (SD=", aic.sd, ")"),
         delta=paste0(delta.mn, " (SD=", delta.sd, ")"),
         wt=paste0(wt.mn, " (SD=", wt.sd, ")")) %>% 
  arrange(delta.mn)

aic.wellcc.boom.tbl <- aic.wellcc.boom %>%   
  dplyr::select(Modnames, k, aic, delta, wt)
View(aic.wellcc.boom.tbl)

ggplot(aic.wellcc.boom) +
  geom_bar(aes(x=Modnames, y=delta.mn), stat="identity") +
  geom_errorbar(aes(x=Modnames, ymin=delta.mn-delta.ci, ymax=delta.mn+delta.ci), width=0.2, position=position_dodge(0.9))

#Coefficients----
mod.wellcc.boom <- rbindlist(occ.wellcc.boom.mod) %>% 
  group_by(var) %>% 
  summarize(est.mn=mean(Estimate),
            est.sd=sd(Estimate),
            se.mn=mean(SE),
            se.sd=sd(SE)) %>% 
  ungroup()
mod.wellcc.boom

#H. Does # of disturbances & # of types matter####

#1. Wrangle----
site.many <- site.dat %>% 
  group_by(StationKey, DepYear, cell, Latitude, Longitude, pine, wetland, types, count) %>% 
  summarize(time=min(time)) %>% 
  ungroup()

#2. Check for vif & covariation----
vif(site.many %>% 
      dplyr::select(time, pine, wetland, types, count) %>% 
      data.frame())

cor(site.many %>% 
      dplyr::select(time, pine, wetland, types, count) %>% 
      data.frame())
#Types & count covary. Pick one?

#3. Spatial thinning----
boot <- 10 

many.boom.aic <- list()
for(i in 1:boot){
  
  dat.many <- site.many %>% 
    group_by(cell) %>% 
    sample_n(1) %>% 
    ungroup() %>% 
    left_join(dat)
  
  #4. Format into unmarked objects----
  hist.many.boom <- dat.many %>% 
    dplyr::select(ID, n, boom) %>% 
    spread(key = n, value = boom) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.many.call <- dat.many%>% 
    dplyr::select(ID, n, call) %>% 
    spread(key = n, value = call) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.many.set <- dat.many %>% 
    dplyr::select(ID, n, sundiff) %>% 
    spread(key = n, value = sundiff) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame() 
  
  hist.many.Doy <- dat.many %>% 
    dplyr::select(ID, n, yday) %>% 
    spread(key=n, value=yday) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.many.psd1 <- dat.many %>% 
    dplyr::select(ID, n, psd.1) %>% 
    spread(key=n, value=psd.1) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.many.s2n1 <- dat.many %>% 
    dplyr::select(ID, n, s2n.1) %>% 
    spread(key=n, value=s2n.1) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  hist.many.s2n2 <- dat.many %>% 
    dplyr::select(ID, n, s2n.2) %>% 
    spread(key=n, value=s2n.2) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  obs.cov.many <- list(Set = hist.many.set, Doy = hist.many.Doy, psd1 = hist.many.psd1, s2n2 = hist.many.s2n2, s2n1 = hist.many.s2n1)
  
  site.cov.many <- dat.many %>% 
    dplyr::select(ID, pine, wetland, count, types, time) %>% 
    unique() %>% 
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    data.frame()
  
  df.boom.many <- unmarkedFrameOccu(hist.many.boom, siteCovs=site.cov.many, obsCovs=obs.cov.many)
  df.call.many <- unmarkedFrameOccu(hist.many.call, siteCovs=site.cov.many, obsCovs=obs.cov.many)
  
  #5. Model----
  cov1 <- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
               ~pine*time + time*types,
               data=df.boom.many)
  cov2<- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
              ~pine*time + types,
              data=df.boom.many)
  cov3<- occu(~ Set + Doy + s2n1 + I(s2n1^2) + s2n2 + I(s2n2^2)
              ~pine*time,
              data=df.boom.many)
  
  many.boom.aic[[i]] <- aictab(list(cov1, cov2, cov3), sort=FALSE) %>% 
    data.frame() %>% 
    mutate(boot=i)
  
  print(paste0("Finished bootstrap ", i))
  
}

#6. Collapse bootstrapped results----
many.boom <- rbindlist(many.boom.aic) %>% 
  group_by(Modnames) %>% 
  summarize(delta.mn=round(mean(Delta_AICc), 2),
            delta.sd=round(sd(Delta_AICc), 2),
            delta.ci=round(sd(Delta_AICc)/sqrt(100)*1.96, 2),
            wt.mn=round(mean(AICcWt), 2),
            wt.sd=round(sd(AICcWt), 2),
            wt.ci=round(sd(AICcWt)/sqrt(100)*1.96, 2),
            k=mean(K),
            aic.mn=round(mean(AICc), 2),
            aic.sd=round(sd(AICc), 2)) %>% 
  ungroup() %>% 
  mutate(aic=paste0(aic.mn, " (SD=", aic.sd, ")"),
         delta=paste0(delta.mn, " (SD=", delta.sd, ")"),
         wt=paste0(wt.mn, " (SD=", wt.sd, ")")) %>% 
  arrange(delta.mn)

many.boom.tbl <- many.boom %>%   
  dplyr::select(Modnames, k, aic, delta, wt)
View(many.boom.tbl)

ggplot(many.boom) +
  geom_bar(aes(x=Modnames, y=delta.mn), stat="identity") +
  geom_errorbar(aes(x=Modnames, ymin=delta.mn-delta.ci, ymax=delta.mn+delta.ci), width=0.2, position=position_dodge(0.9))
#NOPE

summary(cov5)
