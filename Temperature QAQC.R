library(tidyverse)
library(deltareportr)
library(mgcv)
library(lubridate)
library(hms)
library(sf)
library(stars)

Delta<-st_read("Delta Subregions")%>%
  filter(!SubRegion%in%c("South Bay", "San Francisco Bay", "San Pablo Bay", "Upper Yolo Bypass", 
                         "Upper Napa River", "Lower Napa River", "Carquinez Strait"))%>%
  dplyr::select(SubRegion)



Data <- DeltaDater(Start_year = 1900, 
                   WQ_sources = c("EMP", "STN", "FMWT", "EDSM", "DJFMP", "SKT", "20mm", "Suisun", "Baystudy", "USBR", "USGS"), 
                   Variables = "Water quality", 
                   Regions = NULL)%>%
  filter(!is.na(Temperature) & !is.na(Datetime) & !is.na(Latitude) & !is.na(Longitude) & !is.na(Date))%>%
  filter(hour(Datetime)>=5 & hour(Datetime)<=20)%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)%>%
  st_transform(crs=st_crs(Delta))%>%
  st_join(Delta, join=st_intersects)%>%
  filter(!is.na(SubRegion))%>%
  mutate(Datetime = with_tz(Datetime, tz="America/Phoenix"),
         Date = with_tz(Date, tz="America/Phoenix"),
         Julian_day = yday(Date),
         Month_fac=factor(Month))%>%
  mutate(Date_num = as.numeric(Date),
         Time = as_hms(Datetime))%>%
  mutate(Time_num=as.numeric(Time))%>%
  mutate_at(vars(Date_num, Longitude, Latitude, Time_num, Year, Julian_day), list(s=~(.-mean(., na.rm=T))/sd(., na.rm=T)))

# Pull station locations for major monitoring programs
# This will be used to set a boundary for this analysis focused on well-sampled regions.
WQ_stations<-Data%>%
  filter(Source%in%c("FMWT", "STN", "SKT", "20mm", "EMP", "Suisun"))%>%
  group_by(StationID, Source)%>%
  summarise(N=n())%>%
  filter(N>50 & !StationID%in%c("20mm 918", "STN 918"))%>% # These 2 stations are far south of the rest of the well-sampled sites and are not sampled year round, so we're removing them to exclude that far southern region
  st_join(Delta)

Delta <- Delta%>%
  filter(SubRegion%in%unique(WQ_stations$SubRegion) | SubRegion=="Georgiana Slough")
# Visualize sampling regions of major surveys

Data<-Data%>%
  filter(SubRegion%in%unique(Delta$SubRegion))%>%
  st_join(WQ_stations%>%
            st_union()%>%
            st_convex_hull()%>%
            st_as_sf()%>%
            mutate(IN=TRUE),
          join=st_intersects)%>%
  filter(IN)%>%
  dplyr::select(-IN)

ggplot()+
  geom_sf(data=Delta, aes(fill=SubRegion))+
  #geom_sf_label(data=st_centroid(Delta)%>%st_transform(crs=4326), aes(label=SubRegion))+
  geom_sf(data=WQ_stations%>%st_union()%>%st_convex_hull(), alpha=0.1, color="red", size=2)+
  geom_sf(data=WQ_stations)
#Give all datasets the same ending year
#max_date <- Data%>%group_by(Source)%>%summarise(Date=max(Date))%>%pull(Date)%>%min()
#Data <- filter(Data, Year<=year(max_date) & !(Source=="SKT" & Field_coords))

model <- gam(Temperature ~ t2(Date_num_s, Longitude_s, Latitude_s, d=c(1,2)) + t2(Julian_day, bs="cc") + poly(Time_num_s, 2),
             data = Data, method="REML")

modelb <- gam(Temperature ~ s(Longitude_s, Latitude_s) + s(Date_num_s) + ti(Date_num_s, Longitude_s, Latitude_s) + s(Julian_day, bs="cc") + poly(Time_num_s, 2),
              data = Data, method="REML")

modelc <- gam(Temperature ~ s(Longitude_s, Latitude_s) + s(Date_num_s) + ti(Date_num_s, Longitude_s, Latitude_s, d=c(1,2)) + s(Julian_day, bs="cc") + poly(Time_num_s, 2),
              data = Data, method="REML")

modeld <- gam(Temperature ~ s(Longitude_s, Latitude_s) + s(Date_num_s) + s(Julian_day, bs="cc") + poly(Time_num_s, 2),
              data = Data, method="REML")

modele <- gam(Temperature ~ s(Longitude_s, Latitude_s) + s(Year_s) + te(Julian_day_s, Time_num_s, bs="cc"),
              data = Data, method="REML")

modelf <- gam(Temperature ~ s(Longitude_s, Latitude_s) + s(Year_s) + ti(Year_s, Longitude_s, Latitude_s, d=c(1,2)) + te(Julian_day_s, Time_num_s, bs="cc"),
              data = Data, method="REML")

modelg <- gam(Temperature ~ te(Year_s, Longitude_s, Latitude_s, d=c(1,2)) + te(Julian_day_s, Time_num_s, bs="cc"),
              data = Data, method="REML")

modelh <- gam(Temperature ~ te(Year_s, Longitude_s, Latitude_s, d=c(1,2)) + te(Julian_day_s, Time_num_s, bs=c("cc", "cr"), d=c(1,1)),
              data = Data, method="REML")

modeli <- gamm(Temperature ~ te(Year_s, Longitude_s, Latitude_s, d=c(1,2)) + te(Julian_day_s, Time_num_s, bs=c("cc", "cr"), d=c(1,1)), random=list(Source=~1),
              data = Data, method="REML")

modelj <- gam(Temperature ~ te(Year_s, Longitude_s, Latitude_s, by=factor(Season), d=c(1,2)) + te(Julian_day_s, Time_num_s, bs=c("cc", "cr"), d=c(1,1)),
              data = Data, method="REML")

modelk <- gam(Temperature ~ te(Year_s, Longitude_s, Latitude_s, Month_fac, d=c(1,2,1), bs=c("cr", "tp", "fs")) + te(Julian_day_s, Time_num_s, bs=c("cc", "cr"), d=c(1,1)),
              data = Data, method="REML")

modell <- gam(Temperature ~ te(Year_s, Longitude_s, Latitude_s, Julian_day_s, d=c(1,2,1), bs=c("cr", "tp", "cc")) + s(Time_num_s),
              data = Data, method="REML")

modelm <- gamm(Temperature ~ te(Year_s, Longitude_s, Latitude_s, Julian_day_s, d=c(1,2,1), bs=c("cr", "tp", "cc")) + s(Time_num_s), random=list(Source=~1),
              data = Data, method="REML")

modelh <- gamm(Temperature ~ s(Longitude_s, Latitude_s) + s(Date_num_s) + ti(Date_num_s, Longitude_s, Latitude_s) + s(Julian_day, bs="cc") + poly(Time_num_s, 2),
               data = Data, correlation = corCAR1(form = ~ Date_num_s), method="REML")

gam.check(model)
concurvity(model, full=TRUE) # Check for values over 0.8
#If concurvity is bad, run again with full=FALSE

modeld2 <- gam(Temperature ~ s(Longitude_s, Latitude_s) + s(Date_num_s, k=120) + s(Julian_day, bs="cc") + poly(Time_num_s, 2),
               data = Data, method="REML")

modele2 <- gam(Temperature ~ s(Longitude_s, Latitude_s, k=80) + s(Year_s, k=15) + te(Julian_day_s, Time_num_s, bs="cc", k=10),
               data = Data, method="REML")

modelf2 <- gam(Temperature ~ s(Longitude_s, Latitude_s, k=80) + s(Year_s, k=15) + ti(Year_s, Longitude_s, Latitude_s) + te(Julian_day_s, Time_num_s, bs="cc", k=10),
               data = Data, method="REML")

modelb2 <- gam(Temperature ~ s(Longitude_s, Latitude_s, k=60) + s(Date_num_s, k=40) + ti(Date_num_s, Longitude_s, Latitude_s) + s(Julian_day, bs="cc", k=16) + poly(Time_num_s, 2),
               data = Data, method="REML")

modelb3 <- gam(Temperature ~ s(Longitude_s, Latitude_s, k=80) + s(Date_num_s, k=200) + ti(Date_num_s, Longitude_s, Latitude_s) + s(Julian_day, bs="cc", k=16) + poly(Time_num_s, 2),
               data = Data, method="REML")

modelg2 <- gam(Temperature ~ te(Year_s, Longitude_s, Latitude_s, d=c(1,2), k=c(15, 30)) + te(Julian_day_s, Time_num_s, bs="cc", k=10),
               data = Data, method="REML")

modelg3 <- gam(Temperature ~ te(Year_s, Longitude_s, Latitude_s, d=c(1,2), k=c(15, 60)) + te(Julian_day_s, Time_num_s, bs="cc", k=10),
               data = Data, method="REML")

modelh2 <- gam(Temperature ~ te(Year_s, Longitude_s, Latitude_s, d=c(1,2), k=c(15, 60)) + te(Julian_day_s, Time_num_s, Year_s, bs=c("cc", "cr", "cr"), k=c(8, 5, 15)),
               data = Data, method="REML")

modeli2 <- gamm(Temperature ~ te(Year_s, Longitude_s, Latitude_s, d=c(1,2), k=c(15, 30)) + te(Julian_day_s, Time_num_s, bs=c("cc", "cr"), d=c(1,1), k=c(10,10)), random=list(Source=~1),
               data = Data, method="REML")

modell2 <- gam(Temperature ~ te(Year_s, Longitude_s, Latitude_s, Julian_day_s, d=c(1,2,1), bs=c("cr", "tp", "cc"), k=c(15, 30, 10)) + s(Time_num_s),
              data = Data, method="REML")

modelm2 <- gamm(Temperature ~ te(Year_s, Longitude_s, Latitude_s, Julian_day_s, d=c(1,2,1), bs=c("cr", "tp", "cc"), k=c(10, 15, 7)) + s(Time_num_s, k=5), random=list(Source=~1),
               data = Data, method="REML")

modelm2_bottom <- gamm(Temperature_bottom ~ te(Year_s, Longitude_s, Latitude_s, Julian_day_s, d=c(1,2,1), bs=c("cr", "tp", "cc"), k=c(10, 15, 7)) + s(Time_num_s, k=5), random=list(Source=~1),
                data = filter(Data, !is.na(Temperature_bottom)), method="REML")

gam.check(modelm2$gam)
concurvity(modelm2$gam, full=TRUE)

#model 2 seems good
plot(model2, all.terms=TRUE, residuals=TRUE, shade=TRUE)
#vis.gam

Data<-Data%>%
  mutate(Residuals = residuals(modelg3),
         Fitted=fitted(modelg3))%>%
  mutate(Flag=if_else(abs(Residuals)>sd(Residuals)*3, "Bad", "Good"))

Space<-evaluate_smooth(modelg, "s(Longitude_s,Latitude_s)", dist=0.05)%>%
  mutate(Latitude = Latitude_s*sd(Data$Latitude)+mean(Data$Latitude),
         Longitude = Longitude_s*sd(Data$Longitude)+mean(Data$Longitude))

ggplot(Space)+
  geom_raster(aes(x=Longitude, y=Latitude, fill=est))+
  geom_contour(aes(x=Longitude, y=Latitude, z=est))+ 
  geom_sf(data=deltareportr::deltaregions%>%st_transform(crs=4326), aes(color=Stratum))+
  scale_fill_distiller(palette = "RdBu", type = "div")+
  guides(fill = guide_colourbar(title = "Effect",
                                direction = "vertical",
                                barheight = grid::unit(0.25, "npc")))
ageom_ribbon(data=Year, aes(x=Year_s, ymax=est+2*se, ymin=est-2*se), alpha=0.1)

ggplot(data=Data)+
  geom_point(aes(x=Temperature, y=Fitted, fill=Flag), shape=21)


Delta<-st_read("EDSM_Subregions")%>%
  st_transform(crs=26910)%>%
  filter(!SubRegion%in%c("South Bay", "San Francisco Bay", "San Pablo Bay"))
Coords<-Data%>%
  dplyr::select(Latitude, Longitude, StationID)%>%
  distinct()%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)%>%
  st_transform(crs=26910)
DEM<-read_stars("~/dem_bay_delta_10m_20181128/dem_bay_delta_10m_20181128.tif", proxy=TRUE)%>%
  st_crop(Delta)
Coords_joined<-aggregate(DEM, Coords, function(x) mean(x))

## Need to use velox package (when fixed) or rgis package

#DEM<-read_stars("~/dem_bay_delta_10m_20181128/dem_bay_delta_10m_20181128.tif", rasterIO=list(nXSize=18527, nYSize=16660, nBufXSize = 1852, nBufYSize = 1666))%>%
#  st_crop(Delta)
#saveRDS(DEM, file="DEM_cropped.rds", compress="xz")

Delta_water <- spacetools::Delta%>%
  st_transform(crs=st_crs(Delta))%>%
  st_crop(Delta)%>%
  st_rasterize(.,options="ALL_TOUCHED=TRUE")%>%
  st_join(Delta)%>%
  mutate(Include=if_else(is.na(SubRegion), TRUE, FALSE))

n=100

Points<-st_make_grid(Delta, n=n)%>%
  st_as_sf(crs=st_crs(Delta))%>%
  st_join(spacetools::Delta%>%
            dplyr::select(Shape_Area)%>%
            st_transform(crs=st_crs(Delta)))%>%
  filter(!is.na(Shape_Area))%>%
  st_join(WQ_stations%>%
            st_union()%>%
            st_convex_hull()%>%
            st_as_sf()%>%
            mutate(IN=TRUE),
          join=st_intersects)%>%
  filter(IN)%>%
  dplyr::select(-IN)%>%
  st_centroid()%>%
  st_transform(crs=4326)%>%
  st_coordinates()%>%
  as_tibble()%>%
  mutate(Location=1:nrow(.))%>%
  dplyr::select(Longitude=X, Latitude=Y, Location)

Data_effort <- Data%>%
  st_drop_geometry()%>%
  group_by(SubRegion, Season, Year)%>%
  summarise(N=n())%>%
  ungroup()%>%
  left_join(Delta, by="SubRegion")%>%
  dplyr::select(-geometry)

newdata<-expand.grid(Year_s= seq(min(Data$Year_s)+0.2, max(Data$Year_s)-0.2, length.out=9),
                     Location=1:nrow(Points),
                     Julian_day_s=seq(min(Data$Julian_day_s), max(Data$Julian_day_s), length.out=5)[1:4],# min and max are basically the same so excluding the max
                     Time_num_s=0)%>%
  left_join(Points, by="Location")%>%
  mutate(Latitude_s=(Latitude-mean(Data$Latitude, na.rm=T))/sd(Data$Latitude, na.rm=T),
         Longitude_s=(Longitude-mean(Data$Longitude, na.rm=T))/sd(Data$Longitude, na.rm=T),
         Year=round(Year_s*sd(Data$Year)+mean(Data$Year)),
         Julian_day = Julian_day_s*sd(Data$Julian_day, na.rm=T)+mean(Data$Julian_day, na.rm=T),
         Season=case_when(Julian_day<=80 | Julian_day>=356 ~ "Winter",
                          Julian_day>80 & Julian_day<=172 ~ "Spring",
                          Julian_day>173 & Julian_day<=264 ~ "Summer",
                          Julian_day>265 & Julian_day<=355 ~ "Fall"))%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)%>%
  st_transform(crs=st_crs(Delta))%>%
  st_join(Delta, join = st_intersects)%>%
  filter(!is.na(SubRegion))%>%
  left_join(Data_effort, by=c("SubRegion", "Season", "Year"))%>%
  filter(!is.na(N))

pred<-predict(modelm2$gam, newdata=newdata, type="response", se.fit=TRUE)

newdata<-newdata%>%
  mutate(Prediction=pred$fit,
         L95=pred$fit-pred$se.fit*1.96,
         U95=pred$fit+pred$se.fit*1.96)%>%
  mutate(Date=as.Date(Julian_day, origin=as.Date(paste(Year, "01", "01", sep="-"))))

bounds<-st_bbox(Delta%>%st_transform(crs=4326))

pred_plot <- function(data, season){
  ggplot(filter(data, Season==season))+
    geom_point(aes(x=Longitude, y=Latitude, color=Prediction))+
    ggtitle(paste(month(unique(filter(data, Season==season)$Date), label=TRUE), day(unique(filter(data, Season==season)$Date)), ":", season))+
    facet_wrap(~Year)+
    scale_colour_viridis_c()+
    theme_bw()+
    theme(strip.background = element_blank(), plot.title = element_text(hjust=0.5))
}

p<-map(set_names(c("Winter", "Spring", "Summer", "Fall")), ~pred_plot(newdata, .))

#walk(set_names(c("Winter", "Spring", "Summer", "Fall")), ~ggsave(plot=p[[.]], filename=paste0("C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/Prediction", ., " modelm2.png"),
#                                                                                              dpi=300, width=7, height=7, units="in", device="png"))

ggplot(newdata)+
  geom_tile(aes(x=Longitude, y=Latitude, fill=Prediction), width=(bounds["xmax"]-bounds["xmin"])/n, height=(bounds["ymax"]-bounds["ymin"])/n)+
  facet_wrap(~Year)+
  scale_fill_viridis_c()

# Sampling effort
Data_effort <- Data%>%
  st_drop_geometry()%>%
  mutate(Decade=floor(Year/10)*10)%>%
  group_by(SubRegion, Season, Decade)%>%
  summarise(N=n())%>%
  ungroup()%>%
  left_join(Delta, by="SubRegion")%>%
  st_as_sf()

ggplot(Data_effort)+
  geom_sf(aes(fill=N))+
  scale_fill_viridis_c(name="Number of\nsamples")+
  facet_grid(Decade~Season)+
  theme_bw()+
  theme(strip.background=element_blank(), axis.text.x = element_text(angle=45, hjust=1))

test<-st_rasterize(newdata%>%
                     dplyr::select(Year, Julian_day, Prediction, Latitude, Longitude)%>%
                     filter(Year==max(Year) & Julian_day==max(Julian_day))%>%
                     dplyr::select(Prediction), 
                   template=st_as_stars(st_bbox(Delta), dx=diff(st_bbox(Delta)[c(1, 3)])/n, dy=diff(st_bbox(Delta)[c(2, 4)])/n, values = NA_real_))%>%
  st_warp(crs=4326)

ggplot()+
  geom_stars(data=test)+
  scale_fill_viridis_c(na.value="white")+
  coord_equal()+
  ylab("Latitude")+
  xlab("Longitude")+
  theme_bw()
