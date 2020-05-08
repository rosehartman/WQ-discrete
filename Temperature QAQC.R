library(tidyverse)
library(deltareportr)
library(mgcv)
library(lubridate)
library(hms)
library(sf)
library(stars)
require(patchwork)

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



# DEM ---------------------------------------------------------------------


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


# Model predictions -------------------------------------------------------------

WQ_pred<-function(model,
                  Full_data=Data,
                  Delta_subregions=Delta,
                  Delta_water=spacetools::Delta,
                  Stations = WQ_stations,
                  n=100, 
                  Years=round(seq(min(Full_data$Year)+2, max(Full_data$Year)-2, length.out=9)),
                  Julian_days=seq(min(Full_data$Julian_day), max(Full_data$Julian_day), length.out=5)[1:4],
                  Time_num=0){
  
  Points<-st_make_grid(Delta_subregions, n=n)%>%
    st_as_sf(crs=st_crs(Delta_subregions))%>%
    st_join(Delta_water%>%
              dplyr::select(Shape_Area)%>%
              st_transform(crs=st_crs(Delta_subregions)))%>%
    filter(!is.na(Shape_Area))%>%
    st_join(Stations%>%
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
  
  Data_effort <- Full_data%>%
    st_drop_geometry()%>%
    group_by(SubRegion, Season, Year)%>%
    summarise(N=n())%>%
    ungroup()%>%
    left_join(Delta_subregions, by="SubRegion")%>%
    dplyr::select(-geometry)
  
  newdata<-expand.grid(Year= Years,
                       Location=1:nrow(Points),
                       Julian_day=Julian_days,
                       Time_num=0)%>%
    left_join(Points, by="Location")%>%
    mutate(Latitude_s=(Latitude-mean(Full_data$Latitude, na.rm=T))/sd(Full_data$Latitude, na.rm=T),
           Longitude_s=(Longitude-mean(Full_data$Longitude, na.rm=T))/sd(Full_data$Longitude, na.rm=T),
           Year_s=(Year-mean(Full_data$Year, na.rm=T))/sd(Full_data$Year, na.rm=T),
           Julian_day_s = (Julian_day-mean(Full_data$Julian_day, na.rm=T))/sd(Full_data$Julian_day, na.rm=T),
           Season=case_when(Julian_day<=80 | Julian_day>=356 ~ "Winter",
                            Julian_day>80 & Julian_day<=172 ~ "Spring",
                            Julian_day>173 & Julian_day<=264 ~ "Summer",
                            Julian_day>265 & Julian_day<=355 ~ "Fall"),
           Time_num_s=(Time_num-mean(Full_data$Time_num, na.rm=T))/sd(Full_data$Time_num, na.rm=T))%>%
    st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)%>%
    st_transform(crs=st_crs(Delta_subregions))%>%
    st_join(Delta_subregions, join = st_intersects)%>%
    filter(!is.na(SubRegion))%>%
    left_join(Data_effort, by=c("SubRegion", "Season", "Year"))%>%
    filter(!is.na(N))
  
  pred<-predict(model, newdata=newdata, type="response", se.fit=TRUE)
  
  newdata<-newdata%>%
    mutate(Prediction=pred$fit,
           L95=pred$fit-pred$se.fit*1.96,
           U95=pred$fit+pred$se.fit*1.96)%>%
    mutate(Date=as.Date(Julian_day, origin=as.Date(paste(Year, "01", "01", sep="-"))))
  
  return(newdata)
}

newdata <- WQ_pred(modelm2$gam)

# Rasterizing -------------------------------------------------------------


Rasterize_season<-function(season, data, n, out_crs=4326){
  Years <- data%>%
    filter(Season==season)%>%
    pull(Year)%>%
    unique()
  
  preds<-map(Years, function(x) st_rasterize(data%>%
                                               filter(Year==x & Season==season)%>%
                                               dplyr::select(Prediction), 
                                             template=st_as_stars(st_bbox(Delta), dx=diff(st_bbox(Delta)[c(1, 3)])/n, dy=diff(st_bbox(Delta)[c(2, 4)])/n, values = NA_real_))%>%
               st_warp(crs=out_crs))
  
  out <- exec(c, !!!preds, along=list(Year=Years, Season=season))
}

Rasterize_all <- function(data, out_crs=4326){
  
  preds<-map(set_names(unique(data$Season)), function(x) Rasterize_season(season=x, data=data, out_crs=out_crs))
  
  out <- exec(c, !!!preds, along="Season")
  
  return(out)
  
}

raster_plot<-function(data, Years=unique(newdata$Year), labels="All"){
  ggplot()+
    geom_blank(data=tibble(Year=Years, Season=st_get_dimension_values(data, "Season")))+
    geom_stars(data=data)+
    facet_grid(Year~Season)+
    scale_fill_viridis_c(name="Temperature", na.value="white", breaks=seq(6,26,by=2),
                         guide = guide_colorbar(direction="horizontal", title.position = "top", barwidth = 4, ticks.linewidth = 2,
                                                barheight=0.4, title.hjust=0.5, label.position="bottom", label.theme=element_text(size=8), 
                                                title.theme=element_text(size=10)))+
    coord_sf()+
    ylab("Latitude")+
    xlab("Longitude")+
    theme_bw()+
    {if(labels%in%c("None", "Right")){
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())
    }}+
    {if(labels%in%c("None", "Left")){
      theme(strip.text.y=element_blank())
    }}+
    theme(axis.text.x = element_text(angle=45, hjust=1), plot.margin = margin(30,0,0,0), strip.background=element_blank(),
          panel.grid=element_blank(), legend.position = c(0.5,1.06), legend.background = element_rect(color="black"))
}

rastered_preds <- map(set_names(c("Winter", "Spring", "Summer", "Fall")), function(x) Rasterize_season(season=x, data=newdata, n=100))

p<-map2(rastered_preds, c("Left", "None", "None", "Right"), ~raster_plot(data=.x, labels=.y))

p2<-wrap_plots(p)+plot_layout(nrow=1, heights=c(1,1,1,1))

ggsave(plot=p2, filename="C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/Rasterized predictions.png", device=png(), width=7, height=12, units="in")

# Bottom temperature

newdata_bottom <- WQ_pred(modelm2_bottom$gam,
                          Full_data=filter(Data, !is.na(Temperature_bottom)),
                          Delta_subregions=Delta,
                          Delta_water=spacetools::Delta,
                          Stations = WQ_stations,
                          n=100,
                          Time_num=0)

rastered_preds_bottom <- map(set_names(c("Winter", "Spring", "Summer", "Fall")), function(x) Rasterize_season(season=x, data=newdata_bottom, n=100))

p_bottom<-map2(rastered_preds_bottom, c("Left", "None", "None", "Right"), ~raster_plot(data=.x, Years=unique(newdata_bottom$Year), labels=.y))

p2_bottom<-wrap_plots(p_bottom)+plot_layout(nrow=1, heights=c(1,1,1,1))

# Data effort -------------------------------------------------------------

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

