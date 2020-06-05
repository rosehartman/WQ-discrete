library(tidyverse)
library(deltareportr)
library(mgcv)
library(lubridate)
library(hms)
library(sf)
library(stars)
require(patchwork)
require(geofacet)
require(gamm4)
require(dtplyr)

is.even <- function(x) as.integer(x) %% 2 == 0

# Load Delta Shapefile from Brian
Delta<-st_read("Delta Subregions")%>%
  filter(!SubRegion%in%c("South Bay", "San Francisco Bay", "San Pablo Bay", "Upper Yolo Bypass", 
                         "Upper Napa River", "Lower Napa River", "Carquinez Strait"))%>% # Remove regions outside our domain of interest
  dplyr::select(SubRegion)

# Load data
Data <- DeltaDater(Start_year = 1900, 
                   WQ_sources = c("EMP", "STN", "FMWT", "EDSM", "DJFMP", "SKT", "20mm", "Suisun", "Baystudy", "USBR", "USGS"), 
                   Variables = "Water quality", 
                   Regions = NULL)%>%
  filter(!is.na(Temperature) & !is.na(Datetime) & !is.na(Latitude) & !is.na(Longitude) & !is.na(Date))%>% #Remove any rows with NAs in our key variables
  filter(Temperature !=0)%>% #Remove 0 temps
  mutate(Temperature_bottom=if_else(Temperature_bottom>30, NA_real_, Temperature_bottom))%>% #Remove bad bottom temps
  filter(hour(Datetime)>=5 & hour(Datetime)<=20)%>% # Only keep data betwen 5AM and 8PM
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)%>% # Convert to sf object
  st_transform(crs=st_crs(Delta))%>% # Change to crs of Delta
  st_join(Delta, join=st_intersects)%>% # Add subregions
  filter(!is.na(SubRegion))%>% # Remove any data outside our subregions of interest
  mutate(Datetime = with_tz(Datetime, tz="America/Phoenix"), #Convert to a timezone without daylight savings time
         Date = with_tz(Date, tz="America/Phoenix"),
         Julian_day = yday(Date), # Create julian day variable
         Month_fac=factor(Month), # Create month factor variable
         Source_fac=factor(Source),
         Year_fac=factor(Year))%>% 
  mutate(Date_num = as.numeric(Date), # Create numeric version of date for models
         Time = as_hms(Datetime))%>% # Create variable for time-of-day, not date. 
  mutate(Time_num=as.numeric(Time))%>% # Create numeric version of time for models (=seconds since midnight)
  mutate_at(vars(Date_num, Longitude, Latitude, Time_num, Year, Julian_day), list(s=~(.-mean(., na.rm=T))/sd(., na.rm=T))) # Create centered and standardized versions of covariates

# Pull station locations for major monitoring programs
# This will be used to set a boundary for this analysis focused on well-sampled regions.
WQ_stations<-Data%>%
  st_drop_geometry()%>%
  filter(Source%in%c("FMWT", "STN", "SKT", "20mm", "EMP", "Suisun"))%>%
  group_by(StationID, Source, Latitude, Longitude)%>%
  summarise(N=n(), .groups="drop")%>% # Calculate how many times each station was sampled
  filter(N>50 & !StationID%in%c("20mm 918", "STN 918"))%>% # Only keep stations sampled >50 times when deciding which regions to retain. 
  # "20mm 918", "STN 918" are far south of the rest of the well-sampled sites and are not sampled year round, so we're removing them to exclude that far southern region
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)%>% # Convert to sf object
  st_transform(crs=st_crs(Delta))%>%
  st_join(Delta) # Add subregions

# Remove any subregions that do not contain at least one of these >50 samples stations from the major monitoring programs
Delta <- Delta%>%
  filter(SubRegion%in%unique(WQ_stations$SubRegion) | SubRegion=="Georgiana Slough") # Retain Georgiana Slough because it's surrounded by well-sampled regions
# Visualize sampling regions of major surveys

# Now filter data to only include this final set of subregions, and any stations outside the convex hull formed by the >50 samples stations from the major monitoring programs
Data<-Data%>%
  filter(SubRegion%in%unique(Delta$SubRegion))%>%
  st_join(WQ_stations%>%
            st_union()%>%
            st_convex_hull()%>% # Draws a hexagram or pentagram or similar around the outer-most points
            st_as_sf()%>%
            mutate(IN=TRUE),
          join=st_intersects)%>%
  filter(IN)%>%
  dplyr::select(-IN)%>%
  mutate(Group=if_else(is.even(Year), 1, 2))

# Old code to visualize the prior steps

ggplot()+
  geom_sf(data=Delta, aes(fill=SubRegion))+
  #geom_sf_label(data=st_centroid(Delta)%>%st_transform(crs=4326), aes(label=SubRegion))+
  geom_sf(data=WQ_stations%>%st_union()%>%st_convex_hull(), alpha=0.1, color="red", size=2)+
  geom_sf(data=WQ_stations)
#Give all datasets the same ending year
#max_date <- Data%>%group_by(Source)%>%summarise(Date=max(Date))%>%pull(Date)%>%min()
#Data <- filter(Data, Year<=year(max_date) & !(Source=="SKT" & Field_coords))


# Initial models ----------------------------------------------------------


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

# Check models
gam.check(model)
concurvity(model, full=TRUE) # Check for values over 0.8
#If concurvity is bad, run again with full=FALSE


# Models optimizing k values ----------------------------------------------

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

modell2 <- bam(Temperature ~ te(Year_s, Longitude_s, Latitude_s, Julian_day_s, d=c(1,2,1), bs=c("cr", "tp", "cc"), k=c(10, 25, 7)) + s(Time_num_s, k=5),
               data = Data, method="fREML", discrete=T, nthreads=4)
#AIC: 335874.6
#BIC: 344138

modell2b <- bam(Temperature ~ te(Date_num, Longitude_s, Latitude_s, d=c(1,2), bs=c("cr", "tp"), k=c(70, 25)) + s(Time_num_s, k=5),
                data = Data, method="fREML", discrete=T, nthreads=4)

modell2c <- bam(Temperature ~ te(Year_s, Longitude_s, Latitude_s, Julian_day_s, d=c(1,2,1), bs=c("cr", "tp", "cr"), k=c(10, 25, 7)) + s(Time_num_s, k=5),
                data = Data, method="fREML", discrete=T, nthreads=4)
#AIC: 336037.3
#BIC: 344957.1

modell3 <- bam(Temperature ~ te(Year_s, Longitude_s, Latitude_s, Julian_day_s, d=c(1,2,1), bs=c("cr", "tp", "cc"), k=c(15, 35, 7)) + s(Time_num_s, k=5),
               data = Data, method="fREML", discrete=T, nthreads=4)
#AIC: 331733.5
#BIC: 345038.7

modell4 <- bam(Temperature ~ te(Year_s, Longitude_s, Latitude_s, Julian_day_s, d=c(1,2,1), bs=c("cr", "tp", "cc"), k=c(10, 35, 7)) + s(Time_num_s, k=5),
               data = Data, method="fREML", discrete=T, nthreads=4)
#AIC: 335430.2
#BIC: 345330.9

modell5 <- bam(Temperature ~ te(Year_s, Longitude_s, Latitude_s, Julian_day_s, d=c(1,2,1), bs=c("cr", "tp", "cc"), k=c(15, 25, 7)) + s(Time_num_s, k=5),
               data = Data, method="fREML", discrete=T, nthreads=4)
#AIC: 332211.2
#BIC: 343416.9

modell6 <- bam(Temperature ~ te(Year_s, Longitude_s, Latitude_s, Julian_day_s, d=c(1,2,1), bs=c("cr", "tp", "cc"), k=c(15, 30, 7)) + s(Time_num_s, k=5),
               data = Data, method="fREML", discrete=T, nthreads=4)
#AIC: 331920.7
#BIC: 344391.6

modell7 <- bam(Temperature ~ te(Year_s, Longitude_s, Latitude_s, Julian_day_s, d=c(1,2,1), bs=c("cr", "tp", "cc"), k=c(20, 25, 7)) + s(Time_num_s, k=5),
               data = Data, method="fREML", discrete=T, nthreads=4)
#AIC: 326162.7
#BIC: 340341

modell8 <- bam(Temperature ~ te(Year_s, Longitude_s, Latitude_s, Julian_day_s, d=c(1,2,1), bs=c("cr", "tp", "cc"), k=c(25, 25, 7)) + s(Time_num_s, k=5),
               data = Data, method="fREML", discrete=T, nthreads=4)

modell9 <- bam(Temperature ~ te(Year_s, Longitude_s, Latitude_s, Julian_day_s, d=c(1,2,1), bs=c("cr", "tp", "cc"), k=c(35, 25, 7)) + s(Time_num_s, k=5),
               data = Data, method="fREML", discrete=T, nthreads=4)
#AIC: 314035.6
#BIC: 336742

# Now this is best! k-value for Year_fac has no effect. Changing this results in the exact same model
modellb <- bam(Temperature ~ te(Year_fac, Longitude_s, Latitude_s, Julian_day_s, d=c(1,2,1), bs=c("fs", "tp", "cc"), k=c(35, 25, 7)) + s(Time_num_s, k=5),
               data = Data, method="fREML", discrete=T, nthreads=4)

modellc <- bam(Temperature ~ Year_fac + te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 7), by=Year_fac) + s(Time_num_s, k=5),
               data = Data, method="fREML", discrete=T, nthreads=4)

modellc2a <- bam(Temperature ~ Year_fac + te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 15), by=Year_fac) + s(Time_num_s, k=5),
                 data = filter(Data, Group==1)%>%mutate(Year_fac=droplevels(Year_fac)), method="fREML", discrete=T, nthreads=4)
#AIC: 126279.9
#BIC: 158009.3

modellc3a <- bam(Temperature ~ Year_fac + te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(35, 15), by=Year_fac) + s(Time_num_s, k=5),
                 data = filter(Data, Group==1)%>%mutate(Year_fac=droplevels(Year_fac)), method="fREML", discrete=T, nthreads=4)
#AIC: 125299.8
#BIC: 161495.4

modellc4a <- bam(Temperature ~ Year_fac + te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20), by=Year_fac) + s(Time_num_s, k=5),
                 data = filter(Data, Group==1)%>%mutate(Year_fac=droplevels(Year_fac)), method="fREML", discrete=T, nthreads=4)
#AIC: 116576.5
#BIC: 157218.6

modellc5a <- bam(Temperature ~ Year_fac + te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 25), by=Year_fac) + s(Time_num_s, k=5),
                 data = filter(Data, Group==1)%>%mutate(Year_fac=droplevels(Year_fac)), method="fREML", discrete=T, nthreads=4, gc.level=2)

## This is by far the best model by BIC and AIC
modelld <- bam(Temperature ~ Year_fac + te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 7), m=2) + 
                 te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 7), by=Year_fac, m=1) + s(Time_num_s, k=5),
               data = Data, method="fREML", discrete=T, nthreads=4, select=T)

modelld2 <- bam(Temperature ~ Year_fac + t2(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs="tp", k=c(15, 5), m=2, full=T) + 
                  te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs="tp", k=c(25, 15), by=Year_fac, m=1) + s(Time_num_s, k=5),
                data = filter(Data, Group==1)%>%mutate(Year_fac=droplevels(Year_fac)), method="fREML", discrete=F, nthreads=4, select=T)

modelm2 <- gamm(Temperature ~ te(Year_s, Longitude_s, Latitude_s, Julian_day_s, d=c(1,2,1), bs=c("cr", "tp", "cc"), k=c(10, 15, 7)) + s(Time_num_s, k=5), random=list(Source=~1),
                data = Data, method="REML")

modelm2b <- bam(Temperature ~ t2(Year_s, Longitude_s, Latitude_s, Julian_day_s, d=c(1,2,1), bs=c("cr", "tp", "cc"), k=c(10, 15, 7)) + 
                  s(Time_num_s, k=5) + s(Source_fac, bs="re"), 
                data = Data, method="REML")

modelm2b2 <- bam(Temperature ~ t2(Year_s, Longitude_s, Latitude_s, Julian_day_s, d=c(1,2,1), bs=c("cr", "tp", "cc"), k=c(10, 25, 7)) + 
                   s(Time_num_s, k=5) + s(Source_fac, bs="re"), 
                 data = Data, method="REML")

modelm2b <- gamm4(Temperature ~ t2(Year_s, Longitude_s, Latitude_s, Julian_day_s, d=c(1,2,1), bs=c("cr", "tp", "cc"), k=c(10, 15, 7)) + s(Time_num_s, k=5), 
                  random=~(1|Source), data = Data, REML=TRUE, verbose=2)

modelm2_bottom <- gamm(Temperature_bottom ~ te(Year_s, Longitude_s, Latitude_s, Julian_day_s, d=c(1,2,1), bs=c("cr", "tp", "cc"), k=c(10, 15, 7)) + s(Time_num_s, k=5), random=list(Source=~1),
                       data = filter(Data, !is.na(Temperature_bottom)), method="REML")
save.image("~/WQ-discrete/Temperature QAQC.RData")


gam.check(modelm2$gam)
concurvity(modelm2$gam, full=TRUE)

#model 2 seems good
plot(modelm2$gam, all.terms=TRUE, residuals=TRUE, shade=TRUE)
#vis.gam



# DEM ---------------------------------------------------------------------


Delta<-st_read("EDSM_Subregions")%>%
  st_transform(crs=26910)%>%
  filter(!SubRegion%in%c("South Bay", "San Francisco Bay", "San Pablo Bay"))
Coords<-Data%>%
  dplyr::select(Latitude, Longitude, StationID)%>%
  distinct()%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)%>%
  st_transform(crs=26910)
DEM<-read_stars("~/dem_bay_delta_10m_20181128/dem_bay_delta_10m_20181128.tif")%>%
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
                  Julian_days=yday(ymd(paste("2001", c(1,4,7,10), "15", sep="-"))), #Jan, Apr, Jul, and Oct 15 for a non-leap year
                  Time_num=12*60*60, # 12PM x 60 seconds x 60 minutes
                  Source="none",
                  Source_gam_re=TRUE, #Was source included as a random effect in a smoother with bs="re"?
                  Variance="CI" # "CI" or "SE" Anything else won't return any variance
){
  
  # Create point locations on a grid for predictions
  Points<-st_make_grid(Delta_subregions, n=n)%>%
    st_as_sf(crs=st_crs(Delta_subregions))%>%
    st_join(Delta_water%>% # Joining a map of delta waterways (from my spacetools package) to ensure all these points are over water.
              dplyr::select(Shape_Area)%>%
              st_transform(crs=st_crs(Delta_subregions)))%>%
    filter(!is.na(Shape_Area))%>%
    st_join(Stations%>% # Applying the same approach we did to the full data: remove any pounts outside the convex hull formed by major survey stations sampled >50 times
              st_union()%>%
              st_convex_hull()%>%
              st_as_sf()%>%
              mutate(IN=TRUE),
            join=st_intersects)%>%
    filter(IN)%>%
    dplyr::select(-IN)%>%
    st_centroid()%>% # The prior grid was actually a set of polygons, this picks the center point of each
    st_transform(crs=4326)%>%
    st_coordinates()%>%
    as_tibble()%>%
    mutate(Location=1:nrow(.))%>%
    dplyr::select(Longitude=X, Latitude=Y, Location)
  
  # Create dataset for each year and season showing which subregions were sampled
  Data_effort <- Full_data%>%
    st_drop_geometry()%>%
    group_by(SubRegion, Season, Year)%>%
    summarise(N=n())%>%
    ungroup()%>%
    left_join(Delta_subregions, by="SubRegion")%>%
    dplyr::select(-geometry)
  
  
  # Create full dataset for predictions
  newdata<-expand.grid(Year= Years,
                       Location=1:nrow(Points),
                       Julian_day=Julian_days,
                       Time_num=Time_num,
                       Source_fac=Source)%>% # Create all combinations of predictor variables
    left_join(Points, by="Location")%>% #Add Lat/Longs to each location
    mutate(Latitude_s=(Latitude-mean(Full_data$Latitude, na.rm=T))/sd(Full_data$Latitude, na.rm=T), # Standardize each variable based on full dataset for model
           Longitude_s=(Longitude-mean(Full_data$Longitude, na.rm=T))/sd(Full_data$Longitude, na.rm=T),
           Year_s=(Year-mean(Full_data$Year, na.rm=T))/sd(Full_data$Year, na.rm=T),
           Julian_day_s = (Julian_day-mean(Full_data$Julian_day, na.rm=T))/sd(Full_data$Julian_day, na.rm=T),
           Time_num_s=(Time_num-mean(Full_data$Time_num, na.rm=T))/sd(Full_data$Time_num, na.rm=T),
           Year_fac=ordered(Year),
           Season=case_when(Julian_day<=80 | Julian_day>=356 ~ "Winter", # Create a variable for season
                            Julian_day>80 & Julian_day<=172 ~ "Spring",
                            Julian_day>173 & Julian_day<=264 ~ "Summer",
                            Julian_day>265 & Julian_day<=355 ~ "Fall"))%>%
    st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)%>% # Turn into sf object
    st_transform(crs=st_crs(Delta_subregions))%>% # transform to crs of Delta shapefile
    st_join(Delta_subregions, join = st_intersects)%>%
    filter(!is.na(SubRegion))%>% # Make sure all points are within our desired subregions
    left_join(Data_effort, by=c("SubRegion", "Season", "Year"))%>% # Use the Data_effort key created above to remove points in subregions that were not sampled that region, season, and year.
    filter(!is.na(N))
  
  pred<-predict(model, newdata=newdata, type="response", se.fit=TRUE, discrete=T, n.threads=4) # Create predictions
  
  if(Source_gam_re){
    pred$fit<-pred$fit+mean(model$coefficients[grep("Source_fac", names(model$coefficients))]) # Add mean effect of random intercept
  }
  
  # Add predictions to predicter dataset
  newdata<-newdata%>%
    mutate(Prediction=pred$fit)%>%
    {if(Variance=="CI"){
      mutate(., L95=pred$fit-pred$se.fit*1.96,
             U95=pred$fit+pred$se.fit*1.96)
    } else{
      .
    }}%>%
    {if(Variance=="SE"){
      mutate(., SE=pred$se.fit)
    } else{
      .
    }}%>%
    mutate(Date=as.Date(Julian_day, origin=as.Date(paste(Year, "01", "01", sep="-")))) # Create Date variable from Julian Day and Year
  
  return(newdata)
}

# Rasterizing -------------------------------------------------------------

# Function to rasterize season by season. Creates a 4D raster Latitude x Longitude x Year x Season (But Season only has 1 value for the season passed to the function)
Rasterize_season<-function(season, data, n, out_crs=4326){
  
  Years <- data%>%
    filter(Season==season)%>%
    pull(Year)%>%
    unique()
  
  #First rasterize year by year
  preds<-map(Years, function(x) st_rasterize(data%>%
                                               filter(Year==x & Season==season)%>%
                                               dplyr::select(Prediction), 
                                             template=st_as_stars(st_bbox(Delta), dx=diff(st_bbox(Delta)[c(1, 3)])/n, dy=diff(st_bbox(Delta)[c(2, 4)])/n, values = NA_real_))%>%
               st_warp(crs=out_crs))
  
  # Then bind all years together into 1 raster
  out <- exec(c, !!!preds, along=list(Year=Years, Season=season))
}

# Function to rasterize all seasons. Creates a 4D raster Latitude x Longitude x Year x Season (including all 4 seasons). 
# Unfortunately, this requires equal dimension values across all seasons. Since some seasons are missing earlier years, this doesn't work well. 
Rasterize_all <- function(data, out_crs=4326){
  
  preds<-map(set_names(unique(data$Season)), function(x) Rasterize_season(season=x, data=data, out_crs=out_crs))
  
  out <- exec(c, !!!preds, along="Season")
  
  return(out)
  
}

# Plot the rasters. Sizing optimized for 9 different years in predictions. 
raster_plot<-function(data, Years=unique(newdata$Year), labels="All"){
  ggplot()+
    geom_blank(data=tibble(Year=Years, Season=st_get_dimension_values(data, "Season")))+
    geom_stars(data=data)+
    facet_grid(Year~Season)+
    scale_fill_viridis_c(name="Temperature", na.value="white", breaks=seq(6,26,by=1), labels= function(x) ifelse((x/2)==as.integer(x/2), as.character(x), ""),
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
    theme(axis.text.x = element_text(angle=45, hjust=1), plot.margin = margin(40,0,0,0), strip.background=element_blank(),
          panel.grid=element_blank(), legend.position = c(0.5,1.05), legend.background = element_rect(color="black"))
}

# Surface temperature

newdata <- WQ_pred(modellb, Source_gam_re=FALSE) # Run predict function above on modelm2. 

# Rasterize each season
rastered_preds <- map(set_names(c("Winter", "Spring", "Summer", "Fall")), function(x) Rasterize_season(season=x, data=newdata, n=100))

# Plot each season
p<-map2(rastered_preds, c("Left", "None", "None", "Right"), ~raster_plot(data=.x, labels=.y))

# Use Patchwork to bind plots together. Won't look very good until the plot is saved
p2<-wrap_plots(p)+plot_layout(nrow=1, heights=c(1,1,1,1))

# Save plots
ggsave(plot=p2, filename="C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/figures/Rasterized predictions 5.21.20.png", device=png(), width=7, height=12, units="in")

# Do the same for Bottom temperature

newdata_bottom <- WQ_pred(modelm2_bottom$gam,
                          Full_data=filter(Data, !is.na(Temperature_bottom)))

rastered_preds_bottom <- map(set_names(c("Winter", "Spring", "Summer", "Fall")), function(x) Rasterize_season(season=x, data=newdata_bottom, n=100))

p_bottom<-map2(rastered_preds_bottom, c("Left", "None", "None", "Right"), ~raster_plot(data=.x, Years=unique(newdata_bottom$Year), labels=.y))

p2_bottom<-wrap_plots(p_bottom)+plot_layout(nrow=1, heights=c(1,1,1,1))

ggsave(plot=p2_bottom, filename="C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/figures/Rasterized predictions_bottom.png", device=png(), width=7, height=12, units="in")


# Plots by year -----------------------------------------------------------

mygrid <- data.frame(
  name = c("Cache Slough and Lindsey Slough", "Lower Sacramento River Ship Channel", "Liberty Island", "Suisun Marsh", "Middle Sacramento River", "Lower Cache Slough", "Steamboat and Miner Slough", "Upper Mokelumne River", "Lower Mokelumne River", "Georgiana Slough", "Sacramento River near Ryde", "Sacramento River near Rio Vista", "Grizzly Bay", "West Suisun Bay", "Mid Suisun Bay", "Honker Bay", "Confluence", "Lower Sacramento River", "San Joaquin River at Twitchell Island", "San Joaquin River at Prisoners Pt", "Disappointment Slough", "Lower San Joaquin River", "Franks Tract", "Holland Cut", "San Joaquin River near Stockton", "Mildred Island", "Middle River", "Old River", "Upper San Joaquin River", "Grant Line Canal and Old River", "Victoria Canal"),
  row = c(2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 8, 8, 8),
  col = c(4, 6, 5, 2, 8, 6, 7, 9, 9, 8, 7, 6, 2, 1, 2, 3, 4, 5, 6, 8, 9, 5, 6, 7, 9, 8, 8, 7, 9, 8, 7),
  code = c(" 1", " 2", " 3", " 8", " 4", " 5", " 6", " 7", " 9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "30", "29", "31"),
  stringsAsFactors = FALSE
)
#geofacet::grid_preview(mygrid)

newdata_year <- WQ_pred(modellc,
                        Full_data=Data, 
                        Julian_days = yday(ymd(paste("2001", 1:12, "15", sep="-"))),
                        Years=round(min(Data$Year):max(Data$Year)),
                        Variance="SE",
                        Source_gam_re=FALSE) 

Data_year<-Data%>%
  filter(hour(Time)<14 & hour(Time)>10)%>%
  lazy_dt()%>%
  group_by(Year, Month, Season, SubRegion)%>%
  summarize(SD=sd(Temperature), Temperature=mean(Temperature))%>%
  ungroup()%>%
  as_tibble()

newdata_sum<-newdata_year%>%
  mutate(Var=SE^2,
         Month=month(Date))%>%
  lazy_dt()%>%
  group_by(Year, Month, SubRegion)%>%
  summarise(Temperature=mean(Prediction), SE=sqrt(sum(Var)/(n()^2)))%>%
  ungroup()%>%
  as_tibble()%>%
  mutate(L95=Temperature-1.96*SE,
         U95=Temperature+1.96*SE)

# Plot by Season for 1 subregion
ggplot(filter(newdata_sum, SubRegion=="Confluence"))+
  geom_ribbon(aes(x=Year, ymin=L95, ymax=U95), fill="darkorchid4", alpha=0.5)+
  geom_line(aes(x=Year, y=Temperature))+
  geom_pointrange(data=filter(Data_year, SubRegion=="Confluence"), aes(x=Year, y=Temperature, ymin=Temperature-SD, ymax=Temperature+SD))+
  facet_grid(~Month)+
  theme_bw()+
  theme(panel.grid=element_blank())

# Plot by Subregion for 1 season
mapyear<-function(month){
  ggplot(filter(newdata_sum, Month==month))+
    geom_ribbon(aes(x=Year, ymin=L95, ymax=U95), fill="firebrick3", alpha=0.5)+
    geom_line(aes(x=Year, y=Temperature), color="firebrick3")+
    geom_pointrange(data=filter(Data_year, Month==month), aes(x=Year, y=Temperature, ymin=Temperature-SD, ymax=Temperature+SD), size=0.5, alpha=0.4)+
    facet_geo(~SubRegion, grid=mygrid, labeller=label_wrap_gen())+
    theme_bw()+
    theme(panel.grid=element_blank(), axis.text.x = element_text(angle=45, hjust=1))
}

walk(1:12, function(x) ggsave(plot=mapyear(x), filename=paste0("C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/figures/Year predictions month ", x, "5.21.20.png"), device=png(), width=15, height=12, units="in"))
# QAQC by residuals -------------------------------------------------------


Data_qaqc<-Data%>%
  mutate(Residuals = residuals(modellb),
         Fitted=fitted(modellb))%>% # Fitted = model predictipn
  mutate(Flag=if_else(abs(Residuals)>sd(Residuals)*3, "Bad", "Good")) # Anything greater than 3 standard deviations of residuals is "bad"

p<-ggplot(data=Data_qaqc)+
  geom_point(aes(x=Temperature, y=Fitted, fill=Flag), shape=21)+
  geom_abline(intercept=0, slope=1, size=2)+
  ylab("Fitted temperature value from model")+
  xlab("Recorded temperature from surveys")+
  scale_fill_discrete(labels=c("Bad: Residuals > 3 SD", "Good: Residuals < 3 SD"))+
  theme_bw()+
  theme(panel.grid=element_blank(), legend.position=c(0.8,0.1), legend.background = element_rect(color="black"))

ggsave(plot=p, filename="C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/figures/Data qaqc.png", device=png(), width=7, height=7, units="in")


#Refit model with "good" data

modelm2_qaqc <- gamm(Temperature ~ te(Year_s, Longitude_s, Latitude_s, Julian_day_s, d=c(1,2,1), bs=c("cr", "tp", "cc"), k=c(10, 15, 7)) + s(Time_num_s, k=5), random=list(Source=~1),
                     data = filter(Data_qaqc, Flag=="Good"), method="REML")

Data_qaqc2<-Data_qaqc%>%
  filter(Flag=="Good")%>%
  mutate(Residuals = residuals(modelm2_qaqc$gam),
         Fitted=fitted(modelm2_qaqc$gam))%>% # Fitted = model predictipn
  mutate(Flag=if_else(abs(Residuals)>sd(Residuals)*3, "Bad", "Good")) # Anything greater than 3 standard deviations of residuals is "bad"

ggplot(data=Data_qaqc2)+
  geom_point(aes(x=Temperature, y=Fitted, fill=Flag), shape=21)+
  geom_abline(intercept=0, slope=1, size=2)


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
  scale_fill_viridis_c(name="Number of\nsamples per decade")+
  facet_grid(Decade~Season)+
  theme_bw()+
  theme(strip.background=element_blank(), axis.text.x = element_text(angle=45, hjust=1))


# Stratified cross-validation ---------------------------------------------

Data_split<-Data%>%
  mutate(Hour=hour(Time))%>%
  group_by(SubRegion, Year, Season, Hour)%>%
  mutate(Group=sample(1:10, 1, replace=T))%>%
  ungroup()


# Test autocorrelation ----------------------------------------------------

auto<-Data%>%
  filter(Group==1)%>%
  mutate(Resid=resid(modellc4a))%>%
  filter(Source!="EDSM" & !str_detect(Station, "EZ"))%>% # Remove EDSM and EZ stations because they're not fixed
  mutate(Station=paste(Source, Station))%>%
  group_by(Station)%>%
  mutate(N=n())%>%
  filter(N>10)%>%
  summarise(ACF=list(pacf(Resid, plot=F)), N=n(), ci=qnorm((1 + 0.95)/2)/sqrt(n()), .groups="drop")%>%
  rowwise()%>%
  mutate(lag=list(ACF$lag), acf=list(ACF$acf))%>%
  unnest(cols=c(lag, acf))%>%
  arrange(-N)

length(which(abs(auto$acf)>abs(auto$ci)))/nrow(auto)

## Only 6% exceed the CI, very close to the 5% you would expect with our chosen confidence level of 0.95 so I'm taking this as good evidence of no autocorrelation
  