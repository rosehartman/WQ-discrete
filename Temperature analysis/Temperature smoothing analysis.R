require(gstat)
require(sp)
require(spacetime)
library(tidyverse)
library(discretewq)
library(mgcv)
library(lubridate)
library(hms)
library(sf)
library(stars)
require(patchwork)
require(geofacet)
require(gamm4)
require(dtplyr)
require(scales)

### TODO

# 1) Redo all models now that year issue is fixed in the dataset
# 2) Calculate RMSE and pearson's correlation coefficients from CV results


# Data preparation --------------------------------------------------------

is.even <- function(x) as.integer(x) %% 2 == 0

# Load Delta Shapefile from Brian
Delta<-st_read("Delta Subregions")%>%
  filter(!SubRegion%in%c("South Bay", "San Francisco Bay", "San Pablo Bay", "Upper Yolo Bypass", 
                         "Upper Napa River", "Lower Napa River", "Carquinez Strait"))%>% # Remove regions outside our domain of interest
  dplyr::select(SubRegion)

# Load data
Data <- wq()%>%
  filter(!is.na(Temperature) & !is.na(Datetime) & !is.na(Latitude) & !is.na(Longitude) & !is.na(Date))%>% #Remove any rows with NAs in our key variables
  filter(Temperature !=0)%>% #Remove 0 temps
  mutate(Temperature_bottom=if_else(Temperature_bottom>30, NA_real_, Temperature_bottom))%>% #Remove bad bottom temps
  filter(hour(Datetime)>=5 & hour(Datetime)<=20)%>% # Only keep data between 5AM and 8PM
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
  mutate(Time_num=as.numeric(Time)) # Create numeric version of time for models (=seconds since midnight)


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
  mutate(Group=if_else(is.even(Year), 1, 2))%>%
  mutate_at(vars(Date_num, Longitude, Latitude, Time_num, Year, Julian_day), list(s=~(.-mean(., na.rm=T))/sd(., na.rm=T))) # Create centered and standardized versions of covariates

saveRDS(Data, file="Temperature analysis/Discrete Temp Data.Rds")


# Model selection ---------------------------------------------------------


# Chose separate smoothers for each year in order to ensure the most accurate predictions since temperatures fluctuate year-to-year
# Tried including a global smoother for lat, long, & julian_day, but ran into issues with curvilinearity.
# Optimized k-values using BIC comparisons on models fit to the even years of the dataset as follows: 

# Gavin Simpson recommends using AIC, not BIC https://stackoverflow.com/questions/59825442/get-the-aic-or-bic-citerium-from-a-gamm-gam-and-lme-models-how-in-mgcv-and-h

## New best
modellda <- bam(Temperature ~ Year_fac + te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20), by=Year_fac) + 
                  te(Time_num_s, Julian_day_s, bs=c("tp", "cc"), k=c(5, 12)),
                 data = filter(Data, Group==1)%>%mutate(Year_fac=droplevels(Year_fac)), method="fREML", discrete=T, nthreads=8)

#AIC: 120728.7
#BIC: 162087.7
## New best

modellda.2 <- bam(Temperature ~ Year_fac + te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20), by=Year_fac) + 
                  te(Time_num_s, Julian_day_s, bs=c("tp", "cc"), k=c(5, 12)), select=TRUE,
                data = filter(Data, Group==1)%>%mutate(Year_fac=droplevels(Year_fac)), method="fREML", discrete=T, nthreads=8)
#AIC: 120729.9
#BIC: 161957.2

modelld2a <- bam(Temperature ~ Year_fac + te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 10), by=Year_fac) + 
                   te(Time_num_s, Julian_day_s, bs=c("tp", "cc"), k=c(5, 12)),
                 data = filter(Data, Group==1)%>%mutate(Year_fac=droplevels(Year_fac)), method="fREML", discrete=T, nthreads=8)

#AIC: 141331.1
#BIC: 163732.6

modelld3a <- bam(Temperature ~ Year_fac + te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(15, 20), by=Year_fac) + 
                   te(Time_num_s, Julian_day_s, bs=c("tp", "cc"), k=c(5, 12)),
                 data = filter(Data, Group==1)%>%mutate(Year_fac=droplevels(Year_fac)), method="fREML", discrete=T, nthreads=8)

#AIC: 122571.1
#BIC: 155211

modelld4a <- bam(Temperature ~ Year_fac + te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20), by=Year_fac) + 
                  te(Time_num_s, Julian_day_s, bs=c("tp", "cc"), k=c(5, 6)),
                data = filter(Data, Group==1)%>%mutate(Year_fac=droplevels(Year_fac)), method="fREML", discrete=T, nthreads=8)

#AIC: 120732.4
#BIC: 162066.7

modelld5a <- bam(Temperature ~ Year_fac + te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(15, 20), by=Year_fac) + 
                   te(Time_num_s, Julian_day_s, bs=c("tp", "cc"), k=c(5, 6)),
                 data = filter(Data, Group==1)%>%mutate(Year_fac=droplevels(Year_fac)), method="fREML", discrete=T, nthreads=8)
#AIC: 122602.9
#BIC: 155149.3

modelld6a <- bam(Temperature ~ Year_fac + te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(15, 10), by=Year_fac) + 
                   te(Time_num_s, Julian_day_s, bs=c("tp", "cc"), k=c(5, 6)),
                 data = filter(Data, Group==1)%>%mutate(Year_fac=droplevels(Year_fac)), method="fREML", discrete=T, nthreads=3)
#AIC: 142605.2
#BIC: 159811.1

modelld7a <- bam(Temperature ~ Year_fac + te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(10, 20), by=Year_fac) + 
                   te(Time_num_s, Julian_day_s, bs=c("tp", "cc"), k=c(5, 6)),
                 data = filter(Data, Group==1)%>%mutate(Year_fac=droplevels(Year_fac)), method="fREML", discrete=T, nthreads=3)

#AIC: 123986.6
#BIC: 149802.3

modelld8a <- bam(Temperature ~ Year_fac + te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(5, 20), by=Year_fac) + 
                   te(Time_num_s, Julian_day_s, bs=c("tp", "cc"), k=c(5, 6)),
                 data = filter(Data, Group==1)%>%mutate(Year_fac=droplevels(Year_fac)), method="fREML", discrete=T, nthreads=3)
#AIC: 132112.5
#BIC: 147318.7

# Final model -------------------------------------------------------------

require(googleComputeEngineR)
require(googleCloudStorageR)

vm <- gce_vm(template = "rstudio", zone="us-west1-a",
             name = "rstudio-server",
             username = "", password = "",  predefined_type = "e2-highmem-16",
             dynamic_image = "gcr.io/gcer-public/persistent-rstudio",
             disk_size_gb=100)
gce_set_metadata(list(GCS_SESSION_BUCKET = "discretewq"), vm)

#Data used to fit the model are stored in "Simplified data.Rds" as Data_simp


modelld <- bam(Temperature ~ Year_fac + te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20), by=Year_fac) + 
                  te(Time_num_s, Julian_day_s, bs=c("tp", "cc"), k=c(5, 12)),
                data = Data, method="fREML", discrete=T, nthreads=16)


# Data prediction ---------------------------------------------------------

WQ_pred<-function(Full_data=Data,
                  Delta_subregions=Delta,
                  Delta_water=spacetools::Delta,
                  Stations = WQ_stations,
                  n=100, 
                  Years=round(seq(min(Full_data$Year)+2, max(Full_data$Year)-2, length.out=9)),
                  Julian_days=yday(ymd(paste("2001", c(1,4,7,10), "15", sep="-"))), #Jan, Apr, Jul, and Oct 15 for a non-leap year
                  Time_num=12*60*60 # 12PM x 60 seconds x 60 minutes
){
  
  # Create point locations on a grid for predictions
  Points<-st_make_grid(Delta_subregions, n=n)%>%
    st_as_sf(crs=st_crs(Delta_subregions))%>%
    st_join(Delta_water%>% # Joining a map of delta waterways (from my spacetools package) to ensure all these points are over water.
              dplyr::select(Shape_Area)%>%
              st_transform(crs=st_crs(Delta_subregions)))%>%
    filter(!is.na(Shape_Area))%>%
    select(-Shape_Area)%>%
    distinct()%>%
    st_join(Stations%>% # Applying the same approach we did to the full data: remove any points outside the convex hull formed by major survey stations sampled >50 times
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
                       Time_num=Time_num)%>% # Create all combinations of predictor variables
    left_join(Points, by="Location")%>% #Add Lat/Longs to each location
    mutate(Latitude_s=(Latitude-mean(Full_data$Latitude, na.rm=T))/sd(Full_data$Latitude, na.rm=T), # Standardize each variable based on full dataset for model
           Longitude_s=(Longitude-mean(Full_data$Longitude, na.rm=T))/sd(Full_data$Longitude, na.rm=T),
           Year_s=(Year-mean(Full_data$Year, na.rm=T))/sd(Full_data$Year, na.rm=T),
           Julian_day_s = (Julian_day-mean(Full_data$Julian_day, na.rm=T))/sd(Full_data$Julian_day, na.rm=T),
           Time_num_s=(Time_num-mean(Full_data$Time_num, na.rm=T))/sd(Full_data$Time_num, na.rm=T),
           Year_fac=factor(Year),
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
  return(newdata)
}

newdata_year <- WQ_pred(Full_data=Data, 
                        Julian_days = yday(ymd(paste("2001", 1:12, "15", sep="-"))),
                        Years=round(min(Data$Year):max(Data$Year)))

#saveRDS(newdata_year, file="Temperature analysis/Prediction Data.Rds")

# Perform in the cloud

modellc4_predictions<-predict(modellc4, newdata=newdata_year, type="response", se.fit=TRUE, discrete=T, n.threads=16) # Create predictions
modelld_predictions<-predict(modelld, newdata=newdata_year, type="response", se.fit=TRUE, discrete=T, n.threads=16) # Create predictions

# Predictions stored as "modellc4_predictions.Rds"

newdata_year<-readRDS("Temperature analysis/model outputs and validations/Prediction Data.Rds")
modellc4_predictions<-readRDS("Temperature analysis/model outputs and validations/modellc4_predictions.Rds")
modelld_predictions<-readRDS("Temperature analysis/model outputs and validations/modelld_predictions.Rds")

newdata<-newdata_year%>%
  mutate(Prediction=modellc4_predictions$fit)%>%
  mutate(SE=modellc4_predictions$se.fit,
         L95=Prediction-SE*1.96,
         U95=Prediction+SE*1.96)%>%
  mutate(Date=as.Date(Julian_day, origin=as.Date(paste(Year, "01", "01", sep="-")))) # Create Date variable from Julian Day and Year

newdata_d<-newdata_year%>%
  mutate(Prediction=modelld_predictions$fit)%>%
  mutate(SE=modelld_predictions$se.fit,
         L95=Prediction-SE*1.96,
         U95=Prediction+SE*1.96)%>%
  mutate(Date=as.Date(Julian_day, origin=as.Date(paste(Year, "01", "01", sep="-")))) # Create Date variable from Julian Day and Year

# Year predictions --------------------------------------------------------

mygrid <- data.frame(
  name = c("Upper Sacramento River Ship Channel", "Cache Slough and Lindsey Slough", "Lower Sacramento River Ship Channel", "Liberty Island", "Suisun Marsh", "Middle Sacramento River", "Lower Cache Slough", "Steamboat and Miner Slough", "Upper Mokelumne River", "Lower Mokelumne River", "Georgiana Slough", "Sacramento River near Ryde", "Sacramento River near Rio Vista", "Grizzly Bay", "West Suisun Bay", "Mid Suisun Bay", "Honker Bay", "Confluence", "Lower Sacramento River", "San Joaquin River at Twitchell Island", "San Joaquin River at Prisoners Pt", "Disappointment Slough", "Lower San Joaquin River", "Franks Tract", "Holland Cut", "San Joaquin River near Stockton", "Mildred Island", "Middle River", "Old River", "Upper San Joaquin River", "Grant Line Canal and Old River", "Victoria Canal"),
  row = c(2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 8, 8, 8),
  col = c(7, 4, 6, 5, 2, 8, 6, 7, 9, 9, 8, 7, 6, 2, 1, 2, 3, 4, 5, 6, 8, 9, 5, 6, 7, 9, 8, 8, 7, 9, 8, 7),
  code = c(" 1", " 1", " 2", " 3", " 8", " 4", " 5", " 6", " 7", " 9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "30", "29", "31"),
  stringsAsFactors = FALSE
)

Data_effort <- Data%>%
  st_drop_geometry()%>%
  group_by(SubRegion, Month, Year)%>%
  summarise(N=n(), .groups="drop")

newdata_year2<-newdata%>%
  select(-N)%>%
  mutate(Month=month(Date))%>%
  left_join(Data_effort, by=c("SubRegion", "Month", "Year"))%>% 
  filter(!is.na(N))

newdata_year2_d<-newdata_d%>%
  select(-N)%>%
  mutate(Month=month(Date))%>%
  left_join(Data_effort, by=c("SubRegion", "Month", "Year"))%>% 
  filter(!is.na(N))

Data_year<-Data%>%
  filter(hour(Time)<14 & hour(Time)>10)%>%
  lazy_dt()%>%
  group_by(Year, Month, Season, SubRegion)%>%
  summarize(SD=sd(Temperature), Temperature=mean(Temperature))%>%
  ungroup()%>%
  as_tibble()

newdata_sum<-newdata_year2%>%
  mutate(Var=SE^2,
         Month=month(Date))%>%
  lazy_dt()%>%
  group_by(Year, Month, SubRegion)%>%
  summarise(Temperature=mean(Prediction), SE=sqrt(sum(Var)/(n()^2)))%>%
  ungroup()%>%
  as_tibble()%>%
  mutate(L95=Temperature-1.96*SE,
         U95=Temperature+1.96*SE)

newdata_sum_d<-newdata_year2_d%>%
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
  ggplot(filter(newdata_sum_d, Month==month))+
    geom_ribbon(aes(x=Year, ymin=L95, ymax=U95), fill="firebrick3", alpha=0.5)+
    geom_line(aes(x=Year, y=Temperature), color="firebrick3")+
    geom_pointrange(data=filter(Data_year, Month==month), aes(x=Year, y=Temperature, ymin=Temperature-SD, ymax=Temperature+SD), size=0.5, alpha=0.4)+
    facet_geo(~SubRegion, grid=mygrid, labeller=label_wrap_gen())+
    theme_bw()+
    theme(panel.grid=element_blank(), axis.text.x = element_text(angle=45, hjust=1))
}

walk(1:12, function(x) ggsave(plot=mapyear(x), filename=paste0("C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/figures/Year predictions month ", x, " 11.6.20_d.png"), device=png(), width=15, height=12, units="in"))


# Rasterized predictions --------------------------------------------------

# Function to rasterize season by season. Creates a 4D raster Latitude x Longitude x Year x Season (But Season only has 1 value for the season passed to the function)
Rasterize_season<-function(season, data, n, out_crs=4326){
  
  Years <- data %>%
    #filter(Season == season) %>%
    pull(Year) %>%
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

# Function to rasterize all dates. Creates a 3D raster Latitude x Longitude x Date 
Rasterize_all <- function(data, var, out_crs=4326, n=100){
  var<-rlang::enquo(var)
  rlang::as_name(var)
  preds<-map(unique(data$Date), function(x) st_rasterize(data%>%
                                                           filter(Date==x)%>%
                                                           dplyr::select(!!var), 
                                                         template=st_as_stars(st_bbox(Delta), dx=diff(st_bbox(Delta)[c(1, 3)])/n, dy=diff(st_bbox(Delta)[c(2, 4)])/n, values = NA_real_))%>%
               st_warp(crs=out_crs))
  
  # Then bind all dates together into 1 raster
  out <- exec(c, !!!preds, along=list(Date=unique(data$Date)))
  return(out)
}

Data_effort <- Data%>%
  st_drop_geometry()%>%
  group_by(SubRegion, Month, Year)%>%
  summarise(N=n(), .groups="drop")

newdata_rast <- newdata%>%
  mutate(Month=month(Date))%>%
  select(-N)%>%
  left_join(Data_effort, by=c("SubRegion", "Month", "Year"))%>% 
  mutate(across(c(Prediction, SE, L95, U95), ~if_else(is.na(N), NA_real_, .)))

newdata_rast_d <- newdata_d%>%
  mutate(Month=month(Date))%>%
  select(-N)%>%
  left_join(Data_effort, by=c("SubRegion", "Month", "Year"))%>% 
  mutate(across(c(Prediction, SE, L95, U95), ~if_else(is.na(N), NA_real_, .)))

# Create full rasterization of all predictions for interactive visualizations
rastered_preds<-Rasterize_all(newdata_rast, Prediction)
rastered_preds_d<-Rasterize_all(newdata_rast_d, Prediction)
# Same for SE
rastered_SE<-Rasterize_all(newdata_rast, SE)
rastered_SE_d<-Rasterize_all(newdata_rast_d, SE)
# Bind SE and predictions together
rastered_predsSE<-c(rastered_preds, rastered_SE)
rastered_predsSE_d<-c(rastered_preds_d, rastered_SE_d)

#saveRDS(rastered_predsSE, file="Shiny app/Rasterized modellc4 predictions.Rds")
#saveRDS(rastered_predsSE_d, file="Shiny app/Rasterized modelld predictions.Rds")

raster_plot<-function(data, Years=unique(newdata_rast_season$Year), labels="All"){
  ggplot()+
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

newdata_rast_season <- newdata%>%
  mutate(Month=month(Date))%>%
  select(-N)%>%
  filter(Year%in%round(seq(min(Data$Year)+2, max(Data$Year)-2, length.out=9)) & Month%in%c(1,4,7,10))%>%
  left_join(Data_effort, by=c("SubRegion", "Month", "Year"))%>% 
  mutate(across(c(Prediction, SE, L95, U95), ~if_else(is.na(N), NA_real_, .)))

newdata_rast_season_d <- newdata_d%>%
  mutate(Month=month(Date))%>%
  select(-N)%>%
  filter(Year%in%seq(1970, 2018, length.out=9) & Month%in%c(1,4,7,10))%>%
  left_join(Data_effort, by=c("SubRegion", "Month", "Year"))%>% 
  mutate(across(c(Prediction, SE, L95, U95), ~if_else(is.na(N), NA_real_, .)))

# Rasterize each season
rastered_preds_season <- map(set_names(c("Winter", "Spring", "Summer", "Fall")), function(x) Rasterize_season(season=x, data=newdata_rast_season, n=100))
rastered_preds_season_d <- map(set_names(c("Winter", "Spring", "Summer", "Fall")), function(x) Rasterize_season(season=x, data=newdata_rast_season_d, n=100))

# Plot each season
p<-map2(rastered_preds_season, c("Left", "None", "None", "Right"), ~raster_plot(data=.x, labels=.y))
p_d<-map2(rastered_preds_season_d, c("Left", "None", "None", "Right"), ~raster_plot(data=.x, labels=.y))

# Use Patchwork to bind plots together. Won't look very good until the plot is saved
p2<-wrap_plots(p)+plot_layout(nrow=1, heights=c(1,1,1,1))
p2_d<-wrap_plots(p_d)+plot_layout(nrow=1, heights=c(1,1,1,1))

# Save plots
ggsave(plot=p2_d, filename="C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/figures/Rasterized predictions 11.6.20_d.png", device=png(), width=7, height=12, units="in")

# Model error by region ---------------------------------------------------

#modellc4<-readRDS("Temperature analysis/model outputs and validations/modellc4.Rds")
#modellc4_residuals <- modellc4$residuals
#saveRDS(modellc4_residuals, file="Temperature analysis/model outputs and validations/modellc4_residuals.Rds")

#modellc4_fitted <- modellc4$fitted.values
#saveRDS(modellc4_fitted, file="Temperature analysis/model outputs and validations/modellc4_fitted.Rds")

modelld<-readRDS("Temperature analysis/model outputs and validations/modelld.Rds")
modelld_residuals <- modelld$residuals
saveRDS(modelld_residuals, file="Temperature analysis/model outputs and validations/modelld_residuals.Rds")

modelld_fitted <- modelld$fitted.values
saveRDS(modelld_fitted, file="Temperature analysis/model outputs and validations/modelld_fitted.Rds")

# Stored as modellc4_residuals.Rds
modellc4_residuals<-readRDS("Temperature analysis/model outputs and validations/modellc4_residuals.Rds")

Data_resid<-Data%>%
  mutate(Residuals = modellc4_residuals)

Data_resid_d<-Data%>%
  mutate(Residuals = modelld_residuals)

Resid_sum<-Data_resid%>%
  lazy_dt()%>%
  group_by(Year, Month, SubRegion)%>%
  summarise(Resid=mean(Residuals), SD=sd(Residuals))%>%
  ungroup()%>%
  as_tibble()

Resid_sum_d<-Data_resid_d%>%
  lazy_dt()%>%
  group_by(Year, Month, SubRegion)%>%
  summarise(Resid=mean(Residuals), SD=sd(Residuals))%>%
  ungroup()%>%
  as_tibble()

p_resid<-ggplot(Resid_sum)+
  geom_tile(aes(x=Year, y=Month, fill=Resid))+
  scale_fill_gradient2(high = muted("red"),
                       low = muted("blue"),
                       breaks=seq(-3,5.5, by=0.5),
                       guide=guide_colorbar(barheight=40))+
  scale_x_continuous(breaks=unique(Resid_sum$Year), labels = if_else((unique(Resid_sum$Year)/2)%% 2 == 0, as.character(unique(Resid_sum$Year)), ""))+
  scale_y_continuous(breaks=unique(Resid_sum$Month), labels = if_else(unique(Resid_sum$Month)%% 2 == 0, as.character(unique(Resid_sum$Month)), ""))+
  facet_geo(~SubRegion, grid=mygrid, labeller=label_wrap_gen())+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1), panel.grid=element_blank(), panel.background = element_rect(fill="black"))

p_resid_d<-ggplot(Resid_sum_d)+
  geom_tile(aes(x=Year, y=Month, fill=Resid))+
  scale_fill_gradient2(high = muted("red"),
                       low = muted("blue"),
                       breaks=seq(-3,5.5, by=0.5),
                       guide=guide_colorbar(barheight=40))+
  scale_x_continuous(breaks=unique(Resid_sum_d$Year), labels = if_else((unique(Resid_sum_d$Year)/2)%% 2 == 0, as.character(unique(Resid_sum_d$Year)), ""))+
  scale_y_continuous(breaks=unique(Resid_sum_d$Month), labels = if_else(unique(Resid_sum_d$Month)%% 2 == 0, as.character(unique(Resid_sum_d$Month)), ""))+
  facet_geo(~SubRegion, grid=mygrid, labeller=label_wrap_gen())+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1), panel.grid=element_blank(), panel.background = element_rect(fill="black"))

ggsave(plot=p_resid_d, filename="C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/figures/Residuals 11.6.20_d.png", device=png(), width=20, height=12, units="in")


# Plot sampling effort ----------------------------------------------------

p_effort<-ggplot(Data_effort)+
  geom_tile(aes(x=Year, y=Month, fill=N))+
  scale_fill_viridis_c(breaks=seq(0,140, by=10),
                       guide=guide_colorbar(barheight=40))+
  scale_x_continuous(breaks=unique(Data_effort$Year), labels = if_else((unique(Data_effort$Year)/2)%% 2 == 0, as.character(unique(Data_effort$Year)), ""))+
  scale_y_continuous(breaks=unique(Data_effort$Month), labels = if_else(unique(Data_effort$Month)%% 2 == 0, as.character(unique(Data_effort$Month)), ""))+
  facet_geo(~SubRegion, grid=mygrid, labeller=label_wrap_gen())+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1), panel.grid=element_blank())

ggsave(plot=p_effort, filename="C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/figures/Effort 11.6.20.png", device=png(), width=20, height=12, units="in")


# Stratified cross-validation ---------------------------------------------
set.seed(100)
Data_split<-Data%>%
  mutate(Resid=modelld_residuals,
         Fitted=modelld_fitted)%>%
  group_by(SubRegion, Year, Season, Group)%>%
  mutate(Fold=sample(1:10, 1, replace=T))%>%
  ungroup()
set.seed(NULL)

#saveRDS(Data_split, file="Temperature analysis/Split data for cross validation.Rds")

# Saved as "Split data for cross validation.Rds"


CV_fit_1=list()
#~2 hours per model run
for(i in 1:10){
  out<-bam(Temperature ~ Year_fac + te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20), by=Year_fac) + s(Time_num_s, k=5),
           data = filter(Data_split, Group==1 & Fold!=i)%>%mutate(Year_fac=droplevels(Year_fac)), method="fREML", discrete=T, nthreads=4)
  saveRDS(out, file=paste0("Temperature analysis/model outputs and validations/CV_model_1_", i, ".Rds"))
  CV_fit_1[[i]]<-predict(out, newdata=filter(Data_split, Group==1 & Fold==i), type="response", se.fit=TRUE, discrete=T, n.threads=4)
  rm(out)
  gc()
}

saveRDS(CV_fit_1, file="Temperature analysis/model outputs and validations/Group 1 CV predictions.Rds")

CV_fit_2=list()
for(i in 1:10){
  out<-bam(Temperature ~ Year_fac + te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20), by=Year_fac) + s(Time_num_s, k=5),
           data = filter(Data_split, Group==2 & Fold!=i)%>%mutate(Year_fac=droplevels(Year_fac)), method="fREML", discrete=T, nthreads=4)
  saveRDS(out, file=paste0("Temperature analysis/model outputs and validations/CV_model_2_", i, ".Rds"))
  CV_fit_2[[i]]<-predict(out, newdata=filter(Data_split, Group==2 & Fold==i), type="response", se.fit=TRUE, discrete=T, n.threads=4)
  rm(out)
  gc()
  message(paste0("Finished run ", i, "/10"))  
}

saveRDS(CV_fit_2, file="Temperature analysis/model outputs and validations/Group 2 CV predictions.Rds")

CV_bind<-function(group, fold){
  if(group==1){
    fit<-CV_fit_1[[fold]]$fit
  }
  
  if(group==2){
    fit<-CV_fit_2[[fold]]$fit
  }
  
  Out<-Data_split%>%
    filter(Group==group & Fold==fold)%>%
    mutate(Fitted_CV=fit)
}

Data_split_CV<-map2_dfr(rep(c(1,2), each=10), rep(1:10,2), ~CV_bind(group=.x, fold=.y))%>%
  mutate(Resid_CV=Fitted_CV-Temperature,
         Fitted_resid=Fitted_CV-Fitted)

# EMP (first year-round survey) started in 1974 so restricting analysis to those years
Resid_CV_sum<-Data_split_CV%>%
  filter(Year>=1974)%>%
  lazy_dt()%>%
  group_by(Year, Month, SubRegion)%>%
  summarise(SD=sd(Resid_CV), Resid_CV=mean(Resid_CV), Fitted_resid=mean(Fitted_resid))%>%
  ungroup()%>%
  as_tibble()

# First plot deviation of predicted values from true values
p_resid_CV<-ggplot(Resid_CV_sum)+
  geom_tile(aes(x=Year, y=Month, fill=Resid_CV))+
  scale_fill_gradient2(high = muted("red"),
                       low = muted("blue"),
                       breaks=-9:7,
                       guide=guide_colorbar(barheight=40))+
  scale_x_continuous(breaks=unique(Resid_CV_sum$Year), labels = if_else((unique(Resid_CV_sum$Year)/2)%% 2 == 0, as.character(unique(Resid_CV_sum$Year)), ""))+
  scale_y_continuous(breaks=unique(Resid_CV_sum$Month), labels = if_else(unique(Resid_CV_sum$Month)%% 2 == 0, as.character(unique(Resid_CV_sum$Month)), ""))+
  facet_geo(~SubRegion, grid=mygrid, labeller=label_wrap_gen())+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1), panel.grid=element_blank(), panel.background = element_rect(fill="black"))
p_resid_CV
ggsave(plot=p_resid_CV, filename="C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/figures/CV Residuals 8.11.20.png", device=png(), width=20, height=12, units="in")


# Next plot deviation of predicted values from fitted values from full model
p_resid_CV2<-ggplot(Resid_CV_sum)+
  geom_tile(aes(x=Year, y=Month, fill=Fitted_resid))+
  scale_fill_gradient2(high = muted("red"),
                       low = muted("blue"))+
  scale_x_continuous(breaks=unique(Resid_CV_sum$Year), labels = if_else((unique(Resid_CV_sum$Year)/2)%% 2 == 0, as.character(unique(Resid_CV_sum$Year)), ""))+
  scale_y_continuous(breaks=unique(Resid_CV_sum$Month), labels = if_else(unique(Resid_CV_sum$Month)%% 2 == 0, as.character(unique(Resid_CV_sum$Month)), ""))+
  facet_geo(~SubRegion, grid=mygrid, labeller=label_wrap_gen())+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1), panel.grid=element_blank(), panel.background = element_rect(fill="black"))
ggsave(plot=p_resid_CV2, filename="C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/figures/CV Residuals2 6.16.20.png", device=png(), width=20, height=12, units="in")

# model d

CVd_fit_1=list()
#~2 hours per model run
for(i in 1:10){
  out<-bam(Temperature ~ Year_fac + te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20), by=Year_fac) + 
             te(Time_num_s, Julian_day_s, bs=c("tp", "cc"), k=c(5, 12)),
           data = filter(Data_split, Group==1 & Fold!=i)%>%mutate(Year_fac=droplevels(Year_fac)), method="fREML", discrete=T, nthreads=4)
  saveRDS(out, file=paste0("Temperature analysis/model outputs and validations/CVd_model_1_", i, ".Rds"))
  CVd_fit_1[[i]]<-predict(out, newdata=filter(Data_split, Group==1 & Fold==i), type="response", se.fit=TRUE, discrete=T, n.threads=4)
  rm(out)
  gc()
  message(paste0("Finished run ", i, "/10")) 
}

saveRDS(CVd_fit_1, file="Temperature analysis/model outputs and validations/Group 1 CV predictions d.Rds")

CVd_fit_2=list()
for(i in 1:10){
  out<-bam(Temperature ~ Year_fac + te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20), by=Year_fac) + 
             te(Time_num_s, Julian_day_s, bs=c("tp", "cc"), k=c(5, 12)),
           data = filter(Data_split, Group==2 & Fold!=i)%>%mutate(Year_fac=droplevels(Year_fac)), method="fREML", discrete=T, nthreads=4)
  saveRDS(out, file=paste0("Temperature analysis/model outputs and validations/CVd_model_2_", i, ".Rds"))
  CVd_fit_2[[i]]<-predict(out, newdata=filter(Data_split, Group==2 & Fold==i), type="response", se.fit=TRUE, discrete=T, n.threads=4)
  rm(out)
  gc()
  message(paste0("Finished run ", i, "/10"))  
}

saveRDS(CVd_fit_2, file="Temperature analysis/model outputs and validations/Group 2 CV predictions d.Rds")

CV_bind<-function(group, fold){
  if(group==1){
    fit<-CVd_fit_1[[fold]]$fit
  }
  
  if(group==2){
    fit<-CVd_fit_2[[fold]]$fit
  }
  
  Out<-Data_split%>%
    filter(Group==group & Fold==fold)%>%
    mutate(Fitted_CV=fit)
}

Data_split_CV_d<-map2_dfr(rep(c(1,2), each=10), rep(1:10,2), ~CV_bind(group=.x, fold=.y))%>%
  mutate(Resid_CV=Fitted_CV-Temperature,
         Fitted_resid=Fitted_CV-Fitted)

# EMP (first year-round survey) started in 1974 so restricting analysis to those years
Resid_CV_sum_d<-Data_split_CV_d%>%
  filter(Year>=1974)%>%
  lazy_dt()%>%
  group_by(Year, Month, SubRegion)%>%
  summarise(SD=sd(Resid_CV), Resid_CV=mean(Resid_CV), Fitted_resid=mean(Fitted_resid))%>%
  ungroup()%>%
  as_tibble()

RMSE <- function(m, o){
  sqrt(mean((m - o)^2))
}

CV_sum<-Data_split_CV_d%>%
  st_drop_geometry()%>%
  group_by(Group, Fold)%>%
  summarise(RMSE=sqrt(mean(Resid_CV^2)), 
            r=cor(Fitted_CV, Temperature, method="pearson"), .groups="drop")

# First plot deviation of predicted values from true values
p_resid_CV<-ggplot(Resid_CV_sum_d)+
  geom_tile(aes(x=Year, y=Month, fill=Resid_CV))+
  scale_fill_gradient2(high = muted("red"),
                       low = muted("blue"),
                       breaks=-9:7,
                       guide=guide_colorbar(barheight=40))+
  scale_x_continuous(breaks=unique(Resid_CV_sum_d$Year), labels = if_else((unique(Resid_CV_sum_d$Year)/2)%% 2 == 0, as.character(unique(Resid_CV_sum_d$Year)), ""))+
  scale_y_continuous(breaks=unique(Resid_CV_sum_d$Month), labels = if_else(unique(Resid_CV_sum_d$Month)%% 2 == 0, as.character(unique(Resid_CV_sum_d$Month)), ""))+
  facet_geo(~SubRegion, grid=mygrid, labeller=label_wrap_gen())+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1), panel.grid=element_blank(), panel.background = element_rect(fill="black"))
p_resid_CV
ggsave(plot=p_resid_CV, filename="C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/figures/CV Residuals 11.9.20 d.png", device=png(), width=20, height=12, units="in")


# Next plot deviation of predicted values from fitted values from full model
p_resid_CV2<-ggplot(Resid_CV_sum_d)+
  geom_tile(aes(x=Year, y=Month, fill=Fitted_resid))+
  scale_fill_gradient2(high = muted("red"),
                       low = muted("blue"))+
  scale_x_continuous(breaks=unique(Resid_CV_sum_d$Year), labels = if_else((unique(Resid_CV_sum_d$Year)/2)%% 2 == 0, as.character(unique(Resid_CV_sum_d$Year)), ""))+
  scale_y_continuous(breaks=unique(Resid_CV_sum_d$Month), labels = if_else(unique(Resid_CV_sum_d$Month)%% 2 == 0, as.character(unique(Resid_CV_sum_d$Month)), ""))+
  facet_geo(~SubRegion, grid=mygrid, labeller=label_wrap_gen())+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1), panel.grid=element_blank(), panel.background = element_rect(fill="black"))
ggsave(plot=p_resid_CV2, filename="C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/figures/CV Residuals2 6.16.20.png", device=png(), width=20, height=12, units="in")


# Test autocorrelation ----------------------------------------------------

auto<-Data%>%
  mutate(Resid=modellc4_residuals)%>%
  filter(Source!="EDSM" & !str_detect(Station, "EZ"))%>% # Remove EDSM and EZ stations because they're not fixed
  mutate(Station=paste(Source, Station))%>%
  group_by(Station)%>%
  mutate(N=n())%>%
  filter(N>10)%>%
  summarise(ACF=list(pacf(Resid, plot=F)), N=n(), ci=qnorm((1 + 0.95)/2)/sqrt(n()), .groups="drop")%>% # ci formula from https://stackoverflow.com/questions/14266333/extract-confidence-interval-values-from-acf-correlogram
  rowwise()%>%
  mutate(lag=list(ACF$lag), acf=list(ACF$acf))%>%
  unnest(cols=c(lag, acf))%>%
  arrange(-N)%>%
  mutate(Station=factor(Station, levels=unique(Station)))

length(which(abs(auto$acf)>abs(auto$ci)))/nrow(auto)

# Only 6% exceed the CI, very close to the 5% you would expect with our chosen confidence level of 0.95 so I'm taking this as good evidence of no autocorrelation

length(which(abs(filter(auto, lag==1)$acf)>abs(filter(auto, lag==1)$ci)))/nrow(filter(auto, lag==1))
# Around 30% exceed the CI at a lag of 1
# Limitations: Not all surveys sample monthly, so lags aren't constant

ggplot(filter(auto, lag==1))+
  geom_point(aes(x=Station, y=abs(acf)), fill="black", shape=21)+
  geom_point(data=filter(auto, lag==1 & abs(acf)>abs(ci)), aes(x=Station, y=abs(acf)), fill="red", shape=21)+
  geom_point(aes(x=Station, y=abs(ci)), fill="white", shape=21)+
  geom_segment(aes(x=Station, y=abs(acf), xend=Station, yend=abs(ci)), linetype=2)+
  geom_segment(data=filter(auto, lag==1 & abs(acf)>abs(ci)), aes(x=Station, y=abs(acf), xend=Station, yend=abs(ci)), color="red")+
  theme_bw()+
  theme(panel.grid=element_blank(), axis.text.x=element_text(angle=45, hjust=1))

# Using a continuous method
require(gstat)
require(sp)
require(spacetime)
Vario_data <- data.frame(Residuals=modellc4_residuals/sd(Data$Temperature), Longitude=Data$Longitude, Latitude=Data$Latitude)
coordinates(Vario_data)<-c("Longitude", "Latitude")
vario<-variogram(Residuals~1, Vario_data)
plot(vario)

sp <- SpatialPoints(coords=data.frame(Longitude=Data$Longitude, Latitude=Data$Latitude))
sp2<-STIDF(sp, time=Data$Date, data=data.frame(Residuals=modellc4_residuals/sd(Data$Temperature)))
vario2<-variogramST(Residuals~1, data=sp2, tunit="weeks", cores=3)
save(vario, vario2, file="Variograms.Rds")

ggplot(vario2, aes(x=timelag, y=gamma, color=avgDist, group=avgDist))+
  geom_line()+
  geom_point()
