#devtools::install_github("InteragencyEcologicalProgram/discretewq")
require(dplyr)
require(stringr)
require(tidyr)
require(ggplot2)
require(geofacet)
require(lubridate)
require(hms)
require(sf)
require(dtplyr)
require(discretewq)

mygrid <- data.frame(
  name = c("Upper Sacramento River Ship Channel", "Cache Slough and Lindsey Slough", "Lower Sacramento River Ship Channel", "Liberty Island", "Suisun Marsh", "Middle Sacramento River", "Lower Cache Slough", "Steamboat and Miner Slough", "Upper Mokelumne River", "Lower Mokelumne River", "Georgiana Slough", "Sacramento River near Ryde", "Sacramento River near Rio Vista", "Grizzly Bay", "West Suisun Bay", "Mid Suisun Bay", "Honker Bay", "Confluence", "Lower Sacramento River", "San Joaquin River at Twitchell Island", "San Joaquin River at Prisoners Pt", "Disappointment Slough", "Lower San Joaquin River", "Franks Tract", "Holland Cut", "San Joaquin River near Stockton", "Mildred Island", "Middle River", "Old River", "Upper San Joaquin River", "Grant Line Canal and Old River", "Victoria Canal"),
  row = c(2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 8, 8, 8),
  col = c(7, 4, 6, 5, 2, 8, 6, 7, 9, 9, 8, 7, 6, 2, 1, 2, 3, 4, 5, 6, 8, 9, 5, 6, 7, 9, 8, 8, 7, 9, 8, 7),
  code = c(" 1", " 1", " 2", " 3", " 8", " 4", " 5", " 6", " 7", " 9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "30", "29", "31"),
  stringsAsFactors = FALSE
)

# Create overall dataset --------------------------

is.even <- function(x) as.integer(x) %% 2 == 0

## Load Delta Shapefile from Brian
Delta<-st_read("Delta Subregions")%>%
  filter(!SubRegion%in%c("South Bay", "San Francisco Bay", "San Pablo Bay", "Upper Yolo Bypass", 
                         "Upper Napa River", "Lower Napa River", "Carquinez Strait"))%>% # Remove regions outside our domain of interest
  dplyr::select(SubRegion)

## Load data
Data <- wq(Sources=c("EMP", "STN", "FMWT", "EDSM", "DJFMP", 
                     "SDO", "SKT", "SLS", "20mm", "Suisun", "NCRO",
                     "Baystudy", "YBFMP"))%>% # "USBR", "USGS_SFBS", and "USGS_CAWSC" don't have secchi data
  filter(!is.na(Secchi) & !is.na(Datetime) & !is.na(Latitude) & !is.na(Longitude) & !is.na(Date))%>% #Remove any rows with NAs in our key variables
  filter(hour(Datetime)>=5 & hour(Datetime)<=20)%>% # Only keep data between 5AM and 8PM
  mutate(Datetime = with_tz(Datetime, tz="Etc/GMT+7"), #Convert to a timezone without daylight savings time
         Date = with_tz(Date, tz="Etc/GMT+7"),
         Time=as_hms(Datetime), # Create variable for time-of-day, not date. 
         Noon_diff=abs(hms(hours=12)-Time))%>% # Calculate difference from noon for each data point for later filtering
  group_by(Station, Source, Date)%>%
  filter(Noon_diff==min(Noon_diff))%>% # Select only 1 data point per station and date, choose data closest to noon
  filter(Time==min(Time))%>% # When points are equidistant from noon, select earlier point
  ungroup()%>%
  distinct(Date, Station, Source, .keep_all = TRUE)%>% # Finally, fix the 2 straggling datapoints from the same time and station
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)%>% # Convert to sf object
  st_transform(crs=st_crs(Delta))%>% # Change to crs of Delta
  st_join(Delta, join=st_intersects)%>% # Add subregions
  filter(!is.na(SubRegion))%>% # Remove any data outside our subregions of interest
  mutate(Julian_day = yday(Date), # Create julian day variable
         Month_fac=factor(Month), # Create month factor variable
         Source_fac=factor(Source),
         Year_fac=factor(Year))%>% 
  mutate(Date_num = as.numeric(Date))%>%  # Create numeric version of date for models
  mutate(Time_num=as.numeric(Time)) # Create numeric version of time for models (=seconds since midnight)


## Pull station locations for major monitoring programs
### This will be used to set a boundary for this analysis focused on well-sampled regions.
WQ_stations<-Data%>%
  st_drop_geometry()%>%
  filter(Source%in%c("FMWT", "STN", "SKT", "20mm", "EMP", "Suisun", "YBFMP", "DJFMP", "NCRO"))%>%
  group_by(StationID, Source, Latitude, Longitude)%>%
  summarise(N=n(), .groups="drop")%>% # Calculate how many times each station was sampled
  filter(N>50 & !StationID%in%c("20mm 918", "STN 918"))%>% # Only keep stations sampled >50 times when deciding which regions to retain. 
  # "20mm 918", "STN 918" are far south of the rest of the well-sampled sites and are not sampled year round, so we're removing them to exclude that far southern region
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)%>% # Convert to sf object
  st_transform(crs=st_crs(Delta))%>%
  st_join(Delta) # Add subregions

#Which stations will be removed in the next step?
setdiff(Delta$SubRegion, WQ_stations$SubRegion)

## Remove any subregions that do not contain at least one of these >50 samples stations from the major monitoring programs
Delta <- Delta%>%
  filter(SubRegion%in%unique(WQ_stations$SubRegion) | SubRegion=="Georgiana Slough") # Retain Georgiana Slough because it's surrounded by well-sampled regions

ggplot()+
  geom_sf(data=Delta)+
  geom_sf(data=WQ_stations%>%
            st_union()%>%
            st_convex_hull()%>% # Draws a hexagram or pentagram or similar around the outer-most points
            st_as_sf(),
          fill=NA, color="red")+
  geom_sf(data=filter(Data, Source!="EDSM"), alpha=0.1, mapping=aes(color=Source))+
  scale_color_discrete(guide=guide_legend(override.aes = list(alpha=1)))+
  theme_bw()

# Skipping the convex hull step used in the Temp analyses because here it woiuld only remove stations that are very close to other well-sampled stations

## Now filter data to only include this final set of subregions, and any stations outside the convex hull formed by the >50 samples stations from the major monitoring programs
Data<-Data%>%
  filter(SubRegion%in%unique(Delta$SubRegion))%>%
  mutate(Group=if_else(is.even(Year), 1, 2))%>%
  mutate_at(vars(Date_num, Longitude, Latitude, Time_num, Year, Julian_day), list(s=~(.-mean(., na.rm=T))/sd(., na.rm=T))) # Create centered and standardized versions of covariates

saveRDS(Data, file="Secchi analysis/Discrete Secchi Data.Rds")
#Data<-readRDS("Secchi analysis/Discrete Secchi Data.Rds")

Data_analysis<-Data%>%
  filter(Source%in%c("FMWT", "Baystudy", "STN", "Suisun", "SDO", "SKT", "SLS", "20mm", "DJFMP", "EMP", "YBFMP", "DWR_NCRO") & !str_detect(Station, "EZ") & 
           !(Source=="SKT" & Station=="799" & Latitude>38.2) & !(Source=="SKT" & Station=="999"))%>%
  mutate(Station=paste(Source, Station),
         Noon_diff=abs(hms(hours=12)-as_hms(Datetime)),
         mday_15_diff=abs(mday(Date)-15))%>% # Find how far each date is from the 15th of the month
  group_by(Station, Month, Year)%>%
  filter(mday_15_diff==min(mday_15_diff))%>%
  filter(Date==min(Date))%>% # Deal with dates equidistant from the 15th of the month
  filter(Noon_diff==min(Noon_diff))%>%
  ungroup()%>%
  mutate(YearStation=paste(Year, Station),
         Date_num2=as.numeric(Date)/(3600*24*30), # Create a numeric date variable in units of ~ 1 month. 
         Month_fac=factor(Month))%>%
  arrange(Station, Date)%>%
  group_by(Station)%>%
  mutate(Lag=Date-lag(Date, order_by = Date))%>%
  ungroup()%>%
  mutate(Start=if_else(is.na(Lag) | Lag>3600*24*60, TRUE, FALSE),
         Series_ID=1:n(),
         Series_ID=if_else(Start, Series_ID, NA_integer_),
         Series_ID=as.integer(as.factor(Series_ID)))%>%
  fill(Series_ID, .direction="down")

saveRDS(Data_analysis, file="Secchi analysis/Discrete Secchi Analysis Data.Rds")
