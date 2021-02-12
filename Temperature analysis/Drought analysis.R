require(sp)
require(gstat)
require(spacetime)
require(tidybayes)
require(dplyr)
require(tidyr)
require(stringr)
require(dtplyr)
require(broom)
require(brms)
require(mgcv)
require(ggplot2)
require(geofacet)
require(lubridate)
require(hms)
require(scales)
require(sf)
require(stars)
require(itsadug)
require(purrr)
require(readr)
require(slider)
require(colorspace)

mygrid <- data.frame(
  name = c("Upper Sacramento River Ship Channel", "Cache Slough and Lindsey Slough", "Lower Sacramento River Ship Channel", "Liberty Island", "Suisun Marsh", "Middle Sacramento River", "Lower Cache Slough", "Steamboat and Miner Slough", "Upper Mokelumne River", "Lower Mokelumne River", "Georgiana Slough", "Sacramento River near Ryde", "Sacramento River near Rio Vista", "Grizzly Bay", "West Suisun Bay", "Mid Suisun Bay", "Honker Bay", "Confluence", "Lower Sacramento River", "San Joaquin River at Twitchell Island", "San Joaquin River at Prisoners Pt", "Disappointment Slough", "Lower San Joaquin River", "Franks Tract", "Holland Cut", "San Joaquin River near Stockton", "Mildred Island", "Middle River", "Old River", "Upper San Joaquin River", "Grant Line Canal and Old River", "Victoria Canal"),
  row = c(2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 8, 8, 8),
  col = c(7, 4, 6, 5, 2, 8, 6, 7, 9, 9, 8, 7, 6, 2, 1, 2, 3, 4, 5, 6, 8, 9, 5, 6, 7, 9, 8, 8, 7, 9, 8, 7),
  code = c(" 1", " 1", " 2", " 3", " 8", " 4", " 5", " 6", " 7", " 9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "30", "29", "31"),
  stringsAsFactors = FALSE
)

WY<-waterYearType::water_year_indices%>%
  bind_rows(tibble(WY=rep(c(2018, 2019), each=2), location=rep(c("Sacramento Valley", "San Joaquin Valley"), 2), 
                   Index=c(7.14, 3.03, 10.34, 4.94), WYsum=c(12.86, 4.76, 24.77, 9.28)))%>% # Add data for 2018 and 2019
  group_by(WY)%>%
  summarise(Index=sum(Index), WYsum=sum(WYsum), .groups="drop")

Data_D<-readRDS("Temperature analysis/Data_CC4.Rds")%>%
  mutate(WY=Year)%>%
  mutate(WY=if_else(Month%in%10:12, WY+1, WY))%>%
  left_join(WY, by="WY")%>%
  filter(!is.na(WYsum))%>%
  arrange(Station, Date)%>%
  group_by(Station)%>%
  mutate(Lag=Date-lag(Date, order_by = Date))%>%
  ungroup()%>%
  mutate(Start=if_else(is.na(Lag) | Lag>3600*24*60, TRUE, FALSE),
         Series_ID=1:n(),
         Series_ID=if_else(Start, Series_ID, NA_integer_),
         Series_ID=as.integer(as.factor(Series_ID)))%>%
  fill(Series_ID, .direction="down")%>%
  mutate(WYsum_c=WYsum-mean(WYsum),
         WY_s=(WY-mean(unique(WY)))/sd(unique(WY)))

Delta<-st_read("Delta Subregions")%>%
  select(SubRegion)%>%
  filter(SubRegion%in%unique(Data_D$SubRegion))

D_gam_NOAR <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                    te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=WYsum_c) + 
                    s(Time_num_s, k=5), family=scat, data = Data_D, method="fREML", discrete=T, nthreads=3)
r <- start_value_rho(D_gam_NOAR, plot=TRUE)


# Best model
D_gam_AR <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                  te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=WYsum_c) + 
                  s(Time_num_s, k=5), family=scat, rho=r, AR.start=Start, data = Data_D, method="fREML", discrete=T, nthreads=3)
#AIC: 195059.8
#AIC: 195059.8

# Best model

D_gam_AR2 <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(50, 13)) + 
                   te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=WYsum_c) + 
                   s(Time_num_s, k=5), family=scat, rho=r, AR.start=Start, data = Data_D, method="fREML", discrete=T, nthreads=3)
#AIC: 194590.7
#BIC: 199277.9

# D_gam_AR2 has a better AIC value but worse BIC and model predictions are almost identical, so picking D_gam_AR as best model. 


# Now try adding a temperature term
D_gam_NOAR2 <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                     te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=WYsum_c) + 
                     te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=WY_s) + 
                     s(Time_num_s, k=5), family=scat, data = Data_D, method="fREML", discrete=T, nthreads=3)
r2 <- start_value_rho(D_gam_NOAR2, plot=TRUE)

## Best
D_gam_AR3 <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                   te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=WYsum_c) + 
                   te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=WY_s) + 
                   s(Time_num_s, k=5), family=scat, rho=r2, AR.start=Start, data = Data_D, method="fREML", discrete=T, nthreads=3)
##Best

D_gam_AR4 <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(50, 13)) + 
                   te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=WYsum_c) + 
                   te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=WY_s) + 
                   s(Time_num_s, k=5), family=scat, rho=r2, AR.start=Start, data = Data_D, method="fREML", discrete=T, nthreads=3)
# Just as before, adding higher k does not change model predictions

# Does the order of the terms matter?

D_gam_AR5 <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                   te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=WY_s) + 
                   te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=WYsum_c) + 
                   s(Time_num_s, k=5), family=scat, rho=r2, AR.start=Start, data = Data_D, method="fREML", discrete=T, nthreads=3)

# No, predictions still identical
##Just as before, no benefit to increased k. Predictions are identical.

# Predict effect

newdata<-readRDS("Temperature analysis/Prediction Data.Rds")
D_newdata<-newdata%>%
  st_drop_geometry()%>%
  as_tibble()%>%
  select(-Year_fac, -Year, -Year_s, -N,)%>%
  distinct()%>%
  mutate(WYsum_c=10,
         WY_s=2)

D_pred<-predict(D_gam_AR, newdata=D_newdata, type="terms", se.fit=TRUE, discrete=T, n.threads=4)

D_newdata_pred<-D_newdata%>%
  mutate(Slope=D_pred$fit[,"te(Julian_day_s,Latitude_s,Longitude_s):WYsum_c"],
         Slope_se=D_pred$se.fit[,"te(Julian_day_s,Latitude_s,Longitude_s):WYsum_c"],
         Intercept=D_pred$fit[,"te(Julian_day_s,Latitude_s,Longitude_s)"]+D_pred$fit[,"s(Time_num_s)"],
  )%>%
  mutate(across(c(Slope, Slope_se), ~(.x/WYsum_c)))%>%
  mutate(Slope_se=abs(Slope_se))%>%
  mutate(Date=as.Date(Julian_day, origin=as.Date(paste("2000", "01", "01", sep="-"))),
         Month=month(Date, label = T),
         Slope_l95=Slope-Slope_se*qnorm(0.995),
         Slope_u95=Slope+Slope_se*qnorm(0.995))%>%
  mutate(Sig=if_else(Slope_u95>0 & Slope_l95<0, "ns", "*"))%>%
  filter(Sig=="*")%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=F)%>%
  st_transform(crs=st_crs(Delta))


# Function to rasterize all dates. Creates a 3D raster Latitude x Longitude x Date 
Rasterize_all <- function(data, var, region=Delta, out_crs=4326, n=100){
  var<-rlang::enquo(var)
  rlang::as_name(var)
  preds<-map(unique(data$Date), function(x) st_rasterize(data%>%
                                                           filter(Date==x)%>%
                                                           dplyr::select(!!var), 
                                                         template=st_as_stars(st_bbox(region), dx=diff(st_bbox(region)[c(1, 3)])/n, dy=diff(st_bbox(region)[c(2, 4)])/n, values = NA_real_))%>%
               st_warp(crs=out_crs))
  
  # Then bind all dates together into 1 raster
  out <- exec(c, !!!preds, along=list(Date=unique(data$Date)))
  return(out)
}

newdata_D_pred_rast<-Rasterize_all(D_newdata_pred, Slope)

p_D_gam<-ggplot()+
  geom_stars(data=newdata_D_pred_rast)+
  facet_wrap(~month(Date, label=T), drop=F)+
  scale_fill_continuous_diverging(palette="Blue-Red 3", breaks=(-15:7)/100, na.value = "white",
                                  name="Temperature change\nper change in\ntotal runoff (°C/maf)", guide=guide_colorbar(barheight=20))+
  ylab("Latitude")+
  xlab("Longitude")+
  coord_sf()+
  theme_bw()+
  theme(strip.background=element_blank(), axis.text.x = element_text(angle=45, hjust=1), panel.grid=element_blank())

ggsave(p_D_gam, filename="C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/figures/D_gam 12.17.20.png",
       device="png", width=7, height=5, units="in")

# Has runoff changed over time?

runoff_year<-Data_D%>%
  select(WYsum, WY)%>%
  distinct()

p_runoff<-ggplot(runoff_year, aes(x=WY, y=WYsum))+
  geom_point()+
  ylab("Total Sac + SJ runoff (maf)")+
  theme_bw()

ggsave(p_runoff, filename="C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/figures/runoff year.png",
       device="png", width=4, height=3, units="in")

D_CC_pred<-predict(D_gam_AR3, newdata=D_newdata, type="terms", se.fit=TRUE, discrete=T, n.threads=4)

D_CC_newdata_pred<-D_newdata%>%
  mutate(Slope=D_CC_pred$fit[,"te(Julian_day_s,Latitude_s,Longitude_s):WYsum_c"],
         Slope_se=D_CC_pred$se.fit[,"te(Julian_day_s,Latitude_s,Longitude_s):WYsum_c"],
         Intercept=D_CC_pred$fit[,"te(Julian_day_s,Latitude_s,Longitude_s)"]+D_pred$fit[,"s(Time_num_s)"],
  )%>%
  mutate(across(c(Slope, Slope_se), ~(.x/WYsum_c)))%>%
  mutate(Slope_se=abs(Slope_se))%>%
  mutate(Date=as.Date(Julian_day, origin=as.Date(paste("2000", "01", "01", sep="-"))),
         Month=month(Date, label = T),
         Slope_l95=Slope-Slope_se*qnorm(0.995),
         Slope_u95=Slope+Slope_se*qnorm(0.995))%>%
  mutate(Sig=if_else(Slope_u95>0 & Slope_l95<0, "ns", "*"))%>%
  filter(Sig=="*")%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=F)%>%
  st_transform(crs=st_crs(Delta))

newdata_D_CC_pred_rast<-Rasterize_all(D_CC_newdata_pred, Slope)

p_D_CC_gam<-ggplot()+
  geom_stars(data=newdata_D_CC_pred_rast)+
  facet_wrap(~month(Date, label=T), drop=F)+
  scale_fill_continuous_diverging(palette="Blue-Red 3", breaks=(-15:7)/100, na.value = "white",
                                  name="Temperature change\nper change in\ntotal runoff (°C/maf)", guide=guide_colorbar(barheight=20))+
  ylab("Latitude")+
  xlab("Longitude")+
  coord_sf()+
  theme_bw()+
  theme(strip.background=element_blank(), axis.text.x = element_text(angle=45, hjust=1), panel.grid=element_blank())

ggsave(p_D_CC_gam, filename="C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/figures/D_CC_gam WYsum 12.17.20.png",
       device="png", width=7, height=5, units="in")

D_CC_newdata_pred_CC<-D_newdata%>%
  mutate(Slope=D_CC_pred$fit[,"te(Julian_day_s,Latitude_s,Longitude_s):WY_s"],
         Slope_se=D_CC_pred$se.fit[,"te(Julian_day_s,Latitude_s,Longitude_s):WY_s"],
         Intercept=D_CC_pred$fit[,"te(Julian_day_s,Latitude_s,Longitude_s)"]+D_CC_pred$fit[,"s(Time_num_s)"])%>%
  mutate(across(c(Slope, Slope_se), ~(.x/WY_s)/sd(unique(Data_D$WY))))%>%
  mutate(Slope_se=abs(Slope_se))%>%
  mutate(Date=as.Date(Julian_day, origin=as.Date(paste("2000", "01", "01", sep="-"))),
         Month=month(Date, label = T),
         Slope_l95=Slope-Slope_se*qnorm(0.995),
         Slope_u95=Slope+Slope_se*qnorm(0.995))%>%
  mutate(Sig=if_else(Slope_u95>0 & Slope_l95<0, "ns", "*"))%>%
  filter(Sig=="*")%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=F)%>%
  st_transform(crs=st_crs(Delta))

newdata_D_CC_pred_rast_CC<-Rasterize_all(D_CC_newdata_pred_CC, Slope)

p_D_gam_CC<-ggplot()+
  geom_stars(data=newdata_D_CC_pred_rast_CC)+
  facet_wrap(~month(Date, label=T), drop=F)+
  scale_fill_viridis_c(limits=c(-0.025, 0.065), breaks=(-6:7)/100, name="Temperature change\nper year (°C)", guide=guide_colorbar(barheight=20), na.value="white")+
  ylab("Latitude")+
  xlab("Longitude")+
  coord_sf()+
  theme_bw()+
  theme(strip.background=element_blank(), axis.text.x = element_text(angle=45, hjust=1), panel.grid=element_blank())

ggsave(p_D_gam_CC, filename="C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/figures/D_CC_gam WY 12.17.20.png",
       device="png", width=7, height=5, units="in")
# Very similar to climate change model output


# Now using total inflow (QTOT) from Dayflow
dayflow<-read_csv(file.path("Temperature analysis", "data", "Dayflow1997 2019.csv"), col_types = cols_only(Date="c", TOT="d"))%>%
  bind_rows(read_csv(file.path("Temperature analysis", "data", "Dayflow1984 1996.csv"), col_types = cols_only(Date="c", TOT="d")))%>%
  bind_rows(read_csv(file.path("Temperature analysis", "data", "Dayflow1970 1983.csv"), col_types = cols_only(Date="c", TOT="d")))%>%
  bind_rows(read_csv(file.path("Temperature analysis", "data", "Dayflow1956 1969.csv"), col_types = cols_only(Date="c", TOT="d")))%>%
    mutate(Date=parse_date_time(Date, "%m/%d/%Y", tz = "America/Los_Angeles"))%>%
  arrange(Date)

dayflow$TOT_mean30<-slide_index_dbl(dayflow$TOT, dayflow$Date, .before=days(30), .f=mean, .complete = T)

dayflow$TOT_mean7<-slide_index_dbl(dayflow$TOT, dayflow$Date, .before=days(7), .f=mean, .complete = T)

Data_D2<-readRDS("Temperature analysis/Data_CC4.Rds")%>%
  mutate(WY=Year)%>%
  mutate(WY=if_else(Month%in%10:12, WY+1, WY))%>%
  left_join(dayflow, by="Date")%>%
  filter(!is.na(TOT))%>%
  arrange(Station, Date)%>%
  group_by(Station)%>%
  mutate(Lag=Date-lag(Date, order_by = Date))%>%
  ungroup()%>%
  mutate(Start=if_else(is.na(Lag) | Lag>3600*24*60, TRUE, FALSE),
         Series_ID=1:n(),
         Series_ID=if_else(Start, Series_ID, NA_integer_),
         Series_ID=as.integer(as.factor(Series_ID)))%>%
  fill(Series_ID, .direction="down")%>%
  group_by(Month)%>%
  mutate(across(c(TOT, TOT_mean7, TOT_mean30), list(c_month=~.x-mean(.x), s_month=~(.x-mean(.x))/sd(.x), monthmean=~mean(.x), monthsd=~sd(.x))))%>%
  ungroup()%>%
  mutate(across(c(TOT, TOT_mean7, TOT_mean30), list(c=~.x-mean(.x), s=~(.x-mean(.x))/sd(.x))))%>%
  mutate(WY_s=(WY-mean(unique(WY)))/sd(unique(WY)))
  

Delta2<-st_read("Delta Subregions")%>%
  select(SubRegion)%>%
  filter(SubRegion%in%unique(Data_D2$SubRegion))

D2_gam_NOAR <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                    te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=TOT_mean30_c) + 
                    s(Time_num_s, k=5), family=scat, data = Data_D2, method="fREML", discrete=T, nthreads=2)
D2r <- start_value_rho(D2_gam_NOAR, plot=TRUE)

D2_gam_AR <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                  te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=TOT_mean30_c) + 
                  s(Time_num_s, k=5), family=scat, rho=D2r, AR.start=Start, data = Data_D2, method="fREML", discrete=T, nthreads=2)

D2_gam_AR2 <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(50, 13)) + 
                   te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=TOT_mean30_c) + 
                   s(Time_num_s, k=5), family=scat, rho=D2r, AR.start=Start, data = Data_D2, method="fREML", discrete=T, nthreads=2)
# Model predictions are almost identical, so picking D2_gam_AR as best model. 

# Try centered and standardized

D2_gam_NOAR2 <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                     te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=TOT_mean30_s) + 
                     s(Time_num_s, k=5), family=scat, data = Data_D2, method="fREML", discrete=T, nthreads=3)
D2r2 <- start_value_rho(D2_gam_NOAR2, plot=TRUE)

D2_gam_AR3 <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                   te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=TOT_mean30_s) + 
                   s(Time_num_s, k=5), family=scat, rho=D2r2, AR.start=Start, data = Data_D2, method="fREML", discrete=T, nthreads=3)

# Try centered and standardized separately for each month
## These results are much easier to interpret, so going with these moving forward
D2_gam_NOAR3 <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                      te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=TOT_mean30_s_month) + 
                      s(Time_num_s, k=5), family=scat, data = Data_D2, method="fREML", discrete=T, nthreads=3)
D2r3 <- start_value_rho(D2_gam_NOAR3, plot=TRUE)

D2_gam_AR4 <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                    te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=TOT_mean30_s_month) + 
                    s(Time_num_s, k=5), family=scat, rho=D2r3, AR.start=Start, data = Data_D2, method="fREML", discrete=T, nthreads=3)

D2_gam_AR5 <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(50, 13)) + 
                    te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=TOT_mean30_s_month) + 
                    s(Time_num_s, k=5), family=scat, rho=D2r3, AR.start=Start, data = Data_D2, method="fREML", discrete=T, nthreads=3)
# Model predictions are almost identical, so picking D2_gam_AR as best model. 

newdata<-readRDS("Temperature analysis/Prediction Data.Rds")
D2_newdata<-newdata%>%
  st_drop_geometry()%>%
  as_tibble()%>%
  select(-Year_fac, -Year, -Year_s, -N,)%>%
  distinct()%>%
  mutate(TOT_mean30_c=50000,
         TOT_mean30_s=2,
         TOT_mean30_s_month=2,
         WY_s=2)

D2_pred<-predict(D2_gam_AR4, newdata=D2_newdata, type="terms", se.fit=TRUE, discrete=T, n.threads=4)

D2_newdata_pred<-D2_newdata%>%
  mutate(Slope=D2_pred$fit[,"te(Julian_day_s,Latitude_s,Longitude_s):TOT_mean30_s_month"],
         Slope_se=D2_pred$se.fit[,"te(Julian_day_s,Latitude_s,Longitude_s):TOT_mean30_s_month"],
         Intercept=D2_pred$fit[,"te(Julian_day_s,Latitude_s,Longitude_s)"]+D2_pred$fit[,"s(Time_num_s)"],
  )%>%
  mutate(across(c(Slope, Slope_se), ~(.x/TOT_mean30_s)))%>%
  mutate(Slope_se=abs(Slope_se))%>%
  mutate(Slope_l95=Slope-Slope_se*qnorm(0.995),
         Slope_u95=Slope+Slope_se*qnorm(0.995),
         Date=as.Date(Julian_day, origin=as.Date(paste("2000", "01", "01", sep="-"))),
         Month=month(Date, label=T))%>%
  mutate(Sig=if_else(Slope_u95>0 & Slope_l95<0, "ns", "*"))%>%
  filter(Sig=="*")%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=F)%>%
  st_transform(crs=st_crs(Delta2))

newdata_D2_pred_rast<-Rasterize_all(D2_newdata_pred, Slope, region=Delta2)

p_D2_gam<-ggplot()+
  geom_stars(data=newdata_D2_pred_rast)+
  facet_wrap(~month(Date, label=T), drop=F)+
  scale_fill_continuous_diverging(palette="Blue-Red 3", breaks=seq(-12, 6, by=2)/10, na.value = "white",
                                  name="Temperature change\nper change in\ntotal inflow (°C/monthly sd[cfs])", guide=guide_colorbar(barheight=20))+
  ylab("Latitude")+
  xlab("Longitude")+
  coord_sf()+
  theme_bw()+
  theme(strip.background=element_blank(), axis.text.x = element_text(angle=45, hjust=1), panel.grid=element_blank())

ggsave(p_D2_gam, filename="C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/figures/D2_gam 02.01.21.png",
       device="png", width=7, height=5, units="in")

## Magnitude of result too large, try centering and/or scaling inflow separately for each month.