#devtools::install_github("sbashevkin/discretewq", ref="v1.1.0")

require(sp)
require(gstat)
require(spacetime)
require(tidybayes)
require(dplyr)
require(stringr)
require(dtplyr)
require(tidyr)
require(mgcv)
require(ggplot2)
require(geofacet)
require(lubridate)
require(hms)
require(sf)
require(stars)
require(purrr)
require(scales)
require(itsadug)
require(colorspace)
require(patchwork)
require(ggstance)
require(discretewq)
source("Utility_functions.R")

mygrid <- data.frame(
  name = c("Upper Sacramento River Ship Channel", "Cache Slough and Lindsey Slough", "Lower Sacramento River Ship Channel", "Liberty Island", "Suisun Marsh", "Middle Sacramento River", "Lower Cache Slough", "Steamboat and Miner Slough", "Upper Mokelumne River", "Lower Mokelumne River", "Georgiana Slough", "Sacramento River near Ryde", "Sacramento River near Rio Vista", "Grizzly Bay", "West Suisun Bay", "Mid Suisun Bay", "Honker Bay", "Confluence", "Lower Sacramento River", "San Joaquin River at Twitchell Island", "San Joaquin River at Prisoners Pt", "Disappointment Slough", "Lower San Joaquin River", "Franks Tract", "Holland Cut", "San Joaquin River near Stockton", "Mildred Island", "Middle River", "Old River", "Upper San Joaquin River", "Grant Line Canal and Old River", "Victoria Canal"),
  row = c(2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 8, 8, 8),
  col = c(7, 4, 6, 5, 2, 8, 6, 7, 9, 9, 8, 7, 6, 2, 1, 2, 3, 4, 5, 6, 8, 9, 5, 6, 7, 9, 8, 8, 7, 9, 8, 7),
  code = c(" 1", " 1", " 2", " 3", " 8", " 4", " 5", " 6", " 7", " 9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "30", "29", "31"),
  stringsAsFactors = FALSE
)

Data<-readRDS("Temperature analysis/Discrete Temp Data.Rds")

# GAMs --------------------------------------------------------------------

Data_CC4<-readRDS("Temperature analysis/Data_CC4.Rds")

# Split time-series into groupings of contiguous months within each station.
# Split into separate time-series after a gap of 60 days or greater
Data_CC4.3 <- Data_CC4%>%
  arrange(Station, Date)%>%
  group_by(Station)%>%
  mutate(Lag=Date-lag(Date, order_by = Date))%>%
  ungroup()%>%
  mutate(Start=if_else(is.na(Lag) | Lag>3600*24*60, TRUE, FALSE),
         Series_ID=1:n(),
         Series_ID=if_else(Start, Series_ID, NA_integer_),
         Series_ID=as.integer(as.factor(Series_ID)))%>%
  fill(Series_ID, .direction="down")

# First fit model with AR term to find optimal AR rho parameter

CC_gam8d7b_NOAR5 <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                          te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=Year_s) + 
                          s(Time_num_s, k=5), family=scat, data = Data_CC4.3, method="fREML", discrete=T, nthreads=4)
r6 <- start_value_rho(CC_gam8d7b_NOAR5, plot=TRUE)

#########Best Model####################

CC_gam8d7b_AR7 <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                        te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=Year_s) + 
                        s(Time_num_s, k=5), family=scat, rho=r6, AR.start=Start, data = Data_CC4.3, method="fREML", discrete=T, nthreads=4)
#AIC: 199876.5
#BIC: 202706.4

#########Best Model####################

# Try higher Rho
CC_gam8d7b_AR7B <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                        te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=Year_s) + 
                        s(Time_num_s, k=5), family=scat, rho=r6*2, AR.start=Start, data = Data_CC4.3, method="fREML", discrete=T, nthreads=4)
#AIC: 192436.7


# Create spatiotemporal variogram -----------------------------------------

p_variogram_AR7B<-ST_variogram(CC_gam8d7b_AR7, Data_CC4.3, 5)

ggsave(p_variogram, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Temperature change/Figures/Figure S6 variogram.png",
       device="png", width=8, height=5, units="in")

# Trying higher rho

p_variogram<-ST_variogram(CC_gam8d7b_AR7B, Data_CC4.3, 5)

## Now gamma is too high for time points > 1 lag. So sticking with original Rho
  
#timelag is the average of successive specified tlags. It must be computing correlation between tlags. 

gam.check(CC_gam8d7b_AR7)
# s(Time_num_s) is OK (p=0.33)
# te(Julian_day_s,Latitude_s,Longitude_s):Year_s is OK: significant p-value but edf is 82.45 compared to k' of 300.0
# te(Julian_day_s,Latitude_s,Longitude_s) may be able to be improved (edf=222.36 and k'=299.0)


# Find SubRegion and Month combinations with representation in the data
# For filtering the newdata after predictions
Data_CC4_effort<-Data_CC4%>%
  distinct(SubRegion, Month)%>%
  mutate(Keep=TRUE) 

newdata<-readRDS("Temperature analysis/Prediction Data.Rds")
CC_newdata<-newdata%>%
  st_drop_geometry()%>%
  as_tibble()%>%
  select(-Year_fac, -Year, -Year_s, -N,)%>%
  distinct()%>%
  mutate(Year=2000,
         Year_s=(Year-mean(Data$Year))/sd(Data$Year),
         Year_fac="2000",
         Month=as.integer(as.factor(Julian_day)))%>%
  left_join(Data_CC4_effort, by=c("Month", "SubRegion"))%>%
  filter(Keep)%>% # Only retain SubRegion and Month combinations with representation in the data
  select(-Keep)

#saveRDS(CC_newdata, "Temperature analysis/Climate Change Prediction Data.Rds")
CC_newdata<-readRDS("Temperature analysis/Climate Change Prediction Data.Rds")

CC_pred<-predict(CC_gam8d7b_AR7, newdata=CC_newdata, type="terms", se.fit=TRUE, discrete=T, n.threads=4)

Delta<-st_read("Delta Subregions")%>%
  select(SubRegion)%>%
  filter(SubRegion%in%unique(Data_CC4$SubRegion))

newdata_CC_pred<-CC_newdata%>%
  mutate(Slope=CC_pred$fit[,"te(Julian_day_s,Latitude_s,Longitude_s):Year_s"],
         Slope_se=CC_pred$se.fit[,"te(Julian_day_s,Latitude_s,Longitude_s):Year_s"])%>%
  mutate(across(c(Slope, Slope_se), ~(.x/Year_s)/sd(Data$Year)))%>%
  mutate(Slope_se=abs(Slope_se))%>%
  mutate(Date=as.Date(Julian_day, origin=as.Date(paste(Year, "01", "01", sep="-"))),
         Month=month(Date, label = T),
         Slope_l99=Slope-Slope_se*qnorm(0.9995),
         Slope_u99=Slope+Slope_se*qnorm(0.9995),
         Sig=if_else(Slope_u99>0 & Slope_l99<0, "ns", "*"))

newdata_CC_pred_rast<-newdata_CC_pred%>%
  filter(Sig=="*")%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=F)%>%
  st_transform(crs=st_crs(Delta))%>%
  Rasterize_all(., Slope, region=Delta)

# Create background raster of all locations
base<-CC_newdata%>%
  mutate(Date=parse_date_time(paste(Year, Month, "15", sep="-"), "%Y-%m-%d"))%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)%>%
  st_transform(crs=st_crs(Delta))%>%
  Rasterize_all(Location, region=Delta)%>%
  st_as_sf(long=T, connect8=T)%>%
  filter(!is.na(Location))

p_CC_gam<-predict_plot(data=newdata_CC_pred_rast, 
             base=base, 
             scale_fill_continuous_diverging, 
             guide=guide_colorbar(barheight=20), 
             palette="Blue-Red 3",
             limits=c(0, 0.057),
             na.value=NA, 
             breaks=(-6:7)/100, 
             name="Temperature change\nper year (°C)")

ggsave(p_CC_gam, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Temperature change/Figures/Figure 3 Climate change signal.tiff",
       device="tiff", dpi=400, width=7, height=5, units="in")

# Plot all slopes, significant or not

newdata_CC_pred_all_rast<-newdata_CC_pred%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=F)%>%
  st_transform(crs=st_crs(Delta))%>%
  Rasterize_all(., Slope, region=Delta)

p_CC_gam_all<-predict_plot(data=newdata_CC_pred_all_rast, 
                           base=base, 
                           scale_fill_continuous_diverging, 
                           guide=guide_colorbar(barheight=20), 
                           palette="Blue-Red 3",
                           limits=c(-0.015, 0.057),
                           na.value=NA, 
                           breaks=(-6:7)/100, 
                           name="Temperature change\nper year (°C)")

ggsave(p_CC_gam_all, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Temperature change/Figures/Figure S7 Climate change signal all.png",
       device="png", width=7, height=5, units="in")

# Create dataframe of slope deviations for each month and region f --------

Slope_summary<-newdata_CC_pred%>%
  group_by(Month, SubRegion)%>%
  summarise(Slope=mean(Slope), Slope_se=mean(Slope_se), Sig_prop=length(which(Sig=="*"))/n(), .groups="drop")%>%
  mutate(Slope_mean=mean(Slope), 
         Slope_mult=Slope/Slope_mean,
         Slope_se_mean=mean(Slope_se), 
         Slope_se_mult=Slope/Slope_se_mean)
#saveRDS(Slope_summary, "Temperature analysis/Slope summary.Rds")

# Plot slope summary for each region and month
p_slope_sum<-ggplot(Slope_summary)+
  geom_pointrange(aes(y=reorder(Month, desc(Month)), x=Slope, xmax=Slope+Slope_se, xmin=Slope-Slope_se, fill=Sig_prop, color=Sig_prop))+
  facet_geo(~SubRegion, grid=mygrid, labeller=label_wrap_gen(width=15))+
  geom_vline(xintercept=0)+
  scale_y_discrete(breaks=c("Jan", "Mar", "May", "Jul", "Sep", "Nov"))+
  scale_color_viridis_c(name="Proportion\nsignificant slopes", limits=c(0,1), aesthetics = c("fill", "color"),
                        guide=guide_colorbar(barheight=15))+
  ylab("Month")+
  xlab("Temperature change per year (°C)")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1), text=element_text(size=16), legend.position=c(0.4, 0.65), 
        legend.background = element_rect(color="black"), panel.background = element_rect(color="black"), legend.margin=margin(10,10,15,10))

ggsave(p_slope_sum, file="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Temperature change/Figures/Climate change slopes.png",
       device="png", width=15, height=18, units="in")

# Plot all slopes with CIs

SubRegion_levels<-c("Upper Sacramento River Ship Channel", "Middle Sacramento River",
                    "Lower Sacramento River Ship Channel", "Steamboat and Miner Slough",
                    "Cache Slough and Lindsey Slough", "Liberty Island",
                    "Lower Cache Slough", "Sacramento River near Ryde",
                    "Georgiana Slough", "Upper Mokelumne River",
                    "Suisun Marsh", "West Suisun Bay",
                    "Grizzly Bay", "Mid Suisun Bay",
                    "Honker Bay", "Confluence",
                    "Lower Sacramento River", "Sacramento River near Rio Vista",
                    "San Joaquin River at Twitchell Island", "San Joaquin River at Prisoners Pt",
                    "Lower Mokelumne River", "Disappointment Slough",
                    "Lower San Joaquin River", "Franks Tract", 
                    "Holland Cut", "Mildred Island",
                    "San Joaquin River near Stockton", "Old River",
                    "Middle River", "Grant Line Canal and Old River",
                    "Victoria Canal", "Upper San Joaquin River")

Slope_sum2<-newdata_CC_pred%>%
  group_by(Month, SubRegion)%>%
  summarise(Slope_mean=mean(Slope), Slope_u99_max=max(Slope_u99), Slope_u99_min=min(Slope_u99), Slope_l99_max=min(Slope_l99), Slope_l99_min=max(Slope_l99), .groups="drop")%>%
  mutate(SubRegion=factor(SubRegion, levels=SubRegion_levels))

P_slope_sum2<-ggplot(Slope_sum2)+
  geom_vline(xintercept=0)+
  geom_linerange(aes(y=reorder(Month, desc(Month)), xmax=Slope_l99_max, xmin=Slope_u99_max), size=1, color="#d7191c")+
  geom_linerange(aes(y=reorder(Month, desc(Month)), xmax=Slope_l99_min, xmin=Slope_u99_min), size=2, color="#fdae61")+
  geom_point(aes(y=reorder(Month, desc(Month)), x=Slope_mean), size=1, color="#2c7bb6")+
  facet_geo(~SubRegion, grid=mygrid, labeller=label_wrap_gen(width=18))+
  scale_y_discrete(breaks=c("Jan", "Mar", "May", "Jul", "Sep", "Nov"))+
  ylab("Month")+
  xlab("Temperature change per year (°C)")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1), text=element_text(size=16), panel.background = element_rect(color="black"))

ggsave(P_slope_sum2, file="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Temperature change/Figures/Figure S8 Climate change all slopes.png",
       device="png", width=15, height=18, units="in")

# Plot sampling effort for each region, month. and year

Data_effort<-Data_CC4%>%
  group_by(Year, SubRegion, Month)%>%
  summarise(N=n(), .groups="drop")

p_effort<-ggplot(Data_effort)+
  geom_tile(aes(x=Year, y=reorder(month(Month, label=T), desc(month(Month, label=T))), fill=N))+
  scale_fill_viridis_c(breaks=seq(0,140, by=10),
                       guide=guide_colorbar(barheight=15, barwidth = 3))+
  scale_x_continuous(breaks=seq(1970, 2020, by=10))+
  scale_y_discrete(breaks=c("Jan", "Mar", "May", "Jul", "Sep", "Nov"))+
  facet_geo(~SubRegion, grid=mygrid, labeller=label_wrap_gen(width=18))+
  ylab("Month")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1), panel.grid=element_blank(), text=element_text(size=16), legend.position=c(0.4, 0.65), 
        legend.background = element_rect(color="black"), panel.background = element_rect(color="black"), legend.margin=margin(10,10,15,10))

ggsave(p_effort, file="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Temperature change/Figures/Figure S1 Climate change effort.png",
       device="png", width=15, height=18, units="in")

# Now try a model with higher spatial K value -----------------------------

# Check if higher k on te(Julian_day_s,Latitude_s,Longitude_s) would help improve model

CC_gam8d7b_AR8 <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(50, 13)) + 
                        te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=Year_s) + 
                        s(Time_num_s, k=5), family=scat, rho=r6, AR.start=Start, data = Data_CC4.3, method="fREML", discrete=T, nthreads=4)
#AIC: 199510.1
#BIC: 203406.4

gam.check(CC_gam8d7b_AR8)

CC_pred_AR8<-predict(CC_gam8d7b_AR8, newdata=CC_newdata, type="terms", se.fit=TRUE, discrete=T, n.threads=4)

newdata_CC_pred_AR8<-CC_newdata%>%
  mutate(Slope=CC_pred_AR8$fit[,"te(Julian_day_s,Latitude_s,Longitude_s):Year_s"],
         Slope_se=CC_pred_AR8$se.fit[,"te(Julian_day_s,Latitude_s,Longitude_s):Year_s"])%>%
  mutate(across(c(Slope, Slope_se), ~(.x/Year_s)/sd(Data$Year)))%>%
  mutate(Slope_se=abs(Slope_se))%>%
  mutate(Date=as.Date(Julian_day, origin=as.Date(paste(Year, "01", "01", sep="-"))),
         Month=month(Date, label = T),
         Slope_l99=Slope-Slope_se*qnorm(0.9995),
         Slope_u99=Slope+Slope_se*qnorm(0.9995))%>%
  mutate(Sig=if_else(Slope_u99>0 & Slope_l99<0, "ns", "*"))%>%
  filter(Sig=="*")%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=F)%>%
  st_transform(crs=st_crs(Delta))

newdata_CC_pred_rast_AR8<-Rasterize_all(newdata_CC_pred_AR8, Slope, region=Delta)

p_CC_gam_AR8<-predict_plot(data=newdata_CC_pred_rast_AR8, 
                           base=base, 
                           scale_fill_continuous_diverging, 
                           guide=guide_colorbar(barheight=20), 
                           palette="Blue-Red 3",
                           limits=c(0, 0.057),
                           na.value=NA, 
                           breaks=(-6:7)/100, 
                           name="Temperature change\nper year (°C)")

ggsave(p_CC_gam_AR8, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Temperature change/Figures/Figure S5 Climate change signal_higherk.png",
       device="png", width=7, height=5, units="in")

###CC_gam8d7b_AR8 has a lower AIC but higher AIC and  predicted slope values are almost identical to CC_gam8d7b_AR7, so using CC_gam8d7b_AR7 as the final model

# Recreate gam check plots ------------------------------------------------

p_check<-model_validation(CC_gam8d7b_AR7, Data_CC4.3$Temperature)

ggsave(p_check, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Temperature change/Figures/Figure 2 Climate change model validation.tiff",
       device="tiff", dpi=400, width=10, height=7, units="in")

# Fit separate models to each 25 year period ------------------------------

Data<-readRDS("Temperature analysis/Discrete Temp Data.Rds")

Year_mean<-mean(Data$Year)
Year_sd<-sd(Data$Year)
Data_CC4<-readRDS("Temperature analysis/Data_CC4.Rds")


Delta<-st_read("Delta Subregions")%>%
  select(SubRegion)%>%
  filter(SubRegion%in%unique(Data_CC4$SubRegion))

start_years<-set_names(seq(1970, 1995, by=5))
period_length<-25

model_fitter<-function(start_year){
  
  data<-Data_CC4%>%
    filter(Year>=start_year & Year<(start_year+period_length))%>%
    mutate(Year_s=Year_s-mean(Year_s))%>%
    arrange(Station, Date)%>%
    group_by(Station)%>%
    mutate(Lag=Date-lag(Date, order_by = Date))%>%
    ungroup()%>%
    mutate(Start=if_else(is.na(Lag) | Lag>3600*24*60, TRUE, FALSE),
           Series_ID=1:n(),
           Series_ID=if_else(Start, Series_ID, NA_integer_),
           Series_ID=as.integer(as.factor(Series_ID)))%>%
    fill(Series_ID, .direction="down")
  
  noAR <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                     te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=Year_s) + 
                     s(Time_num_s, k=5), family=scat, data = data, method="fREML", discrete=T, nthreads=4)
  r <- start_value_rho(noAR, plot=TRUE)
  
  AR <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                          te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=Year_s) + 
                          s(Time_num_s, k=5), family=scat, rho=r, AR.start=Start, data = data, method="fREML", discrete=T, nthreads=4)
  
  message(paste(start_year, "finished"))
  return(AR)
}

models_period<-map(start_years, model_fitter)

# Find out which regions and months have data in each 25 year period
Data_period_effort<-Data_CC4%>%
  group_by(SubRegion, Month)%>%
  summarise(Sampled_1970=length(which(Year>=1970 & Year<1995)),
            Sampled_1975=length(which(Year>=1975 & Year<2000)),
            Sampled_1980=length(which(Year>=1980 & Year<2005)),
            Sampled_1985=length(which(Year>=1985 & Year<2010)),
            Sampled_1990=length(which(Year>=1990 & Year<2015)),
            Sampled_1995=length(which(Year>=1995 & Year<2020)), .groups="drop")%>%
  mutate(across(all_of(paste("Sampled", start_years, sep="_")), ~if_else(.x>1, TRUE, FALSE)))

newdata<-readRDS("Temperature analysis/Prediction Data.Rds")
CC_newdata_period<-newdata%>%
  st_drop_geometry()%>%
  as_tibble()%>%
  select(-Year_fac, -Year, -Year_s, -N,)%>%
  distinct()%>%
  mutate(Month=as.integer(as.factor(Julian_day)))%>%
  left_join(Data_period_effort, by=c("Month", "SubRegion"))

model_predictor<-function(start_year, sig="filter"){
  CC_newdata<-CC_newdata_period%>%
    filter(across(all_of(paste0("Sampled_", start_year))))%>%
    mutate(Year=start_year,
           Year_s=1)
  
  CC_pred<-predict(models_period[[as.character(start_year)]], newdata=CC_newdata, type="terms", se.fit=TRUE, discrete=T, n.threads=2)
  
  newdata_CC_pred<-CC_newdata%>%
    mutate(Slope=CC_pred$fit[,"te(Julian_day_s,Latitude_s,Longitude_s):Year_s"],
           Slope_se=CC_pred$se.fit[,"te(Julian_day_s,Latitude_s,Longitude_s):Year_s"],
           Intercept=CC_pred$fit[,"te(Julian_day_s,Latitude_s,Longitude_s)"]+CC_pred$fit[,"s(Time_num_s)"])%>%
    mutate(across(c(Slope, Slope_se), ~(.x/Year_s)/Year_sd))%>%
    mutate(Slope_se=abs(Slope_se))%>%
    mutate(Date=as.Date(Julian_day, origin=as.Date(paste(Year, "01", "01", sep="-"))),
           Month=month(Date, label = T),
           Slope_l99=Slope-Slope_se*qnorm(0.9995),
           Slope_u99=Slope+Slope_se*qnorm(0.9995))%>%
    mutate(Sig=if_else(Slope_u99>0 & Slope_l99<0, "ns", "*"))%>%
    {if(sig=="filter"){
      filter(., Sig=="*")
    }else{
      .
    }}%>%
    st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=F)%>%
    st_transform(crs=st_crs(Delta))
  
  return(newdata_CC_pred)
}

preds_period<-map_dfr(start_years, model_predictor)

preds_period_rast<-Rasterize_all(preds_period, Slope, region=Delta)

base_period<-CC_newdata_period%>%
  pivot_longer(all_of(paste0("Sampled_", start_years)), names_to="Period", values_to="Sampled")%>%
  mutate(Year=recode(Period, !!!set_names(start_years, paste0("Sampled_", start_years))),
         Date=as.Date(Julian_day, origin=as.Date(paste(Year, "01", "01", sep="-"))))%>%
  filter(Sampled)%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)%>%
  st_transform(crs=st_crs(Delta))%>%
  Rasterize_all(Location, region=Delta)%>%
  st_as_sf(long=T)%>%
  filter(!is.na(Location))

p_period<-ggplot()+
  geom_sf(data=base_period, color=NA, fill="gray80", lwd = 0)+
  geom_stars(data=preds_period_rast)+
  facet_grid(month(Date)~year(Date), drop=F, labeller=as_labeller(c(set_names(as.character(month(1:12, label=T)), 1:12), set_names(paste0(start_years, "-", start_years+25), start_years))))+
  scale_fill_continuous_diverging(palette="Blue-Red 3", name="Temperature change\nper year (°C)", breaks=seq(-0.05, 0.15, by=0.05), guide=guide_colorbar(barheight=20), na.value=NA)+
  ylab("Latitude")+
  xlab("Longitude")+
  coord_sf()+
  theme_bw()+
  theme(strip.background=element_blank(), text=element_text(size=7), axis.text = element_blank(), axis.ticks=element_blank(), panel.grid=element_blank(),
        panel.spacing=unit(0, "lines"))

ggsave(p_period, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Temperature change/Figures/Figure 4 Climate change signal over time rasters.tiff",
       device="tiff", dpi=400, width=7, height=7, units="in")

###Maybe instead of rasters just plot the slope over time for each region

preds_period_sum<-preds_period%>%
  group_by(Year, SubRegion, Month)%>%
  summarise(Slope=mean(Slope), .groups="drop")%>%
  complete(Month, SubRegion, Year)%>%
  mutate(Years=paste0(Year, "-", Year+25))

ggplot(preds_period_sum)+
  geom_line(aes(x=Years, y=Slope, color=Month, group=Month))+
  geom_point(aes(x=Years, y=Slope, color=Month))+
  geom_hline(yintercept = 0, linetype=2)+
  facet_geo(~SubRegion, grid=mygrid, labeller=label_wrap_gen())+
  #scale_color_brewer(palette="Paired")+
  theme_bw()+
  theme(panel.grid=element_blank(), axis.text.x = element_text(angle=45, hjust=1))

## Plot slope over time for each month

preds_period_month_sum<-map_dfr(start_years, model_predictor, sig="no")%>%
  st_drop_geometry()%>%
  group_by(Year, Month)%>%
  summarise(Slope_mean=mean(Slope), Slope_sd=sd(Slope), Sig_prop=length(which(Sig=="*"))/n(), .groups="drop")%>%
  mutate(Years=factor(paste0(Year, "-", Year+25), levels=paste0(start_years, "-", start_years+25)))

p_period_slope_sum<-ggplot(preds_period_month_sum)+
  geom_hline(yintercept=0)+
  geom_pointrange(aes(x=Years, y=Slope_mean, ymin=Slope_mean-Slope_sd, ymax=Slope_mean+Slope_sd, fill=Sig_prop), shape=21, color="black")+
  facet_wrap(~Month)+
  scale_fill_viridis_c(name="Proportion\nsignificant")+
  ylab("Temperature change per year (°C)")+
  xlab("")+
  theme_bw()+
  theme(strip.background = element_blank(), axis.text.x=element_text(angle=45, hjust=1))

ggsave(p_period_slope_sum, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Temperature change/Figures/Climate change signal over time summary_BDSC.png",
       device="png", width=6, height=6, units="in")

# Plot climate change signal in each priority restoration area from the Delta Plan --------

## First plot a map of those areas ----------------------------------------
PHRA<-st_read("200813_DeMartino_GIS_Export", layer = "ER_P3")%>%
  select(PHRA=Region_)%>%
  mutate(PHRA=recode(PHRA, `Cosumunes-Mokelumne`="Cosumnes - Mokelumne"))%>%
  mutate(PHRA=factor(PHRA, levels=c("Lower San Joaquin River Floodplain", "Western Delta", "Suisun Marsh", "Cosumnes - Mokelumne", "Cache Slough", "Yolo Bypass")))

ggplot()+
  geom_sf(data=Delta)+
  geom_sf(data=PHRA, aes(fill=PHRA))+
  geom_sf(data=deltamapr::WW_Delta%>%st_transform(crs=st_crs(Delta)))

PHRA_grid<-base%>%
  st_transform(crs=st_crs(PHRA))%>%
  st_join(PHRA%>%
            select(PHRA))%>%
  st_transform(crs=4326)

ggplot()+
  geom_sf(data=PHRA_grid, aes(color=PHRA, fill=PHRA))+
  facet_wrap(~month(Date, label=T))

boundary<-Delta%>%
  st_union()%>%
  st_boundary()%>%
  st_polygonize()

boundary<-data.frame(boundary[[1]][[1]][[1]])%>%
  rename(X=X1, Y=X2)
base_PHRA<-deltamapr::WW_Delta%>%
  st_transform(crs=st_crs(Delta))%>%
  st_crop(st_union(st_as_sf(boundary, coords=c("X", "Y"), crs=st_crs(Delta))
                   %>%st_make_valid(), 
                   PHRA%>%
                     st_make_valid()))

colors<-RColorBrewer::brewer.pal(7,"Set1")

p_map<-ggplot()+
  geom_sf(data=base_PHRA, fill="slategray1", color="slategray2")+
  geom_sf(data=PHRA, aes(fill=PHRA, color=PHRA), alpha=0.2)+
  geom_path(data=boundary, color="black", aes(x=X, y=Y))+
  scale_fill_brewer(palette="Dark2", aesthetics = c("color", "fill"), labels = function(x) str_wrap(x, width = 18),
                    name="Priority Habitat\nRestoration Area", guide=guide_legend(reverse = TRUE))+
  ylab("")+
  xlab("")+
  theme_bw()+
  theme(legend.position="none", text=element_text(size=18), legend.background=element_rect(color="black"), plot.margin = margin(r=30))

## Now Plot a more summarized version of the slopes for final figure

PHRA_slopes<-newdata_CC_pred%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)%>%
  st_transform(crs=st_crs(PHRA))%>%
  st_join(PHRA)%>%
  st_drop_geometry()%>%
  filter(!is.na(PHRA))%>%
  group_by(Month, PHRA)%>%
  summarise(Slope_mean=mean(Slope), Slope_u99_max=max(Slope_u99), Slope_u99_min=min(Slope_u99), 
            Slope_l99_max=min(Slope_l99), Slope_l99_min=max(Slope_l99), 
            Slope_sd=sd(Slope), Sig_prop=length(which(Sig=="*"))/n(), .groups="drop")%>%
  mutate(PHRA=factor(PHRA, levels=c("Lower San Joaquin River Floodplain", "Western Delta", "Suisun Marsh", "Cosumnes - Mokelumne", "Cache Slough", "Yolo Bypass")))

p_PHRA_sum<-ggplot(PHRA_slopes%>%mutate(PHRA=factor(PHRA, levels=rev(levels(PHRA_slopes$PHRA)))))+
  geom_vline(xintercept=0)+
  geom_pointrange(aes(y=Month, x=Slope_mean, xmin=Slope_mean-Slope_sd, xmax=Slope_mean+Slope_sd, fill=Sig_prop), shape=21, size=1, color="black")+
  facet_wrap(~PHRA, ncol=1)+
  scale_y_discrete(limits=rev, breaks=month(1:12, label=T)[c(1,3,5,7,9,11)])+
  xlab("Temperature change per year (°C)")+
  ylab("Month")+
  scale_fill_viridis_c(name="Proportion\nsignificant")+
  theme_bw()+
  theme(text=element_text(size=18), strip.text=element_text(color="white"))

g <- ggplot_gtable(ggplot_build(p_PHRA_sum))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- RColorBrewer::brewer.pal(6, "Dark2")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

g_ggplot<-ggplotify::as.ggplot(g)
#ggsave(g_ggplot, file="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Temperature change/Figures/PHRA_sum.png",
#       width=5, height=8, units="in")

p_PHRA<-p_map+g_ggplot+plot_layout(ncol=2, widths = c(1,1))+ plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size=18))

ggsave(p_PHRA, file="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Temperature change/Figures/Figure 5 PHRA.tiff",
       device="tiff", dpi=400, width=14, height=10, units="in")

## Alternative plot of the climate change signal in each area ------------------------------------

P_PHRA_slopes<-ggplot(PHRA_slopes)+
  geom_vline(xintercept=0)+
  geom_linerange(aes(y=PHRA, xmax=Slope_l99_max, xmin=Slope_u99_max), size=1, color="#d7191c")+
  geom_linerange(aes(y=PHRA, xmax=Slope_l99_min, xmin=Slope_u99_min), size=2, color="#fdae61")+
  geom_point(aes(y=PHRA, x=Slope_mean), size=1, color="#2c7bb6")+
  facet_wrap(~month(Month, label=T))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 18))+
  xlab("Temperature change per year (°C)")+
  ylab("Priority Habitat Restoration Area")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1), text=element_text(size=18), panel.background = element_rect(color="black"), strip.background = element_blank())

p_PHRA<-p_map+P_PHRA_slopes+plot_layout(ncol=2, widths = c(1,1))+ plot_annotation(tag_levels = "A")


# Plot our climate change signal relative to other systems ----------------

Slope_stats<-CC_newdata%>%
  mutate(Slope=CC_pred$fit[,"te(Julian_day_s,Latitude_s,Longitude_s):Year_s"],
         Slope_se=CC_pred$se.fit[,"te(Julian_day_s,Latitude_s,Longitude_s):Year_s"])%>%
  mutate(across(c(Slope, Slope_se), ~(.x/Year_s)/sd(Data$Year)))%>%
  mutate(Slope_se=abs(Slope_se))

CC_trends<-tibble(Location=c("San Francisco Estuary\naverage (our study)", "San Francisco Estuary Air", "Mallard Island", "Confluence (our study)", "Hudson River", "Chesapeake", "Narragansett Bay"),
                  Slope=c(mean(Slope_stats$Slope), NA_real_, 0.007, mean(filter(Slope_stats, SubRegion=="Confluence")$Slope), 0.015, rep(NA_real_, 2)),
                  Slope_sd=c(sd(Slope_stats$Slope), rep(NA_real_, 2), sd(filter(Slope_stats, SubRegion=="Confluence")$Slope), rep(NA_real_, 3)),
                  Min=c(NA_real_, 0.009, rep(NA_real_, 3), 0.021, 0.027),
                  Max=c(NA_real_, 0.016, rep(NA_real_, 3), 0.04, 0.032))%>%
  mutate(Location=factor(Location, levels=Location))

p_other_slopes<-ggplot(CC_trends)+
  geom_point(aes(x=Location, y=Slope))+
  geom_linerange(aes(x=Location, ymin=Slope-Slope_sd, ymax=Slope+Slope_sd))+
  geom_point(aes(x=Location, y=Min))+
  geom_point(aes(x=Location, y=Max))+
  geom_segment(aes(x=Location, xend=Location, y=Min, yend=Max))+
  scale_y_continuous(limits=c(0,0.05), expand = expansion(0,0))+
  ylab("Temperature change per year (°C)")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1), panel.grid=element_blank())

ggsave(p_other_slopes, file="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Temperature change/Figures/Slope comparison_BDSC.png",
       device="png", width=4, height=4, units="in")