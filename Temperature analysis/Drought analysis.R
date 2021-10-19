require(sp)
require(gstat)
require(spacetime)
require(dplyr)
require(tidyr)
require(stringr)
require(dtplyr)
require(mgcv)
require(ggplot2)
require(geofacet)
require(lubridate)
require(hms)
require(scales)
require(patchwork)
require(sf)
require(stars)
require(itsadug)
require(purrr)
require(readr)
require(slider)
require(colorspace)
require(ggstance)
require(discretewq)
source("Utility_functions.R")
options(scipen=999)

# setup -------------------------------------------------------------------

mygrid <- data.frame(
  name = c("Upper Sacramento River Ship Channel", "Cache Slough and Lindsey Slough", "Lower Sacramento River Ship Channel", "Liberty Island", "Suisun Marsh", "Middle Sacramento River", "Lower Cache Slough", "Steamboat and Miner Slough", "Upper Mokelumne River", "Lower Mokelumne River", "Georgiana Slough", "Sacramento River near Ryde", "Sacramento River near Rio Vista", "Grizzly Bay", "West Suisun Bay", "Mid Suisun Bay", "Honker Bay", "Confluence", "Lower Sacramento River", "San Joaquin River at Twitchell Island", "San Joaquin River at Prisoners Pt", "Disappointment Slough", "Lower San Joaquin River", "Franks Tract", "Holland Cut", "San Joaquin River near Stockton", "Mildred Island", "Middle River", "Old River", "Upper San Joaquin River", "Grant Line Canal and Old River", "Victoria Canal"),
  row = c(2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 8, 8, 8),
  col = c(7, 4, 6, 5, 2, 8, 6, 7, 9, 9, 8, 7, 6, 2, 1, 2, 3, 4, 5, 6, 8, 9, 5, 6, 7, 9, 8, 8, 7, 9, 8, 7),
  code = c(" 1", " 1", " 2", " 3", " 8", " 4", " 5", " 6", " 7", " 9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "30", "29", "31"),
  stringsAsFactors = FALSE
)

# Now using total inflow (QTOT) from Dayflow ------------------------------

## Create dataset ---------------------------------------------------------
dayflow<-read_csv(file.path("Temperature analysis", "data", "Dayflow1997 2019.csv"), col_types = cols_only(Date="c", TOT="d", PREC="d"))%>%
  bind_rows(read_csv(file.path("Temperature analysis", "data", "Dayflow1984 1996.csv"), col_types = cols_only(Date="c", TOT="d", PREC="d")))%>%
  bind_rows(read_csv(file.path("Temperature analysis", "data", "Dayflow1970 1983.csv"), col_types = cols_only(Date="c", TOT="d", PREC="d")))%>%
  bind_rows(read_csv(file.path("Temperature analysis", "data", "Dayflow1956 1969.csv"), col_types = cols_only(Date="c", TOT="d", PREC="d")))%>%
    mutate(Date=parse_date_time(Date, "%m/%d/%Y", tz = "America/Los_Angeles"),
           PREC_acrefeetperday=(PREC*3600*24)/(43560),
           PREC_feetperday=if_else(Date>=parse_date_time("10/01/1980", "%m/%d/%Y", tz = "America/Los_Angeles"), PREC_acrefeetperday/682230, PREC_acrefeetperday/738000),
           PREC_CFS=(PREC_feetperday*682230*43560)/(3600*24),
           Date = as_date(Date))%>% # Remove times to make sure merging works correctly
  arrange(Date)%>%
  mutate(TOT_mean30=slide_index_dbl(.$TOT, .$Date, .before=days(30), .f=mean, .complete = T),
         PREC_mean30=slide_index_dbl(.$PREC, .$Date, .before=days(30), .f=mean, .complete = T),
         Month=month(Date))%>%
  group_by(Month)%>%
  mutate(across(c(TOT_mean30, PREC_mean30), list(sd= ~sd(.x, na.rm=T), mean=~mean(.x, na.rm=T))))%>%
  ungroup()%>%
  select(-Month)

Data_D2<-readRDS("Temperature analysis/Data_CC4.Rds")%>%
  mutate(WY=Year)%>%
  mutate(WY=if_else(Month%in%10:12, WY+1, WY),
         Date=as_date(Date))%>% # Remove times to make sure merging works correctly
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
  mutate(across(c(TOT_mean30, PREC_mean30), list(s=~(.x-mean(.x))/sd(.x))))%>%
  mutate(WY_s=(WY-mean(unique(WY)))/sd(unique(WY)),
         TOT_mean30_s_month=(TOT_mean30-TOT_mean30_mean)/TOT_mean30_sd,
         PREC_mean30_s_month=(PREC_mean30-PREC_mean30_mean)/PREC_mean30_sd)

#saveRDS(Data_D2, "Temperature analysis/Data_D2.Rds")
  

Delta2<-st_read("Delta Subregions")%>%
  select(SubRegion)%>%
  filter(SubRegion%in%unique(Data_D2$SubRegion))


## First fit non-month adjusted models -------------------------------------

D2_nomonth_gam_NOAR <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                     te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=TOT_mean30_s) + 
                     s(Time_num_s, k=5), family=scat, data = Data_D2, method="fREML", discrete=T, nthreads=3)
D2_nomonthr <- start_value_rho(D2_nomonth_gam_NOAR, plot=TRUE)

D2_nomonth_gam_AR <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                   te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=TOT_mean30_s) + 
                   s(Time_num_s, k=5), family=scat, rho=D2_nomonthr, AR.start=Start, data = Data_D2, method="fREML", discrete=T, nthreads=3)

## Autocorrelation ---------------------------------------------------------

p_D2_nomonth_variogram<-ST_variogram(D2_nomonth_gam_AR, Data_D2, 5)

ggsave(p_D2_nomonth_variogram, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Delta Inflow temperature/Figures/model 2 inflow nomonth model variogram.png",
       device="png", width=8, height=5, units="in")

### Set up model predictions -------------------------------------------------

newdata<-readRDS("Temperature analysis/Prediction Data.Rds")

##### Find SubRegion and Month combinations with representation in the data
Data_D2_effort<-Data_D2%>%
  distinct(SubRegion, Month)%>%
  mutate(Keep=TRUE) 

D2_newdata<-newdata%>%
  st_drop_geometry()%>%
  as_tibble()%>%
  select(-Year_fac, -Year, -Year_s, -N,)%>%
  distinct()%>%
  mutate(TOT_mean30_s=2,
         TOT_mean30_s_month=2,
         PREC_mean30_s_month=2,
         WY_s=2,
         Month=as.integer(as.factor(Julian_day)))%>%
  left_join(Data_D2_effort, by=c("Month", "SubRegion"))%>%
  filter(Keep)%>% # Only retain SubRegion and Month combinations with representation in the data
  select(-Keep)

#### Create background raster of all locations 
base<-D2_newdata%>%
  mutate(Date=parse_date_time(paste("2000", Month, "15", sep="-"), "%Y-%m-%d"))%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)%>%
  st_transform(crs=st_crs(Delta2))%>%
  Rasterize_all(Location, region=Delta2)%>%
  st_as_sf(long=T, connect8=T)%>%
  filter(!is.na(Location))

### Predict for main model D2_nomonth_gam_AR -----------------------------------------

D2_nomonth_pred<-predict(D2_nomonth_gam_AR, newdata=D2_newdata, type="terms", se.fit=TRUE, discrete=T, n.threads=4)

D2_nomonth_newdata_pred<-D2_newdata%>%
  mutate(Slope=D2_nomonth_pred$fit[,"te(Julian_day_s,Latitude_s,Longitude_s):TOT_mean30_s"],
         Slope_se=D2_nomonth_pred$se.fit[,"te(Julian_day_s,Latitude_s,Longitude_s):TOT_mean30_s"])%>%
  mutate(across(c(Slope, Slope_se), ~(.x/TOT_mean30_s)/sd(Data_D2$TOT_mean30)))%>%
  mutate(Slope_se=abs(Slope_se))%>%
  mutate(Slope_l99.9=Slope-Slope_se*qnorm(0.9995),
         Slope_u99.9=Slope+Slope_se*qnorm(0.9995),
         Date=as.Date(Julian_day, origin=as.Date(paste("2000", "01", "01", sep="-"))),
         Month=month(Date, label=T))%>%
  mutate(Sig=if_else(Slope_u99.9>0 & Slope_l99.9<0, "ns", "*"))%>%
  filter(Sig=="*")%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=F)%>%
  st_transform(crs=st_crs(Delta2))

newdata_D2_nomonth_pred_rast<-Rasterize_all(D2_nomonth_newdata_pred, Slope, region=Delta2)

p_D2_nomonth_gam<-predict_plot(data=newdata_D2_nomonth_pred_rast*100000, 
                       base=base, 
                       scale_fill_continuous_diverging, 
                       guide=guide_colorbar(barheight=15), 
                       na.value=NA, 
                       palette="Blue-Red 3", 
                       breaks=seq(-14, 6, by=2),
                       name="Temperature change\nper change in\ntotal inflow\n(°C/100,000cfs)")

ggsave(p_D2_nomonth_gam, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Delta Inflow temperature/Figures/model 2 inflow nomonth temp.png",
       device="png", width=7, height=5, units="in")

### Model validation ----------------

p_D2_nomonth_check<-model_validation(D2_nomonth_gam_AR, Data_D2$Temperature, type="Inflow")

ggsave(p_D2_nomonth_check, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Delta Inflow temperature/Figures/model 2 inflow nomonth model validation.png",
       device="png", width=10, height=7, units="in")

### Slope summary -----------------------------------------------------------

Slope_summary_D2_nomonth<-D2_newdata%>%
  mutate(Slope=D2_nomonth_pred$fit[,"te(Julian_day_s,Latitude_s,Longitude_s):TOT_mean30_s"],
         Slope_se=D2_nomonth_pred$se.fit[,"te(Julian_day_s,Latitude_s,Longitude_s):TOT_mean30_s"])%>%
  mutate(across(c(Slope, Slope_se), ~(.x/TOT_mean30_s)/sd(Data_D2$TOT_mean30)))%>%
  mutate(Slope_se=abs(Slope_se))%>%
  mutate(Slope_l99.9=Slope-Slope_se*qnorm(0.9995),
         Slope_u99.9=Slope+Slope_se*qnorm(0.9995),
         Date=as.Date(Julian_day, origin=as.Date(paste("2000", "01", "01", sep="-"))),
         Sig=if_else(Slope_u99.9>0 & Slope_l99.9<0, "ns", "*"))%>%
  arrange(Month, SubRegion, Slope)%>%
  group_by(Month, SubRegion)%>%
  summarise(Slope_mean=mean(Slope), .groups="drop")

## Fit month-adjusted models ------------------------------------------------

#### Try centered and standardized separately for each month
D2_gam_NOAR <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                      te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=TOT_mean30_s_month) + 
                      s(Time_num_s, k=5), family=scat, data = Data_D2, method="fREML", discrete=T, nthreads=3)
D2r <- start_value_rho(D2_gam_NOAR, plot=TRUE)

D2_gam_AR <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                    te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=TOT_mean30_s_month) + 
                    s(Time_num_s, k=5), family=scat, rho=D2r, AR.start=Start, data = Data_D2, method="fREML", discrete=T, nthreads=3)

D2_gam_AR_higherk <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(50, 13)) + 
                    te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=TOT_mean30_s_month) + 
                    s(Time_num_s, k=5), family=scat, rho=D2r, AR.start=Start, data = Data_D2, method="fREML", discrete=T, nthreads=3)
##### Model predictions are almost identical, so picking D2_gam_AR as best model. 


## Autocorrelation ---------------------------------------------------------

p_D2_variogram<-ST_variogram(D2_gam_AR, Data_D2, 4)

ggsave(p_D2_variogram, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Delta Inflow temperature/Figures/model 1 inflow model variogram.png",
       device="png", width=8, height=5, units="in")


### Predict for main model D2_gam_AR -----------------------------------------

D2_pred<-predict(D2_gam_AR, newdata=D2_newdata, type="terms", se.fit=TRUE, discrete=T, n.threads=4)

D2_newdata_pred<-D2_newdata%>%
  mutate(Slope=D2_pred$fit[,"te(Julian_day_s,Latitude_s,Longitude_s):TOT_mean30_s_month"],
         Slope_se=D2_pred$se.fit[,"te(Julian_day_s,Latitude_s,Longitude_s):TOT_mean30_s_month"])%>%
  mutate(across(c(Slope, Slope_se), ~(.x/TOT_mean30_s_month)))%>%
  mutate(Slope_se=abs(Slope_se))%>%
  mutate(Slope_l99.9=Slope-Slope_se*qnorm(0.9995),
         Slope_u99.9=Slope+Slope_se*qnorm(0.9995),
         Date=as.Date(Julian_day, origin=as.Date(paste("2000", "01", "01", sep="-"))),
         Month=month(Date, label=T))%>%
  mutate(Sig=if_else(Slope_u99.9>0 & Slope_l99.9<0, "ns", "*"))%>%
  filter(Sig=="*")%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=F)%>%
  st_transform(crs=st_crs(Delta2))

newdata_D2_pred_rast<-Rasterize_all(D2_newdata_pred, Slope, region=Delta2)

p_D2_gam<-predict_plot(data=newdata_D2_pred_rast, 
                       base=base, 
                       scale_fill_continuous_diverging, 
                       guide=guide_colorbar(barheight=15), 
                       na.value=NA, 
                       palette="Blue-Red 3", 
                       breaks=seq(-12, 6, by=2)/10,
                       name="Temperature change\nper change in\ntotal inflow\n(°C/monthly sd[cfs])")

ggsave(p_D2_gam, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Delta Inflow temperature/Figures/Figure 4 model 1 inflow temp.png",
       device="png", width=7, height=5, units="in")

### Model validation ----------------

p_D2_check<-model_validation(D2_gam_AR, Data_D2$Temperature, type="Inflow")

ggsave(p_D2_check, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Delta Inflow temperature/Figures/Figure 3 model 1 inflow model validation.png",
       device="png", width=10, height=7, units="in")

#### Create dataframe of slopes for each month and region --------

Slope_summary_D2<-D2_newdata%>%
  mutate(Slope=D2_pred$fit[,"te(Julian_day_s,Latitude_s,Longitude_s):TOT_mean30_s_month"],
         Slope_se=D2_pred$se.fit[,"te(Julian_day_s,Latitude_s,Longitude_s):TOT_mean30_s_month"])%>%
  mutate(across(c(Slope, Slope_se), ~(.x/TOT_mean30_s_month)))%>%
  mutate(Slope_se=abs(Slope_se))%>%
  mutate(Slope_l99.9=Slope-Slope_se*qnorm(0.9995),
         Slope_u99.9=Slope+Slope_se*qnorm(0.9995),
         Date=as.Date(Julian_day, origin=as.Date(paste("2000", "01", "01", sep="-"))),
         Sig=if_else(Slope_u99.9>0 & Slope_l99.9<0, "ns", "*"))%>%
  arrange(Month, SubRegion, Slope)%>%
  group_by(Month, SubRegion)%>%
  summarise(Slope_mean=mean(Slope), Slope_u99.9_max=max(Slope_u99.9), Slope_u99.9_min=min(Slope_u99.9), Slope_l99.9_max=min(Slope_l99.9), Slope_l99.9_min=max(Slope_l99.9), .groups="drop")

P_slope_sum<-ggplot(Slope_summary_D2)+
  geom_vline(xintercept=0)+
  geom_linerange(aes(y=reorder(month(Month, label=T), desc(month(Month, label=T))), xmax=Slope_l99.9_max, xmin=Slope_u99.9_max), size=1, color="#d7191c")+
  geom_linerange(aes(y=reorder(month(Month, label=T), desc(month(Month, label=T))), xmax=Slope_l99.9_min, xmin=Slope_u99.9_min), size=2, color="#fdae61")+
  geom_point(aes(y=reorder(month(Month, label=T), desc(month(Month, label=T))), x=Slope_mean), size=1, color="#2c7bb6")+
  facet_geo(~SubRegion, grid=mygrid, labeller=label_wrap_gen(width=18))+
  scale_y_discrete(breaks=c("Jan", "Mar", "May", "Jul", "Sep", "Nov"))+
  ylab("Month")+
  xlab("Temperature change per change in total inflow (°C/monthly sd[cfs])")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1), text=element_text(size=16), panel.background = element_rect(color="black"))

ggsave(P_slope_sum, file="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Delta Inflow temperature/Figures/model 1 inflow all slopes.png",
       device="png", width=15, height=18, units="in")

### Predict for higher k model D2_gam_AR_higherk -----------------------------------------

D2_pred_higherk<-predict(D2_gam_AR_higherk, newdata=D2_newdata, type="terms", se.fit=TRUE, discrete=T, n.threads=4)

D2_newdata_pred_higherk<-D2_newdata%>%
  mutate(Slope=D2_pred_higherk$fit[,"te(Julian_day_s,Latitude_s,Longitude_s):TOT_mean30_s_month"],
         Slope_se=D2_pred_higherk$se.fit[,"te(Julian_day_s,Latitude_s,Longitude_s):TOT_mean30_s_month"])%>%
  mutate(across(c(Slope, Slope_se), ~(.x/TOT_mean30_s_month)))%>%
  mutate(Slope_se=abs(Slope_se))%>%
  mutate(Slope_l99.9=Slope-Slope_se*qnorm(0.9995),
         Slope_u99.9=Slope+Slope_se*qnorm(0.9995),
         Date=as.Date(Julian_day, origin=as.Date(paste("2000", "01", "01", sep="-"))),
         Month=month(Date, label=T))%>%
  mutate(Sig=if_else(Slope_u99.9>0 & Slope_l99.9<0, "ns", "*"))%>%
  filter(Sig=="*")%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=F)%>%
  st_transform(crs=st_crs(Delta2))

newdata_D2_pred_higherk_rast<-Rasterize_all(D2_newdata_pred_higherk, Slope, region=Delta2)

p_D2_gam_higherk<-predict_plot(data=newdata_D2_pred_higherk_rast, 
                       base=base, 
                       scale_fill_continuous_diverging, 
                       guide=guide_colorbar(barheight=15), 
                       na.value=NA, 
                       palette="Blue-Red 3", 
                       breaks=seq(-12, 6, by=2)/10,
                       name="Temperature change\nper change in\ntotal inflow\n(°C/monthly sd[cfs])")

ggsave(p_D2_gam_higherk, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Delta Inflow temperature/Figures/model 1B inflow temp higher k.png",
       device="png", width=7, height=5, units="in")


# Model climate change signal and inflow effects --------------------------

D2_CC_gam_NOAR <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                      te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=TOT_mean30_s_month) + 
                        te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=WY_s) + 
                      s(Time_num_s, k=5), family=scat, data = Data_D2, method="fREML", discrete=T, nthreads=5)
D2CCr <- start_value_rho(D2_CC_gam_NOAR, plot=TRUE)

D2_CC_gam_AR <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                    te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=TOT_mean30_s_month) + 
                      te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=WY_s) + 
                    s(Time_num_s, k=5), family=scat, rho=D2CCr, AR.start=Start, data = Data_D2, method="fREML", discrete=T, nthreads=5)

## Autocorrelation ---------------------------------------------------------

p_D2_CC_variogram<-ST_variogram(D2_CC_gam_AR, Data_D2, 5)

ggsave(p_D2_CC_variogram, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Delta Inflow temperature/Figures/model 3 inflow_CC model variogram.png",
       device="png", width=8, height=5, units="in")

## Predict ------------------------------------------------------------------

D2_CC_pred<-predict(D2_CC_gam_AR, newdata=D2_newdata, type="terms", se.fit=TRUE, discrete=T, n.threads=4)

D2_CC_newdata_pred<-D2_newdata%>%
  mutate(Slope_D=D2_CC_pred$fit[,"te(Julian_day_s,Latitude_s,Longitude_s):TOT_mean30_s_month"],
         Slope_D_se=D2_CC_pred$se.fit[,"te(Julian_day_s,Latitude_s,Longitude_s):TOT_mean30_s_month"],
         Slope_CC=D2_CC_pred$fit[,"te(Julian_day_s,Latitude_s,Longitude_s):WY_s"],
         Slope_CC_se=D2_CC_pred$se.fit[,"te(Julian_day_s,Latitude_s,Longitude_s):WY_s"])%>%
  mutate(across(c(Slope_D, Slope_D_se), ~(.x/TOT_mean30_s_month)),
         across(c(Slope_CC, Slope_CC_se), ~(.x/WY_s)/sd(unique(Data_D2$WY))))%>%
  mutate(across(c(Slope_D_se, Slope_CC_se), ~abs(.x)))%>%
  mutate(Slope_D_l99=Slope_D-Slope_D_se*qnorm(0.9995),
         Slope_D_u99=Slope_D+Slope_D_se*qnorm(0.9995),
         Slope_CC_l99=Slope_CC-Slope_CC_se*qnorm(0.9995),
         Slope_CC_u99=Slope_CC+Slope_CC_se*qnorm(0.9995),
         Date=as.Date(Julian_day, origin=as.Date(paste("2000", "01", "01", sep="-"))),
         Month=month(Date, label=T))%>%
  mutate(Sig_D=if_else(Slope_D_u99>0 & Slope_D_l99<0, "ns", "*"),
         Sig_CC=if_else(Slope_CC_u99>0 & Slope_CC_l99<0, "ns", "*"))

newdata_D2_CC_pred_rast_D<-D2_CC_newdata_pred%>%
  filter(Sig_D=="*")%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=F)%>%
  st_transform(crs=st_crs(Delta2))%>%
  Rasterize_all(Slope_D, region=Delta2)

newdata_D2_CC_pred_rast_CC<-D2_CC_newdata_pred%>%
  filter(Sig_CC=="*")%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=F)%>%
  st_transform(crs=st_crs(Delta2))%>%
  Rasterize_all(Slope_CC, region=Delta2)

p_D2_CC_gam_D<-predict_plot(data=newdata_D2_CC_pred_rast_D, 
                       base=base, 
                       scale_fill_continuous_diverging, 
                       guide=guide_colorbar(barheight=15), 
                       na.value=NA, 
                       palette="Blue-Red 3", 
                       breaks=seq(-12, 6, by=2)/10,
                       name="Temperature change\nper change in\ntotal inflow\n(°C/monthly sd[cfs])")

ggsave(p_D2_CC_gam_D, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Delta Inflow temperature/Figures/model 3 inflow temp inflow_CC model.png",
       device="png", width=7, height=5, units="in")

p_D2_CC_gam_CC<-predict_plot(data=newdata_D2_CC_pred_rast_CC, 
                            base=base, 
                            scale_fill_viridis_c, 
                            guide=guide_colorbar(barheight=15), 
                            na.value=NA, 
                            breaks=(-6:7)/100, 
                            name="Temperature change\nper year (°C)")

ggsave(p_D2_CC_gam_CC, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Delta Inflow temperature/Figures/model 3 climate change signal inflow_CC model.png",
       device="png", width=7, height=5, units="in")


## Model validation -----------------

p_D2_CC_check<-model_validation(D2_CC_gam_AR, Data_D2$Temperature, type="Inflow")

ggsave(p_D2_CC_check, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Delta Inflow temperature/Figures/model 3 inflow_CC model validation.png",
       device="png", width=10, height=7, units="in")

## Slope summary -----------------------------------------------------------

Slope_summary_D2_CC<-D2_newdata%>%
  mutate(Slope=D2_CC_pred$fit[,"te(Julian_day_s,Latitude_s,Longitude_s):TOT_mean30_s_month"],
         Slope_se=D2_CC_pred$se.fit[,"te(Julian_day_s,Latitude_s,Longitude_s):TOT_mean30_s_month"])%>%
  mutate(across(c(Slope, Slope_se), ~(.x/TOT_mean30_s_month)))%>%
  mutate(Slope_se=abs(Slope_se))%>%
  mutate(Slope_l99.9=Slope-Slope_se*qnorm(0.9995),
         Slope_u99.9=Slope+Slope_se*qnorm(0.9995),
         Date=as.Date(Julian_day, origin=as.Date(paste("2000", "01", "01", sep="-"))),
         Sig=if_else(Slope_u99.9>0 & Slope_l99.9<0, "ns", "*"))%>%
  arrange(Month, SubRegion, Slope)%>%
  group_by(Month, SubRegion)%>%
  summarise(Slope_mean=mean(Slope), .groups="drop")

# Use salinity instead of geography ---------------------------------------


## Prepare data again from scratch -----------------------------------------


### Load Delta Shapefile from Brian ----------------------------------
Delta_D3<-st_read("Delta Subregions")%>%
  filter(!SubRegion%in%c("South Bay", "San Francisco Bay", "San Pablo Bay", "Upper Yolo Bypass", 
                         "Upper Napa River", "Lower Napa River", "Carquinez Strait"))%>% # Remove regions outside our domain of interest
  dplyr::select(SubRegion)

### Load data ----------------------------------
Data_D3_pre <- wq()%>%
  filter(!is.na(Temperature) & !is.na(Salinity) & !is.na(Datetime) & !is.na(Latitude) & !is.na(Longitude) & !is.na(Date))%>% #Remove any rows with NAs in our key variables
  filter(Temperature !=0)%>% #Remove 0 temps
  filter(hour(Datetime)>=5 & hour(Datetime)<=20)%>% # Only keep data between 5AM and 8PM
  mutate(Datetime = with_tz(Datetime, tz="Etc/GMT+7"), #Convert to a timezone without daylight savings time
         Date = with_tz(Date, tz="Etc/GMT+7"),
         Time=as_hms(Datetime), # Create variable for time-of-day, not date. 
         Noon_diff=abs(hms(hours=12)-Time))%>% # Calculate difference from noon for each data point for later filtering
  group_by(Station, Source, Date)%>%
  filter(Noon_diff==min(Noon_diff))%>% # Select only 1 data point per station and date, choose data closest to noon
  filter(Time==min(Time))%>% # When points are equidistant from noon, select earlier point
  ungroup()%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)%>% # Convert to sf object
  st_transform(crs=st_crs(Delta_D3))%>% # Change to crs of Delta
  st_join(Delta_D3, join=st_intersects)%>% # Add subregions
  filter(!is.na(SubRegion))%>% # Remove any data outside our subregions of interest
  mutate(Julian_day = yday(Date), # Create julian day variable
         Month_fac=factor(Month), # Create month factor variable
         Source_fac=factor(Source),
         Year_fac=factor(Year))%>% 
  mutate(Date_num = as.numeric(Date))%>%  # Create numeric version of date for models
  mutate(Time_num=as.numeric(Time)) # Create numeric version of time for models (=seconds since midnight)


### Pull station locations for major monitoring programs ----------------------------------
#### This will be used to set a boundary for this analysis focused on well-sampled regions.
WQ_stations_D3<-Data_D3_pre%>%
  st_drop_geometry()%>%
  filter(Source%in%c("FMWT", "STN", "SKT", "20mm", "EMP", "Suisun"))%>%
  group_by(StationID, Source, Latitude, Longitude)%>%
  summarise(N=n(), .groups="drop")%>% # Calculate how many times each station was sampled
  filter(N>50 & !StationID%in%c("20mm 918", "STN 918"))%>% # Only keep stations sampled >50 times when deciding which regions to retain. 
  # "20mm 918", "STN 918" are far south of the rest of the well-sampled sites and are not sampled year round, so we're removing them to exclude that far southern region
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)%>% # Convert to sf object
  st_transform(crs=st_crs(Delta_D3))%>%
  st_join(Delta_D3) # Add subregions

### Remove any subregions that do not contain at least one of these >50 samples stations from the major monitoring programs ----------------------------------
Delta_D3 <- Delta_D3%>%
  filter(SubRegion%in%unique(WQ_stations_D3$SubRegion) | SubRegion=="Georgiana Slough") # Retain Georgiana Slough because it's surrounded by well-sampled regions

### Now prepare the final dataset  ----------------------------------
#### Using code from the final steps of creating Data and the creation of Data_CC4 and Data_D2
#### Filter data to only include this final set of subregions, and any stations outside the convex hull formed by the >50 samples stations from the major monitoring programs
Data_D3<-Data_D3_pre%>%
  filter(SubRegion%in%unique(Delta_D3$SubRegion))%>%
  st_join(WQ_stations_D3%>%
            st_union()%>%
            st_convex_hull()%>% # Draws a hexagram or pentagram or similar around the outer-most points
            st_as_sf()%>%
            mutate(IN=TRUE),
          join=st_intersects)%>%
  filter(IN)%>%
  dplyr::select(-IN)%>%
  filter(Source%in%c("EMP", "STN", "FMWT", "DJFMP", "SKT", "20mm", "Suisun", "Baystudy", "USGS") & !str_detect(Station, "EZ") & 
           !(Source=="SKT" & Station=="799" & Latitude>38.2) & !(Source=="SKT" & Station=="999"))%>%
  group_by(SubRegion)%>%
  mutate(Sal_sd=sd(Salinity), Sal_ext=(Salinity-mean(Salinity))/Sal_sd)%>%
  ungroup()%>%
  filter(Sal_ext<10)%>%
  mutate(Station=paste(Source, Station),
         Noon_diff=abs(hms(hours=12)-as_hms(Datetime)),
         mday_15_diff=abs(mday(Date)-15))%>% # Find how far each date is from the 15th of the month
  group_by(Station, Month, Year)%>%
  filter(mday_15_diff==min(mday_15_diff))%>%
  filter(Date==min(Date))%>% # Deal with 2 dates equidistant from the 15th of the month
  filter(Noon_diff==min(Noon_diff))%>%
  ungroup()%>%
  lazy_dt()%>%
  group_by(Date, Date_num, Month, SubRegion, Julian_day, Year, Year_fac, Station, Source, Latitude, Longitude)%>%
  summarise(Temperature=mean(Temperature), Salinity=mean(Salinity), Time_num=mean(Time_num))%>%
  as_tibble()%>%
  ungroup()%>%
  mutate(ID=paste(Station, Date_num))%>%
  filter(!(ID%in%ID[which(duplicated(ID))]))%>%
  mutate(WY=Year)%>%
  mutate(WY=if_else(Month%in%10:12, WY+1, WY),
         Date=as_date(Date))%>% # Remove times to make sure merging works correctly
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
  mutate(across(c(TOT_mean30, Salinity, Time_num, Julian_day, Latitude, Longitude), list(s=~(.x-mean(.x))/sd(.x))))%>%
  mutate(WY_s=(WY-mean(unique(WY)))/sd(unique(WY)),
         TOT_mean30_s_month=(TOT_mean30-TOT_mean30_mean)/TOT_mean30_sd)

# Tried log(x+1) transforming salinity prior to use in these analyses, but results were similar to just using salinity. Using a smaller constant
# like log(x+0.01) was producing weird results

D3_gam_NOAR <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                      te(Salinity_s, Julian_day_s, bs=c("tp", "cc"), k=c(15, 13), by=TOT_mean30_s_month) + 
                      s(Time_num_s, k=5), family=scat, data = Data_D3, method="fREML", discrete=T, nthreads=3)

D3r <- start_value_rho(D3_gam_NOAR, plot=TRUE)

D3_gam_AR <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                   te(Salinity_s, Julian_day_s, bs=c("tp", "cc"), k=c(15, 13), by=TOT_mean30_s_month) + 
                    s(Time_num_s, k=5), family=scat, rho=D3r, AR.start=Start, data = Data_D3, method="fREML", discrete=T, nthreads=3)

## Autocorrelation ---------------------------------------------------------

p_D3_variogram<-ST_variogram(D3_gam_AR, Data_D3, 5)

ggsave(p_D3_variogram, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Delta Inflow temperature/Figures/model 4 inflow salinity model variogram.png",
       device="png", width=8, height=5, units="in")

## Predict -----------------------------------------------------------------

D3_effort<-Data_D3%>%
  mutate(Month=month(Month, label=T))%>%
  group_by(Month)%>%
  summarise(Sal_min=quantile(Salinity, probs=0.05), Sal_max=quantile(Salinity, probs=0.95), .groups="drop")

D3_newdata<-expand_grid(Salinity=seq(quantile(Data_D3$Salinity, probs=0.05), quantile(Data_D3$Salinity, probs=0.95), length.out=100),
                        Julian_day=yday(ymd(paste("2001", 1:12, "15", sep="-"))))%>%
  mutate(Julian_day_s=(Julian_day-mean(Data_D3$Julian_day))/sd(Data_D3$Julian_day),
         Salinity_s=(Salinity-mean(Data_D3$Salinity))/sd(Data_D3$Salinity),
         Latitude_s=0,
         Longitude_s=0,
         TOT_mean30_s_month=2,
         Time_num=12*60*60,
         Time_num_s=Time_num*sd(Data_D3$Time_num)+mean(Data_D3$Time_num),
         Date=as.Date(Julian_day, origin=as.Date(paste("2000", "01", "01", sep="-"))),
         Month=month(Date, label=T))%>%
  left_join(D3_effort, by="Month")%>%
  filter(Salinity >= Sal_min & Salinity <= Sal_max)

D3_pred<-predict(D3_gam_AR, newdata=D3_newdata, type="terms", se.fit=TRUE, discrete=T, n.threads=2)

D3_newdata_pred<-D3_newdata%>%
  mutate(Slope=D3_pred$fit[,"te(Julian_day_s,Salinity_s):TOT_mean30_s_month"],
         Slope_se=D3_pred$se.fit[,"te(Julian_day_s,Salinity_s):TOT_mean30_s_month"])%>%
  mutate(across(c(Slope, Slope_se), ~(.x/TOT_mean30_s_month)))%>%
  mutate(Slope_se=abs(Slope_se))%>%
  mutate(Slope_l99.9=Slope-Slope_se*qnorm(0.9995),
         Slope_u99.9=Slope+Slope_se*qnorm(0.9995),
         Slope_l99=Slope-Slope_se*qnorm(0.995),
         Slope_u99=Slope+Slope_se*qnorm(0.995),
         Slope_l95=Slope-Slope_se*qnorm(0.975),
         Slope_u95=Slope+Slope_se*qnorm(0.975))%>%
  mutate(Sig=if_else(Slope_u99.9>0 & Slope_l99.9<0, "ns", "*"))

# Calculate monthly mean salinity for each subregion and month
sal_sum_month<-Data_D3%>%
  group_by(SubRegion, Month)%>%
  summarise(Salinity_mean=mean(Salinity), Salinity_sd=sd(Salinity), .groups="drop")

p_D3<-ggplot(D3_newdata_pred, aes(x=Salinity, y=Slope, ymin=Slope_l99.9, ymax=Slope_u99.9))+
  geom_ribbon(alpha=0.3)+
  geom_ribbon(alpha=0.3, aes(ymin=Slope_l99, ymax=Slope_u99))+
  geom_ribbon(alpha=0.3, aes(ymin=Slope_l95, ymax=Slope_u95))+
  geom_line(alpha=0.3)+
  geom_hline(yintercept=0)+
  geom_vline(data=sal_sum_month%>%filter(SubRegion=="Confluence")%>%mutate(Month=month(Month, label=T)), 
             aes(xintercept=Salinity_mean), color="red", linetype=2)+
  ylab("Temperature change per change in total inflow (°C/monthly sd[cfs])")+
  facet_wrap(~Month)+
  theme_bw()+
  theme(strip.background=element_blank())


ggsave(p_D3, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Delta Inflow temperature/Figures/Figure 5 model 4 inflow temp by salinity.png",
       device="png", width=7, height=5, units="in")

### Model validation

p_D3_check<-model_validation(D3_gam_AR, Data_D3$Temperature, type="Inflow")

ggsave(p_D3_check, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Delta Inflow temperature/Figures/model 4 inflow salinity model validation.png",
       device="png", width=10, height=7, units="in")

### Now plot the salinity field  ----------------------------------

#### Plot overall salinity --------------------------------------------------

sal_sum<-Data_D3%>%
  group_by(SubRegion)%>%
  summarise(Salinity_mean=mean(Salinity), Salinity_sd=sd(Salinity), .groups="drop")

Delta_D3_sal<-Delta_D3%>%
  left_join(sal_sum, by="SubRegion")

water_mean<-st_make_grid(Delta_D3_sal, n=200)%>%
  st_as_sf()%>%
  st_join(deltamapr::WW_Delta%>%
            st_transform(crs=st_crs(Delta_D3_sal))%>%
            mutate(In=TRUE)%>%
            select(In))%>%
  filter(In)%>%
  st_transform(crs=st_crs(Delta2))%>%
  st_join(Delta_D3_sal)%>%
  select(Salinity_mean)%>%
  st_rasterize(template=st_as_stars(st_bbox(Delta_D3_sal), dx=diff(st_bbox(Delta_D3_sal)[c(1, 3)])/200, dy=diff(st_bbox(Delta_D3_sal)[c(2, 4)])/200, values = NA_real_))%>%
  st_warp(crs=4326)

water_sd<-st_make_grid(Delta_D3_sal, n=200)%>%
  st_as_sf()%>%
  st_join(deltamapr::WW_Delta%>%
            st_transform(crs=st_crs(Delta_D3_sal))%>%
            mutate(In=TRUE)%>%
            select(In))%>%
  filter(In)%>%
  st_transform(crs=st_crs(Delta2))%>%
  st_join(Delta_D3_sal)%>%
  select(Salinity_sd)%>%
  st_rasterize(template=st_as_stars(st_bbox(Delta_D3_sal), dx=diff(st_bbox(Delta_D3_sal)[c(1, 3)])/200, dy=diff(st_bbox(Delta_D3_sal)[c(2, 4)])/200, values = NA_real_))%>%
  st_warp(crs=4326)


p_sal_mean<-ggplot()+
  geom_stars(data=water_mean)+
  scale_fill_viridis_c(trans="log", labels=function(x) round(x, 2), name="Mean salinity", na.value = NA, aesthetics=c("fill", "color"), 
                       guide=guide_colorbar(direction = "horizontal", title.position = "top", title.hjust = 0.5, barwidth = 10))+
  coord_sf()+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw()+
  theme(legend.position = c(0.25, 0.8), legend.background=element_rect(color="black"))

p_sal_sd<-ggplot()+
  geom_stars(data=water_sd)+
  scale_fill_viridis_c(trans="log", labels=function(x) round(x, 2), name="SD salinity", na.value = NA, aesthetics=c("fill", "color"),
                       guide=guide_colorbar(direction = "horizontal", title.position = "top", title.hjust = 0.5, barwidth = 10))+
  coord_sf()+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw()+
  theme(legend.position = c(0.25, 0.8), legend.background=element_rect(color="black"))

p_sal<-p_sal_mean+p_sal_sd+plot_annotation(tag_levels="A")

ggsave(p_sal, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Delta Inflow temperature/Figures/Overall salinity map.png",
       device="png", width=10, height=5, units="in")

#### Plot monthly salinity --------------------------------------------------


water_mean2<-base%>%
  mutate(Month=month(Date))%>%
  st_transform(crs=st_crs(Delta_D3))%>%
  st_join(Delta_D3)%>%
  left_join(sal_sum_month, by=c("SubRegion", "Month"))

p_sal_month_mean<-ggplot()+
  geom_sf(data=water_mean2, aes(fill=Salinity_mean, color=Salinity_mean))+
  facet_wrap(~month(Month, label=T))+
  scale_fill_viridis_c(trans="log", labels=function(x) round(x, 2), name="Mean salinity", na.value = NA, aesthetics=c("fill", "color"),
                       guide=guide_colorbar(barwidth = 15))+
  coord_sf()+
  theme_bw()+
  theme(legend.position = "top", strip.background=element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid=element_blank(), legend.margin=margin(0,0,-10,0), panel.spacing.x = unit(0, "lines"))

p_sal_month_sd<-ggplot()+
  geom_sf(data=water_mean2, aes(fill=Salinity_sd, color=Salinity_sd))+
  facet_wrap(~month(Month, label=T))+
  scale_fill_viridis_c(trans="log", labels=function(x) round(x, 2), name="SD salinity", na.value = NA, aesthetics=c("fill", "color"),
                       guide=guide_colorbar(barwidth = 15))+
  coord_sf()+
  theme_bw()+
  theme(legend.position = "top", strip.background=element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid=element_blank(), legend.margin=margin(0,0,-10,0), panel.spacing.x = unit(0, "lines"))

p_sal_month<-p_sal_month_mean+p_sal_month_sd+plot_annotation(tag_levels="A")

ggsave(p_sal_month, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Delta Inflow temperature/Figures/Figure 6 salinity month map.png",
       device="png", width=13, height=7, units="in")


# Model salinity?  --------------------------------------------------


Data_D3_sal<-Data_D3%>%
  mutate(Hour=floor(Time_num/3600),
         Minute=floor((Time_num-3600*Hour)/60),
         Datetime=parse_date_time(paste0(year(Date), "-", month(Date), "-", day(Date), " ", Hour, ":", Minute), "%Y-%m-%d %H:%M", tz="Etc/GMT+7"))%>%
  filter(Salinity!=0)%>%
  mutate(Salinity_l=log(Salinity),
         ID=1:n())

tide<-list()
for(i in unique(Data_D3_sal$SubRegion)){
  data<-filter(Data_D3_sal, SubRegion==i)
  tidal_model<-oce::tidem(t=data$Datetime, x=data$Salinity, latitude=mean(data$Latitude))
  tide[[i]]<-tibble(ID=data$ID, Tide=predict(tidal_model))
}

Data_D3_sal<-Data_D3_sal%>%
  left_join(bind_rows(tide), by="ID")%>%
  mutate(across(c(Tide, Year), list(s=~(.x-mean(.x))/sd(.x))))

D3_sal_gam_NOAR<-bam(Salinity ~ te(Latitude_s, Longitude_s, Julian_day_s, Tide_s, d=c(2, 1, 1), bs=c("tp", "cc", "cc"), k=c(25, 13, 5)),
                     family=scat, data = Data_D3_sal, method="fREML", discrete=T, nthreads=4)
D3_salr <- start_value_rho(D3_sal_gam_NOAR, plot=TRUE)

D3_sal_gam_NOAR2<-bam(Salinity ~ te(Latitude_s, Longitude_s, Tide_s, d=c(2, 1), bs=c("tp", "cc"), k=c(25, 5)),
                     family=scat, data = Data_D3_sal, method="fREML", discrete=T, nthreads=4)

D3_sal_gam_NOAR3<-bam(Salinity_l ~ te(Latitude_s, Longitude_s, Julian_day_s, Tide_s, d=c(2, 1, 1), bs=c("tp", "cc", "cc"), k=c(25, 13, 5)),
                     data = Data_D3_sal, method="fREML", discrete=T, nthreads=4)

D3_sal_gam_NOAR4<-bam(Salinity_l ~ te(Latitude_s, Longitude_s, Julian_day_s, Tide_s, d=c(2, 1, 1), bs=c("tp", "cc", "cc"), k=c(25, 13, 5)),
                      family=scat, data = Data_D3_sal, method="fREML", discrete=T, nthreads=4)

D3_sal_gam_NOAR5<-bam(Salinity ~ te(Latitude_s, Longitude_s, Julian_day_s, Tide_s, d=c(2, 1, 1), bs=c("tp", "cc", "cc"), k=c(25, 13, 5)),
                      family=tw, data = Data_D3_sal, method="fREML", discrete=T, nthreads=4)

D3_sal_gam_NOAR6<-bam(Salinity_l ~ te(Latitude_s, Longitude_s, Julian_day_s, Tide_s, d=c(2, 1, 1), bs=c("tp", "cc", "cc"), k=c(25, 13, 5)),
                      family=scat, data = Data_D3_sal, method="fREML", discrete=T, nthreads=4)

D3_sal_gam_NOAR6B<-bam(Salinity_l ~ te(Latitude_s, Longitude_s, Tide_s, d=c(2, 1), bs=c("tp", "cc"), k=c(25, 5)),
                      family=scat, data = Data_D3_sal, method="fREML", discrete=T, nthreads=4)

D3_sal_gam_NOAR7<-bam(Salinity_l ~ te(Latitude_s, Longitude_s, Tide_s, d=c(2, 1), bs=c("tp", "cc"), k=c(25, 15)),
                      family=scat, data = Data_D3_sal, method="fREML", discrete=T, nthreads=4)

D3_sal_gam_NOAR8<-bam(Salinity_l ~ te(Latitude_s, Longitude_s, Tide_s, d=c(2, 1), bs=c("tp", "tp"), k=c(25, 15)),
                      family=scat, data = Data_D3_sal, method="fREML", discrete=T, nthreads=4)

D3_sal_gam_NOAR9<-bam(Salinity_l ~ te(Latitude_s, Longitude_s, Tide_s, d=c(2, 1), bs=c("tp", "tp"), k=c(50, 15)),
                      family=scat, data = Data_D3_sal, method="fREML", discrete=T, nthreads=4)

D3_sal_gam_NOAR10<-bam(Salinity_l ~ te(Latitude_s, Longitude_s, Tide_s, Year_s, d=c(2, 1, 1), bs=c("tp", "tp", "tp"), k=c(25, 5, 10)),
                      family=scat, data = Data_D3_sal, method="fREML", discrete=T, nthreads=4)

D3_gam_AR <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                   te(Salinity_l_s, Julian_day_s, bs=c("tp", "cc"), k=c(15, 13), by=TOT_mean30_s_month) + 
                   s(Time_num_s, k=5), family=scat, rho=D3_salr, AR.start=Start, data = Data_D3, method="fREML", discrete=T, nthreads=2)



# Plot monthly sd of inflow for each dataset ------------------------------

dayflow_sum<-dayflow%>%
  mutate(Month=month(Date, label=T))%>%
  select(Month, TOT_mean30_mean, TOT_mean30_sd, PREC_mean30_mean, PREC_mean30_sd)%>%
  distinct()

dayflow_sum%>%
  mutate(across(c(TOT_mean30_mean, TOT_mean30_sd, PREC_mean30_mean, PREC_mean30_sd), ~format(round(.x), big.mark = ",")))%>%
  rename(`Mean inflow`=TOT_mean30_mean, `SD inflow`=TOT_mean30_sd, `Mean precip`=PREC_mean30_mean, `SD precip`=PREC_mean30_sd)%>%
  arrange(Month)%>%
  write_csv("C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Delta Inflow temperature/Figures/Table 1.csv")

dayflow_sum2<-dayflow%>%
  mutate(Month=month(Date, label=T),
         Year=year(Date))%>%
  group_by(Month, Year)%>%
  summarise(across(c(TOT_mean30, PREC_mean30), list(mean=mean, sd=sd)), .groups="drop")%>%
  filter(Year>=min(Data_D2$Year))

p_inflow<-ggplot(data=dayflow_sum2, aes(x=Year, y=TOT_mean30_mean, ymin=TOT_mean30_mean-TOT_mean30_sd, ymax=TOT_mean30_mean+TOT_mean30_sd))+
  geom_rect(data=dayflow_sum, aes(xmin=1969, xmax=2019, ymin=TOT_mean30_mean, ymax=TOT_mean30_mean+TOT_mean30_sd, fill="1 SD above mean"), alpha=0.2, inherit.aes = F)+
  geom_rect(data=dayflow_sum, aes(xmin=1969, xmax=2019, ymin=TOT_mean30_mean-TOT_mean30_sd, ymax=TOT_mean30_mean, fill="1 SD below mean"), alpha=0.2, inherit.aes = F)+
  geom_line(color="gray50")+
  geom_point()+
  facet_wrap(~Month, scales = "free_y")+
  scale_fill_manual(values=c("1 SD above mean" = "firebrick3", "1 SD below mean" = "dodgerblue3"), name="")+
  ylab("30-day lagged mean inflow (cfs)")+
  theme_bw()+
  theme(strip.background=element_blank(), legend.position="top", axis.text.x=element_text(angle=45, hjust=1))

ggsave(p_inflow, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Delta Inflow temperature/Figures/inflow by month and year.png",
       device="png", width=9, height=9, units="in")

# Now try precipitation ---------------------------------------------------

D2_gam_PREC_NOAR <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                          te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=PREC_mean30_s_month) + 
                          s(Time_num_s, k=5), family=scat, data = Data_D2, method="fREML", discrete=T, nthreads=5)
D2r_PREC <- start_value_rho(D2_gam_PREC_NOAR, plot=TRUE)

D2_gam_PREC_AR <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                        te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=PREC_mean30_s_month) + 
                        s(Time_num_s, k=5), family=scat, rho=D2r_PREC, AR.start=Start, data = Data_D2, method="fREML", discrete=T, nthreads=5)

## Autocorrelation ---------------------------------------------------------

p_D2_PREC_variogram<-ST_variogram(D2_gam_PREC_AR, Data_D2, 5)

ggsave(p_D2_PREC_variogram, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Delta Inflow temperature/Figures/model 5 precip model variogram.png",
       device="png", width=8, height=5, units="in")

## Model validation ----------------------------------------------------------

p_D2_PREC_check<-model_validation(D2_gam_PREC_AR, Data_D2$Temperature, type="Inflow")

ggsave(p_D2_PREC_check, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Delta Inflow temperature/Figures/model 5 precip model validation.png",
       device="png", width=10, height=7, units="in")

## Predict ------------------------------------------------------------------

D2_PREC_pred<-predict(D2_gam_PREC_AR, newdata=D2_newdata, type="terms", se.fit=TRUE, discrete=T, n.threads=2)

D2_newdata_PREC_pred<-D2_newdata%>%
  mutate(Slope=D2_PREC_pred$fit[,"te(Julian_day_s,Latitude_s,Longitude_s):PREC_mean30_s_month"],
         Slope_se=D2_PREC_pred$se.fit[,"te(Julian_day_s,Latitude_s,Longitude_s):PREC_mean30_s_month"])%>%
  mutate(across(c(Slope, Slope_se), ~(.x/PREC_mean30_s_month)))%>%
  mutate(Slope_se=abs(Slope_se))%>%
  mutate(Slope_l99.9=Slope-Slope_se*qnorm(0.9995),
         Slope_u99.9=Slope+Slope_se*qnorm(0.9995),
         Date=as.Date(Julian_day, origin=as.Date(paste("2000", "01", "01", sep="-"))),
         Month=month(Date, label=T))%>%
  mutate(Sig=if_else(Slope_u99.9>0 & Slope_l99.9<0, "ns", "*"))%>%
  filter(Sig=="*")%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=F)%>%
  st_transform(crs=st_crs(Delta2))

newdata_D2_PREC_pred_rast<-Rasterize_all(D2_newdata_PREC_pred, Slope, region=Delta2)

p_D2_gam_PREC<-predict_plot(data=newdata_D2_PREC_pred_rast, 
                            base=base, 
                            scale_fill_continuous_diverging, 
                            guide=guide_colorbar(barheight=15), 
                            na.value=NA, 
                            palette="Blue-Red 3", 
                            breaks=seq(-12, 6, by=2)/10,
                            name="Temperature change\nper change in\nprecipitation\n(°C/monthly sd[cfs])")

ggsave(p_D2_gam_PREC, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Delta Inflow temperature/Figures/Figure 7 model 5 precip temp.png",
       device="png", width=7, height=5, units="in")

## Slope summary -----------------------------------------------------------

Slope_summary_D2_PREC<-D2_newdata%>%
  mutate(Slope=D2_PREC_pred$fit[,"te(Julian_day_s,Latitude_s,Longitude_s):PREC_mean30_s_month"],
         Slope_se=D2_PREC_pred$se.fit[,"te(Julian_day_s,Latitude_s,Longitude_s):PREC_mean30_s_month"])%>%
  mutate(across(c(Slope, Slope_se), ~(.x/PREC_mean30_s_month)))%>%
  mutate(Slope_se=abs(Slope_se))%>%
  mutate(Slope_l99.9=Slope-Slope_se*qnorm(0.9995),
         Slope_u99.9=Slope+Slope_se*qnorm(0.9995),
         Date=as.Date(Julian_day, origin=as.Date(paste("2000", "01", "01", sep="-"))),
         Sig=if_else(Slope_u99.9>0 & Slope_l99.9<0, "ns", "*"))%>%
  arrange(Month, SubRegion, Slope)%>%
  group_by(Month, SubRegion)%>%
  summarise(Slope_mean=mean(Slope), .groups="drop")


# Combine all slope summaries and plot pairs plot -------------------------

Slope_summary<-Slope_summary_D2_nomonth%>%
  rename(Model_1=Slope_mean)%>%
  left_join(Slope_summary_D2%>%
              select(Month, SubRegion, Model_2=Slope_mean),
            by=c("Month", "SubRegion"))%>%
  left_join(Slope_summary_D2_CC%>%
              rename(Model_3=Slope_mean),
            by=c("Month", "SubRegion"))%>%
  left_join(Slope_summary_D2_PREC%>%
              rename(Model_5=Slope_mean),
            by=c("Month", "SubRegion"))

pairs(formula=~Model_1 + Model_2 + Model_3 + Model_5, data=Slope_summary)

w_ggally_cor <- wrap(ggally_cor, lineheight=5, display_grid=FALSE)

slope_comparison<-ggpairs(Slope_summary, 
                          aes(color=month(Month, label=TRUE)), 
                          columns=3:6, 
                          labeller=as_labeller(function(x) str_replace(x, "_", " ")))+
  scale_color_viridis_d(end=0.9, aesthetics=c("fill", "color"))+
  theme_bw()+
  theme(strip.background = element_blank(), axis.text.x = element_text(angle=45, hjust=1))



ggsave(slope_comparison, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Delta Inflow temperature/Figures/Slope comparison.png",
       device="png", width=8, height=10, units="in")
