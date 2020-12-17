require(sp)
require(gstat)
require(spacetime)
require(tidybayes)
require(dplyr)
require(stringr)
require(dtplyr)
require(tidyr)
require(broom)
require(brms)
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

Data<-readRDS("Temperature analysis/Discrete Temp Data.Rds")

# Bayesian mixed models ---------------------------------------------------


CC_brm<-brm(Temperature~Year_s+(Year_s|Month*SubRegion),
            data=Data, family=gaussian,
            prior=prior(normal(0,5), class="Intercept")+
              prior(normal(0,5), class="b")+
              prior(cauchy(0,5), class="sigma")+
              prior(cauchy(0,5), class="sd"),
            iter=5e3, warmup=1250, cores=1, chains=1)
CC_brm<-add_criterion(CC_brm, criterion=c("waic", "loo"))

CC_brm2<-brm(Temperature~Year_s + s(Time_num_s, k=5) + (Year_s|Month*SubRegion),
             data=Data, family=gaussian,
             prior=prior(normal(0,5), class="Intercept")+
               prior(normal(0,5), class="b")+
               prior(cauchy(0,5), class="sigma")+
               prior(cauchy(0,5), class="sd"),
             iter=5e3, warmup=1250, cores=1, chains=1)
CC_brm2<-add_criterion(CC_brm2, criterion=c("waic", "loo"))
#Much better model with time incorporated

#Without EDSM data to try and reduce autocorrelation

CC_brm3<-brm(Temperature~Year_s + s(Time_num_s, k=5) + (Year_s|Month*SubRegion),
             data=filter(Data, Source!="EDSM"), family=gaussian,
             prior=prior(normal(0,5), class="Intercept")+
               prior(normal(0,5), class="b")+
               prior(cauchy(0,5), class="sigma")+
               prior(cauchy(0,5), class="sd"),
             iter=5e3, warmup=1250, cores=1, chains=1)

CC_brm_EMP<-brm(Temperature~Year_s + s(Time_num_s, k=5) + (Year_s|Month*SubRegion),
                data=filter(Data, Source=="EMP"), family=gaussian,
                prior=prior(normal(0,5), class="Intercept")+
                  prior(normal(0,5), class="b")+
                  prior(cauchy(0,5), class="sigma")+
                  prior(cauchy(0,5), class="sd"),
                iter=5e3, warmup=1250, cores=1, chains=1)

CC_brm4<-brm(Temperature~Year_s + s(Time_num_s, k=5) + (Year_s|Month*SubRegion) + (1|Source),
             data=Data, family=gaussian,
             prior=prior(normal(0,5), class="Intercept")+
               prior(normal(0,5), class="b")+
               prior(cauchy(0,5), class="sigma")+
               prior(cauchy(0,5), class="sd"),
             iter=5e3, warmup=1250, cores=1, chains=1)

CC_predictor<-function(model){
  Post_CC<-model %>%
    recover_types()%>%
    spread_draws(`r_Month:SubRegion`[MonthSubRegion,term])%>%
    ungroup() %>%
    filter(term=="Year_s")%>%
    select(-term, -.draw, -.chain)%>%
    separate(MonthSubRegion, into=c("Month", "SubRegion"), sep="_")%>%
    mutate(Month=as.integer(Month))%>%
    left_join(model %>%
                recover_types()%>%
                spread_draws(b_Year_s, r_Month[Month,term], r_SubRegion[SubRegion,term], )%>%
                ungroup()%>%
                filter(term=="Year_s")%>%
                select(-term, -.draw, -.chain), by=c("Month", "SubRegion", ".iteration"))
  return(Post_CC)
}

Post_CC<-CC_predictor(CC_brm2)

Post_CC_sum<-Post_CC%>%
  mutate(Slope=b_Year_s+r_Month+r_SubRegion+`r_Month:SubRegion`)%>%
  group_by(Month, SubRegion)%>%
  mean_qi(Slope, .width = c(0.99, 0.999))%>%
  ungroup()%>%
  mutate(SubRegion=str_replace_all(SubRegion, fixed("."), " "))

p_regmonth<-ggplot(Post_CC_sum, aes(y = Slope/sd(Data$Year), x = Month, ymin = .lower/sd(Data$Year), ymax = .upper/sd(Data$Year))) +
  geom_pointinterval()+
  facet_geo(~SubRegion, grid=mygrid, labeller=label_wrap_gen())+
  geom_hline(yintercept = 0)+
  ylab("Slope (°C / year)")+
  scale_x_continuous(breaks=seq(1,12, by=2))+
  theme_bw()+
  theme(panel.grid=element_blank())

ggsave(p_regmonth, file="C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/figures/CC regmonth.png",
       device="png", width=15, height=12, units="in")

Post_CC_sum_month<-Post_CC%>%
  select(-SubRegion, -r_SubRegion, -`r_Month:SubRegion`)%>%
  distinct()%>%
  mutate(Slope=b_Year_s+r_Month)%>%
  group_by(Month)%>%
  mean_qi(Slope, .width = c(0.99, 0.999))%>%
  ungroup()

ggplot(Post_CC_sum_month, aes(y = Slope/sd(Data$Year), x = Month, ymin = .lower/sd(Data$Year), ymax = .upper/sd(Data$Year))) +
  geom_pointinterval()+
  geom_hline(yintercept = 0)+
  scale_x_continuous(breaks=seq(1,12, by=2))+
  ylab("Slope (°C/year)")+
  theme_bw()+
  theme(panel.grid=element_blank())

p_month<-ggplot(Post_CC_sum_month, aes(y = Slope/sd(Data$Year), x = Month, ymin = .lower/sd(Data$Year), ymax = .upper/sd(Data$Year))) +
  geom_pointinterval()+
  geom_pointinterval(data=Post_CC_sum_month3, aes(x=Month+0.2), color="dodgerblue3")+
  geom_hline(yintercept = 0)+
  ylab("Slope (°C / year)")+
  scale_x_continuous(breaks=seq(1,12, by=2))+
  theme_bw()+
  theme(panel.grid=element_blank())

ggsave(p_month, file="C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/figures/CC month.png",
       device="png", width=4, height=4, units="in")

Data_CC<-Data%>%
  filter(Source!="EDSM" & !str_detect(Station, "EZ") & !(Source=="SKT" & Station=="799" & Latitude>38.2))%>%
  mutate(Station=paste(Source, Station))%>%
  lazy_dt()%>%
  group_by(Month, SubRegion, Year, Year_s, Station)%>%
  summarise(Temperature=mean(Temperature))%>%
  as_tibble()%>%
  ungroup()%>%
  mutate(Date=parse_date_time(paste(Year, Month, "01", sep="-"), "%Y-%m-%d"))

CC_brm_ar<-brm(Temperature~Year_s+(Year_s|Month*SubRegion) + ar(time = Date, gr=Station),
               data=Data_CC, family=gaussian,
               prior=prior(normal(0,5), class="Intercept")+
                 prior(normal(0,5), class="b")+
                 prior(cauchy(0,5), class="sigma")+
                 prior(cauchy(0,5), class="sd"),
               iter=5e3, warmup=1250, cores=1, chains=1)


# GAMs --------------------------------------------------------------------

#Try reducing temporal resolution to once every month.
# This time, pick days closest to the middle of the month
# And split the filtering to 2 steps to fix error that removed too much data by requiring the earliest day AND the time closest to noon
Data_CC4<-Data%>%
  filter(Source%in%c("EMP", "STN", "FMWT", "DJFMP", "SKT", "20mm", "Suisun", "Baystudy", "USGS") & !str_detect(Station, "EZ") & 
           !(Source=="SKT" & Station=="799" & Latitude>38.2) & !(Source=="SKT" & Station=="999"))%>%
  mutate(Station=paste(Source, Station),
         Noon_diff=abs(hms(hours=12)-as_hms(Datetime)),
         mday_15_diff=abs(mday(Date)-15))%>% # Find how far each date is from the 15th of the month
  group_by(Station, Month, Year)%>%
  filter(mday_15_diff==min(mday_15_diff))%>%
  filter(Date==min(Date))%>% # Deal with 2 dates equidistant from the 15th of the month
  filter(Noon_diff==min(Noon_diff))%>%
  ungroup()%>%
  lazy_dt()%>%
  group_by(Date, Date_num, Date_num_s, Month, SubRegion, Julian_day_s, Julian_day, Year, Year_s, Year_fac, Station, Source, Latitude_s, Longitude_s, Latitude, Longitude)%>%
  summarise(Temperature=mean(Temperature), Time_num=mean(Time_num), Time_num_s=mean(Time_num_s))%>%
  as_tibble()%>%
  ungroup()%>%
  mutate(ID=paste(Station, Date_num))%>%
  filter(!(ID%in%ID[which(duplicated(ID))]))%>%
  mutate(YearStation=paste(Year, Station),
         Date_num2=as.numeric(Date)/(3600*24*30), # Create a numeric date variable in units of ~ 1 month. 
         Month_fac=factor(Month)) 

#saveRDS(Data_CC4, "Temperature analysis/Data_CC4.Rds")

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
                          s(Time_num_s, k=5), family=scat, data = Data_CC4.3, method="fREML", discrete=T, nthreads=2)
r6 <- start_value_rho(CC_gam8d7b_NOAR5, plot=TRUE)

#########Best Model####################

CC_gam8d7b_AR7 <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                        te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=Year_s) + 
                        s(Time_num_s, k=5), family=scat, rho=r6, AR.start=Start, data = Data_CC4.3, method="fREML", discrete=T, nthreads=2)
#AIC: 200051.6
#BIC: 202901.6

#########Best Model####################


resid_norm_CC_gam8d7b_AR7<-resid_gam(CC_gam8d7b_AR7, incl_na=TRUE)
sp <- SpatialPoints(coords=data.frame(Longitude=Data_CC4.3$Longitude, Latitude=Data_CC4.3$Latitude))
sp2<-STIDF(sp, time=Data_CC4.3$Date, data=data.frame(Residuals=resid_norm_CC_gam8d7b_AR7))
CC_gam8d7b_AR7_vario<-variogramST(Residuals~1, data=sp2, tunit="weeks", cores=4, tlags=(30/7)*1:10)

p_time<-ggplot(CC_gam8d7b_AR7_vario, aes(x=timelag, y=gamma, color=avgDist, group=avgDist))+
  geom_line()+
  geom_point()+
  scale_color_viridis_c(name="Distance")+
  xlab("Time difference")+
  theme_bw()+
  theme(legend.justification = "left")

p_space<-ggplot(CC_gam8d7b_AR7_vario, aes(x=dist, y=gamma, color=timelag, group=timelag))+
  geom_line()+
  geom_point()+
  scale_color_viridis_c(name="Time difference")+
  xlab("Distance")+
  theme_bw()+
  theme(legend.justification = "left")

p_variogram<-p_time/p_space

ggsave(p_variogram, filename="C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/figures/CC_gam8d7b_AR7 variogram.png",
       device="png", width=8, height=5, units="in")
  
#timelag is the average of successive specified tlags. It must be computing correlation between tlags. 

gam.check(CC_gam8d7b_AR7)
# s(Time_num_s) is OK (p=0.460)
# te(Julian_day_s,Latitude_s,Longitude_s):Year_s is OK: significant p-value but edf is 82.9 compared ot k' of 300.0
# te(Julian_day_s,Latitude_s,Longitude_s) may be able to be improved (edf=223.3 and k'=299.0)

# Check if higher k on te(Julian_day_s,Latitude_s,Longitude_s) would help improve model

CC_gam8d7b_AR8 <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(50, 13)) + 
                        te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=Year_s) + 
                        s(Time_num_s, k=5), family=scat, rho=r6, AR.start=Start, data = Data_CC4.3, method="fREML", discrete=T, nthreads=2)
#AIC: 199646.3
#BIC: 203599.9

###CC_gam8d7b_AR8 has a lower AIC but predicted slope values are almost identical to CC_gam8d7b_AR7, so using CC_gam8d7b_AR7 as the final model

newdata<-readRDS("Temperature analysis/Prediction Data.Rds")
CC_newdata<-newdata%>%
  st_drop_geometry()%>%
  as_tibble()%>%
  select(-Year_fac, -Year, -Year_s, -N,)%>%
  distinct()%>%
  mutate(Year=2000,
         Year_s=(Year-mean(Data$Year))/sd(Data$Year),
         Year_fac="2000",
         Month=as.integer(as.factor(Julian_day)))

#saveRDS(CC_newdata, "Temperature analysis/Climate Change Prediction Data.Rds")
CC_newdata<-readRDS("Temperature analysis/Climate Change Prediction Data.Rds")
# For filtering the newdata after predictions

CC_pred<-predict(CC_gam8d7b_AR7, newdata=CC_newdata, type="terms", se.fit=TRUE, discrete=T, n.threads=4)


Delta<-st_read("Delta Subregions")%>%
  select(SubRegion)%>%
  filter(SubRegion%in%unique(Data_CC4$SubRegion))

newdata_CC_pred<-CC_newdata%>%
  mutate(Slope=CC_pred$fit[,"te(Julian_day_s,Latitude_s,Longitude_s):Year_s"],
         Slope_se=CC_pred$se.fit[,"te(Julian_day_s,Latitude_s,Longitude_s):Year_s"],
         Intercept=CC_pred$fit[,"te(Julian_day_s,Latitude_s,Longitude_s)"]+CC_pred$fit[,"s(Time_num_s)"])%>%
  mutate(across(c(Slope, Slope_se), ~(.x/Year_s)/sd(Data$Year)))%>%
  mutate(Slope_se=abs(Slope_se))%>%
  mutate(Date=as.Date(Julian_day, origin=as.Date(paste(Year, "01", "01", sep="-"))),
         Month=month(Date, label = T),
         Slope_l95=Slope-Slope_se*qnorm(0.995),
         Slope_u95=Slope+Slope_se*qnorm(0.995))%>%
  mutate(Sig=if_else(Slope_u95>0 & Slope_l95<0, "ns", "*"))%>%
  filter(Sig=="*")%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=F)%>%
  st_transform(crs=st_crs(Delta))


# Create dataframe of slope deviations for each month and region f --------

Slope_summary<-CC_newdata%>%
  mutate(Slope=CC_pred$fit[,"te(Julian_day_s,Latitude_s,Longitude_s):Year_s"],
         Slope_se=CC_pred$se.fit[,"te(Julian_day_s,Latitude_s,Longitude_s):Year_s"])%>%
  mutate(across(c(Slope, Slope_se), ~(.x/Year_s)/sd(Data$Year)))%>%
  mutate(Slope_se=abs(Slope_se))%>%
  mutate(Date=as.Date(Julian_day, origin=as.Date(paste(Year, "01", "01", sep="-"))),
         Month=month(Date),
         Slope_l95=Slope-Slope_se*qnorm(0.995),
         Slope_u95=Slope+Slope_se*qnorm(0.995))%>%
  mutate(Sig=if_else(Slope_u95>0 & Slope_l95<0, "ns", "*"))%>%
  group_by(Month, SubRegion)%>%
  summarise(Slope=mean(Slope), Slope_se=mean(Slope_se), .groups="drop")%>%
  mutate(Slope_mean=mean(Slope), 
         Slope_mult=Slope/Slope_mean,
         Slope_se_mean=mean(Slope_se), 
         Slope_se_mult=Slope/Slope_se_mean)
saveRDS(Slope_summary, "Temperature analysis/Slope summary.Rds")

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

newdata_CC_pred_rast<-Rasterize_all(newdata_CC_pred, Slope)

p_CC_gam<-ggplot()+
  geom_stars(data=newdata_CC_pred_rast)+
  facet_wrap(~month(Date, label=T), drop=F)+
  scale_fill_viridis_c(breaks=(-6:7)/100, name="Temperature change\nper year (°C)", guide=guide_colorbar(barheight=20), na.value="white")+
  ylab("Latitude")+
  xlab("Longitude")+
  coord_sf()+
  theme_bw()+
  theme(strip.background=element_blank(), axis.text.x = element_text(angle=45, hjust=1), panel.grid=element_blank())

ggsave(p_CC_gam, filename="C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/figures/CC_gam 12.14.20b.png",
       device="png", width=7, height=5, units="in")


# Recreate gam check plots ------------------------------------------------

#QQ plot
resids <- resid_gam(CC_gam8d7b_AR7, incl_na=TRUE)
quantiles<-qq.gam(CC_gam8d7b_AR7)
qq_data<-as_tibble(qqplot(test, resids, plot.it=FALSE))

p_qq<-ggplot(qq_data, aes(x=x, y=y))+
  geom_abline(intercept=0, slope=1, color="firebrick3", size=1)+
  geom_point(size=0.5)+
  xlab("Theoretical quantiles")+
  ylab("Residuals")+
  theme_bw()

p_hist<-ggplot(data.frame(Residuals=resids), aes(x=Residuals))+
  geom_histogram()+
  scale_y_continuous(expand=expansion(mult=c(0,0.1)))+
  theme_bw()

linpred<-napredict(CC_gam8d7b_AR7$na.action, CC_gam8d7b_AR7$linear.predictors)

p_pred_resid<-ggplot(tibble(Predictor=linpred, Residuals=resids), aes(x=Predictor, y=Residuals))+
  geom_point(size=0.5)+
  xlab("Linear predictor")+
  theme_bw()

fitted<-predict(CC_gam8d7b_AR7, type = "response")

p_fitted_resid<-ggplot(data=tibble(Fitted=fitted, Response=Data_CC4.3$Temperature), aes(x=Fitted, y=Response))+
  geom_point(size=0.5)+
  xlab("Fitted values")+
  theme_bw()

require(patchwork)

p_check<-(p_qq|p_pred_resid)/(p_hist|p_fitted_resid)

ggsave(p_check, filename="C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/figures/CC_gam8d7b_AR7 model validation.png",
       device="png", width=10, height=7, units="in")

# Visualize raw climate change signal -------------------------------------

Times<-readRDS("Shiny app/Time_correction.Rds")
Data_CC4<-readRDS("Temperature analysis/Data_CC4.Rds")

Data_CC_plot<-Data_CC4%>%
  group_by(Station, Month)%>%
  mutate(N=n(), Anomoly=Temperature-mean(Temperature))%>%
  ungroup()%>%
  mutate(Time=as.character(round(Time_num_s, 1)))%>%
  left_join(Times, by=c("Month", "Time"))%>%
  mutate(Temp2=Temperature+Correction)%>%
  group_by(Station, Month)%>%
  mutate(Anomoly2=Temp2-mean(Temp2))%>%
  ungroup()

CC_data_plot<-function(Month, yrange=range(filter(test, Month%in%c(9,10,11))$Anomoly2)){
  ggplot(filter(test, Month==Month), aes(x=Year, y=Anomoly2))+
    geom_point()+
    facet_wrap(~SubRegion)+
    ylab("Temperature anomaly")+
    scale_x_continuous(breaks=c(1970, 1990, 2010))+
    coord_cartesian(ylim=yrange)+
    ggtitle(month(Month, label=TRUE, abbr = FALSE))+
    theme_bw()+
    theme(text=element_text(size=16))
}