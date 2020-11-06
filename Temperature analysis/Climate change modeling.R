require(sp)
require(gstat)
require(spacetime)
require(tidybayes)
require(dplyr)
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

mygrid <- data.frame(
  name = c("Cache Slough and Lindsey Slough", "Lower Sacramento River Ship Channel", "Liberty Island", "Suisun Marsh", "Middle Sacramento River", "Lower Cache Slough", "Steamboat and Miner Slough", "Upper Mokelumne River", "Lower Mokelumne River", "Georgiana Slough", "Sacramento River near Ryde", "Sacramento River near Rio Vista", "Grizzly Bay", "West Suisun Bay", "Mid Suisun Bay", "Honker Bay", "Confluence", "Lower Sacramento River", "San Joaquin River at Twitchell Island", "San Joaquin River at Prisoners Pt", "Disappointment Slough", "Lower San Joaquin River", "Franks Tract", "Holland Cut", "San Joaquin River near Stockton", "Mildred Island", "Middle River", "Old River", "Upper San Joaquin River", "Grant Line Canal and Old River", "Victoria Canal"),
  row = c(2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 8, 8, 8),
  col = c(4, 6, 5, 2, 8, 6, 7, 9, 9, 8, 7, 6, 2, 1, 2, 3, 4, 5, 6, 8, 9, 5, 6, 7, 9, 8, 8, 7, 9, 8, 7),
  code = c(" 1", " 2", " 3", " 8", " 4", " 5", " 6", " 7", " 9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "30", "29", "31"),
  stringsAsFactors = FALSE
)

Data<-readRDS("Temperature analysis/Discrete Temp Data.Rds")
newdata_year<-readRDS("Temperature analysis/Prediction Data.Rds")
modellc4_predictions<-readRDS("Temperature analysis/modellc4_predictions.Rds")


# Using model-produced data  ----------------------------------------------------

newdata<-newdata_year%>%
  mutate(Prediction=modellc4_predictions$fit)%>%
  mutate(SE=modellc4_predictions$se.fit,
         L95=Prediction-SE*1.96,
         U95=Prediction+SE*1.96)%>%
  mutate(Date=as.Date(Julian_day, origin=as.Date(paste(Year, "01", "01", sep="-")))) # Create Date variable from Julian Day and Year

Data_effort <- Data%>%
  st_drop_geometry()%>%
  group_by(SubRegion, Month, Year)%>%
  summarise(N=n(), .groups="drop")

newdata_year2<-newdata%>%
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

CC_models<-newdata_sum%>%
  nest_by(Month, SubRegion)%>%
  mutate(model=list(lm(Temperature~Year, data=data)))%>%
  summarise(N=nrow(data), broom::tidy(model), .groups="drop")%>%
  filter(term=="Year")
map_CC<-function(month){
  ggplot(filter(newdata_sum, Month==month))+
    geom_point(aes(x=Year, y=Temperature))+
    geom_smooth(aes(x=Year, y=Temperature), method="lm", formula=y ~ x)+
    {if(nrow(filter(CC_models, p.value<0.05 & Month==month))>0){
      geom_text(data=filter(CC_models, p.value<0.05 & Month==month), aes(x=2018, y=max(filter(newdata_sum, Month==month)$Temperature)-1, label=if_else(p.value<0.01, "**", "*")), color="firebrick3", size=6)
    }}+
    facet_geo(~SubRegion, grid=mygrid, labeller=label_wrap_gen())+
    theme_bw()+
    theme(panel.grid=element_blank(), axis.text.x = element_text(angle=45, hjust=1))
}


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

# Remove non-fixed stations, reduce temporal frequency to no more than 1 sample per week
Data_CC<-Data%>%
  filter(Source%in%c("EMP", "STN", "FMWT", "DJFMP", "SKT", "20mm", "Suisun", "Baystudy", "USGS") & !str_detect(Station, "EZ") & !(Source=="SKT" & Station=="799" & Latitude>38.2))%>%
  mutate(Station=paste(Source, Station),
         Week=week(Date),
         Noon_diff=abs(hms(hours=12)-as_hms(Datetime)))%>%
  group_by(Station, Week, Year)%>%
  filter(Date==min(Date) & Noon_diff==min(Noon_diff))%>%
  ungroup()%>%
  lazy_dt()%>%
  group_by(Date, Date_num, Date_num_s, Week, SubRegion, Julian_day_s, Julian_day, Year, Year_s, Year_fac, Station, Source, Latitude_s, Longitude_s, Latitude, Longitude)%>%
  summarise(Temperature=mean(Temperature), Time_num=mean(Time_num), Time_num_s=mean(Time_num_s))%>%
  as_tibble()%>%
  ungroup()%>%
  mutate(ID=paste(Station, Date_num))%>%
  filter(!(ID%in%ID[which(duplicated(ID))]))%>%
  mutate(YearStation=paste(Year, Station))

saveRDS(Data_CC, "Temperature analysis/Discrete Temp Data CC.Rds")

CC_gam <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20)) + 
                te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20), by=Year_s) + 
                s(Time_num_s, k=5),
              data = Data, method="fREML", discrete=T, nthreads=4)
# Adding m=1 to first smoother and m=2 to second smoother resulted in exact same concurvity

CC_gam2 <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20)) + 
                 te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20), by=Year_s) + 
                 s(Time_num_s, k=5),
               data = filter(Data, Source!="EDSM"), method="fREML", discrete=T, nthreads=4)

CC_gam3 <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, Year_fac, d=c(2,1,1), bs=c("tp", "cc", "re"), k=c(15, 10)) + 
                 te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20), by=Year_s) + 
                 s(Time_num_s, k=5),
               data = Data, method="fREML", discrete=T, nthreads=4)
# Weird results, probably too much collinearity

CC_gam4 <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20)) + 
                  te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20), by=Year_s) + 
                  s(Time_num_s, k=5), correlation = corCAR1(form=~Date|Station),
                data = Data_CC, method="REML")

CC_gam5a <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20)) + 
                   te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20), by=Year_s) + 
                   s(Time_num_s, k=5),
                 data = Data_CC, method="REML")

# Find the optimal correlation structure
sp <- SpatialPoints(coords=data.frame(Longitude=Data_CC$Longitude, Latitude=Data_CC$Latitude))
sp2<-STIDF(sp, time=Data_CC$Date, data=data.frame(Residuals=residuals(CC_gam5a$lme, type="normalized")))
CC_gam5a_vario<-variogramST(Residuals~1, data=sp2, tunit="weeks", cores=3)

ggplot(CC_gam5a_vario, aes(x=timelag, y=gamma, color=avgDist, group=avgDist))+
  geom_line()+
  geom_point()

CC_gam5 <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20)) + 
                  te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20), by=Year_s) + 
                  s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num_s|Station),
                data = Data_CC, method="REML")
acf(residuals(CC_gam5$gam),main="raw residual ACF")
resid_norm_CC_gam5<-residuals(CC_gam5$lme,type="normalized")
acf(resid_norm_CC_gam5,main="standardized residual ACF")

sp <- SpatialPoints(coords=data.frame(Longitude=Data_CC$Longitude, Latitude=Data_CC$Latitude))
sp2<-STIDF(sp, time=Data_CC$Date, data=data.frame(Residuals=resid_norm_CC_gam5))
CC_gam5_vario<-variogramST(Residuals~1, data=sp2, tunit="weeks", cores=3)

ggplot(CC_gam5_vario, aes(x=timelag, y=gamma, color=avgDist, group=avgDist))+
  geom_line()+
  geom_point()

CC_gam6 <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20)) + 
                  te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20), by=Year_s) + 
                  s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num_s|YearStation),
                data = Data_CC, method="REML")

CC_gam7 <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20)) + 
                  te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20), by=Year_s) + 
                  s(Time_num_s, k=5), correlation = corExp(form=~Date_num_s|YearStation),
                data = Data_CC, method="REML")
# corExp Won't fit

#Try reducing temporal resolution to once every 2 weeks.
Data_CC2<-Data%>%
  filter(Source%in%c("EMP", "STN", "FMWT", "DJFMP", "SKT", "20mm", "Suisun", "Baystudy", "USGS") & !str_detect(Station, "EZ") & !(Source=="SKT" & Station=="799" & Latitude>38.2))%>%
  mutate(Station=paste(Source, Station),
         Week2=round(week(Date)/2),
         Noon_diff=abs(hms(hours=12)-as_hms(Datetime)))%>%
  group_by(Station, Week2, Year)%>%
  filter(Date==min(Date) & Noon_diff==min(Noon_diff))%>%
  ungroup()%>%
  lazy_dt()%>%
  group_by(Date, Date_num, Date_num_s, Week2, SubRegion, Julian_day_s, Julian_day, Year, Year_s, Year_fac, Station, Source, Latitude_s, Longitude_s, Latitude, Longitude)%>%
  summarise(Temperature=mean(Temperature), Time_num=mean(Time_num), Time_num_s=mean(Time_num_s))%>%
  as_tibble()%>%
  ungroup()%>%
  mutate(ID=paste(Station, Date_num))%>%
  filter(!(ID%in%ID[which(duplicated(ID))]))%>%
  mutate(YearStation=paste(Year, Station),
         Date_num2=as.numeric(Date)/(3600*24*14)) # Create a numeric date variable in units of 2 weeks. 

CC_gam5b <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20)) + 
                   te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20), by=Year_s) + 
                   s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num_s|Station),
                 data = Data_CC2, method="REML")

CC_gam6b <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20)) + 
                   te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20), by=Year_s) + 
                   s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num_s|YearStation),
                 data = Data_CC2, method="REML")

# corExp won't fit

CC_gam7b <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20)) + 
                   te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20), by=Year_s) + 
                   s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num2|Station),
                 data = Data_CC2, method="REML")
#Seems to be exactly the same using the other form of Date_num, so standardization of the date variable does not have an effect on model predictions


sp <- SpatialPoints(coords=data.frame(Longitude=Data_CC2$Longitude, Latitude=Data_CC2$Latitude))
sp2<-STIDF(sp, time=Data_CC2$Date, data=data.frame(Residuals=residuals(CC_gam7b$lme, type="normalized")))
CC_gam7b_vario<-variogramST(Residuals~1, data=sp2, tunit="weeks", cores=3, tlags=seq(2,20,by=2))

#Try reducing temporal resolution to once every month.
Data_CC3<-Data%>%
  filter(Source%in%c("EMP", "STN", "FMWT", "DJFMP", "SKT", "20mm", "Suisun", "Baystudy", "USGS") & !str_detect(Station, "EZ") & !(Source=="SKT" & Station=="799" & Latitude>38.2))%>%
  mutate(Station=paste(Source, Station),
         Noon_diff=abs(hms(hours=12)-as_hms(Datetime)))%>%
  group_by(Station, Month, Year)%>%
  filter(Date==min(Date) & Noon_diff==min(Noon_diff))%>%
  ungroup()%>%
  lazy_dt()%>%
  group_by(Date, Date_num, Date_num_s, Month, SubRegion, Julian_day_s, Julian_day, Year, Year_s, Year_fac, Station, Source, Latitude_s, Longitude_s, Latitude, Longitude)%>%
  summarise(Temperature=mean(Temperature), Time_num=mean(Time_num), Time_num_s=mean(Time_num_s))%>%
  as_tibble()%>%
  ungroup()%>%
  mutate(ID=paste(Station, Date_num))%>%
  filter(!(ID%in%ID[which(duplicated(ID))]))%>%
  mutate(YearStation=paste(Year, Station),
         Date_num2=as.numeric(Date)/(3600*24*30)) # Create a numeric date variable in units of ~ 1 month. 

CC_gam7c <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20)) + 
                   te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20), by=Year_s) + 
                   s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num2|Station),
                 data = Data_CC3, method="REML")
# BIC: 177774.8

sp <- SpatialPoints(coords=data.frame(Longitude=Data_CC3$Longitude, Latitude=Data_CC3$Latitude))
sp2<-STIDF(sp, time=Data_CC3$Date, data=data.frame(Residuals=residuals(CC_gam7c$lme, type="normalized")))
CC_gam7c_vario<-variogramST(Residuals~1, data=sp2, tunit="weeks", tlags=(30/7)*1:10, cores=3)

CC_gam7c2 <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20)) + 
                    te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(30, 20), by=Year_s) + 
                    s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num2|Station),
                  data = Data_CC3, method="REML")
# BIC: 177775.2
# No improvement over CC_gam7c

CC_gam7c3 <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20)) + 
                    te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 25), by=Year_s) + 
                    s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num2|Station),
                  data = Data_CC3, method="REML")
# BIC: 177691.2

CC_gam7c4 <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20)) + 
                    te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 35), by=Year_s) + 
                    s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num2|Station),
                  data = Data_CC3, method="REML")

# BIC: 177545.9

CC_gam7c5 <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20)) + 
                    te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 12), by=Year_s) + 
                    s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num2|Station),
                  data = Data_CC3, method="REML")

# BIC: 178235.2

#Try reducing temporal resolution to once every month.
# This time, pick days closest to the middle of the month
# And split the filtering to 2 steps to fix error that removed too much data by requiring the earliest day AND the time closest to noon
Data_CC4<-Data%>%
  filter(Source%in%c("EMP", "STN", "FMWT", "DJFMP", "SKT", "20mm", "Suisun", "Baystudy", "USGS") & !str_detect(Station, "EZ") & !(Source=="SKT" & Station=="799" & Latitude>38.2))%>%
  mutate(Station=paste(Source, Station),
         Noon_diff=abs(hms(hours=12)-as_hms(Datetime)),
         mday_15_diff=abs(mday(Date)-15))%>% # Find how far each date is from the 15th of the month
  group_by(Station, Month, Year)%>%
  filter(mday_15_diff==min(mday_15_diff))%>%
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

CC_gam8d <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Month, d=c(2,1), bs=c("tp", "cc"), k=c(25, 12)) + 
                    te(Latitude_s, Longitude_s, Month, d=c(2,1), bs=c("tp", "cc"), k=c(25, 12), by=Year_s) + 
                    s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num2|Station), knots=list(Month = c(0.5, seq(1, 12, length = 10), 12.5)),
                  data = Data_CC4, method="REML") # knots from https://fromthebottomoftheheap.net/2015/11/21/climate-change-and-spline-interactions/
# BIC: 207035.5

CC_gam8d2 <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Month, d=c(2,1), bs=c("tp", "cc"), k=c(25, 14)) + 
                   te(Latitude_s, Longitude_s, Month, d=c(2,1), bs=c("tp", "cc"), k=c(25, 14), by=Year_s) + 
                   s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num2|Station), knots=list(Month = c(0.5, 1:12, 12.5)),
                 data = Data_CC4, method="REML")
#BIC: 206880.3

CC_gam8d3 <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Month, d=c(2,1), bs=c("tp", "tp"), k=c(25, 12)) + 
                    te(Latitude_s, Longitude_s, Month, d=c(2,1), bs=c("tp", "tp"), k=c(25, 12), by=Year_s) + 
                    s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num2|Station),
                  data = Data_CC4, method="REML")
#BIC: 206967.2

CC_gam8d4 <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Month, d=c(2,1), bs=c("tp", "cc"), k=c(25, 14)) + 
                    te(Latitude_s, Longitude_s, Month_fac, d=c(2,1), bs=c("tp", "fs"), k=c(25), by=Year_s) + 
                    s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num2|Station), knots=list(Month = c(0.5, 1:12, 12.5)),
                  data = Data_CC4, method="REML")
#BIC: 206975.5

CC_gam8d5 <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Month, d=c(2,1), bs=c("tp", "cc"), k=c(25, 14)) + 
                    te(Latitude_s, Longitude_s, Month_fac, d=c(2,1), bs=c("tp", "re"), k=c(25, 12), by=Year_s) + 
                    s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num2|Station), knots=list(Month = c(0.5, 1:12, 12.5)),
                  data = Data_CC4, method="REML")
#BIC: 206975.5

CC_gam8d6 <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 14)) + 
                    te(Latitude_s, Longitude_s, Month_fac, d=c(2,1), bs=c("tp", "fs"), k=c(25), by=Year_s) + 
                    s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num2|Station),
                  data = Data_CC4, method="REML")
# BIC: 198489.9

CC_gam8d7 <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                    te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 14), by=Year_s) + 
                    s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num2|Station),
                  data = Data_CC4, method="REML")

#BIC: 198500.1

CC_gam8d7b <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                    te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=Year_s) + 
                    s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num2|Station),
                  data = Data_CC4, method="REML")

CC_gam8d8 <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20)) + 
                    te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=Year_s) + 
                    s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num2|Station),
                  data = Data_CC4, method="REML")
#BIC: 198233.8

resid_norm_CC_gam8d8<-residuals(CC_gam8d8$lme,type="normalized")
sp <- SpatialPoints(coords=data.frame(Longitude=Data_CC4$Longitude, Latitude=Data_CC4$Latitude))
sp2<-STIDF(sp, time=Data_CC4$Date, data=data.frame(Residuals=resid_norm_CC_gam8d8))
CC_gam8d8_vario<-variogramST(Residuals~1, data=sp2, tunit="weeks", cores=3, tlags=(30/7)*1:10)

CC_gam8d9 <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20)) + 
                    te(Latitude_s, Longitude_s, Month_fac, d=c(2,1), bs=c("tp", "fs"), k=c(25), by=Year_s) + 
                    s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num2|Station),
                  data = Data_CC4, method="REML")

# BIC: 198314

ggplot(CC_gam8d8_vario, aes(x=avgDist, y=gamma, color=timelag, group=timelag))+
  geom_line()+
  geom_point()

auto<-Data_CC3%>%
  mutate(Resid_raw=residuals(CC_gam7c$gam),
         Resid_norm=residuals(CC_gam7c$lme,type="normalized"))%>%
  filter(Source!="EDSM" & !str_detect(Station, "EZ"))%>% # Remove EDSM and EZ stations because they're not fixed
  mutate(Station=paste(Source, Station))%>%
  group_by(Station)%>%
  mutate(N=n())%>%
  filter(N>10)%>%
  summarise(ACF_raw=list(pacf(Resid_raw, plot=F)), ACF_norm=list(pacf(Resid_norm, plot=F)), N=n(), ci=qnorm((1 + 0.95)/2)/sqrt(n()), .groups="drop")%>% # ci formula from https://stackoverflow.com/questions/14266333/extract-confidence-interval-values-from-acf-correlogram
  rowwise()%>%
  mutate(lag_raw=list(ACF_raw$lag), acf_raw=list(ACF_raw$acf),
         lag_norm=list(ACF_norm$lag), acf_norm=list(ACF_norm$acf))%>%
  unnest(cols=c(lag_raw, acf_raw, lag_norm, acf_norm))%>%
  arrange(-N)%>%
  mutate(Station=factor(Station, levels=unique(Station)))%>%
  select(-ACF_raw, -ACF_norm)%>%
  mutate(across(c(acf_norm, lag_norm, lag_raw, acf_raw), ~as.vector(.x)))

ggplot(filter(auto, lag_norm%in%1:4))+
  geom_point(aes(x=Station, y=abs(acf_norm)), fill="black", shape=21)+
  geom_point(data=filter(auto, lag_norm%in%1:4 & abs(acf_norm)>abs(ci)), aes(x=Station, y=abs(acf_norm)), fill="red", shape=21)+
  geom_point(aes(x=Station, y=abs(ci)), fill="white", shape=21)+
  geom_segment(aes(x=Station, y=abs(acf_norm), xend=Station, yend=abs(ci)), linetype=2)+
  geom_segment(data=filter(auto, lag_norm%in%1:4 & abs(acf_norm)>abs(ci)), aes(x=Station, y=abs(acf_norm), xend=Station, yend=abs(ci)), color="red")+
  facet_wrap(~lag_norm)+
  theme_bw()+
  theme(panel.grid=element_blank(), axis.text.x=element_text(angle=45, hjust=1))

#Try fitting this on top of the prior gam?

#Fit model
# then predict over a range of locations and months using type="terms" or type="iterms" to get value of the Year slope for each location and month

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

saveRDS(CC_newdata, "Temperature analysis/Climate Change Prediction Data.Rds")
# For filtering the newdata after predictions

CC_effort<-newdata%>%
  st_drop_geometry()%>%
  as_tibble()%>%
  mutate(Date=as.Date(Julian_day, origin=as.Date(paste(Year, "01", "01", sep="-"))),
         Month=month(Date))%>%
  select(SubRegion, Year, Month)%>%
  distinct()%>%
  group_by(Month, SubRegion)%>%
  summarise(N=n(), .groups="drop")

CC_pred<-predict(CC_gam8d8$gam, newdata=CC_newdata, type="terms", se.fit=TRUE, discrete=T, n.threads=4)

newdata_CC_pred<-CC_newdata%>%
  mutate(Slope=CC_pred$fit[,"te(Latitude_s,Longitude_s,Julian_day_s):Year_s"],
         Slope_se=CC_pred$se.fit[,"te(Latitude_s,Longitude_s,Julian_day_s):Year_s"],
         Intercept=CC_pred$fit[,"te(Latitude_s,Longitude_s,Julian_day_s)"]+CC_pred$fit[,"s(Time_num_s)"])%>%
  mutate(across(c(Slope, Slope_se), ~(.x/Year_s)/sd(Data$Year)))%>%
  mutate(Slope_se=abs(Slope_se))%>%
  mutate(Date=as.Date(Julian_day, origin=as.Date(paste(Year, "01", "01", sep="-"))),
         Month=month(Date, label = T),
         Slope_l95=Slope-Slope_se*qnorm(0.975),
         Slope_u95=Slope+Slope_se*qnorm(0.975))%>%
  mutate(Sig=if_else(Slope_u95>0 & Slope_l95<0, "ns", "*"))

p_CC_gam<-ggplot(filter(newdata_CC_pred, Sig=="*"), aes(x=Longitude, y=Latitude, color=Slope))+
  geom_point()+
  facet_wrap(~Month)+
  scale_color_gradient2(high = muted("red"), low = muted("blue"), breaks=(-6:7)/100, guide=guide_colorbar(barheight=20))+
  theme_bw()+
  theme(strip.background=element_blank(), axis.text.x = element_text(angle=45, hjust=1))

ggsave(p_CC_gam, file="C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/figures/CC_gam8d2.png",
       device="png", width=7, height=5, units="in")

# Autocorrelation
require(gstat)
require(sp)
require(spacetime)

sp <- SpatialPoints(coords=data.frame(Longitude=Data$Longitude, Latitude=Data$Latitude))
sp2<-STIDF(sp, time=Data$Date, data=data.frame(Residuals=residuals(CC_brm2, type="pearson")[,1]))
CC_vario<-variogramST(Residuals~1, data=sp2, tunit="weeks", cores=3)

sp_EMP <- SpatialPoints(coords=data.frame(Longitude=filter(Data, Source=="EMP")$Longitude, Latitude=filter(Data, Source=="EMP")$Latitude))
sp2_EMP<-STIDF(sp_EMP, time=filter(Data, Source=="EMP")$Date, data=data.frame(Residuals=residuals(CC_brm_EMP, type="pearson")[,1]))
CC_vario_EMP<-variogramST(Residuals~1, data=sp2_EMP, tunit="weeks", cores=3, tlags=1:10*4)

save(CC_vario, CC_vario_EMP, CC_brm, CC_brm2, CC_brm3, CC_brm4, CC_brm_EMP, file="CC models.Rds")

ggplot(CC_vario, aes(x=timelag, y=gamma, color=avgDist, group=avgDist))+
  geom_line()+
  geom_point()

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