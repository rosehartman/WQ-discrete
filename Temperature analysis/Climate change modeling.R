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
require(sf)
require(scales)

mygrid <- data.frame(
  name = c("Upper Sacramento River Ship Channel", "Cache Slough and Lindsey Slough", "Lower Sacramento River Ship Channel", "Liberty Island", "Suisun Marsh", "Middle Sacramento River", "Lower Cache Slough", "Steamboat and Miner Slough", "Upper Mokelumne River", "Lower Mokelumne River", "Georgiana Slough", "Sacramento River near Ryde", "Sacramento River near Rio Vista", "Grizzly Bay", "West Suisun Bay", "Mid Suisun Bay", "Honker Bay", "Confluence", "Lower Sacramento River", "San Joaquin River at Twitchell Island", "San Joaquin River at Prisoners Pt", "Disappointment Slough", "Lower San Joaquin River", "Franks Tract", "Holland Cut", "San Joaquin River near Stockton", "Mildred Island", "Middle River", "Old River", "Upper San Joaquin River", "Grant Line Canal and Old River", "Victoria Canal"),
  row = c(2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 8, 8, 8),
  col = c(7, 4, 6, 5, 2, 8, 6, 7, 9, 9, 8, 7, 6, 2, 1, 2, 3, 4, 5, 6, 8, 9, 5, 6, 7, 9, 8, 8, 7, 9, 8, 7),
  code = c(" 1", " 1", " 2", " 3", " 8", " 4", " 5", " 6", " 7", " 9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "30", "29", "31"),
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

# Fit overall smoothing model to this reduced dataset to identify optimal k parameters

is.even <- function(x) as.integer(x) %% 2 == 0


model_CC4 <- bam(Temperature ~ Year_fac + te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20), by=Year_fac) + 
                   s(Time_num_s, k=5),
                 data = filter(Data_CC4, is.even(Year))%>%mutate(Year_fac=droplevels(Year_fac)), method="fREML", discrete=T, nthreads=4)

model_CC4b <- bam(Temperature ~ Year_fac + te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 20), by=Year_fac) + 
                   s(Time_num_s, k=5), family=scat,
                 data = filter(Data_CC4, is.even(Year))%>%mutate(Year_fac=droplevels(Year_fac)), method="fREML", discrete=T, nthreads=4)



CC_gam8d7b <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                    te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=Year_s) + 
                    s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num2|Station),
                  data = Data_CC4, method="REML")
resid_norm_CC_gam8d7b<-residuals(CC_gam8d7b$lme,type="normalized")
sp <- SpatialPoints(coords=data.frame(Longitude=Data_CC4$Longitude, Latitude=Data_CC4$Latitude))
sp2<-STIDF(sp, time=Data_CC4$Date, data=data.frame(Residuals=resid_norm_CC_gam8d7b))
CC_gam8d7b_vario<-variogramST(Residuals~1, data=sp2, tunit="weeks", cores=4, tlags=(30/7)*1:10)

ggplot(CC_gam8d7b_vario, aes(x=timelag, y=gamma, color=avgDist, group=avgDist))+
  geom_line()+
  geom_point()

#AIC: 201249.1
#BIC: 201375.1

CC_gam8d10<-gamm(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(15, 13)) + 
                   te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(15, 13), by=Year_s) + 
                   s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num2|Station),
                 data = Data_CC4, method="REML")
#AIC: 201429.3
#BIC: 201555.3

CC_gam8d11 <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(50, 13)) + 
                     te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(50, 13), by=Year_s) + 
                     s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num2|Station),
                   data = Data_CC4, method="REML")
#AIC: 201049.7
#BIC: 201175.7

CC_gam8d12<-gamm(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(8, 7)) + 
                   te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(8, 7), by=Year_s) + 
                   s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num2|Station),
                 data = Data_CC4, method="REML")

#Fit model

CC_gam8d13<-gamm(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(50, 13)) + 
                   te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(80, 13), by=Year_s) + 
                   s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num2|Station),
                 data = Data_CC4, method="REML")

# Predicted outputs from this model are indistinguishable from CC_gam8d11

# then predict over a range of locations and months using type="terms" or type="iterms" to get value of the Year slope for each location and month


# Inspect residuals for suitability of k-value ----------------------------

resd_CC_gam8d7b<-residuals(CC_gam8d7b$lme,type="normalized")

CC_gam8d7b_residA <- bam(resd_CC_gam8d7b ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("cs", "cc"), k=c(50, 13)),
                    gamma=1.4, data = Data_CC4, method="fREML", discrete=T, nthreads=3)
summary(CC_gam8d7b_residA)
gam.check(CC_gam8d7b_residA)
#Very low R2 (0.01) and no visible correlation between fitted values and response, so increased k-value here is not improving model fit.

CC_gam8d7b_residB <- bam(resd_CC_gam8d7b ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("cs", "cc"), k=c(50, 13), by=Year_s),
                         gamma=1.4, data = Data_CC4, method="fREML", discrete=T, nthreads=3)
summary(CC_gam8d7b_residB)
gam.check(CC_gam8d7b_residB)
#Very low R2 (0.002) and no visible correlation between fitted values and response, so increased k-value here is not improving model fit.




resd_CC_gam8d10<-residuals(CC_gam8d10$lme,type="normalized")

CC_gam8d10_residA <- bam(resd_CC_gam8d10 ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("cs", "cc"), k=c(25, 13)),
                         gamma=1.4, data = Data_CC4, method="fREML", discrete=T, nthreads=3)
summary(CC_gam8d10_residA)
gam.check(CC_gam8d10_residA)
#Very low R2 (0.01) and no visible correlation between fitted values and response, so increased k-value here is not improving model fit.

CC_gam8d10_residB <- bam(resd_CC_gam8d10 ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("cs", "cc"), k=c(25, 13), by=Year_s),
                         gamma=1.4, data = Data_CC4, method="fREML", discrete=T, nthreads=3)
summary(CC_gam8d10_residB)
gam.check(CC_gam8d10_residB)
#Very low R2 (0.001) and no visible correlation between fitted values and response, so increased k-value here is not improving model fit.



resd_CC_gam8d12<-residuals(CC_gam8d12$lme,type="normalized")

CC_gam8d12_residA <- bam(resd_CC_gam8d12 ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("cs", "cc"), k=c(15, 13)),
                         gamma=1.4, data = Data_CC4, method="fREML", discrete=T, nthreads=3)
summary(CC_gam8d12_residA)
gam.check(CC_gam8d12_residA)
#Very low R2 (0.01) and no visible correlation between fitted values and response, so increased k-value here is not improving model fit.

CC_gam8d12_residB <- bam(resd_CC_gam8d12 ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("cs", "cc"), k=c(15, 13), by=Year_s),
                         gamma=1.4, data = Data_CC4, method="fREML", discrete=T, nthreads=3)
summary(CC_gam8d12_residB)
gam.check(CC_gam8d12_residB)
#Very low R2 (0.001) and no visible correlation between fitted values and response, so increased k-value here is not improving model fit.


CC_gam8d12_residC <- bam(resd_CC_gam8d12 ~ s(Latitude_s, Longitude_s, bs="tp", k=15),
                         gamma=1.4, data = Data_CC4, method="fREML", discrete=T, nthreads=3)
summary(CC_gam8d12_residC)
gam.check(CC_gam8d12_residC)


CC_gam8d12_residD <- bam(resd_CC_gam8d12 ~ s(Julian_day_s, bs="cc", k=13),
                         gamma=1.4, data = Data_CC4, method="fREML", discrete=T, nthreads=3)
summary(CC_gam8d12_residD)
gam.check(CC_gam8d12_residD)


CC_gam8d12_residE <- bam(resd_CC_gam8d12 ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("cs", "cc"), k=c(25, 13))+
                           te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("cs", "cc"), k=c(25, 13), by=Year_s),
                         gamma=1.4, data = Data_CC4, method="fREML", discrete=T, nthreads=3)
summary(CC_gam8d12_residE)
gam.check(CC_gam8d12_residE)


CC_gam8d12_residF <- bam(resd_CC_gam8d12 ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("cs", "cc"), k=c(50, 13))+
                           te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("cs", "cc"), k=c(50, 13), by=Year_s),
                         gamma=1.4, data = Data_CC4, method="fREML", discrete=T, nthreads=3)
summary(CC_gam8d12_residF)
gam.check(CC_gam8d12_residF)


CC_gam8d12_residG <- bam(resd_CC_gam8d12 ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("cs", "cc"), k=c(50, 13))+
                           te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("cs", "cc"), k=c(80, 13), by=Year_s),
                         gamma=1.4, data = Data_CC4, method="fREML", discrete=T, nthreads=3)
summary(CC_gam8d12_residG)
gam.check(CC_gam8d12_residG)

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

CC_pred<-predict(CC_gam8d13$gam, newdata=CC_newdata, type="terms", se.fit=TRUE, discrete=T, n.threads=4)

newdata_CC_pred<-CC_newdata%>%
  mutate(Slope=CC_pred$fit[,"te(Latitude_s,Longitude_s,Julian_day_s):Year_s"],
         Slope_se=CC_pred$se.fit[,"te(Latitude_s,Longitude_s,Julian_day_s):Year_s"],
         Intercept=CC_pred$fit[,"te(Latitude_s,Longitude_s,Julian_day_s)"]+CC_pred$fit[,"s(Time_num_s)"])%>%
  mutate(across(c(Slope, Slope_se), ~(.x/Year_s)/sd(Data$Year)))%>%
  mutate(Slope_se=abs(Slope_se))%>%
  mutate(Date=as.Date(Julian_day, origin=as.Date(paste(Year, "01", "01", sep="-"))),
         Month=month(Date, label = T),
         Slope_l95=Slope-Slope_se*qnorm(0.995),
         Slope_u95=Slope+Slope_se*qnorm(0.995))%>%
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