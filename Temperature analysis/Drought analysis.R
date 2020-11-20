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
require(sf)

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
  mutate(WYsum_c=WYsum-mean(WYsum),
         WY_s=(WY-mean(unique(WY)))/sd(unique(WY)))

D_gam1 <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                     te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=WYsum_c) + 
                     s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num2|Station),
                   data = Data_D, method="REML")
#AIC: 196747.1
#BIC: 196872.9

D_gam2 <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(40, 13)) + 
                 te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(40, 13), by=WYsum_c) + 
                 s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num2|Station),
               data = Data_D, method="REML")
#AIC: 196583.1
#BIC: 196708.9

D_gam3 <- gamm(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                 te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=WYsum_c) + 
                 te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=WY_s) + 
                 s(Time_num_s, k=5), correlation = corCAR1(form=~Date_num2|Station),
               data = Data_D, method="REML")

# Checking k-values

resd<-residuals(D_gam1$lme,type="normalized")

D_gam1_resid <- bam(resd ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("cs", "cc"), k=c(50, 13)),
                    gamma=1.4, data = Data_D, method="fREML", discrete=T, nthreads=3)
gam.check(D_gam1_resid)
summary(D_gam1_resid)
#Very low R2 (0.01) and no visible correlation between fitted values and response, so increased k-value here is not improving model fit.

D_gam1_resid2 <- bam(resd ~  te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("cs", "cc"), k=c(50, 13), by=WYsum_c),
                    gamma=1.4, data = Data_D, method="fREML", discrete=T, nthreads=3)

gam.check(D_gam1_resid2)
summary(D_gam1_resid2)
#Very low R2 (0.0001) and no visible correlation between fitted values and response, so increased k-value here is not improving model fit.

# Predict effect

newdata<-readRDS("Temperature analysis/Prediction Data.Rds")
D_newdata<-newdata%>%
  st_drop_geometry()%>%
  as_tibble()%>%
  select(-Year_fac, -Year, -Year_s, -N,)%>%
  distinct()%>%
  mutate(WYsum_c=10,
         WY_s=2)

D_pred<-predict(D_gam3$gam, newdata=D_newdata, type="terms", se.fit=TRUE, discrete=T, n.threads=4)

D_newdata_pred<-D_newdata%>%
  mutate(Slope=D_pred$fit[,"te(Latitude_s,Longitude_s,Julian_day_s):WYsum_c"],
         Slope_se=D_pred$se.fit[,"te(Latitude_s,Longitude_s,Julian_day_s):WYsum_c"],
         Intercept=D_pred$fit[,"te(Latitude_s,Longitude_s,Julian_day_s)"]+D_pred$fit[,"s(Time_num_s)"],
         Slope2=D_pred$fit[,"te(Latitude_s,Longitude_s,Julian_day_s):WY_s"],
         Slope_se2=D_pred$se.fit[,"te(Latitude_s,Longitude_s,Julian_day_s):WY_s"])%>%
  mutate(across(c(Slope, Slope_se), ~(.x/WYsum_c)),
         across(c(Slope2, Slope_se2), ~(.x/WY_s)/sd(unique(Data_D$WY))))%>%
  mutate(Slope_se=abs(Slope_se),
         Slope_se2=abs(Slope_se2))%>%
  mutate(Date=as.Date(Julian_day, origin=as.Date(paste("2000", "01", "01", sep="-"))),
         Month=month(Date, label = T),
         Slope_l95=Slope-Slope_se*qnorm(0.995),
         Slope_u95=Slope+Slope_se*qnorm(0.995),
         Slope2_l95=Slope2-Slope_se2*qnorm(0.995),
         Slope2_u95=Slope2+Slope_se2*qnorm(0.995))%>%
  mutate(Sig=if_else(Slope_u95>0 & Slope_l95<0, "ns", "*"),
         Sig2=if_else(Slope2_u95>0 & Slope2_l95<0, "ns", "*"))

p_D_gam<-ggplot(filter(D_newdata_pred, Sig=="*"), aes(x=Longitude, y=Latitude, color=Slope))+
  geom_point()+
  facet_wrap(~Month)+
  scale_color_gradient2(high = muted("red"), low = muted("blue"), breaks=(-15:7)/100, guide=guide_colorbar(barheight=20))+
  theme_bw()+
  theme(strip.background=element_blank(), axis.text.x = element_text(angle=45, hjust=1))
