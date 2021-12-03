require(dplyr)
require(lubridate)
require(tidyr)
require(mgcv)
require(readr)

Data<-readRDS("Temperature analysis/Discrete Temp Data.Rds")

modellf<-readRDS("C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/Models/modellf.Rds")
Noon<-12*3600

Times<-expand_grid(Time_num=c(seq(5*60*60, 21*60*60, by=10)),
                   Year_fac="2000", 
                   Longitude_s=0, 
                   Latitude_s=0, 
                   Julian_day=yday(ymd(paste("2001", 1:12, "15", sep="-"))))%>%
  mutate(Julian_day_s=(Julian_day-mean(Data$Julian_day))/sd(Data$Julian_day),
         Time_num_s=(Time_num-mean(Data$Time_num))/sd(Data$Time_num))

Time_factor<-predict(modellf, newdata=Times, type="terms", terms="te(Time_num_s,Julian_day_s)")

Time_correction<-Times%>%
  mutate(Correction=as.vector(Time_factor))%>%
  mutate(Month=as.integer(as.factor(Julian_day)))%>%
  group_by(Month)%>%
  mutate(Noon=Correction[which(Time_num==Noon)])%>%
  ungroup()%>%
  mutate(Correction=Noon-Correction,
         Time_PST=as.integer(Time_num-60*60))%>% # Convert to integer to help with joining and subtract 1 hour to convert from PDT to PST
  select(Time_PST, Month, Correction)

ggplot(Time_correction)+
  geom_point(aes(x=Time_PST, y=Correction))+
  facet_wrap(~Month)

write_csv(Time_correction, "Temperature analysis/Model outputs and validations/Time_correction_PST.csv")
