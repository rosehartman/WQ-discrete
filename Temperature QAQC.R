library(tidyverse)
library(deltareportr)
library(mgcv)
library(lubridate)
library(hms)

Data <- DeltaDater(Start_year = 1900, 
                   WQ_sources = c("EMP", "TNS", "FMWT", "EDSM", "SKT", "20mm", "Suisun"), 
                   Variables = "Water quality", 
                   Regions = NULL)%>%
  filter(!is.na(Temperature) & !is.na(Datetime) & !is.na(Latitude) & !is.na(Longitude) & !is.na(Date))%>%
  mutate(Datetime = with_tz(Datetime, tz="America/Phoenix"),
         Date = with_tz(Date, tz="America/Phoenix"),
         Julian_day = yday(Date))%>%
  mutate(Date_num = as.numeric(Date),
         Time = as_hms(Data$Datetime))%>%
  mutate(Time_num=as.numeric(Time))%>%
  mutate_at(vars(Date_num, Latitude, Longitude, Time_num), list(s=~(.-mean(., na.rm=T))/sd(., na.rm=T)))

model <- gam(Temperature ~ t2(Date_num, Latitude, Longitude, d=c(1,2)) + t2(Julian_day, bs="cc") + poly(Time_num_s, 2),
             data = Data)

model2 <- gam(Temperature ~ t2(Date_num, Latitude, Longitude, d=c(1,2), k=40) + t2(Julian_day, bs="cc", k=24) + poly(Time_num_s, 2),
             data = Data)

Data$Residuals <- residuals(model2)
Data$Fitted<-fitted(model2)
