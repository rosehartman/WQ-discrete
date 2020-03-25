library(tidyverse)
library(deltareportr)
library(mgcv)
library(lubridate)
library(hms)
library(gratia)
library(sf)

Data <- DeltaDater(Start_year = 1900, 
                   WQ_sources = c("EMP", "TNS", "FMWT", "EDSM", "SKT", "20mm", "Suisun"), 
                   Variables = "Water quality", 
                   Regions = NULL)%>%
  filter(!is.na(Temperature) & !is.na(Datetime) & !is.na(Latitude) & !is.na(Longitude) & !is.na(Date))%>%
  mutate(Datetime = with_tz(Datetime, tz="America/Phoenix"),
         Date = with_tz(Date, tz="America/Phoenix"),
         Julian_day = yday(Date))%>%
  mutate(Date_num = as.numeric(Date),
         Time = as_hms(Datetime))%>%
  mutate(Time_num=as.numeric(Time))%>%
  mutate_at(vars(Date_num, Longitude, Latitude, Time_num, Year, Julian_day), list(s=~(.-mean(., na.rm=T))/sd(., na.rm=T)))

model <- gam(Temperature ~ t2(Date_num_s, Longitude_s, Latitude_s, d=c(1,2)) + t2(Julian_day, bs="cc") + poly(Time_num_s, 2),
             data = Data, method="REML")

modelb <- gam(Temperature ~ s(Longitude_s, Latitude_s) + s(Date_num_s) + ti(Date_num_s, Longitude_s, Latitude_s) + s(Julian_day, bs="cc") + poly(Time_num_s, 2),
             data = Data, method="REML")

modelc <- gam(Temperature ~ s(Longitude_s, Latitude_s) + s(Date_num_s) + ti(Date_num_s, Longitude_s, Latitude_s, d=c(1,2)) + s(Julian_day, bs="cc") + poly(Time_num_s, 2),
              data = Data, method="REML")

modeld <- gam(Temperature ~ s(Longitude_s, Latitude_s) + s(Date_num_s) + s(Julian_day, bs="cc") + poly(Time_num_s, 2),
              data = Data, method="REML")

modele <- gam(Temperature ~ s(Longitude_s, Latitude_s) + s(Year_s) + te(Julian_day_s, Time_num_s, bs="cc"),
              data = Data, method="REML")

modelf <- gam(Temperature ~ s(Longitude_s, Latitude_s) + s(Year_s) + ti(Year_s, Longitude_s, Latitude_s, d=c(1,2)) + te(Julian_day_s, Time_num_s, bs="cc"),
              data = Data, method="REML")

modelg <- gam(Temperature ~ te(Year_s, Longitude_s, Latitude_s, d=c(1,2)) + te(Julian_day_s, Time_num_s, bs="cc"),
              data = Data, method="REML")

modelh <- gamm(Temperature ~ s(Longitude_s, Latitude_s) + s(Date_num_s) + ti(Date_num_s, Longitude_s, Latitude_s) + s(Julian_day, bs="cc") + poly(Time_num_s, 2),
              data = Data, correlation = corCAR1(form = ~ Date_num_s), method="REML")

gam.check(model)
concurvity(model, full=TRUE) # Check for values over 0.8
#If concurvity is bad, run again with full=FALSE

modeld2 <- gam(Temperature ~ s(Longitude_s, Latitude_s) + s(Date_num_s, k=120) + s(Julian_day, bs="cc") + poly(Time_num_s, 2),
              data = Data, method="REML")

modele2 <- gam(Temperature ~ s(Longitude_s, Latitude_s, k=80) + s(Year_s, k=15) + te(Julian_day_s, Time_num_s, bs="cc", k=10),
              data = Data, method="REML")

modelf2 <- gam(Temperature ~ s(Longitude_s, Latitude_s, k=80) + s(Year_s, k=15) + ti(Year_s, Longitude_s, Latitude_s) + te(Julian_day_s, Time_num_s, bs="cc", k=10),
               data = Data, method="REML")

modelb2 <- gam(Temperature ~ s(Longitude_s, Latitude_s, k=60) + s(Date_num_s, k=40) + ti(Date_num_s, Longitude_s, Latitude_s) + s(Julian_day, bs="cc", k=16) + poly(Time_num_s, 2),
              data = Data, method="REML")

modelb3 <- gam(Temperature ~ s(Longitude_s, Latitude_s, k=80) + s(Date_num_s, k=200) + ti(Date_num_s, Longitude_s, Latitude_s) + s(Julian_day, bs="cc", k=16) + poly(Time_num_s, 2),
               data = Data, method="REML")

modelg2 <- gam(Temperature ~ te(Year_s, Longitude_s, Latitude_s, d=c(1,2), k=c(15, 40)) + te(Julian_day_s, Time_num_s, bs="cc", k=10),
              data = Data, method="REML")

gam.check(model2)
concurvity(model2, full=TRUE)

#model 2 seems good
plot(model2, all.terms=TRUE, residuals=TRUE, shade=TRUE)
#vis.gam

Data<-Data%>%
  mutate(Residuals = residuals(modele2),
         Fitted=fitted(modele2))%>%
  mutate(Flag=if_else(abs(Data$Residuals)>sd(Data$Residuals)*3, TRUE, FALSE))

Space<-evaluate_smooth(modele2, "s(Longitude_s,Latitude_s)", dist=0.05)%>%
  mutate(Latitude = Latitude_s*sd(Data$Latitude)+mean(Data$Latitude),
         Longitude = Longitude_s*sd(Data$Longitude)+mean(Data$Longitude))

ggplot(Space)+
  geom_raster(aes(x=Longitude, y=Latitude, fill=est))+
  geom_contour(aes(x=Longitude, y=Latitude, z=est))+ 
  geom_sf(data=deltareportr::deltaregions%>%st_transform(crs=4326), aes(color=Stratum))+
  scale_fill_distiller(palette = "RdBu", type = "div")+
  guides(fill = guide_colourbar(title = "Effect",
                                direction = "vertical",
                                barheight = grid::unit(0.25, "npc")))
  ageom_ribbon(data=Year, aes(x=Year_s, ymax=est+2*se, ymin=est-2*se), alpha=0.1)

ggplot(data=Data)+
  geom_point(aes(x=Temperature, y=Fitted, fill=Flag), shape=21)
