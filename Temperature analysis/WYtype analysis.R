require(dplyr)
require(tidyr)
require(mgcv)
require(itsadug)
require(sf)
require(stars)
require(lubridate)
require(colorspace)
require(ggplot2)
source("Utility_functions.R")

WY<-waterYearType::water_year_indices%>%
  bind_rows(tibble(WY=rep(c(2018, 2019, 2020), each=2), location=rep(c("Sacramento Valley", "San Joaquin Valley"), 3), 
                   Index=c(7.14, 3.03, 10.34, 4.94, 6.13, 2.35), WYsum=c(12.86, 4.76, 24.77, 9.28, 9.71, 3.02)))%>% # Add data for 2018 - 2020
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
  mutate(across(c(WY, WYsum), list(s=~(.x-mean(unique(.x)))/sd(unique(.x)))))

Delta<-st_read("Delta Subregions")%>%
  select(SubRegion)%>%
  filter(SubRegion%in%unique(Data_D$SubRegion))

D_gam_NOAR <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                    te(Latitude_s, Longitude_s, Julian_day_s, WYsum_s, d=c(2,1,1), bs=c("tp", "cc", "tp"), k=c(25, 13, 5), by=WY_s) + 
                    s(Time_num_s, k=5), family=scat, data = Data_D, method="fREML", discrete=T, nthreads=5)
D_r <- start_value_rho(D_gam_NOAR, plot=TRUE)

D_gam_AR <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                  te(Latitude_s, Longitude_s, Julian_day_s, WYsum_s, d=c(2,1,1), bs=c("tp", "cc", "tp"), k=c(25, 13, 5), by=WY_s) + 
                  s(Time_num_s, k=5), family=scat, rho=D_r, AR.start=Start, data = Data_D, method="fREML", discrete=T, nthreads=5)

newdata<-readRDS("Temperature analysis/Prediction Data.Rds")

##### Find SubRegion and Month combinations with representation in the data
Data_D_effort<-Data_D%>%
  distinct(SubRegion, Month)%>%
  mutate(Keep=TRUE) 

WYsum_s_vals<-quantile(filter(WY, WY>=min(Data_D$WY))$WYsum, probs=(1:9)/10)

D_newdata<-newdata%>%
  st_drop_geometry()%>%
  as_tibble()%>%
  select(-Year_fac, -Year, -Year_s, -N,)%>%
  distinct()%>%
  mutate(WYsum=WYsum_s_vals[1],
         WY_s=2,
         Month=as.integer(as.factor(Julian_day)))%>%
  left_join(Data_D_effort, by=c("Month", "SubRegion"))%>%
  filter(Keep)%>% # Only retain SubRegion and Month combinations with representation in the data
  select(-Keep, -Season)%>%
  expand(nesting(Location, Julian_day, Time_num, Longitude, Latitude, Latitude_s, Longitude_s, Julian_day_s, Time_num_s, SubRegion, WY_s, Month), WYsum=WYsum_s_vals)%>%
  mutate(WYsum_s=(WYsum-mean(unique(Data_D$WYsum)))/sd(unique(Data_D$WYsum)))

base<-D_newdata%>%
  select(-WYsum, -WYsum_s)%>%
  distinct()%>%
  mutate(Date=parse_date_time(paste("2000", Month, "15", sep="-"), "%Y-%m-%d"))%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)%>%
  st_transform(crs=st_crs(Delta))%>%
  Rasterize_all(Location, region=Delta)%>%
  st_as_sf(long=T, connect8=T)%>%
  filter(!is.na(Location))

D_pred<-predict(D_gam_AR, newdata=D_newdata, type="terms", se.fit=TRUE, discrete=T, n.threads=4)

D_newdata_pred<-D_newdata%>%
  mutate(Slope=D_pred$fit[,"te(WYsum_s,Julian_day_s,Latitude_s,Longitude_s):WY_s"],
         Slope_se=D_pred$se.fit[,"te(WYsum_s,Julian_day_s,Latitude_s,Longitude_s):WY_s"])%>%
  mutate(across(c(Slope, Slope_se), ~(.x/WY_s)/sd(unique(Data_D$WY))))%>%
  mutate(Slope_se=abs(Slope_se))%>%
  mutate(Slope_l99.9=Slope-Slope_se*qnorm(0.9995),
         Slope_u99.9=Slope+Slope_se*qnorm(0.9995),
         Date=as.Date(Julian_day, origin=as.Date(paste("2000", "01", "01", sep="-"))),
         Month=month(Date, label=T))%>%
  mutate(Sig=if_else(Slope_u99.9>0 & Slope_l99.9<0, "ns", "*"))%>%
  filter(Sig=="*")%>%
  complete(nesting(Latitude, Longitude, Date), WYsum_s)%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=F)%>%
  st_transform(crs=st_crs(Delta))

Rasterize_all_extra<-function(Data, Var, var2, Region=Delta, Out_crs=4326, N=100){
  require(stars)
  require(purrr)
  require(sf)
  Var<-rlang::enquo(Var)
  var2vals<-unique(pull(Data, {{var2}}))
  preds<-map(var2vals, function(x) Rasterize_all(data=filter(Data, {{var2}}==x), var=!!Var, region=Region, out_crs=Out_crs, n=N))
  
  # Then bind all dates together into 1 raster
  out <- exec(c, !!!preds, along=list(Var=var2vals))
  return(out)
}

newdata_D_pred_rast<-Rasterize_all_extra(D_newdata_pred, Var=Slope, var2=WYsum_s, Region=Delta)

p_D<-ggplot()+
  geom_sf(data=base, color=NA, fill="gray80", lwd=0)+
  geom_stars(data=newdata_D_pred_rast)+
  facet_grid(month(Date, label=T)~round(Var, 2), drop=F)+
  scale_fill_continuous_diverging(palette="Blue-Red 3", na.value=NA)+
  ylab("Latitude")+
  xlab("Longitude")+
  coord_sf()+
  theme_bw()+
  theme(strip.background=element_blank(), axis.text = element_blank(), axis.ticks=element_blank(), panel.grid=element_blank())

p_D

ggsave(p_D, filename="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/WYtype2.png",
       device="png", width=15, height=15, units="in")

# Plot histogram of WYsum

ggplot(data=filter(WY, WY>=min(Data_D$WY)), aes(x=WYsum))+
  geom_density()+
  geom_point(y=0, alpha=0.1)+
  geom_vline(data=data.frame(WYsum_s_vals), aes(xintercept=WYsum_s_vals), color="red")+
  theme_bw()

ggplot(data=filter(WY, WY>=min(Data_D$WY)), aes(y=WYsum, x=WY))+
  geom_point()+
  geom_line()+
  geom_hline(data=data.frame(WYsum_s_vals), aes(yintercept=WYsum_s_vals), color="red")+
  theme_bw()
