require(dplyr)
require(tidyr)
require(mgcv)
require(sf)
require(stars)
require(purrr)
require(ggplot2)
require(mgcv)
require(itsadug)
require(lubridate)
# New method: Try simulating data (using the gratia simulate or predicted_sampled functions) from the global smoothing model to use for testing the model performance.
# Simulate with a balanced and imbalanced sampling design to examing that impact
# Also try simulating with all years vs. just years with similar average water temps to ensure model can differentiate real from spurious trends.

# First calculate year-to-year variability
# Start in 1975 when sampling became more regular across different regions
Data_CC4<-readRDS("Temperature analysis/Data_CC4.Rds")

Delta<-st_read("Delta Subregions")%>%
  select(SubRegion)%>%
  filter(SubRegion%in%unique(Data_CC4$SubRegion))

Data_effort<-Data_CC4%>%
  group_by(SubRegion, Month, Year)%>%
  summarise(N=n(), .groups="drop")

max_effort<-Data_effort%>%
  group_by(SubRegion)%>%
  summarise(Max=max(N), .groups="drop")

require(brms)

model_var<-brm(Temperature ~ (1|Year_fac) + (1| Month) + (1|Station), data=Data_CC4,
               family=gaussian(),
               prior=prior(normal(0,5), class="Intercept")+
                 prior(cauchy(0,5), class="sigma")+
                 prior(cauchy(0,5), class="sd"),
               control=list(adapt_delta=0.85, max_treedepth=15),
               iter=8e3, warmup=2e3, cores=3, chains=3)
# Pull value of Year_fac sd from variance model to use below.


model_sum<-summary(model)
year_sd<-model_sum$p.coeff[!names(model_sum$p.coeff) %in% c("(Intercept)", "Year_fac1970", "Year_fac1972", "Year_fac1974")]

Slope_summary<-readRDS("Temperature analysis/Slope summary.Rds")

newdata<-readRDS("Temperature analysis/Prediction Data.Rds")%>%
  filter(Year==2018)%>% # Pick year with lots of data (EDSM in effect) and almost all subregions represented (except Grant Line which ended in 1995)
  st_drop_geometry()%>%
  select(-Season, -N, -Time_num)%>%
  mutate(Month=as.integer(as.factor(Julian_day)))

newdata_stations<-Data_CC4%>%
  select(Station, Latitude, Longitude)%>%
  distinct()%>%
  mutate(Julian_day=yday(ymd(paste("2001", 1, "15", sep="-"))))%>%
  complete(Julian_day=yday(ymd(paste("2001", 1:12, "15", sep="-"))), Station)%>%
  group_by(Station)%>%
  mutate(Latitude=mean(Latitude, na.rm=T),
         Longitude=mean(Longitude, na.rm=T))%>%
  ungroup()%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=F)%>%
  st_join(Delta%>%
            st_transform(crs=4326))%>%
  st_drop_geometry()%>%
  mutate(Month=as.integer(as.factor(Julian_day)))%>%
  mutate(across(c(Julian_day, Latitude, Longitude), list(s=~(.-mean(.))/sd(.))))%>%
  mutate(Time_num_s=unique(newdata$Time_num_s),
         Year=unique(newdata$Year),
         Year_s=unique(newdata$Year_s),
         Year_fac=unique(newdata$Year_fac))%>%
  left_join(Slope_summary%>%
              select(Month, SubRegion, Slope, Slope_se), 
            by=c("Month", "SubRegion"))

# Define function to simulate gam data ------------------------------------

gamsim<-function(model, data=newdata, years=1970:2020, year_slope=0.02, year_sd=0.2, nsim=10, seeds=1:nsim){
  mu<-predict(model, newdata=data, type="response", discrete=T, n.threads=4)
  scale <- model[["sig2"]]
  data<-mutate(data, mu=mu)
  data<-map_dfr(years, ~mutate(data, Year=.x))%>%
    mutate(Date=as.Date(Julian_day, origin=as.Date(paste(Year, "01", "01", sep="-"))),
           Year_s=(Year-mean(years))/sd(years))
  
  # Extract distributional function
  family_fun<-fix.family.rd(family(model))[["rd"]]
  
  sim_fun<-function(Data=data, Years=years, Year_slope=year_slope, Year_sd=year_sd, Nsim=nsim, Scale=scale, seed){
    set.seed(seed)
    year_add<-tibble(Year=Years, add=rnorm(n=Year, mean=(Year-min(Year))*Year_slope, sd=Year_sd))
    
    Data<-Data%>%
      complete(Year=Years)%>%
      left_join(year_add, by="Year")%>%
      mutate(pred=mu+add) #############################################Maybe move the rnorm for add to this step so each data point gets some random deviation 
    
    sims <- family_fun(mu=Data$pred, wt=rep(1, length(Data$pred)), scale=Scale)
    set.seed(NULL)
    return(list(sims=sims, year_add=year_add$add))
  }
  
  out<-map(1:nsim, ~sim_fun(seed=seeds[.x]))
  
  sims<-map(out, ~.x[["sims"]])
  names(sims)<-paste("sim", 1:nsim, sep="_")
  sims<-bind_cols(data, sims)
  
  year_add<-map(out, ~.x[["year_add"]])
  names(year_add)<-paste("sim", 1:nsim, sep="_")
  year_add<-bind_cols(list(Year=years), year_add)
  
  #Try returning data with 1 column per simulation
  
  return(list(sims=sims, year_add=year_add))
  
}


# Simulate for 0.02 degree increase per year ------------------------------

#load("~/WQ-discrete/Temperature analysis/Models/Smoothing models 3.0.RData")
sim_data<-gamsim(modelld14a)

#Create two versions of simulated data with balanced and unbalanced sampling regimes

## Shuffle locations in each region first, then pull the number needed for each month, year, region to ensure
# Most are sampled continuously
set.seed(1)
Shuffled_locations<-newdata%>%
  left_join(max_effort, by="SubRegion")%>%
  group_by(SubRegion)%>%
  summarise(Locations=list(sample(unique(Location), max(c(4, Max)))), .groups="drop")%>%
  rename(SubRegion2=SubRegion)
set.seed(NULL)


# For balanced dataset, reduce number of locations to that of the mean number of samples per region
# From the source dataset 

Data_effort%>%
  group_by(SubRegion)%>%
  summarise(Max=max(N), Mean=mean(N), .groups="drop")%>%
  pull(Mean)%>%
  mean()

# mean of about ~4 stations sampled per month, region, and year.
# Use this number to create balanced dataset


sim_data_balanced<-sim_data$sims%>%
  group_by(SubRegion)%>%
  filter(Location%in%filter(Shuffled_locations, SubRegion2==unique(SubRegion))$Locations[[1]][1:4])%>%
  ungroup()%>%
  arrange(Location, Date)%>%
  group_by(Location)%>%
  mutate(Lag=Date-lag(Date, order_by = Date))%>%
  ungroup()%>%
  mutate(Start=if_else(is.na(Lag) | Lag>3600*24*60, TRUE, FALSE),
         Series_ID=1:n(),
         Series_ID=if_else(Start, Series_ID, NA_integer_),
         Series_ID=as.integer(as.factor(Series_ID)))%>%
  fill(Series_ID, .direction="down")

# Create unbalanced version of simulated dataset by selecting the same number of "stations" for 
# each month, year, and region as were actually sampled in the source dataset

sim_data_unbalanced<-sim_data$sims%>%
  left_join(Data_effort, by=c("Month", "Year", "SubRegion"))%>%
  mutate(N=replace_na(N, 0))%>%
  group_by(Year, SubRegion, Month)%>%
  filter(Location%in%filter(Shuffled_locations, SubRegion2==unique(SubRegion))$Locations[[1]][1:unique(N)])%>%
  filter(N!=0)%>%
  ungroup()%>%
  arrange(Location, Date)%>%
  group_by(Location)%>%
  mutate(Lag=Date-lag(Date, order_by = Date))%>%
  ungroup()%>%
  mutate(Start=if_else(is.na(Lag) | Lag>3600*24*60, TRUE, FALSE),
         Series_ID=1:n(),
         Series_ID=if_else(Start, Series_ID, NA_integer_),
         Series_ID=as.integer(as.factor(Series_ID)))%>%
  fill(Series_ID, .direction="down")

year_sims<-sim_data$year_add%>%
  pivot_longer(cols=all_of(paste("sim", 1:10, sep="_")), names_to="Simulation", values_to="add")

ggplot(year_sims, aes(x=Year, y=add, color=Simulation))+
  geom_line()+
  facet_wrap(~Simulation)+
  theme_bw()


# Now do the same thing for the null model (no climate change sign --------

sim_data_NULL<-gamsim(modelld14a, year_slope=0)

## Balanced dataset

sim_data_NULL_balanced<-sim_data_NULL$sims%>%
  group_by(SubRegion)%>%
  filter(Location%in%filter(Shuffled_locations, SubRegion2==unique(SubRegion))$Locations[[1]][1:4])%>%
  ungroup()%>%
  arrange(Location, Date)%>%
  group_by(Location)%>%
  mutate(Lag=Date-lag(Date, order_by = Date))%>%
  ungroup()%>%
  mutate(Start=if_else(is.na(Lag) | Lag>3600*24*60, TRUE, FALSE),
         Series_ID=1:n(),
         Series_ID=if_else(Start, Series_ID, NA_integer_),
         Series_ID=as.integer(as.factor(Series_ID)))%>%
  fill(Series_ID, .direction="down")

## Unbalanced dataset
sim_data_NULL_unbalanced<-sim_data_NULL$sims%>%
  left_join(Data_effort, by=c("Month", "Year", "SubRegion"))%>%
  mutate(N=replace_na(N, 0))%>%
  group_by(Year, SubRegion, Month)%>%
  filter(Location%in%filter(Shuffled_locations, SubRegion2==unique(SubRegion))$Locations[[1]][1:unique(N)])%>%
  filter(N!=0)%>%
  ungroup()%>%
  arrange(Location, Date)%>%
  group_by(Location)%>%
  mutate(Lag=Date-lag(Date, order_by = Date))%>%
  ungroup()%>%
  mutate(Start=if_else(is.na(Lag) | Lag>3600*24*60, TRUE, FALSE),
         Series_ID=1:n(),
         Series_ID=if_else(Start, Series_ID, NA_integer_),
         Series_ID=as.integer(as.factor(Series_ID)))%>%
  fill(Series_ID, .direction="down")

year_sims_NULL<-sim_data_NULL$year_add%>%
  pivot_longer(cols=all_of(paste("sim", 1:10, sep="_")), names_to="Simulation", values_to="add")

ggplot(year_sims_NULL, aes(x=Year, y=add, color=Simulation))+
  geom_line()+
  facet_wrap(~Simulation)+
  theme_bw()

# Simulate for actual station locations

station_effort<-Data_CC4%>%
  select(Station, Month, Year)%>%
  mutate(Keep=TRUE)

sim_stations<-gamsim(modellea, data=newdata_stations)

sim_data_stations<-sim_stations$sims%>%
  left_join(station_effort, by=c("Station", "Month", "Year"))%>%
  mutate(Keep=replace_na(Keep, FALSE))%>%
  filter(Keep)%>%
  arrange(Station, Date)%>%
  group_by(Station)%>%
  mutate(Lag=Date-lag(Date, order_by = Date))%>%
  ungroup()%>%
  mutate(Start=if_else(is.na(Lag) | Lag>3600*24*60, TRUE, FALSE),
         Series_ID=1:n(),
         Series_ID=if_else(Start, Series_ID, NA_integer_),
         Series_ID=as.integer(as.factor(Series_ID)))%>%
  fill(Series_ID, .direction="down")

sim_stations_NULL<-gamsim(modelld14a, data=newdata_stations, year_slope=0, seeds=2^(1:10))

sim_data_stations_NULL<-sim_stations_NULL$sims%>%
  left_join(station_effort, by=c("Station", "Month", "Year"))%>%
  mutate(Keep=replace_na(Keep, FALSE))%>%
  filter(Keep)%>%
  arrange(Station, Date)%>%
  group_by(Station)%>%
  mutate(Lag=Date-lag(Date, order_by = Date))%>%
  ungroup()%>%
  mutate(Start=if_else(is.na(Lag) | Lag>3600*24*60, TRUE, FALSE),
         Series_ID=1:n(),
         Series_ID=if_else(Start, Series_ID, NA_integer_),
         Series_ID=as.integer(as.factor(Series_ID)))%>%
  fill(Series_ID, .direction="down")

# Modeling simulated data

test_noAR<-bam(sim_1 ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(10, 12)) + 
        te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(10, 12), by=Year_s), 
        family=scat, data = sim_data_balanced, method="fREML", discrete=T, nthreads=2)

r <- start_value_rho(test_noAR, plot=TRUE)

test_AR<-bam(sim_1 ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(10, 12)) + 
                 te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(10, 12), by=Year_s), 
               family=scat, rho=r, AR.start=Start, data = sim_data_balanced, method="fREML", discrete=T, nthreads=4)

test_noAR2<-bam(sim_2 ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(10, 12)) + 
                 te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(10, 12), by=Year_s), 
               family=scat, data = sim_data_balanced, method="fREML", discrete=T, nthreads=4)

r2 <- start_value_rho(test_noAR2, plot=TRUE)

test_AR2<-bam(sim_2 ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(10, 12)) + 
               te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(10, 12), by=Year_s), 
             family=scat, rho=r2, AR.start=Start, data = sim_data_balanced, method="fREML", discrete=T, nthreads=4)

test_noAR_NULL<-bam(sim_2 ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(10, 12)) + 
                  te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(10, 12), by=Year_s), 
                family=scat, data = sim_data_NULL_balanced, method="fREML", discrete=T, nthreads=4)

r_NULL <- start_value_rho(test_noAR_NULL, plot=TRUE)

test_AR_NULL<-bam(sim_2 ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(10, 12)) + 
                te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(10, 12), by=Year_s), 
              family=scat, rho=r_NULL, AR.start=Start, data = sim_data_NULL_balanced, method="fREML", discrete=T)

CC_newdata<-readRDS("Temperature analysis/Climate Change Prediction Data.Rds")%>%
  mutate(Date=as.Date(Julian_day, origin=as.Date(paste(Year, "01", "01", sep="-"))),
         Month=month(Date, label = T))
# For filtering the newdata after predictions

CC_pred<-predict(test_AR_NULL, newdata=CC_newdata2, type="terms", se.fit=TRUE, discrete=T, n.threads=4)

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

sim_tester<-function(data, sim){
  
  withCallingHandlers({
  
  model_noAR<-bam(data[[sim]] ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("cr", "cc"), k=c(25, 12)) + 
                   te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("cr", "cc"), k=c(25, 12), by=Year_s), 
                 family=scat, data = data, method="fREML", discrete=T, nthreads=4)
  
  r <- start_value_rho(model_noAR, plot=TRUE)
  
  model_AR<-bam(data[[sim]] ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("cr", "cc"), k=c(25, 12)) + 
                 te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("cr", "cc"), k=c(25, 12), by=Year_s), 
                family=scat, rho=r, AR.start=Start, data = data, method="fREML", discrete=T, nthreads=4)
  }, warning = function(w) {
    if (startsWith(conditionMessage(w), "step failure in theta estimation"))
      invokeRestart("muffleWarning")
  })
  
  CC_pred<-predict(model_AR, newdata=CC_newdata, type="terms", se.fit=TRUE, discrete=T, n.threads=4)
  
  newdata_CC_pred<-CC_newdata%>%
    mutate(Slope=CC_pred$fit[,"te(Julian_day_s,Latitude_s,Longitude_s):Year_s"],
           Slope_se=CC_pred$se.fit[,"te(Julian_day_s,Latitude_s,Longitude_s):Year_s"])%>%
    mutate(across(c(Slope, Slope_se), ~(.x/Year_s)/sd(1970:2020)))%>%
    mutate(Slope_se=abs(Slope_se))%>%
    mutate(Slope_l95=Slope-Slope_se*qnorm(0.9995),
           Slope_u95=Slope+Slope_se*qnorm(0.9995))%>%
    mutate(Sig=if_else(Slope_u95>0 & Slope_l95<0, "ns", "*"))
  out<-list()
  
  if(all(newdata_CC_pred$Sig=="ns")){
    message("No significant climate change signals detected")
  } else{
    newdata_CC_pred<-newdata_CC_pred%>%
      filter(Sig=="*")%>%
      st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=F)%>%
      st_transform(crs=st_crs(Delta))
    
    newdata_CC_pred_rast<-Rasterize_all(newdata_CC_pred, Slope)
    
    Slope_mean<-round(mean(newdata_CC_pred$Slope, na.rm=T),4)
    Slope_median<-round(median(newdata_CC_pred$Slope, na.rm=T),4)
    
    p<-ggplot()+
      geom_stars(data=newdata_CC_pred_rast)+
      facet_wrap(~month(Date, label=T))+
      scale_fill_viridis_c(limits=c(0, 0.03), guide=guide_colorbar(barheight=20), na.value="white")+
      ylab("Latitude")+
      xlab("Longitude")+
      ggtitle(paste(sim, "Mean:", Slope_mean, "Median:", Slope_median))+
      coord_sf()+
      theme_bw()+
      theme(strip.background=element_blank(), axis.text.x = element_text(angle=45, hjust=1), panel.grid=element_blank())
    out$plot<-p
  }
  out$slopes<-newdata_CC_pred$Slope
  return(out)
}
sims<-set_names(paste("sim", 1:10, sep="_"))

# Slope = 0.02
## Balanced
sim_results_balanced<-map(sims, ~sim_tester(sim_data_balanced, .x))
sim_results_balanced_extracted<-map_dbl(sims, ~median(sim_results_balanced[[.x]]$slopes))

## Unbalanced
sim_results_unbalanced<-map(sims, ~sim_tester(sim_data_unbalanced, .x))
sim_results_unbalanced_extracted<-map_dbl(sims, ~median(sim_results_unbalanced[[.x]]$slopes))

# Slope = 0
## Balanced
sim_results_balanced_NULL<-map(sims, ~sim_tester(sim_data_NULL_balanced, .x))
sim_results_balanced_NULL_extracted<-map_dbl(sims, ~median(sim_results_balanced_NULL[[.x]]$slopes))

## Unbalanced
sim_results_unbalanced_NULL<-map(sims, ~sim_tester(sim_data_NULL_unbalanced, .x))
sim_results_unbalanced_NULL_extracted<-map_dbl(sims, ~median(sim_results_unbalanced_NULL[[.x]]$slopes))

# Stations

## Slope = 0.02
sim_results_stationsb<-map(sims, ~sim_tester(sim_data_stations, .x))
sim_results_stations_extractedb<-map_dbl(sims, ~median(sim_results_stationsb[[.x]]$slopes))

## Slope = 0
sim_results_stations_NULLb<-map(sims, ~sim_tester(sim_data_stations_NULL, .x))
sim_results_stations_NULL_extractedb<-map_dbl(sims, ~median(sim_results_stations_NULLb[[.x]]$slopes))
