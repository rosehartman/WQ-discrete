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
require(patchwork)
source("Utility_functions.R")
# Stop plotting scientific notation
options(scipen=999)

# Simulate data from a spatio-seasonal smoothing model to use for testing the model performance.
# Simulate with a balanced and imbalanced sampling design to examing that impact
# Also try simulating with all years vs. just years with similar average water temps to ensure model can differentiate real from spurious trends.

# First calculate year-to-year variability
# Start in 1975 when sampling became more regular across different regions
Data_CC4<-readRDS("Temperature analysis/Data_CC4.Rds")
Data<-readRDS("Temperature analysis/Discrete Temp Data.Rds")

Delta<-st_read("Delta Subregions")%>%
  select(SubRegion)%>%
  filter(SubRegion%in%unique(Data_CC4$SubRegion))

Data_effort<-Data_CC4%>%
  group_by(SubRegion, Month, Year)%>%
  summarise(N=n(), .groups="drop")

max_effort<-Data_effort%>%
  group_by(SubRegion)%>%
  summarise(Max=max(N), .groups="drop")

Slope_summary<-readRDS("Temperature analysis/Slope summary.Rds")%>%
  mutate(Month=as.integer(Month))

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
  mutate(Time_num_s=(12*60*60-mean(Data$Time_num))/sd(Data$Time_num),
         Year=2018,
         Year_s=(Year-mean(Data$Year))/sd(Data$Year),
         Year_fac=factor(Year))%>%
  left_join(Slope_summary%>%
              select(Month, SubRegion, Slope, Slope_se), 
            by=c("Month", "SubRegion"))%>%
  filter(!is.na(Slope))

# Fit model to use for spatial seasonal patterns
model_2018 <- bam(Temperature ~ te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs=c("cr", "cc"), k=c(30, 13)) + 
                    s(Time_num_s, bs="cr", k=5), family=scat,
                  data = filter(Data, Year==2018), method="fREML", discrete=T, nthreads=4)

# Define function to simulate gam data ------------------------------------

gamsim<-function(model, data=newdata_stations, null=FALSE, years=1970:2020, nsim=10, seeds=1:nsim){
  mu<-predict(model, newdata=data, type="response", discrete=T, n.threads=4)
  scale <- model[["sig2"]]
  data<-mutate(data, mu=mu)
  data<-map_dfr(years, ~mutate(data, Year=.x))%>%
    mutate(Date=as.Date(Julian_day, origin=as.Date(paste(Year, "01", "01", sep="-"))),
           Year_s=(Year-mean(years))/sd(years),
           Slope_year=Year-min(Year))
  
  # Extract distributional function
  family_fun<-fix.family.rd(family(model))[["rd"]]
  
  sim_fun<-function(Data=data, Null=null, Nsim=nsim, Scale=scale, seed){
    set.seed(seed)
    
    Data<-Data%>%
      mutate(add=rnorm(n=n(), mean=ifelse(rep(Null, n()), rep(0, n()), Slope_year*Slope), sd=Slope_se))%>%
      mutate(pred=mu+add)
    
    sims <- family_fun(mu=Data$pred, wt=rep(1, length(Data$pred)), scale=Scale)
    set.seed(NULL)
    return(tibble(sims=sims, add=Data$add, pred=Data$pred))
  }
  
  out<-map(1:nsim, ~sim_fun(seed=seeds[.x]))
  
  sims<-map(out, ~.x[["sims"]])
  adds<-map(out, ~.x[["add"]])
  preds<-map(out, ~.x[["pred"]])
  names(sims)<-paste("sim", 1:nsim, sep="_")
  names(adds)<-paste("add", 1:nsim, sep="_")
  names(preds)<-paste("pred", 1:nsim, sep="_")
  sims<-bind_cols(data, sims, adds, preds)
  
  return(sims)
  
}


# Simulate for actual station locations ------------------------------

## Unbalanced design reflecting actual sampling desing ---------------

### First calculate sampling effort for each station
station_effort<-Data_CC4%>%
  select(Station, Month, Year)%>%
  mutate(Keep=TRUE)



### Climate change signal
sim_data_stations<-gamsim(model_2018, null=FALSE)%>%
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

### Null model (no climate change)
sim_data_stations_NULL<-gamsim(model_2018, null=TRUE)%>%
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

## Balanced design ---------------------------------------------------

### For balanced dataset, reduce number of locations to that of the mean number of samples per region from the source dataset 

Data_effort%>%
  group_by(SubRegion)%>%
  summarise(Max=max(N), Mean=mean(N), .groups="drop")%>%
  pull(Mean)%>%
  mean()

### mean of about ~4 stations sampled per month, region, and year. Use this number to create balanced dataset

### Select 4 stations from each subregion. If Subregion has <4 unique station locations, just use all stations
set.seed(1)
Shuffled_stations<-newdata_stations%>%
  group_by(SubRegion)%>%
  mutate(N_stations=n_distinct(Station))%>%
  summarise(Stations=list(sample(unique(Station), min(c(4, N_stations)))), .groups="drop")%>%
  rename(SubRegion2=SubRegion)
set.seed(NULL)


### Climate change signal
sim_data_stations_balanced<-gamsim(model_2018, null=FALSE)%>%
  group_by(SubRegion)%>%
  filter(Station%in%filter(Shuffled_stations, SubRegion2==unique(SubRegion))$Stations[[1]][1:4])%>%
  ungroup()%>%
  arrange(Station, Date)%>%
  group_by(Station)%>%
  mutate(Lag=Date-lag(Date, order_by = Date))%>%
  ungroup()%>%
  mutate(Start=if_else(is.na(Lag) | Lag>3600*24*60, TRUE, FALSE),
         Series_ID=1:n(),
         Series_ID=if_else(Start, Series_ID, NA_integer_),
         Series_ID=as.integer(as.factor(Series_ID)))%>%
  fill(Series_ID, .direction="down")

### Null model (no climate change)
sim_data_stations_NULL_balanced<-gamsim(model_2018, null=TRUE)%>%
  group_by(SubRegion)%>%
  filter(Station%in%filter(Shuffled_stations, SubRegion2==unique(SubRegion))$Stations[[1]][1:4])%>%
  ungroup()%>%
  arrange(Station, Date)%>%
  group_by(Station)%>%
  mutate(Lag=Date-lag(Date, order_by = Date))%>%
  ungroup()%>%
  mutate(Start=if_else(is.na(Lag) | Lag>3600*24*60, TRUE, FALSE),
         Series_ID=1:n(),
         Series_ID=if_else(Start, Series_ID, NA_integer_),
         Series_ID=as.integer(as.factor(Series_ID)))%>%
  fill(Series_ID, .direction="down")

# Modeling simulated data -------------------------------

## First load prediction data for climate change data 
CC_newdata<-readRDS("Temperature analysis/Climate Change Prediction Data.Rds")%>%
  mutate(Date=as.Date(Julian_day, origin=as.Date(paste(Year, "01", "01", sep="-"))),
         Month=month(Date, label = T))

## Function to fit models on simulated data
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
    mutate(Slope_l99=Slope-Slope_se*qnorm(0.9995),
           Slope_u99=Slope+Slope_se*qnorm(0.9995))%>%
    mutate(Sig=if_else(Slope_u99>0 & Slope_l99<0, "ns", "*"))
  out<-list()
  
  if(all(newdata_CC_pred$Sig=="ns")){
    message("No significant climate change signals detected")
  } else{
    newdata_CC_pred2<-newdata_CC_pred%>%
      filter(Sig=="*")%>%
      st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=F)%>%
      st_transform(crs=st_crs(Delta))
    
    newdata_CC_pred_rast<-Rasterize_all(newdata_CC_pred2, Slope, region=Delta)
    
    Slope_mean<-round(mean(newdata_CC_pred2$Slope, na.rm=T),4)
    Slope_median<-round(median(newdata_CC_pred2$Slope, na.rm=T),4)
    
    p<-ggplot()+
      geom_stars(data=newdata_CC_pred_rast)+
      facet_wrap(~month(Date, label=T))+
      scale_fill_viridis_c(guide=guide_colorbar(barheight=20), na.value="white")+
      ylab("Latitude")+
      xlab("Longitude")+
      ggtitle(paste(sim, "Mean:", Slope_mean, "Median:", Slope_median))+
      coord_sf()+
      theme_bw()+
      theme(strip.background=element_blank(), axis.text.x = element_text(angle=45, hjust=1), panel.grid=element_blank())
    out$plot<-p
  }
  out$slopes<-select(newdata_CC_pred, Slope, Sig)
  return(out)
}

## Create function to process model results (if any slopes are significant)
sim_test_processer<-function(sim_results, data_type, model_type){
  
  Slope<-map(sim_results, ~.x[["slopes"]][["Slope"]])
  
  Sig<-map(sim_results, ~.x[["slopes"]][["Sig"]])
  Sig<-bind_cols(tibble(ID=1:length(Sig[[1]])), Sig)%>%
    pivot_longer(cols=all_of(sims), names_to="Simulation", values_to="Sig")
  
  title<-if_else(data_type=="unbalanced", "Unbalanced (realistic) sampling design", "Balanced sampling design")
  sim_results_stations_slopes<-bind_cols(CC_newdata, Slope)%>%
    mutate(ID=1:n())%>%
    left_join(Slope_summary%>%
                select(Month, SubRegion, Slope, Slope_se)%>%
                mutate(Month=month(Month, label=T)), 
              by=c("Month", "SubRegion"))%>%
    pivot_longer(cols=all_of(sims), names_to="Simulation", values_to="Slope_estimated")%>%
    left_join(Sig, by=c("ID", "Simulation"))%>%
    rename(Slope_true=Slope)%>%
    mutate(Simulation=factor(recode(Simulation, !!!set_names(paste("Simulation", 1:10), sims)), levels=paste("Simulation", 1:10)))
  
  p1<-ggplot(data=sim_results_stations_slopes)+
    {if(model_type=="NULL"){
      list(geom_point(aes(x=Slope_se, y=Slope_estimated, color=Sig), alpha=0.2),
        geom_hline(yintercept=0))
    }else{
      list(geom_point(aes(x=Slope_true, y=Slope_estimated, color=Sig), alpha=0.2),
        geom_abline(slope=1, intercept=0))
    }}+
    geom_blank(data=tibble(Sig=c("*", "ns")), aes(color=Sig))+
    scale_color_viridis_d(direction=-1, name="Significance",
                          guide=guide_legend(override.aes=list(alpha=1)))+
    facet_wrap(~Simulation, nrow=2)+
    ggtitle(title)+
    {if(model_type=="NULL"){
      xlab("'True' Slope (Temperature change per year [°C]) standard error")
    }else{
      xlab("'True' Slope (Temperature change per year [°C])")
    }}+
    ylab("Model-estimated slope (Temperature change per year [°C])")+
    theme_bw()+
    {if(model_type=="NULL"){
      theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=45, hjust=1))
    }else{
      theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5))
    }}
    
  
  if(model_type=="NULL"){
    return(p1)
  }
  
  Sig_sum<-sim_results_stations_slopes%>%
    group_by(Month, SubRegion, Simulation, Slope_se, Slope_true)%>%
    summarise(Sig=length(which(Sig=="*"))/n(), .groups="drop")
  
  p2<-ggplot()+
    geom_point(data=Sig_sum, aes(x=Slope_se, y=Slope_true, color=Sig))+
    scale_color_viridis_c(name="Proportion\nsignificant slopes", limits=c(0,1))+
    facet_wrap(~Simulation, nrow=2)+
    xlab("'True' Slope (Temperature change per year [°C]) standard error")+
    ylab("'True' slope (Temperature change per year [°C])")+
    ggtitle(title)+
    theme_bw()+
    theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=45, hjust=1))
  
  return(list(Estimation=p1, Significance=p2))
}

## Create list of all simulation scenarios
sims<-set_names(paste("sim", 1:10, sep="_"))

## Fit models to simulated data -------------------------------

### Unbalanced (realistic) sampling design -------------------

#### Climate change
sim_results_stations<-map(sims, ~sim_tester(sim_data_stations, .x))

p_sim_unbalanced_CC<-sim_test_processer(sim_results_stations, "unbalanced", "CC")

#### Null model
sim_results_stations_NULL<-map(sims, ~sim_tester(sim_data_stations_NULL, .x))

p_sim_unbalanced_NULL<-sim_test_processer(sim_results_stations_NULL, "unbalanced", "NULL")

### Unbalanced sampling design -------------------

#### Climate change
sim_results_stations_balanced<-map(sims, ~sim_tester(sim_data_stations_balanced, .x))

p_sim_balanced_CC<-sim_test_processer(sim_results_stations_balanced, "balanced", "CC")

## Slope = 0
sim_results_stations_NULL_balanced<-map(sims, ~sim_tester(sim_data_stations_NULL_balanced, .x))

p_sim_balanced_NULL<-sim_test_processer(sim_results_stations_NULL_balanced, "balanced", "NULL")

## Create final graphs ----------------------------------------
p_sim_CC_estimation<-p_sim_unbalanced_CC$Estimation/p_sim_balanced_CC$Estimation+plot_annotation(tag_levels="A")+plot_layout(guides="collect")

ggsave(p_sim_CC_estimation, file="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Temperature change/Figures/Figure S2 Slope_sim_CC.png",
       device="png", width=10, height=10, units="in")

p_sim_NULL_estimation<-p_sim_unbalanced_NULL/p_sim_balanced_NULL+plot_annotation(tag_levels="A")+plot_layout(guides="collect")

ggsave(p_sim_NULL_estimation, file="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Temperature change/Figures/Figure S4 Slope_sim_NULL.png",
       device="png", width=10, height=10, units="in")

p_sim_CC_significance<-p_sim_unbalanced_CC$Significance/p_sim_balanced_CC$Significance+plot_annotation(tag_levels="A")+plot_layout(guides="collect")

ggsave(p_sim_CC_significance, file="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Temperature change/Figures/Figure S3 Slope_sim_CC_Sig.png",
       device="png", width=10, height=10, units="in")
