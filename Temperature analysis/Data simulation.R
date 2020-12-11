require(dplyr)
require(tidyr)
require(mgcv)
require(sf)
require(purrr)
require(ggplot2)
# New method: Try simulating data (using the gratia simulate or predicted_sampled functions) from the global smoothing model to use for testing the model performance.
# Simulate with a balanced and imbalanced sampling design to examing that impact
# Also try simulating with all years vs. just years with similar average water temps to ensure model can differentiate real from spurious trends.

# First calculate year-to-year variability
# Start in 1975 when sampling became more regular across different regions
Data_CC4<-readRDS("Temperature analysis/Data_CC4.Rds")

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

newdata<-readRDS("Temperature analysis/Prediction Data.Rds")%>%
  filter(Year==2018)%>% # Pick year with lots of data (EDSM in effect) and almost all subregions represented (except Grant Line which ended in 1995)
  st_drop_geometry()%>%
  select(-Season, -N, -Time_num)%>%
  mutate(Month=as.integer(as.factor(Julian_day)))


# Define function to simulate gam data ------------------------------------

gamsim<-function(model, data=newdata, years=1970:2020, year_slope=0.02, year_sd=0.56, nsim=10, seeds=1:nsim){
  mu<-predict(model, newdata=data, type="response", discrete=T, n.threads=4)
  scale <- model[["sig2"]]
  data<-mutate(data, mu=mu)
  data<-map_dfr(years, ~mutate(data, Year=.x))
  
  # Extract distributional function
  family_fun<-fix.family.rd(family(model))[["rd"]]
  
  sim_fun<-function(Data=data, Years=years, Year_slope=year_slope, Year_sd=year_sd, Nsim=nsim, Scale=scale, seed){
    set.seed(seed)
    year_add<-tibble(Year=Years, add=rnorm(n=Year, mean=(Year-min(Year))*Year_slope, sd=Year_sd))
    
    Data<-Data%>%
      complete(Year=Years)%>%
      left_join(year_add, by="Year")%>%
      mutate(pred=mu+add)
    
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
  ungroup()

# Create unbalanced version of simulated dataset by selecting the same number of "stations" for 
# each month, year, and region as were actually sampled in the source dataset

sim_data_unbalanced<-sim_data$sims%>%
  left_join(Data_effort, by=c("Month", "Year", "SubRegion"))%>%
  mutate(N=replace_na(N, 0))%>%
  group_by(Year, SubRegion, Month)%>%
  filter(Location%in%filter(Shuffled_locations, SubRegion2==unique(SubRegion))$Locations[[1]][1:unique(N)])%>%
  filter(N!=0)%>%
  ungroup()

year_sims<-sim_data$year_add%>%
  pivot_longer(cols=all_of(paste("sim", 1:10, sep="_")), names_to="Simulation", values_to="add")

ggplot(year_sims, aes(x=Year, y=add, color=Simulation))+
  geom_line()+
  facet_wrap(~Simulation)+
  theme_bw()


# Now do the same thing for the null model (no climate change sign --------

sim_data_NULL<-gamsim(modelld14a, year_slope=0)

# Balanced dataset

sim_data_NULL_balanced<-sim_data_NULL$sims%>%
  group_by(SubRegion)%>%
  filter(Location%in%filter(Shuffled_locations, SubRegion2==unique(SubRegion))$Locations[[1]][1:4])%>%
  ungroup()

# Unbalanced dataset
sim_data_NULL_unbalanced<-sim_data_NULL$sims%>%
  left_join(Data_effort, by=c("Month", "Year", "SubRegion"))%>%
  mutate(N=replace_na(N, 0))%>%
  group_by(Year, SubRegion, Month)%>%
  filter(Location%in%filter(Shuffled_locations, SubRegion2==unique(SubRegion))$Locations[[1]][1:unique(N)])%>%
  filter(N!=0)%>%
  ungroup()

year_sims_NULL<-sim_data_NULL$year_add%>%
  pivot_longer(cols=all_of(paste("sim", 1:10, sep="_")), names_to="Simulation", values_to="add")

ggplot(year_sims_NULL, aes(x=Year, y=add, color=Simulation))+
  geom_line()+
  facet_wrap(~Simulation)+
  theme_bw()
