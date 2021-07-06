require(dplyr)
require(sf)
require(purrr)
require(stars)
require(lubridate)
require(dtplyr)
require(tibble)
require(mgcv)
source("Utility_functions.R")

Data<-readRDS("Temperature analysis/Discrete Temp Data.Rds")

#Region data for shiny app

Delta<-st_read("Delta Subregions")%>%
  select(SubRegion)%>%
  filter(SubRegion%in%unique(Data$SubRegion))

Data_effort <- Data%>%
  st_drop_geometry()%>%
  group_by(SubRegion, Month, Year)%>%
  summarise(N=n(), .groups="drop")


# Raw data ----------------------------------------------------------------

Data_shiny<-Data%>%
  select(Source, Station, Date, Datetime, Temperature, Year, StationID, Latitude, Longitude, Month, Season, SubRegion, Julian_day, Time_num_s)
saveRDS(Data_shiny, file="Shiny app/Raw temp data.Rds")

# Delta Regions -----------------------------------------------------------

Delta_regions<-Delta%>%
  left_join(Data%>%
              st_drop_geometry()%>%
              group_by(SubRegion)%>%
              summarise(N_data=n(), .groups="drop"),
            by="SubRegion")
saveRDS(Delta_regions, file="Shiny app/Delta subregions.Rds")


# Rasterized model predictions --------------------------------------------

newdata_year<-readRDS("Temperature analysis/Prediction Data.Rds")
modellf_predictions<-readRDS("Temperature analysis/Model outputs and validations/modellf_predictions.Rds")

newdata<-newdata_year%>%
  mutate(Prediction=modellf_predictions$fit)%>%
  mutate(SE=modellf_predictions$se.fit,
         L95=Prediction-SE*1.96,
         U95=Prediction+SE*1.96)%>%
  mutate(Date=as.Date(Julian_day, origin=as.Date(paste(Year, "01", "01", sep="-")))) # Create Date variable from Julian Day and Year

newdata_rast <- newdata%>%
  mutate(Month=month(Date))%>%
  select(-N)%>%
  left_join(Data_effort, by=c("SubRegion", "Month", "Year"))%>% 
  mutate(across(c(Prediction, SE, L95, U95), ~if_else(is.na(N), NA_real_, .)))

# Create full rasterization of all predictions for interactive visualizations
rastered_preds<-Rasterize_all(newdata_rast, Prediction, region=Delta)
# Same for SE
rastered_SE<-Rasterize_all(newdata_rast, SE, region=Delta)
# Bind SE and predictions together
rastered_predsSE<-c(rastered_preds, rastered_SE)

saveRDS(rastered_predsSE, file="Shiny app/Rasterized modellf predictions.Rds")


# Model evaluation data ---------------------------------------------------
# Stored as modelld_residuals.Rds
modellf_residuals<-readRDS("Temperature analysis/Model outputs and validations/modellf_residuals.Rds")

Data_resid<-Data%>%
  mutate(Residuals = modellf_residuals)

Resid_sum<-Data_resid%>%
  lazy_dt()%>%
  group_by(Year, Month, SubRegion)%>%
  summarise(Resid=mean(Residuals), Resid_mag=mean(abs(Residuals)), SD=sd(Residuals))%>%
  ungroup()%>%
  as_tibble()

CV_fit_1<-readRDS("Temperature analysis/Model outputs and validations/Group 1 CV predictions f.Rds")
CV_fit_2<-readRDS("Temperature analysis/Model outputs and validations/Group 2 CV predictions f.Rds")
Data_split<-readRDS("Temperature analysis/Split data for cross validation.Rds")


Data_split_CV<-map2_dfr(rep(c(1,2), each=10), rep(1:10,2), ~CV_bind(group=.x, fold=.y))%>%
  mutate(Resid_CV=Fitted_CV-Temperature,
         Fitted_resid=Fitted_CV-Fitted)

# EMP (first year-round survey) started in 1974 so restricting analysis to those years
Resid_CV_sum<-Data_split_CV%>%
  filter(Year>=1974)%>%
  lazy_dt()%>%
  group_by(Year, Month, SubRegion)%>%
  summarise(SD_CV=sd(Resid_CV), Resid_mag_CV=mean(abs(Resid_CV)), Resid_CV=mean(Resid_CV))%>%
  ungroup()%>%
  as_tibble()

Model_eval<-Resid_sum%>%
  left_join(Resid_CV_sum,
            by=c("Year", "Month", "SubRegion"))%>%
  left_join(Data_effort,
            by=c("Year", "Month", "SubRegion"))

saveRDS(Model_eval, file="Shiny app/Model_eval.Rds")

CV_sum<-Data_split_CV%>%
  st_drop_geometry()%>%
  group_by(Group, Fold)%>%
  summarise(RMSE=sqrt(mean(Resid_CV^2)), 
            r=cor(Fitted_CV, Temperature, method="pearson"), .groups="drop")

saveRDS(CV_sum, file="Shiny app/CV_sum.Rds")


# Time correction ---------------------------------------------------------

modellf<-readRDS("C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/Models/modellf.Rds")
Noon<-((12*3600)-mean(Data$Time_num))/sd(Data$Time_num)

Times<-expand_grid(Time_num_s=c(seq(floor(min(Data$Time_num_s)*10)/10, ceiling(max(Data$Time_num_s)*10)/10, by=0.1), Noon),
              Year_fac="2000", 
              Longitude_s=0, 
              Latitude_s=0, 
              Julian_day=yday(ymd(paste("2001", 1:12, "15", sep="-"))))%>%
  mutate(Julian_day_s=(Julian_day-mean(Data$Julian_day))/sd(Data$Julian_day))
### Need to install earlier version of mgcv for this to work
## devtools::install_version("mgcv", version = "1.8-31", repos = "http://cran.us.r-project.org")

Time_correction<-predict(modellf, newdata=Times, type="terms", terms="te(Time_num_s,Julian_day_s)")

Times<-Times%>%
  mutate(Correction=as.vector(Time_correction))%>%
  mutate(Month=as.integer(as.factor(Julian_day)))%>%
  group_by(Month)%>%
  mutate(Noon=Correction[which(Time_num_s==Noon)])%>%
  ungroup()%>%
  mutate(Correction=Noon-Correction,
         Time=as.character(Time_num_s))%>%
  filter(Time_num_s!=Noon)%>%
  select(Time, Month, Correction)

saveRDS(Times, file="Shiny app/Time_correction.Rds")
