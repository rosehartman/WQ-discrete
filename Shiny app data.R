require(dplyr)
require(sf)
require(purrr)
require(stars)
require(lubridate)
require(dtplyr)
require(tibble)
require(mgcv)

Data<-readRDS("Temperature smoothing model/Discrete Temp Data.Rds")



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
            by="SubRegion")%>%
  mutate(SubRegion=droplevels(SubRegion))
saveRDS(Delta_regions, file="Shiny app/Delta subregions.Rds")


# Rasterized model predictions --------------------------------------------

newdata_year<-readRDS("Temperature smoothing model/Prediction Data.Rds")
modellc4_predictions<-readRDS("Temperature smoothing model/modellc4_predictions.Rds")

newdata<-newdata_year%>%
  mutate(Prediction=modellc4_predictions$fit)%>%
  mutate(SE=modellc4_predictions$se.fit,
         L95=Prediction-SE*1.96,
         U95=Prediction+SE*1.96)%>%
  mutate(Date=as.Date(Julian_day, origin=as.Date(paste(Year, "01", "01", sep="-")))) # Create Date variable from Julian Day and Year

Rasterize_all <- function(data, var, out_crs=4326, n=100){
  var<-rlang::enquo(var)
  rlang::as_name(var)
  preds<-map(unique(data$Date), function(x) st_rasterize(data%>%
                                                           filter(Date==x)%>%
                                                           dplyr::select(!!var), 
                                                         template=st_as_stars(st_bbox(Delta), dx=diff(st_bbox(Delta)[c(1, 3)])/n, dy=diff(st_bbox(Delta)[c(2, 4)])/n, values = NA_real_))%>%
               st_warp(crs=out_crs))
  
  # Then bind all years together into 1 raster
  out <- exec(c, !!!preds, along=list(Date=unique(data$Date)))
  return(out)
}

newdata_rast <- newdata%>%
  mutate(Month=month(Date))%>%
  select(-N)%>%
  left_join(Data_effort, by=c("SubRegion", "Month", "Year"))%>% 
  mutate(across(c(Prediction, SE, L95, U95), ~if_else(is.na(N), NA_real_, .)))

# Create full rasterization of all predictions for interactive visualizations
rastered_preds<-Rasterize_all(newdata_rast, Prediction)
# Same for SE
rastered_SE<-Rasterize_all(newdata_rast, SE)
# Bind SE and predictions together
rastered_predsSE<-c(rastered_preds, rastered_SE)

#saveRDS(rastered_predsSE, file="Shiny app/Rasterized modellc4 predictions.Rds")


# Model evaluation data ---------------------------------------------------
# Stored as modellc4_residuals.Rds
modellc4_residuals<-readRDS("Temperature smoothing model/modellc4_residuals.Rds")

Data_resid<-Data%>%
  mutate(Residuals = modellc4_residuals)

Resid_sum<-Data_resid%>%
  lazy_dt()%>%
  group_by(Year, Month, SubRegion)%>%
  summarise(Resid=mean(Residuals), Resid_mag=mean(abs(Residuals)), SD=sd(Residuals))%>%
  ungroup()%>%
  as_tibble()

CV_fit_1<-readRDS("Temperature smoothing model/Group 1 CV predictions.Rds")
CV_fit_2<-readRDS("Temperature smoothing model/Group 2 CV predictions.Rds")
Data_split<-readRDS("Temperature smoothing model/Split data for cross validation.Rds")

CV_bind<-function(group, fold){
  if(group==1){
    fit<-CV_fit_1[[fold]]$fit
  }
  
  if(group==2){
    fit<-CV_fit_2[[fold]]$fit
  }
  
  Out<-Data_split%>%
    filter(Group==group & Fold==fold)%>%
    mutate(Fitted_CV=fit)
}

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


# Time correction ---------------------------------------------------------

modellc4<-readRDS("C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/Models/modellc4.Rds")
Noon<-((12*3600)-mean(Data$Time_num))/sd(Data$Time_num)

Times<-tibble(Time_num_s=c(seq(floor(min(Data$Time_num_s)*10)/10, ceiling(max(Data$Time_num_s)*10)/10, by=0.1), Noon),
              Year_fac="2000", Longitude_s=0, Latitude_s=0, Julian_day_s=0)
### Need to install earlier version of mgcv for this to work
## devtools::install_version("mgcv", version = "1.8-31", repos = "http://cran.us.r-project.org")

Time_correction<-predict(modellc4, newdata=Times, type="terms", terms="s(Time_num_s)")

Times<-Times%>%
  mutate(Correction=as.vector(Time_correction))

Times<-Times%>%
  filter(Time_num_s!=Noon)%>%
  mutate(Correction=filter(Times, Time_num_s==Noon)%>%pull(Correction)-Correction,
         Time=as.character(Time_num_s))%>%
  select(Time, Correction)

saveRDS(Times, file="Shiny app/Time_correction.Rds")
