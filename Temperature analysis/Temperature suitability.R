require(dplyr)
require(mgcv)
require(lubridate)
require(ggplot2)
require(sf)
require(tidyr)
require(geofacet)

mygrid <- data.frame(
  name = c("Upper Sacramento River Ship Channel", "Cache Slough and Lindsey Slough", "Lower Sacramento River Ship Channel", "Liberty Island", "Suisun Marsh", "Middle Sacramento River", "Lower Cache Slough", "Steamboat and Miner Slough", "Upper Mokelumne River", "Lower Mokelumne River", "Georgiana Slough", "Sacramento River near Ryde", "Sacramento River near Rio Vista", "Grizzly Bay", "West Suisun Bay", "Mid Suisun Bay", "Honker Bay", "Confluence", "Lower Sacramento River", "San Joaquin River at Twitchell Island", "San Joaquin River at Prisoners Pt", "Disappointment Slough", "Lower San Joaquin River", "Franks Tract", "Holland Cut", "San Joaquin River near Stockton", "Mildred Island", "Middle River", "Old River", "Upper San Joaquin River", "Grant Line Canal and Old River", "Victoria Canal"),
  row = c(2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 8, 8, 8),
  col = c(7, 4, 6, 5, 2, 8, 6, 7, 9, 9, 8, 7, 6, 2, 1, 2, 3, 4, 5, 6, 8, 9, 5, 6, 7, 9, 8, 8, 7, 9, 8, 7),
  code = c(" 1", " 1", " 2", " 3", " 8", " 4", " 5", " 6", " 7", " 9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "30", "29", "31"),
  stringsAsFactors = FALSE
)

Data<-readRDS("Temperature analysis/Discrete Temp Data.Rds")
newdata_all_sum<-readRDS("Temperature analysis/Model outputs and validations/newdata_all_sum.Rds")
Delta<-st_read("Delta Subregions")%>%
  filter(SubRegion%in%unique(Data$SubRegion))%>% # Remove regions outside our domain of interest
  dplyr::select(SubRegion)

Data_effort <- Data%>%
  st_drop_geometry()%>%
  group_by(SubRegion, Month, Year)%>%
  summarise(N=n(), .groups="drop")

newdata_all_sum_filt <- newdata_all_sum%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)%>%
  st_transform(crs=st_crs(Delta))%>%
  st_join(Delta)%>%
  st_drop_geometry()%>%
  left_join(Data_effort, by=c("SubRegion", "Month", "Year"))%>% 
  filter(!is.na(N))



suitability_data<-function(lower, higher){
  newdata_all_sum_filt%>%
  group_by(SubRegion, Month, Year)%>%
  summarise(N=n(), Good=length(which(Monthly_mean<lower))/N, 
            OK=length(which(Monthly_mean>=lower & Monthly_mean<higher))/N, 
            Bad=length(which(Monthly_mean>=higher))/N,
            .groups="drop")%>%
  complete(SubRegion, Month, Year)%>%
  pivot_longer(cols=c(Good, OK, Bad), names_to = "Suitability", values_to="Freq")%>%
  mutate(Suitability=factor(Suitability, levels=c("Bad", "OK", "Good")))%>%
  arrange(Suitability)%>%
  mutate(Suitability2=case_when(Suitability=="Bad" ~ paste0(Suitability, " (>=", higher, " °C)"),
                                Suitability=="OK" ~ paste0(Suitability, " (>=", lower, " & < ", higher, " °C)"),
                                Suitability=="Good" ~ paste0(Suitability, " (<", lower, " °C)")))%>%
  mutate(Suitability2=reorder(factor(Suitability2), as.integer(Suitability)))
}

suitability_plot<-function(data, month){
  ggplot(filter(data, Month==month), aes(x=Year, y=Freq, fill=Suitability2, color=Suitability2, order=Suitability2))+
  geom_col()+
  scale_fill_viridis_d(direction=-1, drop=F, na.translate = FALSE, 
                    name="Suitability", aesthetics=c("color", "fill"),
                    guide = guide_legend(direction="horizontal", title.position = "top",
                                         title.hjust=0.5, label.position="bottom", reverse=T))+
  facet_geo(~SubRegion, grid=mygrid, labeller=label_wrap_gen())+
  coord_cartesian(expand = FALSE)+
  ylab("Proportion of raster cells")+
  theme_bw()+
  theme(panel.grid=element_blank(), axis.text.x=element_text(angle=45, hjust=1),
        legend.position=c(0.2,0.2), legend.background=element_rect(color="black"))
}

# Delta Smelt
DS<-suitability_data(20,22)

suitability_plot(DS, 7)
suitability_plot(DS, 8)
suitability_plot(DS, 9)

# Chinook salmon
CS<-suitability_data(18,20)
suitability_plot(CS, 6)

newdata_all_sum_filt_av<-newdata_all_sum_filt%>%
  group_by(SubRegion, Month, Year)%>%
  summarise(N=n(), 
            Temperature=mean(Monthly_mean),
            .groups="drop")%>%
  complete(SubRegion, Month, Year)

av_plot<-function(month){
  ggplot(filter(newdata_all_sum_filt_av, Month==month), aes(x=Year, y=Temperature))+
    geom_line()+
    geom_point(aes(color=Temperature))+
    facet_geo(~SubRegion, grid=mygrid, labeller=label_wrap_gen(), scales="free_y")+
    scale_color_viridis_c()+
    coord_cartesian(expand = FALSE)+
    theme_bw()+
    theme(panel.grid=element_blank(), axis.text.x=element_text(angle=45, hjust=1),
          legend.position=c(0.2,0.2), legend.background=element_rect(color="black"))
}

# Delta Smelt
av_plot(8)

# Chinook salmon
av_plot(4)