require(dplyr)
require(sf)
require(ggplot2)
require(mgcv)
require(itsadug)
require(lubridate)
require(patchwork)
require(geofacet)

mygrid <- data.frame(
  name = c("Upper Sacramento River Ship Channel", "Cache Slough and Lindsey Slough", "Lower Sacramento River Ship Channel", "Liberty Island", "Suisun Marsh", "Middle Sacramento River", "Lower Cache Slough", "Steamboat and Miner Slough", "Upper Mokelumne River", "Lower Mokelumne River", "Georgiana Slough", "Sacramento River near Ryde", "Sacramento River near Rio Vista", "Grizzly Bay", "West Suisun Bay", "Mid Suisun Bay", "Honker Bay", "Confluence", "Lower Sacramento River", "San Joaquin River at Twitchell Island", "San Joaquin River at Prisoners Pt", "Disappointment Slough", "Lower San Joaquin River", "Franks Tract", "Holland Cut", "San Joaquin River near Stockton", "Mildred Island", "Middle River", "Old River", "Upper San Joaquin River", "Grant Line Canal and Old River", "Victoria Canal"),
  row = c(2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 8, 8, 8),
  col = c(7, 4, 6, 5, 2, 8, 6, 7, 9, 9, 8, 7, 6, 2, 1, 2, 3, 4, 5, 6, 8, 9, 5, 6, 7, 9, 8, 8, 7, 9, 8, 7),
  code = c(" 1", " 1", " 2", " 3", " 8", " 4", " 5", " 6", " 7", " 9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "30", "29", "31"),
  stringsAsFactors = FALSE
)

Data_analysis<-readRDS("Salinity analysis/Discrete Salinity Analysis Data.Rds")

Data_effort<-Data_analysis%>%
  st_drop_geometry()%>%
  group_by(Year, SubRegion, Month)%>%
  summarise(N=n(), .groups="drop")

p_effort<-ggplot(Data_effort)+
  geom_tile(aes(x=Year, y=reorder(month(Month, label=T), desc(month(Month, label=T))), fill=N))+
  scale_fill_viridis_c(breaks=seq(0,140, by=10),
                       guide=guide_colorbar(barheight=15, barwidth = 3))+
  scale_x_continuous(breaks=seq(1970, 2020, by=10))+
  scale_y_discrete(breaks=c("Jan", "Mar", "May", "Jul", "Sep", "Nov"))+
  facet_geo(~SubRegion, grid=mygrid, labeller=label_wrap_gen(width=18))+
  ylab("Month")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1), panel.grid=element_blank(), text=element_text(size=16), legend.position=c(0.4, 0.65), 
        legend.background = element_rect(color="black"), panel.background = element_rect(color="black"), legend.margin=margin(10,10,15,10))

ggsave(p_effort, file="C:/Users/sbashevkin/deltacouncil/Science Extranet - Discrete water quality synthesis/Salinity change/Figures/Sampling effort.png",
       device="png", width=15, height=18, units="in")

# First fit model with AR term to find optimal AR rho parameter

SC_gam_NOAR <- bam(Salinity ~ te(Latitude_s, Longitude_s, Julian_day_s, Tide, d=c(2,1, 1), bs=c("tp", "cc", "fs"), k=c(25, 13, 5)) + 
                          te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=Year_s),
                          family=scat, data = Data_analysis, method="fREML", discrete=T, nthreads=4)

SC_gam_NOAR2 <- bam(Salinity ~ te(Latitude_s, Longitude_s, Julian_day_s, Tide, d=c(2,1, 1), bs=c("tp", "cc", "fs"), k=c(25, 13, 5)) + 
                     te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=Year_s),
                   family=gaussian, data = Data_analysis, method="fREML", discrete=T, nthreads=4)

SC_gam_NOAR3 <- bam(log(Salinity) ~ te(Latitude_s, Longitude_s, Julian_day_s, Tide, d=c(2,1, 1), bs=c("tp", "cc", "fs"), k=c(25, 13, 5)) + 
                      te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=Year_s),
                    family=gaussian, data = Data_analysis, method="fREML", discrete=T, nthreads=4)

SC_gam_NOAR4 <- bam(log(Salinity) ~ te(Latitude_s, Longitude_s, Julian_day_s, Tide, d=c(2,1, 1), bs=c("tp", "cc", "fs"), k=c(25, 13, 5)) + 
                      te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=Year_s),
                    family=scat, data = Data_analysis, method="fREML", discrete=T, nthreads=4)

SC_gam_NOAR5 <- bam(log(Salinity) ~ te(Latitude_s, Longitude_s, Julian_day_s, Tide, d=c(2,1, 1), bs=c("tp", "cc", "fs"), k=c(25, 13, 5)) + 
                      te(Latitude_s, Longitude_s, Julian_day_s, Year_s, d=c(2,1,1), bs=c("tp", "cc", "tp"), k=c(25, 13, 5)),
                    family=scat, data = Data_analysis, method="fREML", discrete=T, nthreads=4)
r <- start_value_rho(SC_gam_NOAR, plot=TRUE)

#########Best Model####################

CC_gam8d7b_AR7 <- bam(Temperature ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                        te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=Year_s) + 
                        s(Time_num_s, k=5), family=scat, rho=r6, AR.start=Start, data = Data_CC4.3, method="fREML", discrete=T, nthreads=4)