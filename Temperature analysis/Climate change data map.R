require(dplyr)
require(sf)
require(ggplot2)
require(maps)
require(ggspatial)

Data_CC4<-readRDS("Temperature analysis/Data_CC4.Rds")

SubRegions<-deltamapr::R_EDSM_Subregions_Mahardja%>%
  filter(SubRegion%in%unique(Data_CC4$SubRegion))

base<-deltamapr::WW_Delta%>%
  st_transform(crs=st_crs(SubRegions))%>%
  st_crop(SubRegions)

Data<-Data_CC4%>%
  group_by(Station, Latitude, Longitude)%>%
  summarise(N=n(), .groups="drop")%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)%>%
  st_transform(crs=st_crs(SubRegions))


states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))%>%
  st_transform(crs=st_crs(SubRegions))
california<-filter(states, ID=="california")

base2<-base%>%
  st_make_valid()%>%
  st_crop(Data)

station_lims<-st_bbox(Data)

labels<-tibble(label=c("Suisun Bay", "Suisun Marsh", "Confluence", "Cache Slough", "Sacramento River", 
                       "Sacramento\nShipping Channel", "San Joaquin River", "Cosumnes River", "Mokelumne\nRiver"), 
               Y=c(4212000, 4226700, 4211490, 4232164, 4262276, 4262276, 4183000,4247000,4225000), 
               X=c(583318, 590000, 597000, 615970, 625568, 623600,649500,645000,648500),
               label_Y=c(4208500, 4240000, 4200000, 4228686, 4262058, 4262058, 4180000, 4255000,4220000), 
               label_X=c(585000, 590000, 610000, 607572, 640000, 610743, 642000, 642000, 647000))

plot(select(base, geometry),reset=F, col="slategray1", border="slategray2")
plot(select(SubRegions, geometry), add=T, lwd=2)
points(as_tibble(st_coordinates(Data)), col="black", pch=16)
Letter_locs<-locator()

Letters<-tibble(Label=c(letters, paste0("a", letters)[1:6]), 
                SubRegion=c("Upper Sacramento River Ship Channel", "Middle Sacramento River",
                            "Lower Sacramento River Ship Channel", "Steamboat and Miner Slough",
                            "Cache Slough and Lindsey Slough", "Liberty Island",
                            "Lower Cache Slough", "Sacramento River near Ryde",
                            "Georgiana Slough", "Upper Mokelumne River",
                            "Suisun Marsh", "West Suisun Bay",
                            "Grizzly Bay", "Mid Suisun Bay",
                            "Honker Bay", "Confluence",
                            "Lower Sacramento River", "Sacramento River near Rio Vista",
                            "San Joaquin River at Twitchell Island", "San Joaquin River at Prisoners Pt",
                            "Lower Mokelumne River", "Disappointment Slough",
                            "Lower San Joaquin River", "Franks Tract", 
                            "Holland Cut", "Mildred Island",
                            "San Joaquin River near Stockton", "Old River",
                            "Middle River", "Grant Line Canal and Old River",
                            "Victoria Canal", "Upper San Joaquin River"),
                X=c(627233.1, 629416.7, 620376.8, 622167.3, 606969.7, 613040.0, 616315.3, 
                    618236.9, 624700.2, 633347.1, 594610.7, 576661.9, 586793.6, 585920.2, 
                    594916.4, 600593.7, 604960.8, 613258.3, 618586.2, 627407.8, 628718.0, 
                    639854.1, 614393.8, 622254.6, 625398.9, 631949.6, 640334.5, 626621.7, 
                    632735.7, 640683.9, 623914.1, 647234.5),
                Y=c(4271603, 4257104, 4249768, 4248719, 4246711, 4245837, 4228543, 
                    4223128, 4224002, 4232212, 4232648, 4214307, 4219984, 4211599, 
                    4214831, 4214569, 4215180, 4215180, 4218848, 4214743, 4222167, 
                    4219547, 4210289, 4215005, 4212298, 4204961, 4204350, 4203127, 
                    4202253, 4200244, 4193432, 4197886))
paste(Letters$Label, Letters$SubRegion, sep=" - ", collapse=", ")

pout<-ggplot(states)+
  geom_sf(color="slategray1", fill="gray70")+
  geom_sf(data=base2, color="slategray1", fill="slategray1")+
  geom_rect(xmin = station_lims["xmin"]-0.2, xmax = station_lims["xmax"]+0.2, ymin = station_lims["ymin"]-0.2, ymax = station_lims["ymax"]+0.2, 
            fill = NA, colour = "black", size = 1)+
  coord_sf(xlim=c(st_bbox(california)["xmin"], st_bbox(california)["xmax"]), ylim=c(st_bbox(california)["ymin"], st_bbox(california)["ymax"]))+
  theme_bw()+
  theme(panel.background = element_rect(fill = "slategray1"), axis.text.x=element_text(angle=45, hjust=1))
pout

p<-ggplot()+
  geom_sf(data=base, fill="slategray1", color="slategray2")+
  geom_sf(data=SubRegions, alpha=0.1)+
  geom_segment(data=labels, aes(x=label_X, y=label_Y, xend=X, yend=Y), size=1)+
  geom_label(data=labels, aes(label=label, x=label_X, y=label_Y))+
  geom_sf(data=Data, aes(color=N))+
  geom_segment(data=tibble(x=580000, y=4205000, xend=575000, yend=4205000), aes(x=x, y=y, xend=xend, yend=yend), arrow=arrow(length = unit(0.03, "npc")), size=1)+
  geom_label(data=tibble(x=590000, y=4205000, label="to San Francisco Bay"), aes(x=x, y=y, label=label))+
  geom_text(data=Letters, aes(x=X, y=Y, label=Label))+
  ylab("")+
  xlab("")+
  #coord_sf(datum=st_crs(SubRegions))+
  scale_color_viridis_c()+
  theme_bw()+
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "bl", pad_y=unit(0.05, "npc"), which_north = "true")+
  annotation_custom(
    grob = ggplotGrob(pout),
    xmin = -Inf,
    xmax = 599971,
    ymin = 4242767,
    ymax = Inf
  )

ggsave("C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/Manuscripts/Climate change/Figures/Figure 1 map.png", plot=p, device="png", width=8, height=8, units = "in")
