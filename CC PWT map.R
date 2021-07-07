require(sf)
require(ggplot2)
require(dplyr)

SubRegions<-deltamapr::R_EDSM_Subregions_Mahardja%>%
  filter(!SubRegion%in%c("Upper Napa River", "Lower Napa River", "San Pablo Bay", "San Francisco Bay", "South Bay", "Upper Yolo Bypass"))

yolo<-sf::st_read("Yolo Bypass Extent")%>%
  st_transform(crs=st_crs(SubRegions))%>%
  st_union()

add<-tibble(Latitude=c(37.6, 37.6), Longitude=c(-121.7, -121.2))%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)%>%
  st_transform(crs=st_crs(SubRegions))

base<-deltamapr::WW_Watershed%>%
  st_transform(crs=st_crs(SubRegions))%>%
  st_crop(st_union(st_union(st_buffer(yolo, units::set_units(3000, "m")), SubRegions), add))%>%
  filter(!HNAME%in%c("PUTAH CREEK", "LAKE BERRYESSA", "SOUTH FORK PUTAH CREEK", "LAKE CURRY", "UPPER SAN LEANDRO RESERVOIR", "BETHANY RESERVOIR", "LAKE CHABOT", "SAN FRANCISCO BAY"))%>%
  st_union()

stations<-tibble(Station=c("Liberty Island", "Sacramento Deep Water Shipping Channel", "Decker Island", "Hood", "Jersey Point", "Burns Cut", "Mallard Island"),
                 Region=c(rep("North Delta", 2), rep("Sacramento River", 2), rep("San Joaquin River", 2), "Confluence"),
                 Latitude=c(38.24210, 38.25611, 38.09340, 38.36798, 38.05200, 37.98027, 38.04280),
                 Longitude=c(-121.68490, 	-121.66667, -121.73600, -121.51996, -121.68900, -121.38572, -121.92009))%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)%>%
  mutate(Region=factor(Region, levels=c("North Delta", "Sacramento River", "San Joaquin River", "Confluence")))

locations_points<-tibble(Location=c("Stockton", "Freeport", "Antioch", "Rio Vista", "", ""),
                  Latitude=c(37.956984, 38.461562, 38.004094, 38.155604, 37.800913, 37.796647),
                  Longitude=c(-121.290812, -121.499371, -121.805606, -121.691347, -121.620782, -121.585454))%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)

locations_text<-tibble(Location=c("Stockton", "Freeport", "Antioch", "Rio Vista", "Cache Slough", "Yolo Bypass", "Sacramento River", "Old River", "Middle\nRiver", "San Joaquin River"),
                         Latitude=c(37.956984, 38.481562, 37.99, 38.155604, 38.36, 38.5, 38.51, 37.88, 37.85, 37.65),
                         Longitude=c(-121.210812, -121.419371, -121.805606, -121.761347, -121.8, -121.7, -121.4, -121.65, -121.43, -121.4))%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)
  

p<-ggplot()+
  geom_sf(data=yolo, color=NA, fill="gray80", alpha=0.5)+
  geom_sf(data=base, fill="slategray3", color="slategray4")+
  geom_sf(data=locations_points)+
  geom_sf(data=stations, aes(fill=Region), shape=21, color="black", size=3)+
  geom_sf_text(data=locations_text, aes(label=Location), lineheight = 1)+
  scale_fill_manual(values=c("#FFB17A", "#886176", "#69995D", "#E34A6F"), guide="none")+
  theme_void()

ggsave("CC PWT map.png", plot=p, device="png", width=8, height=8, units = "in")
