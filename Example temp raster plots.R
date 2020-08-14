library(shiny)
library(shinyWidgets)
library(dplyr)
library(ggplot2)
library(sf)
library(stars)
library(lubridate)
load("Shiny app/Rasterized modellc4 predictions.Rds")

TempData<-filter(rastered_predsSE, year(Date)<=2018)%>%
  filter(year(Date)>=2010)%>%
  filter(month(Date)==7)

p<-ggplot()+
  geom_stars(data=TempData%>%select(Prediction))+
  scale_fill_viridis_c(expand=expansion(0,0), name="Temperature (°C)", na.value="white", breaks=seq(6,30,by=0.5), labels= function(x) ifelse(x==as.integer(x), as.character(x), ""),
                       guide = guide_colorbar(direction="horizontal", title.position = "top", ticks.linewidth = 2,
                                              title.hjust=0.5, label.position="bottom", barheight = 0.5))+
  coord_sf()+
  ylab("Latitude")+
  xlab("Longitude")+
  facet_wrap(~year(Date))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1), strip.background=element_blank(),
        panel.grid=element_blank(), #legend.background = element_rect(color="black"),
        legend.position=c(0.5, 1.15), legend.key.width = unit(150, "native"), legend.justification="center", plot.margin = margin(50,0,0,0))
ggsave(p, file="~/Temperature raster.png", device="png", units="in", width=5, height=5)