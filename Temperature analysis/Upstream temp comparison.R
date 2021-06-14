require(dplyr)
require(ggplot2)
require(readr)
require(lubridate)
require(tidyr)

download.file("https://portal.edirepository.org/nis/dataviewer?packageid=edi.591.2&entityid=fb147771773b8354667a0b43e3f8e0e4", file.path(tempdir(), "Temp_filtered.csv"))

data<-read_csv(file.path(tempdir(), "Temp_filtered.csv"))

Data_sum<-data%>%
  filter(Station%in%c("FPT", "RCS", "SRV", "ANH", "PTS", "SWE", "GES"))%>%
  group_by(Station, Date)%>%
  summarise(Temp=mean(Temp), .groups="drop")%>%
  mutate(JDay=yday(Date),
         Year=year(Date),
         Position=recode(Station, RCS=4, FPT=3, SRV=2, ANH=1, PTS=0, SWE=NA_real_, GES=NA_real_))

ggplot(filter(Data_sum, Year>=2012 & Year<=2018 & !is.na(Position)), aes(x=JDay, y=Temp, color=Position, group=Position))+
  geom_line()+
  facet_wrap(~Year)+
  scale_color_viridis_c(breaks=4:0, labels=c("Knights landing", "Freeport", "Rio Vista", "Antioch", "Pittsburg"))+
  scale_x_continuous(breaks=yday(parse_date_time(paste("2001", 1:12, 1, sep="-"), orders="%Y-%m-%d")), labels=month(1:12, label=T), 
                     minor_breaks = F, expand=expansion(0,0))+
  theme_bw()

# Freeport vs Rio Vista

ggplot(filter(Data_sum, Station%in%c("FPT", "SRV") & Year>=2009 & Year<=2018), aes(x=JDay, y=Temp, color=Station, group=Station))+
  geom_line()+
  facet_wrap(~Year)+
  scale_color_viridis_d()+
  scale_x_continuous(breaks=yday(parse_date_time(paste("2001", 1:12, 1, sep="-"), orders="%Y-%m-%d")), labels=month(1:12, label=T), 
                     minor_breaks = F, expand=expansion(0,0))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1))

Data_sum_diff<-Data_sum%>%
  filter(Station%in%c("FPT", "SRV"))%>%
  select(-Position)%>%
  pivot_wider(names_from = Station, values_from=Temp)%>%
  filter(!is.na(FPT) & !is.na(SRV))%>%
  mutate(Diff=FPT-SRV,
         Month=month(Date, label=T))

ggplot(Data_sum_diff, aes(x=Month, y=Diff, fill=Year, group=interaction(Year, Month)))+
  geom_boxplot()+
  geom_hline(yintercept=0, color="red")+
  scale_fill_viridis_c()+
  ylab("Temperature difference FPT-SRV")+
  theme_bw()

# Freeport vs antioch

ggplot(filter(Data_sum, Station%in%c("FPT", "ANH") & Year>=2009 & Year<=2018), aes(x=JDay, y=Temp, color=Station, group=Station))+
  geom_line()+
  facet_wrap(~Year)+
  scale_color_viridis_d()+
  scale_x_continuous(breaks=yday(parse_date_time(paste("2001", 1:12, 1, sep="-"), orders="%Y-%m-%d")), labels=month(1:12, label=T), 
                     minor_breaks = F, expand=expansion(0,0))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1))

Data_sum_diff<-Data_sum%>%
  filter(Station%in%c("FPT", "ANH"))%>%
  select(-Position)%>%
  pivot_wider(names_from = Station, values_from=Temp)%>%
  filter(!is.na(FPT) & !is.na(ANH))%>%
  mutate(Diff=FPT-ANH,
         Month=month(Date, label=T))

ggplot(Data_sum_diff, aes(x=Month, y=Diff, fill=Year, group=interaction(Year, Month)))+
  geom_boxplot()+
  geom_hline(yintercept=0, color="red")+
  scale_fill_viridis_c()+
  ylab("Temperature difference FPT-ANH")+
  theme_bw()

# Freeport vs Walnut Grove
Data_diff_FPT_WG<-Data_sum%>%
  filter(Station%in%c("FPT", "SWE", "GES"))%>%
  select(-Position)%>%
  pivot_wider(names_from = Station, values_from=Temp)%>%
  mutate(GES=if_else(is.na(GES), SWE, GES))%>%
  select(-SWE)%>%
  rename(Walnut_grove=GES, Freeport=FPT)%>%
  #filter(!is.na(Freeport) & !is.na(Walnut_grove))%>%
  mutate(Diff=Freeport-Walnut_grove,
         Month=month(Date, label=T))
  
ggplot(Data_diff_FPT_WG, aes(x=JDay))+
  geom_line(aes(y=Freeport, color="Freeport"))+
  geom_line(aes(y=Walnut_grove, color="Walnut grove"))+
  facet_wrap(~Year)+
  scale_color_viridis_d()+
  scale_x_continuous(breaks=yday(parse_date_time(paste("2001", 1:12, 1, sep="-"), orders="%Y-%m-%d")), labels=month(1:12, label=T), 
                     minor_breaks = F, expand=expansion(0,0))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1))

ggplot(Data_diff_FPT_WG, aes(x=Month, y=Diff, fill=Year, group=interaction(Year, Month)))+
  geom_boxplot()+
  geom_hline(yintercept=0, color="red")+
  scale_fill_viridis_c()+
  ylab("Temperature difference Freeport-Walnut grove")+
  theme_bw()
