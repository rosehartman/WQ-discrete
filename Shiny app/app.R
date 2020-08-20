# TODO
# 1) Use geofacet
# 2) Facet by month and region

library(shiny)
library(shinyWidgets)
library(dplyr)
library(ggplot2)
require(tidyr)
library(sf)
library(stars)
library(lubridate)
library(shinythemes)
require(leaflet)
require(mapedit)
require(leafem)
require(leaflet.extras)

rastered_predsSE<-readRDS("Rasterized modellc4 predictions.Rds")
Dates<-st_get_dimension_values(rastered_predsSE, "Date")

all_points_static<-select(rastered_predsSE, Prediction)%>%
    aggregate(by=c(as.Date("1960-01-01"), Sys.Date()), function(x) length(which(!is.na(x))))%>%
    filter(time==as.Date("1960-01-01"))%>%
    mutate(across(Prediction, ~na_if(.x, 0)))%>%
    st_as_sf()%>%
    rename(N=`1960-01-01`)%>%
    mutate(ID=as.character(1:n()),)%>%
    st_transform(crs=4326)

pal_N_static<-colorNumeric("viridis", domain=range(all_points_static$N, na.rm=T), na.color="#00000000")

pal_N_rev_static<-colorNumeric("viridis", domain=range(all_points_static$N, na.rm=T), reverse=T, na.color="#00000000")

Delta_regions_static<-readRDS("Delta subregions.Rds")
Delta_regions_N<-all_points_static%>%
    st_transform(crs=26910)%>%
    st_centroid()%>%
    st_join(Delta_regions_static)%>%
    group_by(SubRegion)%>%
    summarise(N=sum(N), .groups="drop")%>%
    st_drop_geometry()

Delta_regions_static<-Delta_regions_static%>%
    left_join(Delta_regions_N, by="SubRegion")

pal_N2_static<-colorNumeric("viridis", domain=range(Delta_regions_static$N, na.rm=T), na.color="#00000000")

pal_N2_rev_static<-colorNumeric("viridis", domain=range(Delta_regions_static$N, na.rm=T), reverse=T, na.color="#00000000")

# Define UI for application that draws a histogram
ui <- fluidPage(
    theme = shinytheme("cerulean"),
    
    # Application title
    titlePanel("Delta Temperature"),
    
    sidebarLayout(
        # Sidebar with a slider input for number of bins 
        sidebarPanel(width=3,
                     h1("Filters"),
                     dateRangeInput("Date_range", label = "Date range", 
                                    start = min(Dates), end = max(Dates), startview = "year"),
                     pickerInput("Months", "Months:", choices=c("January" = 1, "February" = 2, "March" = 3, "April" = 4, 
                                                                "May" = 5, "June" = 6, "July" = 7, "August" = 8, 
                                                                "September" = 9, "October" = 10, "November" = 11, 
                                                                "December" = 12), 
                                 selected = 1:12, multiple = T, options=list(`actions-box`=TRUE, `selected-text-format` = "count > 3")),
                     h1("Plot options"),
                     radioGroupButtons("variable",
                                       "Plot:",
                                       choices = c("Predicted temperature"="Prediction", "Standard Error"="SE"), selected="Prediction", individual = TRUE, 
                                       checkIcon = list( yes = tags$i(class = "fa fa-circle", style = "color: steelblue"), no = tags$i(class = "fa fa-circle-o", style = "color: steelblue"))),
                     uiOutput("facet_options"),
                     conditionalPanel(condition="input.Tab=='Rasters' && input.Facet!='Year x Month' && input.Facet!='Month x Year'", 
                                      uiOutput("scale_options")),
                     conditionalPanel(condition="input.Tab=='Time series'", 
                                      h5("Select preset regions or draw your own?"),
                                      switchInput("Regions", value = TRUE , offLabel="Preset", onLabel="Draw", onStatus = "success", offStatus = "primary"),
                                      conditionalPanel(condition="input.Regions",
                                                       prettySwitch("Regions_all","Select all regions?", status = "success", fill = TRUE)))
        ),
        mainPanel(width=9,
                  tabsetPanel(type="tabs",
                              id="Tab",
                              tabPanel("Rasters",
                                       conditionalPanel(condition="input.Facet!='Year' && input.Facet!='Year x Month' && input.Facet!='Month x Year'",
                                                        uiOutput("select_Year")),
                                       conditionalPanel(condition="input.Facet!='Month' && input.Facet!='Year x Month' && input.Facet!='Month x Year'",
                                                        uiOutput("select_Month")),
                                       plotOutput("TempPlot")
                              ),
                              tabPanel("Time series",
                                       fluidRow(conditionalPanel(condition="input.Regions", editModUI("Pointplot", height = "80vh", width="100%")),
                                                conditionalPanel(condition="!input.Regions", editModUI("Regionplot", height = "80vh", width="100%"))),
                                       fluidRow(plotOutput("Time_plot"))
                              )
                  ))
    ),
    tags$head(tags$style("#TempPlot{height:80vh !important;}")),
    tags$head(tags$style("#Pointplot{height:80vh !important;}")),
    tags$head(tags$style("#Regionplot{height:80vh !important;}")),
    tags$head(tags$style("#Time_plot{height:80vh !important;}"))
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    TempData<-reactive({
        Data<-rastered_predsSE%>%
            filter(Date>min(input$Date_range))%>%
            filter(Date<max(input$Date_range))%>%
            filter(month(Date)%in%input$Months)
        return(Data)
    }
    )
    
    Data_dates<-reactive({
        req(TempData())
        st_get_dimension_values(TempData(), "Date")
    })
    
    output$select_Year <- renderUI({
        req(!is.null(input$Months))
        sliderInput("Year",
                    "Select year:", min=min(year(Data_dates())), max=max(year(Data_dates())), value =  min(c(2019, max(year(Data_dates())))), step=1, round=T, sep="", 
                    animate=animationOptions(interval=1000), width="100%")
    })
    
    output$select_Month <- renderUI({
        req(!is.null(input$Months))
        sliderInput("Month",
                    "Select month:", min=min(month(Data_dates())), max=max(month(Data_dates())), value =  min(month(Data_dates())), step=1, round=T, sep="", 
                    animate=animationOptions(interval=1000, loop=T), width="100%")
    })
    
    output$facet_options<-renderUI({
        if(input$Tab=='Rasters'){
            Choices<-c("None", "Month", "Year", "Month x Year", "Year x Month")
        }else{
            Choices<-c("None", "Month", "Region")
        }
        pickerInput("Facet",
                    "Facet plot by:",
                    choices = Choices, selected="None")
    })
    
    
    
    output$scale_options <- renderUI({
        if(is.null(input$Facet)){
            Facet<-"None"
        }else{
            Facet<-input$Facet
        }
        
        if(Facet=="None"){
            scale_choices<-c("No", "Yes: Across all months and years", "Yes: Across all years", "Yes: Across all months")
        } else{
            scale_choices<-c("No", "Yes: Across all months and years")
        }
        
        pickerInput("Scale",
                    "Fix scale:",
                    choices = scale_choices, selected="No")
    })
    
    Scale <- reactive({
        req(TempData(), input$Scale, input$Facet)
        if(input$Facet=="None"){
            if(input$Scale=="No"){
                out<-c(NA, NA)
            } else{
                if(input$Scale=="Yes: Across all months and years"){
                    out<-TempData()%>%
                        pull(all_of(input$variable))%>%
                        range(na.rm=T)
                } else{
                    if(input$Scale=="Yes: Across all years"){
                        out<-filter(TempData(), month(Date)==input$Month)%>%
                            pull(all_of(input$variable))%>%
                            range(na.rm=T)
                    } else{
                        out<-filter(TempData(), year(Date)==input$Year)%>%
                            pull(all_of(input$variable))%>%
                            range(na.rm=T)
                    }
                }
            }
            
        } else{
            if(input$Scale=="No"){
                out<-c(NA, NA)
            } else{
                out<-TempData()%>%
                    pull(Prediction)%>%
                    range(na.rm=T)
            }
            
        }
    })
    
    PlotData<-reactive({
        req(TempData(), input$Facet)
        Data<-TempData()%>%
            {if(input$Facet!='Month' & input$Facet!='Year x Month' & input$Facet!='Month x Year'){
                filter(., month(Date)==input$Month)
            } else{
                .
            }}%>%
            {if(input$Facet!='Year' & input$Facet!='Year x Month' & input$Facet!='Month x Year'){
                filter(., year(Date)==input$Year)
            } else{
                .
            }}
        return(Data)
    }
    )
    
    Plot<-reactive({
        req(input$variable)
        if(is.null(input$Months)){
            p<-ggplot()+
                geom_label(data=tibble(x=0.5, y=0.5, label="You must select at least one month."), aes(x=x, y=y, label=label))+
                theme_void()+
                theme(text=element_text(size=50))
            return(p)
        }
        if(input$variable=="Prediction"){
            Breaks<-seq(6,30,by=0.5)
            Labels<- function(x) ifelse(x==as.integer(x), as.character(x), "")
            
        }else{
            Breaks<-seq(0,5,by=0.05)
            Labels<- function(x) ifelse(near(x*10, as.integer(x*10)), as.character(x), "")
        }
        ggplot()+
            geom_stars(data=PlotData()%>%select(all_of(input$variable)))+
            scale_fill_viridis_c(limits=Scale(), expand=expansion(0,0), name=if_else(input$variable=="Prediction", "Temperature (°C)", "Standard error"), na.value="white", breaks=Breaks, labels= Labels,
                                 guide = guide_colorbar(direction="horizontal", title.position = "top", ticks.linewidth = 2,
                                                        title.hjust=0.5, label.position="bottom"))+
            {if(input$Facet=="None"){
                #ggtitle(par("din")[1])
                annotate("label", x=-122, y=38.45, label=paste(month(input$Month, label=T), input$Year), size=10)
            }}+
            {if(input$Facet=="Month"){
                annotate("label", x=-122, y=38.45, label=input$Year, size=10)
            }}+
            coord_sf()+
            {if(input$Facet=="Month"){
                facet_wrap(~month(Date, label=T))
            }}+
            {if(input$Facet=="Year"){
                facet_wrap(~year(Date))
            }}+
            {if(input$Facet=="Year x Month"){
                facet_grid(year(Date)~month(Date, label=T))
            }}+
            {if(input$Facet=="Month x Year"){
                facet_grid(month(Date, label=T)~year(Date))
            }}+
            ylab("Latitude")+
            xlab("Longitude")+
            theme_bw()+
            theme(axis.text.x = element_text(angle=45, hjust=1), strip.background=element_blank(),
                  panel.grid=element_blank(), #legend.background = element_rect(color="black"),
                  legend.position="top", legend.key.width = unit(150, "native"), legend.justification="center",
                  text=element_text(size=18), plot.title = element_text(hjust=0.5))
    })
    
    output$TempPlot <- renderPlot({
        Plot()
    })
    
    
    # Timeseries plots --------------------------------------------------------
    
    # Draw your own regions/areas of interest -------------------------------------------------
    
    all_points<-reactive({
        select(TempData(), all_of(input$variable))%>%
            aggregate(by=c(as.Date("1960-01-01"), Sys.Date()), function(x) length(which(!is.na(x))))%>%
            filter(time==as.Date("1960-01-01"))%>%
            mutate(across(all_of(input$variable), ~na_if(.x, 0)))%>%
            st_as_sf()%>%
            rename(N=`1960-01-01`)%>%
            mutate(ID=as.character(1:n()),)%>%
            st_transform(crs=4326)
        
    })
    
    all_points_plot<-reactive({
        all_points()%>%
            {if(is.null(Points())){
                mutate(., opac=0.2)
            }else{
                mutate(., opac=if_else(ID%in%Points()$ID, 0.9, 0.2))
            }}
    })
    
    pal_N<-reactive({
        colorNumeric("viridis", domain=range(all_points()$N, na.rm=T), na.color="#00000000")
    })
    
    pal_N_rev<-reactive({
        colorNumeric("viridis", domain=range(all_points()$N, na.rm=T), reverse=T, na.color="#00000000")
    })
    
    
    
    #set the namespace for the map
    ns <- shiny::NS("Pointplot")
    
    Point_plot<-leaflet()%>%
        addProviderTiles("Esri.WorldGrayCanvas")%>%
        addFeatures(data=all_points_static, fillColor=~pal_N_static(N), color="black", fillOpacity = 0.2, label=~N, layerId = ~ID, weight=0.4)%>%
        addLegend(data=all_points_static, position="topright", pal = pal_N_rev_static, values = ~N, opacity=0.5, 
                  labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)))
    
    #call the editMod function from 'mapedit' to use in the leaflet map.
    
    edits <- callModule(editMod, "Pointplot", leafmap = Point_plot, 
                        editorOptions=list(polylineOptions=F, circleMarkerOptions=F))
    
    observeEvent(all_points_plot(), {  
        proxy.points <- leafletProxy(ns("map"))
        
        proxy.points %>%
            clearShapes()%>%
            clearControls()%>%
            addFeatures(data=all_points_plot(), fillColor=~pal_N()(N), color="black", fillOpacity = ~opac, label=~N, layerId = ~ID, weight=~opac*2)%>%
            addLegend(data=all_points_plot(), position="topright", pal = pal_N_rev(), values = ~N, opacity=1, 
                      labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)))
    })
    
    
    
    
    # Use preset regions ------------------------------------------------------
    
    ns2 <- shiny::NS("Regionplot")
    
    Delta_regions<-reactive({
        Delta_regions_N<-all_points()%>%
            st_transform(crs=26910)%>%
            st_centroid()%>%
            st_join(Delta_regions_static%>%select(-N))%>%
            group_by(SubRegion)%>%
            summarise(N=sum(N), .groups="drop")%>%
            st_drop_geometry()
        
        Delta_regions<-Delta_regions_static%>%
            select(-N)%>%
            left_join(Delta_regions_N, by="SubRegion")
        
        return(Delta_regions)
    })
    
    Delta_regions_plot<-reactive({
        Delta_regions()%>%
            st_transform(crs=4326)%>%
            {if(is.null(Points())){
                mutate(., opac=0.2)
            }else{
                mutate(., opac=if_else(SubRegion%in%Points()$Region, 0.9, 0.2))
            }}
    })
    
    pal_N2<-reactive({
        colorNumeric("viridis", domain=range(Delta_regions()$N, na.rm=T), na.color="#00000000")
    })
    
    pal_N2_rev<-reactive({
        colorNumeric("viridis", domain=range(Delta_regions()$N, na.rm=T), reverse=T, na.color="#00000000")
    })
    
    Region_plot<-leaflet()%>%
        addProviderTiles("Esri.WorldGrayCanvas")%>%
        addFeatures(data=st_transform(Delta_regions_static, 4326), fillColor=~pal_N2_static(N), color="black", fillOpacity = 0.2, label=~SubRegion, layerId = ~SubRegion, weight=0.4)%>%
        addLegend(data=st_transform(Delta_regions_static, 4326), position="topright", pal = pal_N2_rev_static, values = ~N, opacity=0.5, 
                  labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)))
    
    #call the editMod function from 'mapedit' to use in the leaflet map.
    
    edits_regions <- callModule(editMod, "Regionplot", leafmap = Region_plot, 
                                editorOptions=list(polylineOptions=F, circleMarkerOptions=F))
    
    observeEvent(Delta_regions_plot(), {  
        proxy.points <- leafletProxy(ns2("map"))
        
        proxy.points %>%
            clearShapes()%>%
            clearControls()%>%
            addFeatures(data=Delta_regions_plot(), fillColor=~pal_N2()(N), color="black", fillOpacity = ~opac, label=~SubRegion, layerId = ~SubRegion, weight=~opac*2)%>%
            addLegend(data=Delta_regions_plot(), position="topright", pal = pal_N2_rev(), values = ~N, opacity=1, 
                      labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)))
    })
    

# Extract and plot timeseries data ----------------------------------------

    Points<-reactive({
        req(!(is.null(edits()) & is.null(edits_regions())))
        if(input$Regions){
            if(is.null(edits()$all)){
                return(NULL)
            }
            Area<-edits()$all%>%
                {if(input$Facet!="Region"){
                    st_union(.)
                } else{
                    .
                }}
            
            Points<-all_points()%>%
                st_filter(Area)%>%
                st_drop_geometry()%>%
                pull(ID)%>%
                unique()
            
            if(length(Points)>0){
                if(input$Facet=="Region"){
                    out<-filter(all_points(), ID%in%Points)%>%
                        st_join(Area)%>%
                        rename(Region=X_leaflet_id)
                } else{
                    out<-filter(all_points(), ID%in%Points)%>%
                        st_transform(crs=26910)%>%
                        st_centroid()%>%
                        st_transform(crs=4326)
                }
            }else(
                out<-NULL
            )
        }else{
            if(is.null(edits_regions()$all)){
                return(NULL)
            }
            
            Area<-edits_regions()$all%>%
                st_transform(crs=26910)%>%
                st_union()
            
            Regions<-Delta_regions()%>%
                st_filter(Area)%>%
                st_drop_geometry()%>%
                pull(SubRegion)%>%
                unique()
            
            if(length(Regions)>0){
                out<-filter(Delta_regions(), SubRegion%in%Regions)%>%
                    rename(Region=SubRegion)%>%
                    st_transform(crs=4326)
            }else(
                out<-NULL
            )
        }
        
        return(out)
    })  
    
    timeseries_data<-reactive({
        req(Points())
        if(input$Facet=="Region"){
            TempData()%>%
                select(all_of(input$variable))%>%
                aggregate(by=Points(), mean, na.rm=T)%>%
                st_as_sf(long=F)%>%
                st_drop_geometry()%>%
                mutate(Region=Points()$Region)%>%
                pivot_longer(cols=c(-Region), names_to="Date", values_to="out")%>%
                mutate(Date=parse_date_time(Date, "%Y-%m-%d"))
        } else{
            TempData()%>%
                select(all_of(input$variable))%>%
                {if(unique(st_geometry_type(Points()))=="POLYGON"){
                    aggregate(., by=Points(), mean, na.rm=T)%>%
                    st_as_sf(as_points=T, long=T)
                } else{
                    st_extract(., Points())%>%
                    st_as_sf(as_points=T, long=T)
                }}%>%
                as_tibble()%>%
                mutate(across(all_of(input$variable), ~if_else(is.nan(.x), NA_real_, .x)))%>%
                group_by(Date)%>%
                summarise(out=mean(.data[[input$variable]], na.rm=T), .groups="drop")%>%
                mutate(Month=month(Date),
                       Year=year(Date))
        }
        
    })
    
    time_series_plot<-reactive({
        req(timeseries_data())
        ggplot(timeseries_data(), aes(x=Date, y=out))+
            geom_line()+
            {if(input$Facet=="Month"){
                facet_wrap(~month(Date, label=T))
            }}+
            {if(input$Facet=="Region"){
                facet_wrap(~Region)
            }}+
            ylab(if_else(input$variable=="Prediction", "Temperature (°C)", "Standard error"))+
            theme_bw()+
            theme(strip.background=element_blank(), text=element_text(size=18))
        
    })
    
    output$Time_plot<-renderPlot({
        time_series_plot()
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
