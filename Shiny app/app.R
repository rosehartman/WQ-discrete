#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyWidgets)
library(dplyr)
library(ggplot2)
library(sf)
library(stars)
library(lubridate)
library(shinythemes)
require(ggiraph)
rastered_predsSE<-readRDS("Rasterized modellc4 predictions.Rds")
Dates<-st_get_dimension_values(rastered_predsSE, "Date")

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
                     conditionalPanel(condition="input.Tab=='Rasters'", 
                                      pickerInput("Facet",
                                                  "Facet plot by:",
                                                  choices = c("None", "Month", "Year", "Month x Year", "Year x Month"), selected="None")),
                     conditionalPanel(condition="input.Tab=='Rasters' && input.Facet!='Year x Month' && input.Facet!='Month x Year'", 
                                      uiOutput("scale_options"))
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
                                       girafeOutput("Pointplot", height="100vh", width="120vh")
                              )
                  ))
    ),
    tags$head(tags$style("#TempPlot{height:80vh !important;}")),
    tags$head(tags$style("#Pointplot{height:80vh !important;}"))
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

    Point_plot<-reactive({
        Data<-select(TempData(), Prediction)%>%
            aggregate(by=c(as.Date("1960-01-01"), Sys.Date()), function(x) length(which(!is.na(x))))%>%
            filter(time==as.Date("1960-01-01"))%>%
            mutate(Prediction=na_if(Prediction, 0))%>%
            as_tibble()%>%
            select(Longitude=x, Latitude=y, N=Prediction)%>%
            mutate(ID=as.character(1:n()))%>%
            filter(!is.na(N))
        p<-ggplot()+
            geom_sf(data=spacetools::Delta)+
            geom_point_interactive(data=Data, aes(x=Longitude, y=Latitude, color=N, tooltip=N, data_id=ID), alpha=0.5, size=0.5)+
            scale_color_viridis_c(name="N", na.value="white")+
            ylab("Latitude")+
            xlab("Longitude")+
            theme_bw()+
            theme(panel.grid=element_blank(), text=element_text(size=18))
        return(p)
    })
    
    output$Pointplot <- renderGirafe({
        girafe(ggobj=Point_plot())
    })
    
Points<-reactive({
    
})    
    
    timeseries_data<-reactive({
        TempData()%>%
            select(all_of(input$variable))%>%
            st_extract(Points())%>%
            st_as_sf(as_points=T, long=T)%>%
            as_tibble()
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
