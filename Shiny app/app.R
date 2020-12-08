# Add raw data with same features as time series plots
require(cubelyr)
require(stringr)
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
require(geofacet)
require(ggiraph)
require(scales)

rastered_predsSE<-readRDS("Rasterized modelld predictions.Rds")
Temp_min<-floor(min(rastered_predsSE$Prediction, na.rm=T)*10)/10
Temp_max<-floor(max(rastered_predsSE$Prediction, na.rm=T)*10)/10

# To fix error with different versions of Proj on shinyapps.io vs. the saved object and since the function `st_crs<-` is missing for stars objects
d <- st_dimensions(rastered_predsSE)
xy <- attr(rastered_predsSE, "raster")$dimensions
if (!all(is.na(xy))) { # has x/y spatial dimensions:
    d[[ xy[1] ]]$refsys <- st_crs(9122)
    d[[ xy[2] ]]$refsys <- st_crs(9122)
}
rastered_predsSE<-structure(rastered_predsSE, dimensions = d)

Dates<-st_get_dimension_values(rastered_predsSE, "Date")

all_points_static<-select(rastered_predsSE, Prediction)%>%
    aggregate(by=c(as.Date("1960-01-01"), Sys.Date()), function(x) length(which(!is.na(x))))%>%
    filter(time==as.Date("1960-01-01"))%>%
    mutate(across(Prediction, ~na_if(.x, 0)))%>%
    st_as_sf()%>%
    rename(N=`1960-01-01`)%>%
    mutate(ID=as.character(1:n()))%>%
    st_transform(crs=4326)

pal_N_static<-colorNumeric("viridis", domain=range(all_points_static$N, na.rm=T), na.color="#00000000")

pal_N_rev_static<-colorNumeric("viridis", domain=range(all_points_static$N, na.rm=T), reverse=T, na.color="#00000000")

Point_plot<-leaflet()%>%
    addProviderTiles("Esri.WorldGrayCanvas")%>%
    addFeatures(data=all_points_static, fillColor=~pal_N_static(N), color="black", fillOpacity = 0.2, label=~N, layerId = ~ID, weight=0.4)%>%
    addLegend(data=all_points_static, position="topright", pal = pal_N_rev_static, values = ~N, opacity=0.5, 
              labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)), title="Number of time</br>points with</br>model-predicted data")

Delta_regions_static<-readRDS("Delta subregions.Rds")%>%
    mutate(SubRegion=as.character(SubRegion))

# To fix error with different versions of Proj on shinyapps.io vs. the saved object
st_crs(Delta_regions_static)<-26910

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

Region_plot<-leaflet()%>%
    addProviderTiles("Esri.WorldGrayCanvas")%>%
    addFeatures(data=st_transform(Delta_regions_static, 4326), fillColor=~pal_N2_static(N), 
                color="black", fillOpacity = 0.2, label=lapply(paste0("<h5 align='center'>", Delta_regions_static$SubRegion, "</h5>", "<h6 align='center' style='color:red'>N: ", 
                                                                      format(Delta_regions_static$N, big.mark   = ","), "</h6>"), htmltools::HTML), 
                layerId = ~SubRegion, weight=0.4)%>%
    addLegend(data=st_transform(Delta_regions_static, 4326), position="topright", pal = pal_N2_rev_static, values = ~N, opacity=0.5, 
              labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)), title="Number of raster cells</br>and time points with</br>model-predicted data")

pal_N3_static<-colorNumeric("viridis", domain=range(Delta_regions_static$N_data, na.rm=T), na.color="#00000000")

pal_N3_rev_static<-colorNumeric("viridis", domain=range(Delta_regions_static$N_data, na.rm=T), reverse=T, na.color="#00000000")

Region_plot2<-leaflet()%>%
    addProviderTiles("Esri.WorldGrayCanvas")%>%
    addFeatures(data=st_transform(Delta_regions_static, 4326), fillColor=~pal_N3_static(N_data), 
                color="black", fillOpacity = 0.2, label= lapply(paste0("<h5 align='center'>", Delta_regions_static$SubRegion, "</h5>", "<h6 align='center' style='color:red'>N: ", 
                                                                       format(Delta_regions_static$N_data, big.mark   = ","), "</h6>"), htmltools::HTML),
                layerId = ~SubRegion, weight=0.4)%>%
    addLegend(data=st_transform(Delta_regions_static, 4326), position="topright", pal = pal_N3_rev_static, values = ~N_data, opacity=0.5, 
              labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)), title="Total number of data points</br>(actual temperature records)")

mygrid <- data.frame(
    name = c("Upper Sacramento River Ship Channel", "Cache Slough and Lindsey Slough", "Lower Sacramento River Ship Channel", "Liberty Island", "Suisun Marsh", 
             "Middle Sacramento River", "Lower Cache Slough", "Steamboat and Miner Slough", "Upper Mokelumne River", "Lower Mokelumne River", "Georgiana Slough", 
             "Sacramento River near Ryde", "Sacramento River near Rio Vista", "Grizzly Bay", "West Suisun Bay", "Mid Suisun Bay", "Honker Bay", "Confluence", 
             "Lower Sacramento River", "San Joaquin River at Twitchell Island", "San Joaquin River at Prisoners Pt", "Disappointment Slough", "Lower San Joaquin River", 
             "Franks Tract", "Holland Cut", "San Joaquin River near Stockton", "Mildred Island", "Middle River", "Old River", "Upper San Joaquin River", 
             "Grant Line Canal and Old River", "Victoria Canal"),
    row = c(2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 8, 8, 8),
    col = c(7, 4, 6, 5, 2, 8, 6, 7, 9, 9, 8, 7, 6, 2, 1, 2, 3, 4, 5, 6, 8, 9, 5, 6, 7, 9, 8, 8, 7, 9, 8, 7),
    code = c(" 1", " 1", " 2", " 3", " 8", " 4", " 5", " 6", " 7", " 9", "10", "11", "12", "13", "14", "15", 
             "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "30", "29", "31"),
    stringsAsFactors = FALSE
)

Model_eval<-readRDS("Model_eval.Rds")

CV_sum<-readRDS("CV_sum.Rds")%>%
    mutate(Year=recode(Group, `1`="Even years", `2`="Odd years"))


# Raw data ----------------------------------------------------------------


Data<-readRDS("Raw temp data.Rds")%>%
    mutate(Station=case_when(
        str_detect(Station, "EZ") ~ paste("EMP", str_extract(Station, "([^\\s]+)")),
        Source=="EDSM" ~ paste("EDSM", SubRegion),
        TRUE ~ paste(Source, Station)
    ))
# To fix error with different versions of Proj on shinyapps.io vs. the saved object
st_crs(Data)<-26910

all_stations_static<-Data%>%
    filter(Source!="EDSM" & !str_detect(Station, "EZ") & !StationID%in%c("SKT 799", "SKT 999", "SKT 699"))%>%
    st_drop_geometry()%>%
    group_by(StationID, Latitude, Longitude)%>%
    summarise(N=n(), .groups="drop")%>%
    distinct()%>%
    st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)

pal_N3_static<-colorNumeric("viridis", domain=range(log(all_stations_static$N), na.rm=T), na.color="#00000000")

pal_N3_rev_static<-colorNumeric("viridis", domain=range(log(all_stations_static$N), na.rm=T), reverse=T, na.color="#00000000")

Station_plot<-leaflet()%>%
    addProviderTiles("Esri.WorldGrayCanvas")%>%
    addFeatures(data=all_stations_static, fillColor=~pal_N3_static(log(N)), 
                color="black", fillOpacity = 0.2, label=lapply(paste0("<h5 align='center'>", all_stations_static$StationID, "</h5>", "<h6 align='center' style='color:red'>N: ", 
                                                                      format(all_stations_static$N, big.mark   = ","), "</h6>"), htmltools::HTML), 
                layerId = ~StationID, weight=0.4)%>%
    addLegend(data=all_stations_static, position="topright", pal = pal_N3_rev_static, values = ~log(N), opacity=0.5,
              labFormat = labelFormat(transform = function(x) sort(round(exp(x)), decreasing = TRUE)), title="Number of raw</br>temperature records</br>(log scale)")


# Time correction

Time_correction<-readRDS("Time_correction.Rds")


#Settings for the "data crunching" message. 
info_loading <- "Crunching data"
progress_color <- "black"
progress_background <- "#c5c5c9"



# User interface ----------------------------------------------------------


ui <- fluidPage(
    theme = shinytheme("cerulean"),
    a(shiny::icon("reply"), "Delta Science shinyapps homepage", href="https://deltascience.shinyapps.io/Home/"),
    
    # Application title
    titlePanel(title=div(h1("Delta Discrete Temperature Model", style="display: inline-block"), 
                         a(img(src="DSP_Logo_Horizontal.png", height = 100, align="right", style="display: inline-block"), href="https://deltacouncil.ca.gov/delta-science-program/"),
                         h5("If you encounter any issues, please email ", a("sam.bashevkin@deltacouncil.ca.gov.", 
                                                                            href="mailto:sam.bashevkin@deltacouncil.ca.gov?subject=Discrete%20Temperature%20Shiny%20app"))), 
               windowTitle = "Delta Discrete Temperature"),
    
    sidebarLayout(
        # Sidebar with a slider input for number of bins 
        sidebarPanel(width=3,
                     #Information
                     actionBttn("Information", "Information", style="simple", color="primary", icon=icon("question-circle")),
                     h1("Filters"),
                     dateRangeInput("Date_range", label = "Date range", 
                                    start = min(Dates), end = max(Dates), startview = "year"),
                     pickerInput("Months", "Months:", choices=c("January" = 1, "February" = 2, "March" = 3, "April" = 4, 
                                                                "May" = 5, "June" = 6, "July" = 7, "August" = 8, 
                                                                "September" = 9, "October" = 10, "November" = 11, 
                                                                "December" = 12), 
                                 selected = 1:12, multiple = T, options=list(`actions-box`=TRUE, `selected-text-format` = "count > 3")),
                     #Buttons to download data
                     conditionalPanel(condition="input.Tab!='Raw data'",
                                      actionBttn("Download", "Download data", style="simple", color="royal", icon=icon("file-download"))),
                     h1("Plot options"),
                     conditionalPanel(condition="input.Tab=='Rasters' && !input.Habitat", 
                                      radioGroupButtons("variable",
                                                        "Plot:",
                                                        choices = c("Predicted temperature"="Prediction", "Standard Error"="SE"), selected="Prediction", individual = TRUE, 
                                                        checkIcon = list( yes = tags$i(class = "fa fa-circle", style = "color: steelblue"), 
                                                                          no = tags$i(class = "fa fa-circle-o", style = "color: steelblue")))),
                     conditionalPanel(condition="input.Tab=='Rasters' || input.Tab=='Time series'", 
                                      prettySwitch("Habitat",HTML("<b>Plot habitat suitability?</b>"), status = "success", fill = TRUE, bigger=TRUE),
                                      conditionalPanel(condition="input.Habitat",
                                                       sliderInput("Habitat_range", "Use the two sliders to choose 3 temperature regions of differing suitability for your focal species",
                                                                   min=Temp_min, max=Temp_max, value=c(20, 22), step=0.1, round=-1, post=" °C"),
                                                       wellPanel(style = "background: #cccccc",
                                                                 h3("Define labels for your 3 temperature regions"),
                                                                 textInput("Habitat_range_high", "High", value="Bad"),
                                                                 textInput("Habitat_range_med", "Medium", value="Tolerable"),
                                                                 textInput("Habitat_range_low", "Low", value="Good"))
                                      )),
                     conditionalPanel(condition="input.Tab=='Rasters' || input.Tab=='Time series' || input.Tab=='Raw data'", 
                                      uiOutput("facet_options")),
                     conditionalPanel(condition="input.Tab=='Rasters' && input.Facet!='Year x Month' && input.Facet!='Month x Year' && !input.Habitat", 
                                      uiOutput("scale_options")),
                     conditionalPanel(condition="input.Tab=='Time series'", 
                                      radioGroupButtons("Regions",
                                                        "Draw your own regions or select preset regions?",
                                                        choices = c(`<i class='fa fa-pen'> Draw</i>`="Draw", `<i class='fa fa-shapes'> Preset</i>`="Preset"), 
                                                        selected="Draw", individual = TRUE),
                                      conditionalPanel(condition="input.Regions=='Preset'",
                                                       prettySwitch("Regions_all","Select all regions?", status = "success", fill = TRUE))),
                     conditionalPanel(condition="input.Tab=='Model evaluation'",
                                      radioGroupButtons("Uncertainty",
                                                        "Metric to plot:",
                                                        choices = c("Model residuals", "Cross-validation", "Sample size of measured values"), 
                                                        selected="Model residuals", individual = TRUE, 
                                                        checkIcon = list( yes = tags$i(class = "fa fa-circle", style = "color: steelblue"), 
                                                                          no = tags$i(class = "fa fa-circle-o", style = "color: steelblue"))),
                                      conditionalPanel(condition="input.Uncertainty!='Sample size of measured values'", 
                                                       uiOutput("summary_stats"))),
                     conditionalPanel(condition="input.Tab=='Raw data'",
                                      h5("Include data from non-fixed stations (EDSM, EMP EZ stations, and SKT stations 699, 799, and 999) that are not plotted on map?"), 
                                      materialSwitch(
                                          inputId = "Nonfixed",
                                          value = TRUE, inline=F,
                                          status = "primary"),
                                      h5("Apply time-of-day correction to normalize all data to noon?"), 
                                      materialSwitch(
                                          inputId = "Time_correction",
                                          value = TRUE, inline=F,
                                          status = "primary")),
                     actionBttn("Plot_info", label = "Plot description", style="simple", 
                                color="success", icon=icon("question-circle"))
        ),
        mainPanel(width=9,
                  tabsetPanel(type="tabs",
                              id="Tab",
                              tabPanel("Model evaluation",
                                       fluidRow(girafeOutput("Uncertainty_plot", height="100vh", width="130vh")),
                                       conditionalPanel(condition="input.Uncertainty=='Sample size of measured values' || input.Uncertainty_value=='Mean residuals' 
                                                        || input.Uncertainty_value=='Mean magnitude of residuals' || input.Uncertainty_value=='SD of residuals'", 
                                                        h2("Region map", align = "center"),
                                                        fluidRow(leafletOutput("Regions_plot2", width = "100%", height = "100%")))),
                              tabPanel("Rasters",
                                       conditionalPanel(condition="input.Facet!='Year' && input.Facet!='Year x Month' && input.Facet!='Month x Year'",
                                                        uiOutput("select_Year")),
                                       conditionalPanel(condition="input.Facet!='Month' && input.Facet!='Year x Month' && input.Facet!='Month x Year'",
                                                        uiOutput("select_Month")),
                                       plotOutput("TempPlot")
                              ),
                              tabPanel("Time series",
                                       fluidRow(conditionalPanel(condition="input.Regions=='Draw'", editModUI("Pointplot", height = "80vh", width="100%")),
                                                conditionalPanel(condition="input.Regions=='Preset'", editModUI("Regionplot", height = "80vh", width="100%"))),
                                       fluidRow(conditionalPanel(condition="input.Facet!='Month'",
                                                                 column(2, h5("Plot all months?"), materialSwitch(
                                                                     inputId = "Month2_slider",
                                                                     value = FALSE, inline=F,
                                                                     status = "primary")), 
                                                                 column(10,conditionalPanel(condition="!input.Month2_slider",uiOutput("select_Month2"))))),
                                       fluidRow(girafeOutput("Time_plot", height="100vh", width="130vh"))
                              ),
                              tabPanel("Raw data",
                                       fluidRow(editModUI("Rawmap", height = "80vh", width="100%")),
                                       fluidRow(conditionalPanel(condition="input.Facet!='Month'",
                                                                 column(2, h5("Plot all months?"), materialSwitch(
                                                                     inputId = "Month3_slider",
                                                                     value = FALSE, inline=F,
                                                                     status = "primary")), 
                                                                 column(10,conditionalPanel(condition="!input.Month3_slider",uiOutput("select_Month3"))))),
                                       fluidRow(girafeOutput("Raw_time_plot", height="100vh", width="130vh")))
                  ))
    ),
    tags$head(tags$style("#TempPlot{height:80vh !important;}")),
    tags$head(tags$style("#Pointplot{height:80vh !important;}")),
    tags$head(tags$style("#Regionplot{height:80vh !important;}")),
    tags$head(tags$style("#Time_plot{height:80vh !important;}")),
    tags$head(tags$style("#Regions_plot2{height:80vh !important;}")),
    tags$head(tags$style("#Rawmap{height:80vh !important;}")),
    tags$head(tags$style("#Raw_time_plot{height:80vh !important;}")),
    
    # Display the "data crunching" message. 
    tags$head(tags$style(type="text/css",
                         paste0("
                                             #loadmessage {
                                             position: fixed;
                                             top: 25%;
                                             left: 25%;
                                             width: 50%;
                                             padding: 30px 0px 30px 0px;
                                             text-align: center;
                                             font-weight: bold;
                                             font-size: 200%;
                                             color: ", progress_color,";
                                             background-color: ", progress_background,";
                                             z-index: 105;
                                             border: 2px solid black;
                                             border-radius: 50px;
                                             }
                                             "))),
    conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                     tags$div(tags$i(class = "fa fa-spinner", style = "color: black"), info_loading, tags$i(class = "fa fa-spinner", style = "color: black"), id="loadmessage")),
    tags$style(type="text/css", ".recalculating {opacity: 1.0;}"),
    tags$style(HTML(" [for=Habitat_range]+span>.irs-bar {background: linear-gradient(rgb(221, 221, 221) -50%, rgb(255, 255, 255) 150%); border:1px solid #cccccc}"))
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    
    # Informational popups ----------------------------------------------------
    
    #Popup for app Information
    observeEvent(input$Information, {
        sendSweetAlert(session, title = "Information", 
                       text = tags$span(tags$p("This app presents model-smoothed water temperatures across the Delta 
                                               and Suisun regions. Data are predicted from a generalized additive model 
                                               and represent an approximation of surface water temperatures, not the 
                                               actual measured values. The generalized additive model was fit to an integrated
                                               discrete water quality dataset."),
                                        tags$p("The model was fit with the R package", tags$code("mgcv"), "with the following model structure:"), 
                                        tags$p(align="left", tags$code("bam(Temperature ~ Year_fac + te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs=c('tp', 'cc'), k=c(25, 20), by=Year_fac) + 
                                        te(Time_num_s, Julian_day_s, bs=c('tp', 'cc'), k=c(5, 12)),
    data = Data, method='fREML', discrete=T)")),
                                        tags$p("The first tab 'Model evaluation' presents the uncertainty in model predictions 
                                               and should be explored to understand the limitations of the model."),
                                        tags$p("Use the data filters to control the range of dates and months included in the plots."),
                                        tags$p("The raw data can be plotted in a time-series within your region of interest on the 'Raw data' tab."),
                                        "------------------------------------------",
                                        tags$p(tags$b("Dataset, model, and app created and maintained by Sam Bashevkin (Delta Science Program) in collaboration with Brian Mahardja (USFWS) and Larry Brown (USGS). 
                                                      Please email", a("Sam", href="mailto:sam.bashevkin@deltacouncil.ca.gov?subject=Discrete%20Temperature%20Shiny%20app"), "with any questions.")),
                                        tags$p(tags$i("The integrated dataset used to fit this model included data from the", 
                                                      tags$a("Environmental Monitoring Program, ", href = "https://portal.edirepository.org/nis/mapbrowse?packageid=edi.458.2", target="_blank"), 
                                                      tags$a("Summer Townet, ", href = "https://wildlife.ca.gov/Conservation/Delta/Townet-Survey", target="_blank"),
                                                      tags$a("Fall Midwater Trawl, ", href = "https://wildlife.ca.gov/Conservation/Delta/Fall-Midwater-Trawl", target="_blank"),
                                                      tags$a("Spring Kodiak Trawl, ", href = "https://portal.edirepository.org/nis/mapbrowse?packageid=edi.527.1", target="_blank"),
                                                      tags$a("20-mm Survey, ", href = "https://portal.edirepository.org/nis/mapbrowse?packageid=edi.535.1", target="_blank"),
                                                      tags$a("Bay Study, ", href = "https://wildlife.ca.gov/Conservation/Delta/Bay-Study", target="_blank"),
                                                      tags$a("Enhanced Delta Smelt Monitoring Program, ", href = "https://portal.edirepository.org/nis/mapbrowse?packageid=edi.415.3", target="_blank"),
                                                      tags$a("Delta Juvenile Fish Monitoring Program, ", href = "https://portal.edirepository.org/nis/mapbrowse?packageid=edi.244.4", target="_blank"),
                                                      tags$a("UC Davis Suisun Marsh Fish Study, ", href = "https://watershed.ucdavis.edu/project/suisun-marsh-fish-study", target="_blank"),
                                                      "USBR Sacramento Deepwater Shipping Channel surveys, and the",
                                                      tags$a("USGS San Francisco Bay Water Quality Surveys.", href = "https://www.nature.com/articles/sdata201798", target="_blank")))),
                       type = "info",
                       btn_labels = "Ok", html = F, closeOnClickOutside = TRUE)
    })    
    
    #Plot info
    Plot_info<-reactive({
        if(input$Tab=="Model evaluation"){
            out<-tags$span(h2("Visualizations of model error and sample size"),
                           h4("Cross validation methods"),
                           p("A stratified cross-validation approach was used, in which
                             the dataset was split into 20 parts. First, the dataset was split into even and odd years, then within each year grouping the data
                             were further divided into 10 parts evenly split across regions, years, and seasons. Twenty separate datasets were created
                             each missing 1 of those 20 parts and predictions were generated for the missing parts. CV residuals were not calculated for years
                             before 1974 due to the very small number of data points in those years."),
                           h4("Metric descriptions"),
                           p("Residuals represent the deviations of model predictions from true measured values.
                              The residuals were either calculated from the full model (model residuals) or from the cross-validation models described above (cross-validation residuals)."),
                           p("The 'Mean magnitude' option for visualizing residuals allows you to explore the average magnitude of error. 
                             Since the residuals can be positive or negative, the regular mean may result in large positive or negative values canceling one another out."),
                           p("The standard deviation of the residual values should help inform the range of error magnitudes. Small mean errors could still be problematic if individual
                             errors can be large."),
                           p("The cross-validation results can be further explored with the RMSE (root-mean-square error) or Pearson's correlation (r) options. 
                              Each metric is calculated separately for each of the 20 cross-validation models."),
                           p("The 'Sample size of measured values' options allows you to explore the number of temperature measurements in each month, year, and region.
                             This should give you an idea of the data available to the model during fitting."),
                           p("These plots are interactive, so hover your mouse over the graph to see the data values."),
                           p("Scroll down to see a map of the regions used to produce the plot."))
        } else{
            if(input$Tab=="Rasters"){
                out<-tags$span(h2("Raterized plots of model-predicted surface water temperatures"),
                               p("Model predictions have been spatially and temporally trimmed to the extent of the integrated dataset. Predictions were only generated for
                                 regions sampled in each month and year. The spatial extent of these raster plots changes over time due to variability in this sampling effort."),
                               p("Plots can be facetted by year, month, or both. Facetted plots will be slower to load."),
                               p("Year and month sliders can be used to animate through time-series. This works best with fewer facets."),
                               p(tags$b("It is highly recommended to also 'fix' the scale when scrolling or animating through time so the color scale does not
                                 change with every new time point.")),
                               p("You can also visualize the extent of suitable habitat for a species of interest by clicking the 'Habitat suitability' slider. Select the same
                                  option in the 'Time series' tab to explore changes in the spatial extent of suitable habitat over time."))
            }else{
                if(input$Tab=="Time series"){
                    out<-tags$span(h2("Time series plots"),
                                   p("Select spatial regions of interest and plot a timeseries for those regions."),
                                   p("Model predictions have been spatially and temporally trimmed to the extent of the integrated dataset. Predictions were only generated for
                                 regions sampled in each month and year. Explore the Raster plots, or the sample size plots in the 'Model evaluation' tab to see how the 
                                     spatial extent of sampling has changed over time."),
                                   p("You can visualize changes in the extent of suitable habitat for a species of interest by clicking the 'Habitat suitability' slider. Select the same option in the 
                                     'Rasters' tab to explore the spatial distribution of suitable habitat."),
                                   p("You have the choice to draw your own spatial regions, or select from a set of pre-set regions based on the EDSM sub-regions (controlled with the slider)."),
                                   p("In either case, you can draw rectangles or polygons on the map to select any spatial areas they touch. You can also drop a marker to select invididual cells or regions.
                                 There are also buttons for editing or deleting your drawn areas. Hover over the buttons to see their use. The map will update to highlight the areas you've selected.
                                 It may take a few seconds to load."),
                                   p("If you want to average across all the areas you select to produce one time-series, leave the facets option on 'None'. 
                                   But if you want to plot each region separately, switch it to 'Region'. You can also choose to facet plots by month."),
                                   p("Use the month slider to control which month is represented on the plot. To plot the full time-series across all months, switch the slider labeled 'Plot all months?'"),
                                   p("These plots are interactive, so hover to reveal the data values. For the time-series plot, you need to hover over the points (i.e., vertices in the line)."),
                                   tags$em("NOTE: If you select areas that span multiple regions with different sample sizes, different subsets of spatial locations will be included in the averaged
                                       temperature values for each date. For example, the Confluence is heavily represented, while the Upper San Joaquin River is often missing dates in the timeseries.
                                       Averaging across both these areas would be biased by some time points including the warm Upper San Joaquin River temperatures, and others only including the
                                       cooler Confluence temperatures."))
                }else{
                    out<-tags$span(h2("Raw data time-series"),
                                   p("This section will display the raw temperature data, rather than the model-predicted data shown in the other tabs"),
                                   p("Select sampling stations of interest and plot a timeseries of raw data for those regions."),
                                   p("Only fixed stations are displayed on the map but data from unfixed stations (EDSM, EMP EZ stations, and SKT stations 699, 799, and 999)
                                     within the selected region are also included in plots by default unless the slider on the left is clicked."),
                                   p("By default, time-of-day effects will be eliminated by correcting temperature data to noon. 
                                     This correction is applied using the value of the time-of-day smoother in the GAM presented in the other 3 tabs.
                                     You can disable this correction using the slider to the left."),
                                   tags$em("Plots facetted by month will be slow to load, especially if a large number of stations are selected."))
                }
            }
        }
        return(out)
    })
    
    observeEvent(input$Plot_info, {
        sendSweetAlert(session, title = "Plot info", text = Plot_info(),
                       type = "info",
                       btn_labels = "Ok", html = FALSE, closeOnClickOutside = TRUE)
    })
    
    #Show modal dialog to save data when Download button is clicked
    observeEvent(input$Download, {
        showModal(ModalDownloadData())
    })
    
    #Modal dialog (popup window) to download data
    ModalDownloadData<-function(){
        Data_info<-if(input$Tab=="Model evaluation"){
            span(p("This will download the model evaluation dataset. To download the temperature data, navigate to one of the other tabs."),
                 p("The variable", code('Resid'), "represents the model residuals,", code('SD'), "represents the standard deviation of those model residuals.",
                   code('Resid_CV'), "and", code('SD_CV'), "represent the same for the cross-validation residuals.", code('N'), "represents the number of measured data points in
                   the integrated dataset used to fit these models."))
        } else{
            if(input$Tab=="Rasters"){
                p("This will download the rasterized model predictions in tif format. Only the variable (Predicted Temperature or Standard Error) 
                  selected under 'Plot options' will be downloaded. To download both, you must download 2 separate files selecting a different variable each time.")
            }else{
                span(p("This will download the model predictions converted to a time-series format. Only predictions from selected areas on the map will be included."),
                     p(code("Prediction"), "represents the model-predicted temperature in °C, while", code("SE"), "represents the standard error."),
                     p("Data will be averaged for each date across the region of interest."),
                     h5("If you select", tags$b("Region"), "in the", tags$b("Facet"), "selections under", tags$b('Plot options,'), "your
                       data will come divided by regions. Otherwise, the data will be averaged across all selected regions", style="color:DarkOrchid"))
            }
        }
        
        modalDialog(
            h1("Data info"),
            h5("Data will be clipped to the date and month filters before downloading.", style="color:red"),
            Data_info,
            footer = tagList(modalButton("Cancel"),
                             downloadBttn("Downloaddata", "Download data", style="bordered", color = "primary", size="sm")),
            easyClose=TRUE
        )
    }
    
    #Download handler for data downloading
    output$Downloaddata <- downloadHandler(
        filename = function() {
            if(input$Tab=="Rasters"){
                ext<-".tif"
            }else{
                ext<-".csv"
            }
            
            paste0("Delta temperature ", tolower(input$Tab), ext)
            
        },
        content = function(file) {
            
            if(input$Tab=="Model evaluation"){
                Model_eval%>%
                    filter(Year>year(min(input$Date_range)) & Year<year(max(input$Date_range)) & Month%in%input$Months)%>%
                    write_csv(file)
            }else{
                if(input$Tab=="Rasters"){
                    write_stars(TempData(), file, layer=input$variable)
                }else{
                    timeseries_data()%>%
                        mutate(Date=as.character(Date, format="%Y-%m-%d"))%>%
                        write_csv(file)
                }
            }
            
        })
    
    # Uncertainty plot --------------------------------------------------------
    
    output$summary_stats <- renderUI({
        if(input$Uncertainty=="Cross-validation"){
            choices <- c("Mean residuals", "Mean magnitude of residuals", "SD of residuals", "RMSE", "Pearson's correlation")
        } else{
            choices <- c("Mean residuals", "Mean magnitude of residuals", "SD of residuals")
        }
        
        radioGroupButtons("Uncertainty_value",
                          "Summary statistic:",
                          choices = choices, selected="Mean residuals", individual = TRUE, 
                          checkIcon = list( yes = tags$i(class = "fa fa-circle", style = "color: steelblue"), no = tags$i(class = "fa fa-circle-o", style = "color: steelblue")))
    })
    
    Uncertainty_plot<-reactive({
        req(input$Uncertainty, input$Date_range, any(input$Months%in%1:12), input$Uncertainty_value)
        
        if(input$Uncertainty=="Model residuals"){
            Data<-Model_eval%>%
                rename(Fill=Resid, Fill_mag=Resid_mag, Fill_SD=SD)
            
        }else{
            if(input$Uncertainty=="Sample size of measured values"){
                Data<-Model_eval%>%
                    rename(Fill=N)
            }else{
                if(input$Uncertainty_value%in%c("RMSE", "Pearson's correlation")){
                    Data<-CV_sum%>%
                        {if(input$Uncertainty_value=="RMSE"){
                            rename(., Metric=RMSE)
                        }else{
                            rename(., Metric=r)
                        }}%>%
                        group_by(Group)%>%
                        mutate(Median=median(Metric))%>%
                        ungroup()%>%
                        select(Group, Year, Fold, Metric, Median)
                } else{
                    Data<-Model_eval%>%
                        rename(Fill=Resid_CV, Fill_mag=Resid_mag_CV, Fill_SD=SD_CV)%>%
                        filter(Year>=1974)
                }
            }
            
        }
        
        if(input$Uncertainty=="Sample size of measured values"){
            str_model <- paste0("<tr><td>N: &nbsp</td><td>%s</td></tr>",
                                "<tr><td>Year: &nbsp</td><td>%s</td></tr>", 
                                "<tr><td>Month: &nbsp</td><td>%s</td></tr>")
            Data<-Data%>%
                filter(Year>year(min(input$Date_range)) & Year<year(max(input$Date_range)) & Month%in%input$Months)%>%
                mutate(Month=month(Month, label=T))%>%
                mutate(tooltip=sprintf(str_model, round(Fill, 2), Year, Month),
                       tooltip=paste0( "<table>", tooltip, "</table>" ),
                       ID=1:n())
        } else{
            if(input$Uncertainty_value%in%c("RMSE", "Pearson's correlation")){
                str_model1 <- paste0("<tr><td>Data grouping: &nbsp</td><td>%s</td></tr>",
                                     "<tr><td>Median: &nbsp</td><td>%s</td></tr>")
                
                str_model2 <- paste0("<tr><td>Data grouping: &nbsp</td><td>%s</td></tr>",
                                     "<tr><td>Fold:  &nbsp</td><td>%s</td></tr>",
                                     "<tr><td>Value: &nbsp</td><td>%s</td></tr>")
                
                Data<-Data%>%
                    mutate(tooltip1=sprintf(str_model1, Year, round(Median, 4)),
                           tooltip1=paste0( "<table>", tooltip1, "</table>" ),
                           tooltip2=sprintf(str_model2, Year, Fold, round(Metric, 4)),
                           tooltip2=paste0( "<table>", tooltip2, "</table>" ))
            }else{
                
                str_model <- paste0("<tr><td>Mean resid: &nbsp</td><td>%s</td></tr>",
                                    "<tr><td>Mean resid magnitude: &nbsp</td><td>%s</td></tr>",
                                    "<tr><td>SD: &nbsp</td><td>%s</td></tr>",
                                    "<tr><td>Year: &nbsp</td><td>%s</td></tr>", 
                                    "<tr><td>Month: &nbsp</td><td>%s</td></tr>")
                Data<-Data%>%
                    filter(Year>year(min(input$Date_range)) & Year<year(max(input$Date_range)) & Month%in%input$Months)%>%
                    mutate(Month=month(Month, label=T))%>%
                    mutate(tooltip=sprintf(str_model, round(Fill, 2), round(Fill_mag, 2), round(Fill_SD, 2), Year, Month),
                           tooltip=paste0( "<table>", tooltip, "</table>" ),
                           ID=1:n())
            }
        }
        
        if(input$Uncertainty=="Cross-validation" & input$Uncertainty_value%in%c("RMSE", "Pearson's correlation")){
            
            p_resid<-ggplot(Data, aes(x=as.factor(Group), y=Metric, group=as.factor(Group), fill=as.factor(Group)))+
                geom_boxplot_interactive(aes(tooltip=tooltip1, data_id=Group), show.legend=F, outlier.shape = NA)+
                geom_point_interactive(aes(x=Group+Fold/100, tooltip=tooltip2, data_id=paste(Group, Fold)), show.legend=F)+
                scale_x_discrete(labels=c("Even years", "Odd years"))+
                ylab(if_else(input$Uncertainty_value=="RMSE", "RMSE (°C)", "r"))+
                xlab("Data grouping")+
                theme_bw()
            
            
        }else{
            
            p_resid<-ggplot(Data)+
                {if(input$Uncertainty_value=="Mean" | input$Uncertainty=="Sample size of measured values"){
                    geom_tile_interactive(aes(x=Year, y=Month, fill=Fill, color=Fill, tooltip=tooltip, data_id=ID))
                }else{ 
                    if(input$Uncertainty_value=="Mean magnitude"){
                        geom_tile_interactive(aes(x=Year, y=Month, fill=Fill_mag, color=Fill_mag, tooltip=tooltip, data_id=ID))
                    }else{
                        geom_tile_interactive(aes(x=Year, y=Month, fill=Fill_SD, color=Fill_SD, tooltip=tooltip, data_id=ID))
                    }
                    
                }}+
                {if(input$Uncertainty_value=="Mean" & !input$Uncertainty=="Sample size of measured values"){
                    scale_fill_gradient2(name=paste("Mean", tolower(input$Uncertainty), "(°C)"),
                                         high = muted("red"),
                                         low = muted("blue"),
                                         breaks=seq(-9,6.5, by=0.5),
                                         guide=guide_colorbar(direction="horizontal", title.position = "top", ticks.linewidth = 1,
                                                              title.hjust=0.5, label.position="bottom"),
                                         aesthetics = c("color", "fill"))
                }else{
                    if(input$Uncertainty=="Sample size of measured values"){
                        scale_fill_viridis_c(name="Sample size",
                                             guide=guide_colorbar(direction="horizontal", title.position = "top", ticks.linewidth = 1,
                                                                  title.hjust=0.5, label.position="bottom"),
                                             aesthetics = c("color", "fill"))
                    }else{
                        if(input$Uncertainty_value=="Mean magnitude"){
                            scale_fill_viridis_c(name=paste("Mean magnitude of", tolower(input$Uncertainty), "(°C)"),
                                                 breaks=seq(0,9, by=0.5),
                                                 guide=guide_colorbar(direction="horizontal", title.position = "top", ticks.linewidth = 1,
                                                                      title.hjust=0.5, label.position="bottom"),
                                                 aesthetics = c("color", "fill")) 
                        }else{
                            scale_fill_viridis_c(name=paste("Standard deviation in", tolower(input$Uncertainty), "(°C)"),
                                                 breaks=seq(0,4.5, by=0.5),
                                                 guide=guide_colorbar(direction="horizontal", title.position = "top", ticks.linewidth = 1,
                                                                      title.hjust=0.5, label.position="bottom"),
                                                 aesthetics = c("color", "fill")) 
                        }
                    }
                    
                }}+
                scale_x_continuous(breaks=seq(1970, 2020, by=10))+
                {if(setequal(1:12, input$Months)){
                    scale_y_discrete(breaks=levels(Data$Month), labels = if_else(sort(as.integer(unique(month(sample(1:12, 100, replace=T), label=T)))+1)%% 2 == 0, levels(Data$Month), ""))
                }}+
                facet_geo(~SubRegion, grid=mygrid, labeller=label_wrap_gen())+
                theme_bw()+
                theme(text=element_text(size=4), axis.text.x=element_text(angle=45, hjust=1), 
                      panel.grid=element_blank(),
                      strip.background = element_blank(), panel.spacing = unit(0.2, units="lines"),
                      axis.ticks=element_line(size=0.1), axis.ticks.length = unit(0.2, units="lines"),
                      legend.position="top", legend.key.width = unit(50, "native"), legend.justification="center",
                      legend.box.spacing = unit(0.2, "lines"),legend.box.margin=margin(b=-50))+
                {if(input$Uncertainty_value=="Mean" & !input$Uncertainty=="Sample size of measured values"){
                    theme(panel.background = element_rect(fill="black"))
                }else{
                    theme(panel.background = element_rect(fill="white"))
                }}
        }
        return(p_resid)
    })
    
    output$Uncertainty_plot<-renderGirafe({
        p<-girafe(ggobj=Uncertainty_plot(), width_svg = 5.5,  pointsize=6, options = list(
            opts_sizing(rescale = T, width = 1)))
        girafe_options(p, opts_toolbar(saveaspng = FALSE), opts_selection(type="none"))
    })
    
    output$Regions_plot2<-renderLeaflet({
        Region_plot2
    })
    
    
    # Rasters -----------------------------------------------------------------
    
    TempData<-reactive({
        req(input$Date_range, any(input$Months%in%1:12))
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
    
    output$select_Month2 <- renderUI({
        req(!is.null(input$Months))
        sliderInput("Month2",
                    "Select month:", min=min(month(Data_dates())), max=max(month(Data_dates())), value =  min(month(Data_dates())), step=1, round=T, sep="", 
                    animate=animationOptions(interval=1000, loop=T), width="100%")
    })
    
    output$select_Month3 <- renderUI({
        req(!is.null(input$Months))
        sliderInput("Month3",
                    "Select month:", min=min(Data_filtered()$Month), max=max(Data_filtered()$Month), value =  min(Data_filtered()$Month), step=1, round=T, sep="", 
                    animate=animationOptions(interval=1000, loop=T), width="100%")
    })
    
    output$facet_options<-renderUI({
        if(input$Tab=='Rasters'){
            Choices<-c("None", "Month", "Year", "Month x Year", "Year x Month")
        }else{
            if(input$Tab=="Time series"){
                Choices<-c("None", "Month", "Region")
                
            }else{
                Choices<-c("None", "Month")
            }
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
            }}%>%
            {if(input$Habitat){
                mutate(., Suitability=case_when(
                    Prediction<min(input$Habitat_range) ~ input$Habitat_range_low,
                    Prediction>=min(input$Habitat_range) & Prediction <=max(input$Habitat_range) ~ input$Habitat_range_med,
                    Prediction>max(input$Habitat_range) ~ input$Habitat_range_high
                ))%>%
                    mutate(Suitability=factor(Suitability, levels=c(input$Habitat_range_low, input$Habitat_range_med, input$Habitat_range_high)))
            }else{
                .
            }}
        return(Data)
    }
    )
    
    Plot<-reactive({
        req(input$variable, input$Facet)
        if(input$Facet%in%c("None", "Year")){
            req(input$Month)
        }
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
            {if(input$Habitat){
                geom_stars(data=PlotData()%>%select(Suitability))
                
            }else{
                geom_stars(data=PlotData()%>%select(all_of(input$variable)))
                
            }}+
            {if(input$Habitat){
                scale_fill_manual(values=c("#111D4A", "#FFE74C", "#CE4257"), drop=F, na.translate = FALSE, 
                                  guide = guide_legend(direction="horizontal", title.position = "top",
                                                       title.hjust=0.5, label.position="bottom"))
            }else{
                scale_fill_viridis_c(limits=Scale(), expand=expansion(0,0), name=if_else(input$variable=="Prediction", "Temperature (°C)", "Standard error"), na.value="white", breaks=Breaks, labels= Labels,
                                     guide = guide_colorbar(direction="horizontal", title.position = "top", ticks.linewidth = 2,
                                                            title.hjust=0.5, label.position="bottom"))
            }}+
            coord_sf()+
            {if(input$Facet=="None"){
                ggtitle(paste(month(input$Month, label=T, abbr = FALSE), input$Year))
            }}+
            {if(input$Facet=="Month"){
                facet_wrap(~month(Date, label=T))
            }}+
            {if(input$Facet=="Month"){
                ggtitle(paste("Year:", input$Year))
            }}+
            {if(input$Facet=="Year"){
                facet_wrap(~year(Date))
            }}+
            {if(input$Facet=="Year"){
                ggtitle(paste("Month:", month(input$Month, label=T, abbr = FALSE)))
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
                  panel.grid=element_blank(), 
                  legend.position="top", legend.key.width = unit(150, "native"), legend.justification="center",
                  text=element_text(size=18), plot.title = element_text(size=24, face="bold", hjust=0.5))
    })
    
    output$TempPlot <- renderPlot({
        Plot()
    })
    
    
    # Timeseries plots --------------------------------------------------------
    
    # Draw your own regions/areas of interest -------------------------------------------------
    
    all_points<-reactive({
        select(TempData(), Prediction)%>%
            aggregate(by=c(as.Date("1960-01-01"), Sys.Date()), function(x) length(which(!is.na(x))))%>%
            filter(time==as.Date("1960-01-01"))%>%
            mutate(Prediction=na_if(Prediction, 0))%>%
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
                      labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)))%>%
            {if(!is.null(Points()$Region)){
                addFeatures(., Points()%>%group_by(Region)%>%summarise(geometry=st_union(geometry), .groups="drop"), 
                            label=~Region, color=NULL, fillColor=NULL, weight=0, fillOpacity=0, labelOptions=labelOptions(permanent = TRUE))
            }else{
                .
            }}
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
    
    
    
    #call the editMod function from 'mapedit' to use in the leaflet map.
    
    edits_regions <- callModule(editMod, "Regionplot", leafmap = Region_plot, 
                                editorOptions=list(polylineOptions=F, circleMarkerOptions=F))
    
    observeEvent(Delta_regions_plot(), {  
        proxy.points <- leafletProxy(ns2("map"))
        
        proxy.points %>%
            clearShapes()%>%
            clearControls()%>%
            addFeatures(data=Delta_regions_plot(), fillColor=~pal_N2()(N), color="black", fillOpacity = ~opac, 
                        label=~lapply(paste0("<h5 align='center'>", Delta_regions_plot()$SubRegion, "</h5>", "<h6 align='center' style='color:red'>N: ", 
                                             format(Delta_regions_plot()$N, big.mark   = ","), "</h6>"), htmltools::HTML), 
                        layerId = ~SubRegion, weight=~opac*2)%>%
            addLegend(data=Delta_regions_plot(), position="topright", pal = pal_N2_rev(), values = ~N, opacity=1, 
                      labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)))
    })
    
    
    # Extract and plot timeseries data ----------------------------------------
    
    Points<-reactive({
        req(!(is.null(edits()) & is.null(edits_regions())))
        if(input$Regions=="Draw"){
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
                        rename(Region=X_leaflet_id)%>%
                        mutate(Region=as.integer(as.factor(Region)))
                } else{
                    out<-filter(all_points(), ID%in%Points)
                }
            }else(
                out<-NULL
            )
        }else{
            if(input$Regions_all){
                out<-Delta_regions()%>%
                    rename(Region=SubRegion)%>%
                    st_transform(crs=4326)
                return(out)
            }
            
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
        if(input$Habitat){
            suitability_ag<-function(x){
                case_when(x<min(input$Habitat_range) ~ input$Habitat_range_low,
                          x>=min(input$Habitat_range) & x <=max(input$Habitat_range) ~ input$Habitat_range_med,
                          x>max(input$Habitat_range) ~ input$Habitat_range_high)
            }
            
            if(input$Facet=="Region"){
                Points<-Points()%>%
                    group_by(Region)%>%
                    summarise(geometry=st_union(geometry), .groups="drop")
                
                Low<-TempData()%>%
                    select(Prediction)%>%
                    aggregate(by=Points, function(x) length(which(!is.na(x) & x<min(input$Habitat_range))))%>%
                    st_as_sf(long=F)%>%
                    st_drop_geometry()%>%
                    mutate(Region=Points$Region)%>%
                    pivot_longer(cols=c(-Region), names_to="Date", values_to="Low")
                
                Med<-TempData()%>%
                    select(Prediction)%>%
                    aggregate(by=Points, function(x) length(which(!is.na(x) & x>=min(input$Habitat_range) & x <=max(input$Habitat_range))))%>%
                    st_as_sf(long=F)%>%
                    st_drop_geometry()%>%
                    mutate(Region=Points$Region)%>%
                    pivot_longer(cols=c(-Region), names_to="Date", values_to="Med")
                
                High<-TempData()%>%
                    select(Prediction)%>%
                    aggregate(by=Points, function(x) length(which(!is.na(x) & x>max(input$Habitat_range))))%>%
                    st_as_sf(long=F)%>%
                    st_drop_geometry()%>%
                    mutate(Region=Points$Region)%>%
                    pivot_longer(cols=c(-Region), names_to="Date", values_to="High")
                
                out<-full_join(Low, Med, by=c("Region", "Date"))%>%
                    full_join(High, by=c("Region", "Date"))%>%
                    mutate(Date=parse_date_time(Date, "%Y-%m-%d"))%>%
                    complete(Date=parse_date_time(Data_dates(), "%Y-%m-%d"), Region=unique(Points()$Region))%>%
                    mutate(Month=month(Date),
                           Year=year(Date))%>%
                    pivot_longer(cols = all_of(c("Low", "Med", "High")), names_to="Suitability", values_to="N_cells")%>%
                    mutate(Suitability=recode(Suitability, Low=input$Habitat_range_low, Med=input$Habitat_range_med, High=input$Habitat_range_high),
                           Suitability=factor(Suitability, levels=c(input$Habitat_range_high, input$Habitat_range_med, input$Habitat_range_low)))%>%
                    group_by(Date, Region)%>%
                    mutate(Total_cells=sum(N_cells))%>%
                    ungroup()%>%
                    mutate(Freq=if_else(Total_cells==0 | is.na(Total_cells), NA_real_, N_cells/Total_cells))
                return(out)
                
            }else{
                Low<-TempData()%>%
                    select(Prediction)%>%
                    aggregate(by=st_union(Points()), function(x) length(which(!is.na(x) & x<min(input$Habitat_range))))%>%
                    st_as_sf(as_points=T, long=T)%>%
                    st_drop_geometry()%>%
                    as_tibble()%>%
                    rename(Low=Prediction)
                Med<-TempData()%>%
                    select(Prediction)%>%
                    aggregate(by=st_union(Points()), function(x) length(which(!is.na(x) & x>=min(input$Habitat_range) & x <=max(input$Habitat_range))))%>%
                    st_as_sf(as_points=T, long=T)%>%
                    st_drop_geometry()%>%
                    as_tibble()%>%
                    rename(Med=Prediction)
                High<-TempData()%>%
                    select(Prediction)%>%
                    aggregate(by=st_union(Points()), function(x) length(which(!is.na(x) & x>max(input$Habitat_range))))%>%
                    st_as_sf(as_points=T, long=T)%>%
                    st_drop_geometry()%>%
                    as_tibble()%>%
                    rename(High=Prediction)
                
                out<-full_join(Low, Med, by="Date")%>%
                    full_join(High, by="Date")%>%
                    complete(Date=parse_date_time(Data_dates(), "%Y-%m-%d"))%>%
                    mutate(Month=month(Date),
                           Year=year(Date))%>%
                    pivot_longer(cols = all_of(c("Low", "Med", "High")), names_to="Suitability", values_to="N_cells")%>%
                    mutate(Suitability=recode(Suitability, Low=input$Habitat_range_low, Med=input$Habitat_range_med, High=input$Habitat_range_high),
                           Suitability=factor(Suitability, levels=c(input$Habitat_range_high, input$Habitat_range_med, input$Habitat_range_low)))%>%
                    group_by(Date)%>%
                    mutate(Total_cells=sum(N_cells))%>%
                    ungroup()%>%
                    mutate(Freq=if_else(Total_cells==0 | is.na(Total_cells), NA_real_, N_cells/Total_cells))
                return(out)
            }
            
        }else{
            
            if(input$Facet=="Region"){
                Points<-Points()%>%
                    group_by(Region)%>%
                    summarise(geometry=st_union(geometry), .groups="drop")
                TempData()%>%
                    select(Prediction)%>%
                    aggregate(by=Points, mean, na.rm=T)%>%
                    st_as_sf(long=F)%>%
                    st_drop_geometry()%>%
                    mutate(Region=Points$Region)%>%
                    pivot_longer(cols=c(-Region), names_to="Date", values_to="Prediction")%>%
                    left_join(TempData()%>%
                                  select(SE)%>%
                                  mutate(SE=SE^2)%>%
                                  aggregate(by=Points, function(x) sqrt(sum(x, na.rm=T)/(length(x)^2)))%>%
                                  st_as_sf(long=F)%>%
                                  st_drop_geometry()%>%
                                  mutate(Region=Points$Region)%>%
                                  pivot_longer(cols=c(-Region), names_to="Date", values_to="SE"), by=c("Region", "Date"))%>%
                    mutate(Prediction=if_else(is.nan(Prediction), NA_real_, Prediction),
                           SE=na_if(SE, 0)
                    )%>%
                    mutate(Date=parse_date_time(Date, "%Y-%m-%d"))%>%
                    complete(Date=parse_date_time(Data_dates(), "%Y-%m-%d"), Region=unique(Points()$Region))%>%
                    mutate(Month=month(Date),
                           Year=year(Date))
            } else{
                TempData()%>%
                    aggregate(by=st_union(Points()), mean, na.rm=T)%>%
                    st_as_sf(as_points=T, long=T)%>%
                    st_drop_geometry()%>%
                    as_tibble()%>%
                    mutate(across(c(Prediction, SE), ~if_else(is.nan(.x), NA_real_, .x)))%>%
                    mutate(Var=SE^2)%>%
                    complete(Date=parse_date_time(Data_dates(), "%Y-%m-%d"))%>%
                    mutate(Month=month(Date),
                           Year=year(Date),
                           SE=na_if(SE, 0))
                
            }
        }
    })
    
    timeseries_data_month<-reactive({
        req(input$Facet)
        if(input$Facet!="Month" & !input$Month2_slider){
            timeseries_data()%>%
                filter(Month%in%input$Month2)
        }else{
            timeseries_data()
        }
    })
    
    time_series_plot<-reactive({
        req(timeseries_data_month())
        
        if(input$Habitat){
            
            str_model <- paste0("<tr><td>Year: &nbsp</td><td>%s</td></tr>",
                                "<tr><td>Suitability: &nbsp</td><td>%s</td></tr>",
                                "<tr><td>Proportion of raster cells: &nbsp</td><td>%s</td></tr>")
            
            Data<-timeseries_data_month()%>%
                rowwise()%>%
                mutate(tooltip=sprintf(str_model, Year, Suitability, round(Freq, 2)))%>%
                ungroup()%>%
                mutate(ID=1:n(),
                       tooltip=paste0( "<table>", tooltip, "</table>" ))
            
            ggplot(Data, aes(y=Freq, fill=Suitability, order=Suitability, tooltip=tooltip, data_id=ID))+
                {if(input$Facet=="Month"){
                    geom_col_interactive(aes(x=Year))
                }else{
                    geom_col_interactive(aes(x=Date))
                }}+
                scale_fill_manual(values=c("#CE4257", "#FFE74C", "#111D4A"), drop=F, na.translate = FALSE, 
                                  guide = guide_legend(direction="horizontal", title.position = "top",
                                                       title.hjust=0.5, label.position="bottom", reverse=T))+
                {if(input$Facet=="Month"){
                    facet_wrap(~month(Date, label=T))
                }}+
                {if(input$Facet=="Region"){
                    if(input$Regions_all){
                        facet_geo(~Region, grid=mygrid, labeller=label_wrap_gen())
                    } else{
                        facet_wrap(~Region)
                    }
                }}+
                coord_cartesian(expand = FALSE)+
                ylab("Proportion of raster cells")+
                theme_bw()+
                theme(strip.background=element_blank(), panel.grid=element_blank(),
                      legend.position="top", legend.key.width = unit(50, "native"), 
                      legend.justification="center", axis.text.x=element_text(angle=45, hjust=1))+
                {if(input$Facet=="Region" & input$Regions_all){
                    theme(text=element_text(size=4), panel.spacing = unit(0.2, units="lines"),
                          axis.ticks=element_line(size=0.1), axis.ticks.length = unit(0.2, units="lines"))
                } else{
                    theme(text=element_text(size=14))
                }}
            
        }else{
            
            str_model <- paste0("<tr><td>Mean: &nbsp</td><td>%s</td></tr>",
                                "<tr><td>Lower SE: &nbsp</td><td>%s</td></tr>", 
                                "<tr><td>Upper SE: &nbsp</td><td>%s</td></tr>")
            
            Data<-timeseries_data_month()%>%
                rowwise()%>%
                mutate(tooltip=sprintf(str_model, round(Prediction, 2), round(Prediction-SE, 2), round(Prediction+SE, 2)))%>%
                ungroup()%>%
                mutate(ID=1:n(),
                       tooltip=paste0( "<table>", tooltip, "</table>" ))%>%
                {if(input$Facet=="Month"){
                    mutate(., Group=Month)
                } else{
                    if(input$Facet=="Region"){
                        mutate(., Group=Region)
                    }else{
                        mutate(., Group=1)
                    }
                }}
            
            ggplot(Data, aes(x=Date, y=Prediction, ymin=Prediction-SE, ymax=Prediction+SE, group=Group))+
                geom_ribbon(alpha=0.4, fill="firebrick3")+
                geom_line(color="firebrick3", 
                          size=if_else(input$Facet=="Region" & input$Regions_all, 0.3, 0.5))+
                geom_pointrange_interactive(aes(tooltip=tooltip, data_id=ID), color="firebrick3", alpha=0.5, 
                                            size=if_else(input$Facet=="Region" & input$Regions_all, 0.001, 0.2), 
                                            shape=ifelse(input$Facet=="Region" & input$Regions_all, ".", 16))+
                {if(input$Facet=="Month"){
                    facet_wrap(~month(Date, label=T))
                }}+
                {if(input$Facet=="Region"){
                    if(input$Regions_all){
                        facet_geo(~Region, grid=mygrid, labeller=label_wrap_gen())
                    } else{
                        facet_wrap(~Region)
                    }
                }}+
                ylab("Temperature ± SE (°C)")+
                theme_bw()+
                theme(strip.background=element_blank())+
                {if(input$Facet=="Region" & input$Regions_all){
                    theme(text=element_text(size=4), panel.spacing = unit(0.2, units="lines"),
                          axis.ticks=element_line(size=0.1), axis.ticks.length = unit(0.2, units="lines"))
                } else{
                    theme(text=element_text(size=18))
                }}
        }
        
    })
    
    output$Time_plot<-renderGirafe({
        p<-girafe(ggobj=time_series_plot(), width_svg = if_else(input$Facet=="Region" & input$Regions_all, 5.5, 10),  
                  pointsize=if_else(input$Facet=="Region" & input$Regions_all, 6, 12), 
                  options = list(
                      opts_sizing(rescale = T, width = 1)))
        girafe_options(p, opts_toolbar(saveaspng = FALSE), opts_selection(type="none"))
    })
    
    
    # Raw data ----------------------------------------------------------------
    
    Data_filtered<-reactive({
        req(input$Date_range, any(input$Months%in%1:12))
        Data%>%
            filter(Date>min(input$Date_range) & Date<max(input$Date_range) & month(Date)%in%input$Months)
    })
    
    
    all_stations<-reactive({
        Data_filtered()%>%
            filter(Source!="EDSM" & !str_detect(Station, "EZ") & !StationID%in%c("SKT 799", "SKT 999", "SKT 699"))%>%
            st_drop_geometry()%>%
            group_by(StationID, Latitude, Longitude)%>%
            summarise(N=n(), .groups="drop")%>%
            distinct()%>%
            st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)
        
    })
    
    all_stations_plot<-reactive({
        all_stations()%>%
            {if(is.null(Stations())){
                mutate(., opac=0.2)
            }else{
                mutate(., opac=if_else(StationID%in%Stations()$StationID, 0.9, 0.2))
            }}
    })
    
    pal_N3<-reactive({
        colorNumeric("viridis", domain=range(log(all_stations()$N), na.rm=T), na.color="#00000000")
    })
    
    pal_N3_rev<-reactive({
        colorNumeric("viridis", domain=range(log(all_stations()$N), na.rm=T), reverse=T, na.color="#00000000")
    })
    
    #set the namespace for the map
    ns_Rawmap <- shiny::NS("Rawmap")
    
    
    
    #call the editMod function from 'mapedit' to use in the leaflet map.
    
    edits_Rawmap <- callModule(editMod, "Rawmap", leafmap = Station_plot, 
                               editorOptions=list(polylineOptions=F, circleMarkerOptions=F, markerOptions=F))
    
    observeEvent(all_stations_plot(), {  
        proxy.points <- leafletProxy(ns_Rawmap("map"))
        
        proxy.points %>%
            clearShapes()%>%
            clearControls()%>%
            addFeatures(data=all_stations_plot(), fillColor=~pal_N3()(log(N)), color="black", fillOpacity = ~opac, 
                        label=lapply(paste0("<h5 align='center'>", all_stations_plot()$StationID, "</h5>", "<h6 align='center' style='color:red'>N: ", 
                                            format(all_stations_plot()$N, big.mark   = ","), "</h6>"), htmltools::HTML), 
                        layerId = ~StationID, weight=~opac*2)%>%
            addLegend(data=all_stations_plot(), position="topright", pal = pal_N3_rev(), values = ~log(N), opacity=1, 
                      labFormat = labelFormat(transform = function(x) sort(round(exp(x)), decreasing = TRUE)), title="Number of raw</br>temperature records</br>(log scale)")
    }) 
    
    # Extract and plot raw timeseries data ----------------------------------------
    
    Stations<-reactive({
        req(!(is.null(edits_Rawmap())))
        if(is.null(edits_Rawmap()$all)){
            return(NULL)
        }
        
        Area<-edits_Rawmap()$all%>%
            st_transform(crs=26910)%>%
            st_union()
        
        Selected_stations<-Data_filtered()%>%
            st_transform(crs=26910)%>%
            st_filter(Area)%>%
            st_drop_geometry()%>%
            pull(StationID)%>%
            unique()
        
        if(length(Selected_stations)>0){
            out<-filter(Data_filtered(), StationID%in%Selected_stations)%>%
                st_transform(crs=4326)
        }else(
            out<-NULL
        )
        
        
        return(out)
    })  
    
    raw_data_plot<-reactive({
        req(input$Facet, input$Month3, Stations())
        
        if(input$Nonfixed){
            Data<-Stations()
        }else{
            Data<-Stations()%>%
                filter(Source!="EDSM" & !str_detect(Station, "EZ"))
        }
        
        if(input$Facet!="Month" & !input$Month3_slider){
            Data<-Data%>%
                filter(Month%in%input$Month3)
        }
        
        if(input$Time_correction){
            Data<-Data%>%
                mutate(Time=as.character(round(Time_num_s, 1)))%>%
                left_join(Time_correction, by=c("Time", "Month"))%>%
                mutate(Temperature=Temperature+Correction)
        }
        return(Data)
    })
    
    raw_time_series_plot<-reactive({
        req(raw_data_plot())
        
        str_model <- paste0("<tr><td>Station: &nbsp</td><td>%s</td></tr>",
                            "<tr><td>Date: &nbsp</td><td>%s</td></tr>",
                            "<tr><td>Time: &nbsp</td><td>%s</td></tr>", 
                            "<tr><td>Temperature: &nbsp</td><td>%s</td></tr>")
        
        Data<-raw_data_plot()%>%
            rowwise()%>%
            mutate(tooltip=sprintf(str_model, Station, format(Date, format="%b %d, %Y"), format(Datetime, format="%I:%M %p"), round(Temperature, 2)))%>%
            ungroup()%>%
            mutate(ID=1:n(),
                   tooltip=paste0( "<table>", tooltip, "</table>" ))
        
        ggplot(Data, aes(x=Date, y=Temperature, group=Station, fill=Station, color=Station))+
            geom_line_interactive(size=0.5, aes(data_id=Station, tooltip=Station,
                                                hover_css = "fill:none; stroke:gold; stroke-width:3px;"))+
            geom_point_interactive(aes(tooltip=tooltip, data_id=ID,
                                       hover_css = "stroke-width:3px; fill:gold; stroke:gold;"), alpha=0.5, 
                                   shape=21)+
            {if(input$Facet=="Month"){
                facet_wrap(~month(Month, label=T))
            }}+
            ylab("Temperature (°C)")+
            theme_bw()+
            theme(strip.background=element_blank(), text=element_text(size=18), axis.text.x=element_text(angle=45, hjust=1))
        
    })
    
    output$Raw_time_plot<-renderGirafe({
        p<-girafe(ggobj=raw_time_series_plot(), width_svg = 10, pointsize=12, 
                  options = list(
                      opts_sizing(rescale = T, width = 1)))
        girafe_options(p, opts_toolbar(saveaspng = FALSE), opts_selection(type="none"))
    })
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
