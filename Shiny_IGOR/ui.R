library(shiny)
library(shinycssloaders)
library(shinyjs)
library(ggplot2)
library(RColorBrewer)
library(plotly)
library(DT)
library(reshape2)
library(TMB)

appCSS <- "
#loading-content {
position: absolute;
background: #c5e0d7;
opacity: 0.9;
z-index: 100;
left: 0;
right: 0;
height: 100%;
text-align: center;
color: #000000;
}
"

shinyUI(fluidPage(
  useShinyjs(),
  inlineCSS(appCSS),
  
  # Loading message
  div(
    id = "loading-content",
    h2("Loading IGOR...This may take several minutes")
  ),
  
  # The main app code goes here
  hidden(
    div(
      id = "app-content",
      titlePanel("IGOR+: Fitting Growth Curves with Random Effects"),
      
      fluidRow(
        column(1),
        column(10,
               tabsetPanel(id = "menu",
                 tabPanel("Input Data", icon = icon("upload", lib = "glyphicon"),
                          fluidRow(
                            column(4, fileInput("file", "Choose CSV File for your data",
                                                multiple = FALSE,
                                                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")
                            )),
                            column(4, downloadButton("downloadData", "Download Sample Data File", style = "margin-top:25px"))
                          ),
                          fluidRow(
                            tabsetPanel(id = "data",
                              tabPanel("Input Data",
                                uiOutput("input_control"),
                                DTOutput("input_data")
                              ),
                              tabPanel("Selected Data",
                                p("When you are ready, go to the", 
                                  actionButton("go_analyze", "Analyze", icon = icon("hourglass", lib = "glyphicon")), "tab", 
                                  style = "margin-top:25px; font-size:15px"),
                                withSpinner(plotlyOutput("selected_data_plot")),
                                DTOutput("selected_data")
                              )
                            )
                          )
                 ),
                 tabPanel("Analyze", icon = icon("hourglass", lib = "glyphicon"),
                          fluidRow(
                            p("Fit your data with either a standard nonlinear fit or a random effects model. 
                              Check your results in the", actionButton("go_summaries", "Summaries", icon("edit", lib = "glyphicon")), "tab",
                              style = "margin-top:25px; font-size:15px"),
                            tabsetPanel(id = "fit_choice",
                                        tabPanel("Standard Non-Linear Fit",
                                                 column(4, 
                                                        h3("Standard Non-Linear Fit"),
                                                        wellPanel(
                                                          selectInput("model", "Growth Model", choices = c("Von Bertalanffy")),
                                                          fluidRow(
                                                            column(4, textInput("linf_nls", "Linf", 912)),
                                                            column(4, textInput("kappa_nls", "K", 0.04)),
                                                            column(4, textInput("t0_nls", "t0", 0))
                                                          ),
                                                          selectInput("fit_data", "Fit to", c("1 - Selected Read", "2 - Average Age", "3 - Median Age", 
                                                                                              "4 - constant Length CV", "5 - Lengths with Multiple Ages")),
                                                          uiOutput("selected_read"),
                                                          conditionalPanel(
                                                            condition = "input.fit_data == '4 - constant Length CV'",
                                                            textInput("CV_const", "CV", 0.1)
                                                          ),
                                                          textInput("runname_nls", "Run Name", "Default"),
                                                          actionButton("run_nls", "Run A Standard Non-Linear Estimation",
                                                                       icon("stats", lib = "glyphicon")))
                                                 ),
                                                 column(8, withSpinner(plotlyOutput("ModelNLS_Plot")))
                                        ),
                                        
                                        tabPanel("Random Effects Model",
                                                 column(4, 
                                                        h3("Random Effects Model"),
                                                        wellPanel(
                                                          selectInput("model", "Growth Model", choices = c("Von Bertalanffy")),
                                                          fluidRow(
                                                            column(4, textInput("linf_re", "Linf", 912)),
                                                            column(4, textInput("kappa_re", "K", 0.04)),
                                                            column(4, textInput("t0_re", "t0", 0))
                                                          ),
                                                          fluidRow(
                                                            column(6, textInput("CV_e", "CV of Random Effect", -1),
                                                                   helpText("If negative, internally estimates")),
                                                            column(6, textInput("CV_Lt", "Length Standard", 0.1))
                                                          ),
                                                          h5("Age Likelihood Type for Random Effects"),
                                                          tabsetPanel(id = "likelihoods",
                                                                      tabPanel("Normal", 
                                                                               textInput("alpha", "Mean Population Age", 5),
                                                                               textInput("sigma_age", "StDev of Population", 1)
                                                                      ),
                                                                      tabPanel("Exponential", 
                                                                               textInput("beta", "Rate of Population", 0.2)
                                                                      ),
                                                                      tabPanel("Gamma",
                                                                               textInput("gam_shape", "Shape", 5),
                                                                               textInput("gam_scale", "Scale", 3)
                                                                      )
                                                          ),
                                                          textInput("runname_re", "Run Name", "Default"),
                                                          actionButton("run_re", "Run A Random Effects Estimation",
                                                                       icon("stats", lib = "glyphicon"))
                                                        )
                                                 ),
                                                 column(8,
                                                        fluidRow(style = "margin-top:25px",
                                                                 tabsetPanel(
                                                                   tabPanel("Plot", withSpinner(plotlyOutput("ModelRE_Plot"))),
                                                                   tabPanel("Estimate Z",
                                                                            p("Estimate the Z value based on the age random effects. Typically, choose the
                                                                              coordinates of the tallest bar in the histogram as one of the end points. Click
                                                                              Choose Start/Choose End to select coordinates of the histogram bars or manually
                                                                              enter values for the end points."),
                                                                            wellPanel(
                                                                              fluidRow(
                                                                                column(9, 
                                                                                       fluidRow(
                                                                                         column(3, actionButton("get_start", "Choose Start", style = "margin-top:25px")),
                                                                                         column(4, textInput("start_x", "Start X Coordinate")),
                                                                                         column(4, textInput("start_y", "Start Y Coordinate"))
                                                                                       ),
                                                                                       fluidRow(
                                                                                         column(3, actionButton("get_end", "Choose End", style = "margin-top:25px")),
                                                                                         column(4, textInput("end_x", "End X Coordinate")),
                                                                                         column(4, textInput("end_y", "End Y Coordinate"))
                                                                                       )
                                                                                ),
                                                                                column(3, 
                                                                                       actionButton("get_z", "Get Z Value", style = "margin-top:25px"),
                                                                                       textOutput("z_value")
                                                                                )
                                                                              ),
                                                                              sliderInput(inputId = "bins",
                                                                                          label = "Number of bins:",
                                                                                          min = 1,
                                                                                          max = 50,
                                                                                          value = 30)
                                                                            ),
                                                                            withSpinner(plotlyOutput("hist_RE"))  
                                                                   )
                                                                 )
                                                        )
                                                        
                                                )
                                        )
                            )
                          )
                 ),
                 tabPanel("Summaries", icon = icon("edit", lib = "glyphicon"),
                   DTOutput("summaries")   
                 )
               )
        ),
        column(1)
      )
    )
  )
))