library(shiny)
library(shinycssloaders)
library(shinyjs)
library(ggplot2)
library(RColorBrewer)
library(plotly)
library(DT)
library(reshape2)
library(TMB)

compile("vbre_Gamma.cpp")
compile("vbre_Exponential.cpp")
compile("vbre_Normal.cpp")
compile("vb_likelihood.cpp")
compile("linear_likelihood.cpp")
compile("schnute_likelihood.cpp")
compile("logistic_likelihood.cpp")
compile("gompertz_likelihood.cpp")
compile("linear_Normal.cpp")
compile("linear_Exponential.cpp")
compile("linear_Gamma.cpp")
compile("gompertz_Normal.cpp")
compile("gompertz_Exponential.cpp")
compile("gompertz_Gamma.cpp")
compile("schnute_Normal.cpp")
compile("schnute_Exponential.cpp")
compile("schnute_Gamma.cpp")
compile("logistic_Normal.cpp")
compile("logistic_Exponential.cpp")
compile("logistic_Gamma.cpp")
dyn.load(dynlib("vbre_Gamma"))
dyn.load(dynlib("vbre_Exponential"))
dyn.load(dynlib("vbre_Normal"))
dyn.load(dynlib("vb_likelihood"))
dyn.load(dynlib("linear_likelihood"))
dyn.load(dynlib("schnute_likelihood"))
dyn.load(dynlib("logistic_likelihood"))
dyn.load(dynlib("gompertz_likelihood"))
dyn.load(dynlib("logistic_Normal"))
dyn.load(dynlib("logistic_Exponential"))
dyn.load(dynlib("logistic_Gamma"))
dyn.load(dynlib("schnute_Normal"))
dyn.load(dynlib("schnute_Exponential"))
dyn.load(dynlib("schnute_Gamma"))
dyn.load(dynlib("linear_Normal"))
dyn.load(dynlib("linear_Exponential"))
dyn.load(dynlib("linear_Gamma"))
dyn.load(dynlib("gompertz_Normal"))
dyn.load(dynlib("gompertz_Exponential"))
dyn.load(dynlib("gompertz_Gamma"))


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
  withMathJax(),
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
                            p("Fit your data with either a standard fit or a random effects model. 
                              Check your results in the", actionButton("go_summaries", "Summaries", icon("edit", lib = "glyphicon")), "tab",
                              style = "margin-top:25px; font-size:15px"),
                            tabsetPanel(id = "fit_choice",
                                        tabPanel("Standard Fit",
                                                 column(4, 
                                                        h3("Standard Fit"),
                                                        wellPanel(
                                                          selectInput("model_nls", "Growth Model", 
                                                                      choices = c("Linear", "Logistic", "Von Bertalanffy", "Schnute", "Gompertz")),
                                                          conditionalPanel(
                                                            condition = "input.model_nls == 'Von Bertalanffy'",
                                                            helpText('$$y=linf\\cdot[1-e^{-K\\cdot(x-t0)}]$$'),
                                                            fluidRow(
                                                              column(4, textInput("linf_nls", "Linf (Upper Asymptote)", 912)),
                                                              column(4, textInput("kappa_nls", "K (Growth Rate)", 0.04)),
                                                              column(4, textInput("t0_nls", "t0", 0))
                                                            )
                                                          ),
                                                          conditionalPanel(
                                                            condition = "input.model_nls == 'Gompertz'",
                                                            helpText('$$y=a\\cdot e^{-b\\cdot e^{-kx}}$$'),
                                                            fluidRow(
                                                              column(4, textInput("gom_a_nls", "a (Upper Asymptote)")),
                                                              column(4, textInput("gom_b_nls", "b (Growth Displacement)")),
                                                              column(4, textInput("gom_k_nls", "k (Growth Rate)"))
                                                            )
                                                          ),
                                                          conditionalPanel(
                                                            condition = "input.model_nls == 'Schnute'",
                                                            helpText('$$y=(r0+b\\cdot e^{kx})^m$$'),
                                                            fluidRow(
                                                              column(3, textInput("sch_r0_nls", "r0 (Reference Value")),
                                                              column(3, textInput("sch_b_nls", "b (Growth Range)")),
                                                              column(3, textInput("sch_k_nls", "k (Growth Rate)")),
                                                              column(3, textInput("sch_m_nls", "m (Slope of Growth)"))
                                                            )
                                                          ),
                                                          conditionalPanel(
                                                            condition = "input.model_nls == 'Logistic'",
                                                            helpText('$$y=\\dfrac{a}{1+b\\cdot e^{-kx}}$$'),
                                                            fluidRow(
                                                              column(4, textInput("log_a_nls", "a (Upper Asymptote)")),
                                                              column(4, textInput("log_b_nls", "b (Growth Range)")),
                                                              column(4, textInput("log_k_nls", "k (Growth Rate)"))
                                                            )
                                                          ),
                                                          conditionalPanel(
                                                            condition = "input.model_nls == 'Linear'",
                                                            helpText('$$y=Intercept + Slope\\cdot x$$')
                                                          ),
                                                          conditionalPanel(
                                                            condition = "input.model_nls == 'Linear' && input.fit_type == 'Constant Length CV'",
                                                            fluidRow(
                                                              column(4, textInput("intercept_nls", "Intercept (Length at Age 0)")),
                                                              column(4, textInput("slope_nls", "Slope (Growth Rate)"))
                                                            )
                                                          ),
                                                          radioButtons("fit_type", "Fit type", c("Standard Deviation", "Constant Length CV")),
                                                          selectInput("fit_data", "Fit to", c("1 - Selected Read", "2 - Average Age", "3 - Median Age", 
                                                                                              "4 - Lengths with Multiple Ages")),
                                                          uiOutput("selected_read"),
                                                          conditionalPanel(
                                                            condition = "input.fit_type == 'Constant Length CV'",
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
                                                          selectInput("model_re", "Growth Model", 
                                                                      choices = c("Linear", "Logistic", "Von Bertalanffy", "Schnute", "Gompertz")),
                                                          conditionalPanel(
                                                            condition = "input.model_re == 'Von Bertalanffy'",
                                                            helpText('$$y=linf\\cdot[1-e^{-K\\cdot(x-t0)}]$$'),
                                                            fluidRow(
                                                              column(4, textInput("linf_re", "Linf (Upper Asymptote)", 912)),
                                                              column(4, textInput("kappa_re", "K (Growth Rate)", 0.04)),
                                                              column(4, textInput("t0_re", "t0", 0))
                                                            )
                                                          ),
                                                          conditionalPanel(
                                                            condition = "input.model_re == 'Gompertz'",
                                                            helpText('$$y=a\\cdot e^{-b\\cdot e^{-kx}}$$'),
                                                            fluidRow(
                                                              column(4, textInput("gom_a_re", "a (Upper Asymptote)")),
                                                              column(4, textInput("gom_b_re", "b (Growth Displacement)")),
                                                              column(4, textInput("gom_k_re", "k (Growth Rate)"))
                                                            )
                                                          ),
                                                          conditionalPanel(
                                                            condition = "input.model_re == 'Schnute'",
                                                            helpText('$$y=(r0+b\\cdot e^{kx})^m$$'),
                                                            fluidRow(
                                                              column(3, textInput("sch_r0_re", "r0 (Reference Value")),
                                                              column(3, textInput("sch_b_re", "b (Growth Displacement)")),
                                                              column(3, textInput("sch_k_re", "k (Growth Rate)")),
                                                              column(3, textInput("sch_m_re", "m (Slope of Growth)"))
                                                            )
                                                          ),
                                                          conditionalPanel(
                                                            condition = "input.model_re == 'Logistic'",
                                                            helpText('$$y=\\dfrac{a}{1+b\\cdot e^{-kx}}$$'),
                                                            fluidRow(
                                                              column(4, textInput("log_a_re", "a (Upper Asymptote)")),
                                                              column(4, textInput("log_b_re", "b (Growth Range)")),
                                                              column(4, textInput("log_k_re", "k (Growth Rate)"))
                                                            )
                                                          ),
                                                          conditionalPanel(
                                                            condition = "input.model_re == 'Linear'",
                                                            helpText('$$y=Intercept + Slope\\cdot x$$'),
                                                            fluidRow(
                                                              column(4, textInput("intercept_re", "Intercept (Length at Age 0)")),
                                                              column(4, textInput("slope_re", "Slope (Growth Rate)"))
                                                            )
                                                          ),
                                                          fluidRow(
                                                            column(6, textInput("CV_e", "CV of Random Effect", -1),
                                                                   helpText("If negative, internally calculated")),
                                                            column(6, textInput("CV_Lt", "Length Standard", 0.1))
                                                          ),
                                                          h5("Age Likelihood Type for Random Effects"),
                                                          tabsetPanel(id = "likelihoods", selected = "Gamma",
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
                                                                                       helpText("Get Z Value and update the most", em("recent RE model"), "summary"),
                                                                                       actionButton("get_z", "Get Z"),
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
                   h2("Von Bertlanffy Model Summaries"),
                   DTOutput("vb_summaries"),
                   h2("Linear Model Summaries"),
                   DTOutput("linear_summaries"),
                   h2("Gompertz Model Summaries"),
                   DTOutput("gompertz_summaries"),
                   h2("Logistic Model Summaries"),
                   DTOutput("logistic_summaries"),
                   h2("Schnute Model Summaries"),
                   DTOutput("schnute_summaries")
                 )
               )
        ),
        column(1)
      )
    )
  )
))