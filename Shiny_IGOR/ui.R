# shiny ui
# set spinner color to green
options(spinner.color = "#5cb935")
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
  # Include IntroJS styling
  includeCSS("www/introjs.min.css"),
  
  # Include IntroJS library
  includeScript("www/intro.min.js"),
  
  # Include javascript code to make shiny communicate with introJS
  includeScript("www/app.js"),
  
  includeCSS("www/app.css"),
  
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
      
      
      fluidRow(
        column(1),
        column(10,
               titlePanel("IGOR+: Fitting Growth Curves with Random Effects"),
               tabsetPanel(id = "menu",
                           tabPanel("File Upload", icon = icon("upload", lib = "glyphicon"), value = "tab1",
                                    h3("Welcome!", actionButton(inputId = "startHelp", label = "Start A Guide", class = "btn-success")),
                                    fluidRow(
                                      column(6,
                                             wellPanel(id = "oto_age_upload", style = "height:450px",
                                                       h4(strong("Predicting Ages From Otolith Weights")),
                                                       fileInput("oto_age_file", "Choose CSV File for your data",
                                                                 multiple = FALSE,
                                                                 accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                                                       downloadButton("download_oto_age_data", "Download Sample Data File"),
                                                       helpText("If you want to predict/create", em("one column"), "of age data from a regression relation between
                                 otolith weight and age (you may want to do this because many of the age data is not available),
                                 upload a file here. We will help you generate a data table that is later used for growth analysis.
                                 Refer to the sample file for the data format."),
                                                       p("When you are ready, go to the",
                                                         actionButton("go_oto_age", "Age vs. Otolith Weight Data",
                                                                      icon = icon("list", lib = "glyphicon")), "tab",
                                                         style = "margin-top:25px; font-size:15px")
                                             )
                                      ),
                                      column(6,
                                             wellPanel(id = "age_len_upload", style = "height:450px",
                                                       h4(strong("Fitting Growth Models")),
                                                       fileInput("age_length_file", "Choose CSV File for your data",
                                                                 multiple = FALSE,
                                                                 accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                                                       downloadButton("download_age_length_data", "Download Sample Data File"),
                                                       helpText("If you have age and length data ready for growth analysis, upload a file here.
                                        Refer to the sample file for the data format."),
                                                       p("When you are ready, go to the",
                                                         actionButton("go_age_length", "Length vs. Age Data",
                                                                      icon = icon("list", lib = "glyphicon")), "tab",
                                                         style = "margin-top:25px; font-size:15px")
                                             )
                                      )
                                    ),
                                    p("In order for IGOR to work, you must have all the required columns (not necessarily ordered) as in the sample files")
                           ),
                           tabPanel("Age vs. Otolith Weight Data", icon = icon("list", lib = "glyphicon"), value = "tab2",
                                    p("This analysis is not necessary if you have the age data prepared.",
                                      style = "margin-top:25px; font-size:15px"),
                                    tabsetPanel(id = "tab2_subpanel",
                                                tabPanel("Run A New Otolith Weight vs. Age Model", value = "tab2-1",
                                                         fluidRow(
                                                           column(4, 
                                                                  wellPanel(id = "oto_age_controls", style = "margin-top:25px",
                                                                            selectInput("oto_age_model_type", "Model Type", 
                                                                                        choices = c("Linear", "Quadratic", "1-Breakpoint Piecewise", "2-Breakpoints Piecewise")),
                                                                            conditionalPanel(
                                                                              condition = "input.oto_age_model_type == '1-Breakpoint Piecewise'",
                                                                              textInput("wt_bkpt1", "Otolith Weight Breakpoint 1")
                                                                            ),
                                                                            conditionalPanel(
                                                                              condition = "input.oto_age_model_type == '2-Breakpoints Piecewise'",
                                                                              fluidRow(
                                                                                column(6, textInput("wt_bkpt1_2bp", "Otolith Weight Breakpoint 1")),
                                                                                column(6, textInput("wt_bkpt2_2bp", "Otolith Weight Breakpoint 2"))
                                                                              )
                                                                            ),
                                                                            conditionalPanel(
                                                                              condition = "input.oto_age_model_type == '2-Breakpoints Piecewise' || 
                                                  input.oto_age_model_type == '1-Breakpoint Piecewise'",
                                                                              checkboxInput("use_quantile", "Use quantile to set breakpoints"),
                                                                              helpText("When quantile is checked, enter number(s) between 0 to 1 to set weight breakpoint(s)"),
                                                                              checkboxInput("jitter", "Jitter"),
                                                                              conditionalPanel(
                                                                                condition = "input.jitter == true",
                                                                                numericInput("num_jitter", "Number of times", 5, min = 1, max = 100)
                                                                              )
                                                                            ),
                                                                            fluidRow(
                                                                              column(6, textInput("wt_unit", "Otolith Weight Unit", "mg")),
                                                                              column(6, textInput("age_unit_oto", "Age Unit", "yr"))
                                                                            ),
                                                                            actionButton("run_oto_age", "Predict Ages", icon("stats", lib = "glyphicon"))
                                                                  )
                                                           ),
                                                           column(8, 
                                                                  fluidRow(style = "margin-top: 25px",
                                                                           tabsetPanel(id = "oto_age_menu",
                                                                                       tabPanel("Input Data", value = "tab2-1-1",
                                                                                                DTOutput("oto_age_data")
                                                                                       ),
                                                                                       tabPanel("Summary & Plot", value = "tab2-1-2",
                                                                                                p("You can", downloadButton("download_oto_age_model", "Download This Model"), "for future use.",
                                                                                                  style = "margin-top:25px"),
                                                                                                tableOutput("oto_age_summary"),
                                                                                                tableOutput("brkpt_summary"),
                                                                                                fluidRow(
                                                                                                  column(6, tableOutput("intercept_summary")),
                                                                                                  column(6, tableOutput("slope_summary"))
                                                                                                ),
                                                                                                helpText("The plot legend label \"Selected\" represents the filtered data points and \"Not Selected\"
                                                                                    represents the data left out."),
                                                                                                withSpinner(uiOutput("oto_age_plot")),
                                                                                                # these 6 plots only available for 1-breakpoint jitter model
                                                                                                fluidRow(
                                                                                                  column(6, withSpinner(plotlyOutput("jitter_intercept_plot"))),
                                                                                                  column(6, withSpinner(plotlyOutput("jitter_breakpoint_plot")))
                                                                                                ),
                                                                                                fluidRow(
                                                                                                  column(6, withSpinner(plotlyOutput("jitter_slope1_plot"))),
                                                                                                  column(6, withSpinner(plotlyOutput("jitter_slope2_plot")))
                                                                                                ),
                                                                                                fluidRow(
                                                                                                  column(6, withSpinner(plotlyOutput("jitter_AIC_plot"))),
                                                                                                  column(6, withSpinner(plotlyOutput("jitter_comp_plot")))
                                                                                                )
                                                                                       ),
                                                                                       tabPanel("Fitted Values", value = "tab2-1-3",
                                                                                                div(id = "step1-3",
                                                                                                    p("You can view the predicted ages from the chosen model. If this table looks good, 
                                                                                 you may want to", actionButton("set_age_length_data", 
                                                                                                                "Use Current Model to Predict Ages for Growth Fitting",
                                                                                                                icon = icon("wrench", lib = "glyphicon")), ".",
                                                                                                      style = "margin-top:25px"
                                                                                                    ),
                                                                                                    DTOutput("predicted_age_preview")
                                                                                                )
                                                                                       )
                                                                           )
                                                                  )
                                                           )
                                                         )       
                                                ),
                                                tabPanel("Use An Existing Otolith Weight vs. Age Model", value = "tab2-2",
                                                         fluidRow(
                                                           column(4,
                                                                  wellPanel(id = "oto_model_upload",style = "margin-top:25px",
                                                                            h4(strong("Using Existing Otolith-Age Relationaship To Predict Ages")),
                                                                            fileInput("oto_file", "Choose CSV File for your data",
                                                                                      multiple = FALSE,
                                                                                      accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                                                                            downloadButton("download_oto_data", "Download Sample Data File", style = "margin-bottom:25px"),
                                                                            fileInput("model_oto_file", "Choose RDS File for your model",
                                                                                      multiple = FALSE,
                                                                                      accept = c(".rds")),
                                                                            downloadButton("download_oto_model", "Download Sample Model File"),
                                                                            helpText("If you have an rds object that models the age and otolith weight relationship and want to predict
                                                              ages for fish with known otolith weight, upload two files here.
                                                              Refer to the sample files for the data format."),
                                                                            actionButton("run_existing_oto_model", "Predict Ages with Existing Model", icon("stats", lib = "glyphicon"))
                                                                  )
                                                           ),
                                                           column(8,
                                                                  fluidRow(style = "margin-top: 25px",
                                                                           tabsetPanel(
                                                                             tabPanel("Fitted Values",
                                                                                      div(id = "step1-5",
                                                                                          p("You can view the predicted ages from your uploaded model. If this table looks good, 
                                                                 you may want to", actionButton("set_age_length_data_from_existing_model", 
                                                                                                "Use Current Model to Predict Ages for Growth Fitting",
                                                                                                icon = icon("wrench", lib = "glyphicon")), ".",
                                                                                            style = "margin-top:25px"
                                                                                          ),
                                                                                          DTOutput("existing_model_predicted_data")
                                                                                      )
                                                                             )
                                                                           )
                                                                  )
                                                           )
                                                         )     
                                                )
                                    )
                           ),
                           tabPanel("Length vs. Age Data", icon = icon("list", lib = "glyphicon"), value = "tab3",
                                    fluidRow(
                                      p("Before you start a growth analysis, you can filter and preview your data here.",
                                        style = "margin-top:25px; font-size:15px"),
                                      tabsetPanel(id = "data",
                                                  tabPanel("Input Data", value = "tab3-1",
                                                           div(id = "step2-2",
                                                               uiOutput("input_control"),
                                                               DTOutput("input_data")
                                                           )
                                                  ),
                                                  tabPanel("Selected Data", value = "tab3-2",
                                                           div(id = "step2-3",
                                                               p("When you are ready, go to the", 
                                                                 actionButton("go_analyze", "Growth Curve Analysis", icon = icon("hourglass", lib = "glyphicon")), "tab", 
                                                                 style = "margin-top:25px; font-size:15px")  
                                                           ),
                                                           withSpinner(plotlyOutput("selected_data_plot")),
                                                           DTOutput("selected_data")
                                                  )
                                      )
                                    )
                           ),
                           tabPanel("Growth Curve Analysis", icon = icon("hourglass", lib = "glyphicon"), value = "tab4", 
                                    fluidRow(
                                      p("Fit your data with either a standard fit or a random effects model. 
                              Check your results in the", actionButton("go_summaries", "Growth Curve Summaries", icon("edit", lib = "glyphicon")), "tab",
                                        style = "margin-top:25px; font-size:15px"),
                                      tabsetPanel(id = "fit_choice",
                                                  tabPanel("Standard Fit", value = "tab4-1",
                                                           column(4, 
                                                                  h3("Standard Fit"),
                                                                  wellPanel(id = "age_len_std_control",
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
                                                                            fluidRow(
                                                                              column(6, textInput("len_unit_nls", "Length Unit", "mm")),
                                                                              column(6, textInput("age_unit_nls", "Age Unit", "yr"))
                                                                            ),
                                                                            textInput("runname_nls", "Run Name", "Default"),
                                                                            actionButton("run_nls", "Run A Standard Non-Linear Estimation",
                                                                                         icon("stats", lib = "glyphicon")))
                                                           ),
                                                           column(8, withSpinner(plotlyOutput("ModelNLS_Plot")))
                                                  ),
                                                  
                                                  tabPanel("Random Effects Model", value = "tab4-2",
                                                           column(4, 
                                                                  h3("Random Effects Model"),
                                                                  wellPanel(id = "age_len_re_control",
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
                                                                            fluidRow(
                                                                              column(6, textInput("len_unit_re", "Length Unit", "mm")),
                                                                              column(6, textInput("age_unit_re", "Age Unit", "yr"))
                                                                            ),
                                                                            textInput("runname_re", "Run Name", "Default"),
                                                                            actionButton("run_re", "Run A Random Effects Estimation",
                                                                                         icon("stats", lib = "glyphicon"))
                                                                  )
                                                           ),
                                                           column(8,
                                                                  fluidRow(style = "margin-top:25px",
                                                                           tabsetPanel(id = "re_results",
                                                                                       tabPanel("Plot", value = "tab4-2-1", withSpinner(plotlyOutput("ModelRE_Plot"))),
                                                                                       tabPanel("Estimate Z", value = "tab4-2-2",
                                                                                                p("Estimate the Z value based on the age random effects. Typically, choose the
                                                                              coordinates of the tallest bar in the histogram as one of the end points. Click
                                                                              Choose Start/Choose End to select coordinates of the histogram bars or manually
                                                                              enter values for the end points."),
                                                                                                wellPanel(id = "z_value_control",
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
                           tabPanel("Growth Curve Summaries", icon = icon("edit", lib = "glyphicon"), value = "tab5",
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
