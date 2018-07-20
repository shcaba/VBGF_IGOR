

# computes coefficient of variation for a vector
coeff_var <- function(v) {sd(v, na.rm = TRUE) / mean(v, na.rm = TRUE)}

# report convergence failed message
convg_err <- function() {showModal(modalDialog(
  title = "Important message",
  "The model didn't converge. Check your starting values.",
  easyClose = TRUE
))}

# returns an optimization upon success, otherwise error message alert
# and returns NULL
opt <- function(obj, lower, upper) {
  opt = tryCatch({
    nlminb(obj$par, obj$fn, obj$gr, lower = lower, upper = upper)
  },
  error = function(cond) {
    convg_err()
    return(NULL)
  })
  return(opt)
}

# plot data points and the fitting curve
plot <- function(rv, rep, mode, model) {
  age_max = max(rv$selected_data[c(-1, -2, -3, -4)])
  # get current selected reads (the ones in the filtered data table)
  reads_choices = names(rv$selected_data[c(-1, -2, -3, -4)])
  colors = c(brewer.pal(length(reads_choices), "Greens"), "black")
  names(colors) = c(reads_choices, "Expected")
  
  p = ggplot(data = rv$selected_data) +
    labs(title = "Age and Growth Fit", x = "Age (yr)", y = "Length (mm)")
  for (i in reads_choices) {
    p = p + geom_point(mapping = aes_string(Area = "Area", Sex = "Sex", x = i, y= "Length", color = paste("'", i, "'", sep = "")))
  }
  if (mode == "analyze" ) {
    if (model == "nls") {
      model_name = rv$last_run_type
    } else {
      model_name = rv$last_re_run_type
    }
    if (!is.null(model_name)) {
      if (model_name == "vb") {
        p = p + stat_function(fun = function(x) rep[1,1] * (1 - exp(-rep[2,1] * (x - rep[3,1]))),
                              aes(x = x, color = "Expected"), data = data.frame(x = c(0, age_max))) 
      } else if (model_name == "gompertz") {
        p = p + stat_function(fun = function(x) rep[1,1] * exp(-rep[2,1] * exp(-rep[3,1] * x)),
                              aes(x = x, color = "Expected"), data = data.frame(x = c(0, age_max))) 
      } else if (model_name == "logistic") {
        p = p + stat_function(fun = function(x) rep[1,1] / (1 + rep[2,1] * exp(-rep[3,1] * x)),
                              aes(x = x, color = "Expected"), data = data.frame(x = c(0, age_max))) 
      } else if (model_name == "linear") {
        p = p + stat_function(fun = function(x) rep[1,1] + x * rep[2,1],
                              aes(x = x, color = "Expected"), data = data.frame(x = c(0, age_max))) 
      } else { # model_name == "schnute"
        p = p + stat_function(fun = function(x) (rep[1,1] + rep[2,1] * exp(rep[3,1] * x)) ** rep[4,1],
                              aes(x = x, color = "Expected"), data = data.frame(x = c(0, age_max))) 
      }
    } 
  }
  p = p + scale_colour_manual(name = "", values = colors)
  return(ggplotly(p, tooltip = c("x", "y", "Area", "Sex")) %>% plotly::config(modeBarButtonsToRemove = 
    list("sendDataToCloud", "zoom2d", "pan2d", "select2d", "lasso2d", "zoomIn2d", "zoomOut2d", "autoScale2d",
         "resetScale2d", "hoverCompareCartesian", "hoverClosestCartesian"), displaylogo = FALSE, collaborate = FALSE))
}

shinyServer(function(input, output, session) {

  hide(id = "loading-content", anim = TRUE, animType = "fade")    
  show("app-content")
  
  rv <- reactiveValues(
    # summary table for estimation runs
    vb_summaries = data.frame("Run name" = character(), "Starting Linf" = double(), "Starting K" = double(),
                           "Starting t0" = double(), "Linf mean" = double(), "Linf sd" = double(), "K mean" = double(), "K sd" = double(),
                           "t0 mean" = double(), "t0 sd" = double(), "CV Length Error" = double(), "CV Age Error" = double(),
                           "alpha" = double(), "sigma" = double(), "beta" = double(), "shape" = double(), "scale" = double(), "z" = double(),
                           check.names = FALSE),
    schnute_summaries = data.frame("Run name" = character(), "Starting r0" = double(), "Starting b" = double(), "Starting k" = double(), 
                                   "Starting m" = double(), "r0 mean" = double(), "r0 Std" = double(), "b mean" = double(), "b Std" = double(), 
                                   "k mean" = double(), "k Std" = double(), "m mean" = double(), "m Std" = double(), "CV Length Error" = double(), 
                                   "CV Age Error" = double(), "alpha" = double(), "sigma" = double(), "beta" = double(), "shape" = double(), "scale" = double(),
                                   "z" = double(), check.names = FALSE),
    linear_summaries = data.frame("Run name" = character(), "Starting intercept" = double(), "Starting slope" = double(),
                                  "Intercept mean" = double(), "Intercept Std" = double(), "Slope mean" = double(), "Slope Std" = double(),
                                  "CV Length Error" = double(), "CV Age Error" = double(), "alpha" = double(), "sigma" = double(), 
                                  "beta" = double(), "shape" = double(), "scale" = double(), "z" = double(), check.names = FALSE),
    logistic_summaries = data.frame("Run name" = character(), "Starting a" = double(), "Starting b" = double(), "Starting k" = double(), 
                                    "a mean" = double(), "a Std" = double(), "b mean" = double(), "b Std" = double(),
                                    "k mean" = double(), "k Std" = double(), "CV Length Error" = double(), "CV Age Error" = double(), "alpha" = double(), "sigma" = double(),
                                    "beta" = double(), "shape" = double(), "scale" = double(), "z" = double(), check.names = FALSE),
    gompertz_summaries = data.frame("Run name" = character(), "Starting a" = double(), "Starting b" = double(),
                                    "Starting k" = double(), "a mean" = double(), "a Std" = double(), "b mean" = double(), "b Std" = double(),
                                    "k mean" = double(), "k Std" = double(), "CV Length Error" = double(), "CV Age Error" = double(), "alpha" = double(), "sigma" = double(),
                                    "beta" = double(), "shape" = double(), "scale" = double(), "z" = double(), check.names = FALSE),
    # filtered data table as estimation runs input
    selected_data = NULL,
    type = NULL, # only for random effects model: normal, exponential, or gamma
    get_start = FALSE,
    get_end = FALSE,
    last_run_type = NULL, # vb, gompertz, linear, etc.
    last_re_run_type = NULL, # same as last_run_type
    last_re_run = 0,
    z_plot_layer = NULL
  )
  
  # download a sample data table
  output$downloadData <- downloadHandler(
    filename = "sample_data.csv",
    content <- function(file) {
      file.copy("sample.csv", file)
    }
  )
  
  # original data table uploaded by user
  fish <- reactive({
    req(input$file)
    return(read.csv(input$file$datapath))
  })
  
  # renders the original data table uploaded by user
  output$input_data <- renderDT(
    fish(), rownames = FALSE
  )
  
  # creates control panel to select reads, areas, and sex
  output$input_control <- renderUI({
    reads_choices = paste0("Read", seq_len(ncol(fish()) - 4)) # only include the ReadX columns
    area_choices = unique(fish()$Area)
    tagList(
      fluidRow(
        column(3, 
          selectizeInput("reads_selected", "Reads to use", choices = reads_choices, selected = reads_choices, multiple = TRUE)
        ),
        column(3,
          selectizeInput("areas_selected", "Areas to use", choices = area_choices, selected = area_choices, multiple = TRUE)
        ),
        column(3, selectInput("sex", "Sex", choices = c("Any", "F", "M"))),
        column(3, actionButton("preview_data", "Preview Selected Data", 
                               icon = icon("wrench", lib = "glyphicon"), style = "margin-top:25px"))
      )
    )
  })
  
  # show the Analyze tab when user clicks the analyze button
  observeEvent(input$go_analyze, updateTabsetPanel(session, "menu", selected = "Analyze"))
  
  # renders the selected/filtered data table (downloadable)
  output$selected_data <- renderDT(
    rv$selected_data,
    extensions = 'Buttons',
    options = list(
      dom = 'Bfrtip',
      scrollX = TRUE,
      buttons = 
        list('copy', 'print', list(
          extend = 'collection',
          buttons = c('csv', 'excel'),
          text = 'Download'
        ))
    ),
    rownames = FALSE
  )
  
  # computes the filtered data table when user presses preview data button
  # IMPORTANT: filtered data table and plot only update when preview button is pressed!!!
  observeEvent(input$preview_data, {
    if (length(input$reads_selected) == 0 || length(input$areas_selected) == 0) {
      showModal(modalDialog(
        title = "Important message",
        "You have to select at least one read and one area!",
        easyClose = TRUE
      ))
    } else {
      cols_to_use = c("Sample", "Area", "Sex", "Length", input$reads_selected)
      if (input$sex == "Any") {
        rv$selected_data = subset(fish(), Area %in% input$areas_selected)
      } else {
        rv$selected_data = subset(fish(), Sex == input$sex & Area %in% input$areas_selected)
      }
      rv$selected_data = rv$selected_data[cols_to_use]
      rv$last_run_type = NULL
      rv$last_re_run_type = NULL
      # swith to "Selected Data" tab
      updateTabsetPanel(session, "data", selected = "Selected Data") 
    }
  })
  
  # creates a scatter plot for selected data
  output$selected_data_plot <- renderPlotly({
    # do not plot when the filtered data table is not ready
    if (is.null(rv$selected_data)) return()
    plot(rv, NULL, "input", NULL)
  })
  
  # show the Summaries tab when user clicks the summaries button
  observeEvent(input$go_summaries, updateTabsetPanel(session, "menu", selected = "Summaries"))
  
  # creates a select menu for the select read choice in nonlinear fit model 
  # based on the available reads in the filtered data table
  output$selected_read <- renderUI({
    if (!is.null(rv$selected_data) && input$fit_data == "1 - Selected Read") {
      reads_choices = names(rv$selected_data[c(-1, -2, -3, -4)])
      tagList(
        helpText("Select a single read from the", em("filtered"), "data"),
        selectInput("read_selected", "Read to use", choices = reads_choices)
      )
    }
  })
  
  # creates a downloadable summary table for estimation runs
  output$vb_summaries <- renderDT(
    rv$vb_summaries, extensions = 'Buttons',
    options = list(
      dom = 'Bfrtip',
      scrollX = TRUE,
      buttons = 
        list(I('colvis'), 'copy', 'print', list(
          extend = 'collection',
          buttons = c('csv', 'excel'),
          text = 'Download'
        ))
    ),
    rownames = FALSE
  )
  
  output$linear_summaries <- renderDT(
    rv$linear_summaries, extensions = 'Buttons',
    options = list(
      dom = 'Bfrtip',
      scrollX = TRUE,
      buttons =
        list(I('colvis'), 'copy', 'print', list(
          extend = 'collection',
          buttons = c('csv', 'excel'),
          text = 'Download'
        ))
    ),
    rownames = FALSE
  )

  output$gompertz_summaries <- renderDT(
    rv$gompertz_summaries, extensions = 'Buttons',
    options = list(
      dom = 'Bfrtip',
      scrollX = TRUE,
      buttons =
        list(I('colvis'), 'copy', 'print', list(
          extend = 'collection',
          buttons = c('csv', 'excel'),
          text = 'Download'
        ))
    ),
    rownames = FALSE
  )

  output$logistic_summaries <- renderDT(
    rv$logistic_summaries, extensions = 'Buttons',
    options = list(
      dom = 'Bfrtip',
      scrollX = TRUE,
      buttons =
        list(I('colvis'), 'copy', 'print', list(
          extend = 'collection',
          buttons = c('csv', 'excel'),
          text = 'Download'
        ))
    ),
    rownames = FALSE
  )

  output$schnute_summaries <- renderDT(
    rv$schnute_summaries, extensions = 'Buttons',
    options = list(
      dom = 'Bfrtip',
      scrollX = TRUE,
      buttons =
        list(I('colvis'), 'copy', 'print', list(
          extend = 'collection',
          buttons = c('csv', 'excel'),
          text = 'Download'
        ))
    ),
    rownames = FALSE
  )
  
  observe(rep_nls())
  
  # run estimates for a nonlinear model and returns a report as a dataframe
  rep_nls <- eventReactive(input$run_nls, {
    rv$test = rv$test + 1
    if (!is.null(rv$selected_data)) {
      num_reads = ncol(rv$selected_data) - 4
      len = rv$selected_data$Length
      num = length(len)
      # age matrix: nrow = #reads, ncol = num
      age = t(rv$selected_data[c(-1, -2, -3, -4)])
      
      fit_data = strtoi(substring(input$fit_data, 1, 1))
      
      len_use = len
      if (fit_data == 1) {
        # fit to selected read
        age_use = age[input$read_selected,]
      } else if (fit_data == 2) {
        # fit to average age
        # create an age vector that is the mean across all reads
        age_use = apply(age, 2, mean)
      } else if (fit_data == 3) {
        # fit to median age
        # create an age vector that is the median across all reads
        age_use = apply(age, 2, median)
      } else {
        # fit to multiple ages
        data = melt(rv$selected_data, id = c("Sample", "Area", "Sex", "Length"))
        age_use = data$value
        len_use = data$Length
      }
      
      runname = input$runname_nls
      
      if (input$model_nls == "Linear") {
        rv$last_run_type = "linear"
        if (input$fit_type == "Standard Deviation") {
          model = tryCatch({
            lm(len_use ~ age_use)
          },
          error = function(cond) {
            convg_err()
            return(NULL)
          })
          # returs NULL if model didn't converge
          if (is.null(model)) return(NULL)
          rep = summary(model)[["coefficients"]]
          newRow = format(data.frame("Run name" = runname, "Starting intercept" = NA, "Starting slope" = NA,
                     "Intercept mean" = rep[1,1], "Intercept Std" = rep[1,2], "Slope mean" = rep[2,1], "Slope Std" = rep[2,2],
                     "CV Length Error" = NA, "CV Age Error" = NA, "alpha" = NA, "sigma" = NA, "beta" = NA, "shape" = NA,
                     "scale" = NA, "z" = NA, check.names = FALSE), digits = 4)
        } else {
          intercept = as.numeric(input$intercept_nls)
          slope = as.numeric(input$slope_nls)
          
          data = list(len = len_use, age = age_use)
          CV_Lt = as.numeric(input$CV_const)
          parameters = list(intercept = intercept, slope = slope, CV_Lt = CV_Lt)
          obj = MakeADFun(data = data, parameters = parameters, DLL = "linear_likelihood")
          lower = c(-15, 0, exp(1) ^ (-10))
          upper = c(3 * max(len_use), 1000, exp(1) ^ (2))
          opt = opt(obj, lower, upper)
          
          # return NULL if model didn't converge
          if (is.null(opt)) return(NULL)
          
          rep = summary(sdreport(obj))
          newRow = format(data.frame("Run name" = runname, "Starting intercept" = intercept, "Starting slope" = slope,
                    "Intercept mean" = rep[1,1], "Intercept Std" = rep[1,2], "Slope mean" = rep[2,1], "Slope Std" = rep[2,2],
                    "CV Length Error" = rep[3,1], "CV Age Error" = NA, "alpha" = NA, "sigma" = NA, "beta" = NA, "shape" = NA,
                    "scale" = NA, "z" = NA, check.names = FALSE), digits = 4)
        }
        rv$linear_summaries = rbind(rv$linear_summaries, newRow)
      } else if (input$model_nls == "Gompertz") {
        rv$last_run_type = "gompertz"
        a = as.numeric(input$gom_a_nls)
        b = as.numeric(input$gom_b_nls)
        k = as.numeric(input$gom_k_nls)
        if (input$fit_type == "Standard Deviation") {
          model = tryCatch({
            nls(len ~ a * exp(-b * exp(-k * age)),
                data = list(age = age_use, len = len_use),
                start = list(a = a, b = b, k = k),
                control = list(reltol=0.00000000001))
          },
          error = function(cond) {
            convg_err()
            return(NULL)
          })
          # returs NULL if model didn't converge
          if (is.null(model)) return(NULL)
          rep = summary(model)[["parameters"]]
          newRow = format(data.frame("Run name" = runname, "Starting a" = a, "Starting b" = b,
                    "Starting k" = k, "a mean" = rep[1,1], "a Std" = rep[1,2], "b mean" = rep[2,1], "b Std" = rep[2,2],
                    "k mean" = rep[3,1], "k Std" = rep[3,2], "CV Length Error" = NA, "CV Age Error" = NA, "alpha" = NA, "sigma" = NA,
                    "beta" = NA, "shape" = NA, "scale" = NA, "z" = NA, check.names = FALSE), digits = 4)
        } else {
          data = list(len = len_use, age = age_use)
          CV_Lt = as.numeric(input$CV_const)
          parameters = list(a = a, b = b, k = k, CV_Lt = CV_Lt)
          obj = MakeADFun(data = data, parameters = parameters, DLL = "gompertz_likelihood")
          lower = c(-Inf, -Inf, -Inf, exp(1) ^ (-10))
          upper = c(Inf, Inf, Inf, exp(1) ^ (2))
          opt = opt(obj, lower, upper)
          
          # return NULL if model didn't converge
          if (is.null(opt)) return(NULL)
          
          rep = summary(sdreport(obj))
          newRow = format(data.frame("Run name" = runname, "Starting a" = a, "Starting b" = b,
                     "Starting k" = k, "a mean" = rep[1,1], "a Std" = rep[1,2], "b mean" = rep[2,1], "b Std" = rep[2,2],
                     "k mean" = rep[3,1], "k Std" = rep[3,2], "CV Length Error" = rep[4,1], "CV Age Error" = NA, "alpha" = NA, "sigma" = NA,
                     "beta" = NA, "shape" = NA, "scale" = NA, "z" = NA, check.names = FALSE), digits = 4)
        }
        rv$gompertz_summaries = rbind(rv$gompertz_summaries, newRow)
      } else if (input$model_nls == "Logistic") {
        rv$last_run_type = "logistic"
        a = as.numeric(input$log_a_nls)
        b = as.numeric(input$log_b_nls)
        k = as.numeric(input$log_k_nls)
        if (input$fit_type == "Standard Deviation") {
          model = tryCatch({
            nls(len ~ a / (1 + b * exp(-k * age)),
                data = list(age = age_use, len = len_use),
                start = list(a = a, b = b, k = k),
                control = list(reltol=0.00000000001))
          },
          error = function(cond) {
            convg_err()
            return(NULL)
          })
          # returs NULL if model didn't converge
          if (is.null(model)) return(NULL)
          rep = summary(model)[["parameters"]]
          newRow = format(data.frame("Run name" = runname, "Starting a" = a, "Starting b" = b,
                    "Starting k" = k, "a mean" = rep[1,1], "a Std" = rep[1,2], "b mean" = rep[2,1], "b Std" = rep[2,2],
                    "k mean" = rep[3,1], "k Std" = rep[3,2], "CV Length Error" = NA, "CV Age Error" = NA, "alpha" = NA, "sigma" = NA,
                    "beta" = NA, "shape" = NA, "scale" = NA, "z" = NA, check.names = FALSE), digits = 4)
        } else {
          data = list(len = len_use, age = age_use)
          CV_Lt = as.numeric(input$CV_const)
          parameters = list(a = a, b = b, k = k, CV_Lt = CV_Lt)
          obj = MakeADFun(data = data, parameters = parameters, DLL = "logistic_likelihood")
          lower = c(-Inf, -Inf, -Inf, exp(1) ^ (-10))
          upper = c(Inf, Inf, Inf, exp(1) ^ (2))
          opt = opt(obj, lower, upper)
          
          # return NULL if model didn't converge
          if (is.null(opt)) return(NULL)
          
          rep = summary(sdreport(obj))
          newRow = format(data.frame("Run name" = runname, "Starting a" = a, "Starting b" = b, "Starting k" = k, 
                    "a mean" = rep[1,1], "a Std" = rep[1,2], "b mean" = rep[2,1], "b Std" = rep[2,2],
                     "k mean" = rep[3,1], "k Std" = rep[3,2], "CV Length Error" = rep[4,1], "CV Age Error" = NA, "alpha" = NA, "sigma" = NA,
                     "beta" = NA, "shape" = NA, "scale" = NA, "z" = NA, check.names = FALSE), digits = 4)
        }
        rv$logistic_summaries = rbind(rv$logistic_summaries, newRow)
      } else if (input$model_nls == "Schnute") {
        rv$last_run_type = "schnute"
        r0 = as.numeric(input$sch_r0_nls)
        b = as.numeric(input$sch_b_nls)
        k = as.numeric(input$sch_k_nls)
        m = as.numeric(input$sch_m_nls)
        if (input$fit_type == "Standard Deviation") {
          model = tryCatch({
            nls(len ~ (r0 + b * exp(k * age)) ** m,
                data = list(age = age_use, len = len_use),
                start = list(r0 = r0, b = b, k = k, m = m),
                control = list(reltol=0.00000000001))
          },
          error = function(cond) {
            convg_err()
            return(NULL)
          })
          # returs NULL if model didn't converge
          if (is.null(model)) return(NULL)
          rep = summary(model)[["parameters"]]
          newRow = format(data.frame("Run name" = runname, "Starting r0" = r0, "Starting b" = b,
                     "Starting k" = k, "Starting m" = m, "r0 mean" = rep[1,1], "r0 Std" = rep[1,2],
                     "b mean" = rep[2,1], "b Std" = rep[2,2], "k mean" = rep[3,1], "k Std" = rep[3,2], "m mean" = rep[4,1], "m Std" = rep[4,2],
                      "CV Length Error" = NA, "CV Age Error" = NA, "alpha" = NA, "sigma" = NA, "beta" = NA, "shape" = NA, "scale" = NA,
                      "z" = NA, check.names = FALSE), digits = 4)
        } else {
          data = list(len = len_use, age = age_use)
          CV_Lt = as.numeric(input$CV_const)
          parameters = list(r0 = r0, b = b, k = k, m = m, CV_Lt = CV_Lt)
          obj = MakeADFun(data = data, parameters = parameters, DLL = "schnute_likelihood")
          lower = c(-Inf, -Inf, -Inf, -Inf, exp(1) ^ (-10))
          upper = c(Inf, Inf, Inf, Inf, exp(1) ^ (2))
          opt = opt(obj, lower, upper)
          
          # return NULL if model didn't converge
          if (is.null(opt)) return(NULL)
          
          rep = summary(sdreport(obj))
          newRow = format(data.frame("Run name" = runname, "Starting r0" = r0, "Starting b" = b, "Starting k" = k, 
                    "Starting m" = m, "r0 mean" = rep[1,1], "r0 Std" = rep[1,2], "b mean" = rep[2,1], "b Std" = rep[2,2], 
                    "k mean" = rep[3,1], "k Std" = rep[3,2], "m mean" = rep[4,1], "m Std" = rep[4,2], "CV Length Error" = rep[5,1], 
                    "CV Age Error" = NA, "alpha" = NA, "sigma" = NA, "beta" = NA, "shape" = NA, "scale" = NA,
                     "z" = NA, check.names = FALSE), digits = 4)
        }
        rv$schnute_summaries = rbind(rv$schnute_summaries, newRow)
      } else { # von bertlanffy
        rv$last_run_type = "vb"
        linf = as.numeric(input$linf_nls)
        kappa = as.numeric(input$kappa_nls)
        t0 = as.numeric(input$t0_nls)
        
        if (input$fit_type == "Standard Deviation") {
          model = tryCatch({
            nls(len ~ linf * (1 - exp(-kappa * (age - t0))),
                data = list(age = age_use, len = len_use),
                start = list(linf = linf, kappa = kappa, t0 = t0),
                control = list(reltol=0.00000000001))
          },
          error = function(cond) {
            convg_err()
            return(NULL)
          })
          
          # returs NULL if model didn't converge
          if (is.null(model)) return(NULL)
          
          rep = summary(model)[["parameters"]]
          newRow = format(data.frame("Run name" = runname, "Starting Linf" = linf, "Starting K" = kappa,
                     "Starting t0" = t0, "Linf mean" = rep[1,1], "Linf sd" = rep[1,2], "K mean" = rep[2,1], "K sd" = rep[2,2],
                     "t0 mean" = rep[3,1], "t0 sd" = rep[3,2], "CV Length Error" = NA, "CV Age Error" = NA, "alpha" = NA, "sigma" = NA,
                     "beta" = NA, "shape" = NA, "scale" = NA, "z" = NA, check.names = FALSE), digits = 4)
        } else {
          data = list(len = len_use, age = age_use)
          CV_Lt = as.numeric(input$CV_const)
          parameters = list(linf = linf, kappa = kappa, t0 = t0, CV_Lt = CV_Lt)
          obj = MakeADFun(data = data, parameters = parameters, DLL = "vb_likelihood")
          lower = c(0.75 * max(len_use), 0.0001, -15, exp(1) ^ (-10))
          upper = c(3 * max(len_use), 1, 1, exp(1) ^ (2))
          opt = opt(obj, lower, upper)
          
          # return NULL if model didn't converge
          if (is.null(opt)) return(NULL)
          
          rep = summary(sdreport(obj))
          newRow = format(data.frame("Run name" = runname, "Starting Linf" = linf, "Starting K" = kappa,
                     "Starting t0" = t0, "Linf mean" = rep[1,1], "Linf sd" = rep[1,2], "K mean" = rep[2,1], "K sd" = rep[2,2],
                     "t0 mean" = rep[3,1], "t0 sd" = rep[3,2], "CV Length Error" = rep[4,1], "CV Age Error" = NA, "alpha" = NA, "sigma" = NA,
                     "beta" = NA, "shape" = NA, "scale" = NA, "z" = NA, check.names = FALSE), digits = 4)
        }
        rv$vb_summaries = rbind(rv$vb_summaries, newRow)
      }
      return(rep)
    }
  })
  
  observe(rep_re())
  # run estimates for a random effects model and returns a report as a dataframe
  rep_re <- eventReactive(input$run_re, {
    if (!is.null(rv$selected_data)) {
      
      updateTextInput(session, "start_x", value = "")
      updateTextInput(session, "start_y", value = "")
      updateTextInput(session, "end_x", value = "")
      updateTextInput(session, "end_y", value = "") 
      rv$z_plot_layer = NULL
      rv$get_start = FALSE
      rv$get_end = FALSE
      
      num_reads = ncol(rv$selected_data) - 4
      len = rv$selected_data$Length
      num = length(len)
      # age matrix: nrow = #reads, ncol = num
      age = t(rv$selected_data[c(-1, -2, -3, -4)])
      
      
      if (as.numeric(input$CV_e) < 0) {
        # estimates age error
        cv = apply(age, 2, coeff_var)
        CV_e = sqrt(sum(cv * cv) / num)
      } else {
        CV_e = as.numeric(input$CV_e)
      }
      
      data = list(
        age = age,
        len = len,
        CV_e = CV_e,
        num_reads = num_reads
      )
      
      CV_Lt = as.numeric(input$CV_Lt)
      alpha = as.numeric(input$alpha)
      sigma_age = as.numeric(input$sigma_age)
      beta = as.numeric(input$beta)
      gam_shape = as.numeric(input$gam_shape)
      gam_scale = as.numeric(input$gam_scale)
      runname =  input$runname_re
      
      if (input$model_re == "Linear") {
        rv$last_re_run_type = "linear"
        dll = paste("linear", input$likelihoods, sep = "_")
        intercept = as.numeric(input$intercept_re)
        slope = as.numeric(input$slope_re)
        if (input$likelihoods == "Normal") {
          rv$type = "norm"
          parameters = list(intercept = intercept, slope = slope, CV_Lt = CV_Lt, 
                            alpha = alpha, sigma_age = sigma_age, age_re = rep(mean(age, na.rm = TRUE), num))
          
          lower = c(-Inf, -Inf, exp(1) ^ (-10), 0, exp(1) ^ (-100), rep(min(age), num))
          # upper bounds
          upper = c(Inf, Inf, exp(1) ^ (2), 500, exp(1) ^ (100), rep(max(age), num))
          obj = MakeADFun(data = data, parameters = parameters, random = "age_re", DLL = dll)
          opt = opt(obj, lower, upper)
          
          # return NULL if model didn't converge
          if (is.null(opt)) return(NULL)
          
          rep = summary(sdreport(obj))
          newRow = format(data.frame("Run name" = runname, "Starting intercept" = intercept, "Starting slope" = slope,
                   "Intercept mean" = rep[1,1], "Intercept Std" = rep[1,2], "Slope mean" = rep[2,1], "Slope Std" = rep[2,2],
                   "CV Length Error" = rep[3,1], "CV Age Error" = CV_e, "alpha" = rep[4,1], "sigma" = rep[5,1], "beta" = NA, "shape" = NA,
                   "scale" = NA, "z" = NA, check.names = FALSE), digits = 4)
        } else if (input$likelihoods == "Exponential") {
          rv$type = "exp"
          parameters = list(intercept = intercept, slope = slope, CV_Lt = CV_Lt, beta = beta, age_re = rep(mean(age, na.rm = TRUE), num))
          
          lower = c(-Inf, -Inf, exp(1) ^ (-10), exp(1) ^ (-100), rep(min(age), num))
          # upper bounds
          upper = c(Inf, Inf, exp(1) ^ (2), exp(1) ^ (100), rep(max(age), num) )
          obj = MakeADFun(data = data, parameters = parameters, random = "age_re", DLL = dll)
          opt = opt(obj, lower, upper)
          
          # return NULL if model didn't converge
          if (is.null(opt)) return(NULL)
          
          rep = summary(sdreport(obj))
          newRow = format(data.frame("Run name" = runname, "Starting intercept" = intercept, "Starting slope" = slope,
                   "Intercept mean" = rep[1,1], "Intercept Std" = rep[1,2], "Slope mean" = rep[2,1], "Slope Std" = rep[2,2],
                   "CV Length Error" = rep[3,1], "CV Age Error" = CV_e, "alpha" = NA, "sigma" = NA, "beta" = rep[4,1], "shape" = NA,
                   "scale" = NA, "z" = NA, check.names = FALSE), digits = 4)
        } else {
          rv$type = "gam"
          parameters = list(intercept = intercept, slope = slope, CV_Lt = CV_Lt,
                            gam_shape = gam_shape, gam_scale = gam_scale, age_re = rep(mean(age, na.rm = TRUE), num))
          
          lower = c(-Inf, -Inf, exp(1) ^ (-10), 0, 0, rep(min(age), num))
          # upper bounds
          upper = c(Inf, Inf, exp(1) ^ (2), 100, 100, rep(max(age), num))
          obj = MakeADFun(data = data, parameters = parameters, random = "age_re", DLL = dll)
          opt = opt(obj, lower, upper)
          
          # return NULL if model didn't converge
          if (is.null(opt)) return(NULL)
          
          rep = summary(sdreport(obj))
          newRow = format(data.frame("Run name" = runname, "Starting intercept" = intercept, "Starting slope" = slope,
                   "Intercept mean" = rep[1,1], "Intercept Std" = rep[1,2], "Slope mean" = rep[2,1], "Slope Std" = rep[2,2],
                   "CV Length Error" = rep[3,1], "CV Age Error" = CV_e, "alpha" = NA, "sigma" = NA, "beta" = NA, "shape" = rep[4,1],
                   "scale" = rep[5,1], "z" = NA, check.names = FALSE), digits = 4)
        }
        rv$linear_summaries = rbind(rv$linear_summaries, newRow)
        rv$last_re_run = nrow(rv$linear_summaries)
      } else if (input$model_re == "Gompertz") {
        rv$last_re_run_type = "gompertz"
        dll = paste("gompertz", input$likelihoods, sep = "_")
        a = as.numeric(input$gom_a_re)
        b = as.numeric(input$gom_b_re)
        k = as.numeric(input$gom_k_re)
        if (input$likelihoods == "Normal") {
          rv$type = "norm"
          parameters = list(a = a, b = b, k = k, CV_Lt = CV_Lt, 
                            alpha = alpha, sigma_age = sigma_age, age_re = rep(mean(age, na.rm = TRUE), num))
          
          lower = c(-Inf, -Inf, -Inf, exp(1) ^ (-10), 0, exp(1) ^ (-100), rep(min(age), num))
          # upper bounds
          upper = c(Inf, Inf, Inf, exp(1) ^ (2), 500, exp(1) ^ (100), rep(max(age), num))
          obj = MakeADFun(data = data, parameters = parameters, random = "age_re", DLL = dll)
          opt = opt(obj, lower, upper)
          
          # return NULL if model didn't converge
          if (is.null(opt)) return(NULL)
          
          rep = summary(sdreport(obj))
          newRow = format(data.frame("Run name" = runname, "Starting a" = a, "Starting b" = b,
                   "Starting k" = k, "a mean" = rep[1,1], "a Std" = rep[1,2], "b mean" = rep[2,1], "b Std" = rep[2,2],
                   "k mean" = rep[3,1], "k Std" = rep[3,2], "CV Length Error" = rep[4,1], "CV Age Error" = CV_e, "alpha" = rep[5,1], "sigma" = rep[6,1],
                   "beta" = NA, "shape" = NA, "scale" = NA, "z" = NA, check.names = FALSE), digits = 4)
        } else if (input$likelihoods == "Exponential") {
          rv$type = "exp"
          parameters = list(a = a, b = b, k = k, CV_Lt = CV_Lt, beta = beta, age_re = rep(mean(age, na.rm = TRUE), num))
          
          lower = c(-Inf, -Inf, -Inf, exp(1) ^ (-10), exp(1) ^ (-100), rep(min(age), num))
          # upper bounds
          upper = c(Inf, Inf, Inf, exp(1) ^ (2), exp(1) ^ (100), rep(max(age), num) )
          obj = MakeADFun(data = data, parameters = parameters, random = "age_re", DLL = dll)
          opt = opt(obj, lower, upper)
          
          # return NULL if model didn't converge
          if (is.null(opt)) return(NULL)
          
          rep = summary(sdreport(obj))
          newRow =format(data.frame("Run name" = runname, "Starting a" = a, "Starting b" = b,
                  "Starting k" = k, "a mean" = rep[1,1], "a Std" = rep[1,2], "b mean" = rep[2,1], "b Std" = rep[2,2],
                  "k mean" = rep[3,1], "k Std" = rep[3,2], "CV Length Error" = rep[4,1], "CV Age Error" = CV_e, "alpha" = NA, "sigma" = NA,
                  "beta" = rep[5,1], "shape" = NA, "scale" = NA, "z" = NA, check.names = FALSE), digits = 4)
        } else {
          rv$type = "gam"
          parameters = list(a = a, b = b, k = k, CV_Lt = CV_Lt,
                            gam_shape = gam_shape, gam_scale = gam_scale, age_re = rep(mean(age, na.rm = TRUE), num))
          
          lower = c(-Inf, -Inf, -Inf, exp(1) ^ (-10), 0, 0, rep(min(age), num))
          # upper bounds
          upper = c(Inf, Inf, Inf, exp(1) ^ (2), 100, 100, rep(max(age), num))
          obj = MakeADFun(data = data, parameters = parameters, random = "age_re", DLL = dll)
          opt = opt(obj, lower, upper)
          
          # return NULL if model didn't converge
          if (is.null(opt)) return(NULL)
          
          rep = summary(sdreport(obj))
          newRow = format(data.frame("Run name" = runname, "Starting a" = a, "Starting b" = b,
                   "Starting k" = k, "a mean" = rep[1,1], "a Std" = rep[1,2], "b mean" = rep[2,1], "b Std" = rep[2,2],
                   "k mean" = rep[3,1], "k Std" = rep[3,2], "CV Length Error" = rep[4,1], "CV Age Error" = CV_e, "alpha" = NA, "sigma" = NA,
                   "beta" = NA, "shape" = rep[5,1], "scale" = rep[6,1], "z" = NA, check.names = FALSE), digits = 4)
        }
        rv$gompertz_summaries = rbind(rv$gompertz_summaries, newRow)
        rv$last_re_run = nrow(rv$gompertz_summaries)
      } else if (input$model_re == "Logistic") {
        rv$last_re_run_type = "logistic"
        dll = paste("logistic", input$likelihoods, sep = "_")
        a = as.numeric(input$log_a_re)
        b = as.numeric(input$log_b_re)
        k = as.numeric(input$log_k_re)
        if (input$likelihoods == "Normal") {
          rv$type = "norm"
          parameters = list(a = a, b = b, k = k, CV_Lt = CV_Lt, 
                            alpha = alpha, sigma_age = sigma_age, age_re = rep(mean(age, na.rm = TRUE), num))
          
          lower = c(-Inf, -Inf, -Inf, exp(1) ^ (-10), 0, exp(1) ^ (-100), rep(min(age), num))
          # upper bounds
          upper = c(Inf, Inf, Inf, exp(1) ^ (2), 500, exp(1) ^ (100), rep(max(age), num))
          obj = MakeADFun(data = data, parameters = parameters, random = "age_re", DLL = dll)
          opt = opt(obj, lower, upper)
          
          # return NULL if model didn't converge
          if (is.null(opt)) return(NULL)
          
          rep = summary(sdreport(obj))
          newRow = format(data.frame("Run name" = runname, "Starting a" = a, "Starting b" = b,
                   "Starting k" = k, "a mean" = rep[1,1], "a Std" = rep[1,2], "b mean" = rep[2,1], "b Std" = rep[2,2],
                   "k mean" = rep[3,1], "k Std" = rep[3,2], "CV Length Error" = rep[4,1], "CV Age Error" = CV_e, "alpha" = rep[5,1], "sigma" = rep[6,1],
                   "beta" = NA, "shape" = NA, "scale" = NA, "z" = NA, check.names = FALSE), digits = 4)
        } else if (input$likelihoods == "Exponential") {
          rv$type = "exp"
          parameters = list(a = a, b = b, k = k, CV_Lt = CV_Lt, beta = beta, age_re = rep(mean(age, na.rm = TRUE), num))
          
          lower = c(-Inf, -Inf, -Inf, exp(1) ^ (-10), exp(1) ^ (-100), rep(min(age), num))
          # upper bounds
          upper = c(Inf, Inf, Inf, exp(1) ^ (2), exp(1) ^ (100), rep(max(age), num) )
          obj = MakeADFun(data = data, parameters = parameters, random = "age_re", DLL = dll)
          opt = opt(obj, lower, upper)
          
          # return NULL if model didn't converge
          if (is.null(opt)) return(NULL)
          
          rep = summary(sdreport(obj))
          newRow =format(data.frame("Run name" = runname, "Starting a" = a, "Starting b" = b,
                  "Starting k" = k, "a mean" = rep[1,1], "a Std" = rep[1,2], "b mean" = rep[2,1], "b Std" = rep[2,2],
                  "k mean" = rep[3,1], "k Std" = rep[3,2], "CV Length Error" = rep[4,1], "CV Age Error" = CV_e, "alpha" = NA, "sigma" = NA,
                  "beta" = rep[5,1], "shape" = NA, "scale" = NA, "z" = NA, check.names = FALSE), digits = 4)
        } else {
          rv$type = "gam"
          parameters = list(a = a, b = b, k = k, CV_Lt = CV_Lt,
                            gam_shape = gam_shape, gam_scale = gam_scale, age_re = rep(mean(age, na.rm = TRUE), num))
          
          lower = c(-Inf, -Inf, -Inf, exp(1) ^ (-10), 0, 0, rep(min(age), num))
          # upper bounds
          upper = c(Inf, Inf, Inf, exp(1) ^ (2), 100, 100, rep(max(age), num))
          obj = MakeADFun(data = data, parameters = parameters, random = "age_re", DLL = dll)
          opt = opt(obj, lower, upper)
          
          # return NULL if model didn't converge
          if (is.null(opt)) return(NULL)
          
          rep = summary(sdreport(obj))
          newRow = format(data.frame("Run name" = runname, "Starting a" = a, "Starting b" = b,
                   "Starting k" = k, "a mean" = rep[1,1], "a Std" = rep[1,2], "b mean" = rep[2,1], "b Std" = rep[2,2],
                   "k mean" = rep[3,1], "k Std" = rep[3,2], "CV Length Error" = rep[4,1], "CV Age Error" = CV_e, "alpha" = NA, "sigma" = NA,
                   "beta" = NA, "shape" = rep[5,1], "scale" = rep[6,1], "z" = NA, check.names = FALSE), digits = 4)
        }
        rv$logistic_summaries = rbind(rv$logistic_summaries, newRow)
        rv$last_re_run = nrow(rv$logistic_summaries)
      } else if (input$model_re == "Schnute") {
        rv$last_re_run_type = "schnute"
        dll = paste("Schnute", input$likelihoods, sep = "_")
        r0 = as.numeric(input$sch_r0_re)
        b = as.numeric(input$sch_b_re)
        k = as.numeric(input$sch_k_re)
        m = as.numeric(input$sch_m_re)
        if (input$likelihoods == "Normal") {
          rv$type = "norm"
          parameters = list(r0 = r0, b = b, k = k, m = m, CV_Lt = CV_Lt, 
                   alpha = alpha, sigma_age = sigma_age, age_re = rep(mean(age, na.rm = TRUE), num))
          
          lower = c(-Inf, -Inf, -Inf, -Inf, exp(1) ^ (-10), 0, exp(1) ^ (-100), rep(min(age), num))
          # upper bounds
          upper = c(Inf, Inf, Inf, Inf, exp(1) ^ (2), 500, exp(1) ^ (100), rep(max(age), num))
          obj = MakeADFun(data = data, parameters = parameters, random = "age_re", DLL = dll)
          opt = opt(obj, lower, upper)
          
          # return NULL if model didn't converge
          if (is.null(opt)) return(NULL)
          
          rep = summary(sdreport(obj))
          newRow = format(data.frame("Run name" = runname, "Starting r0" = r0, "Starting b" = b, "Starting k" = k, 
                   "Starting m" = m, "r0 mean" = rep[1,1], "r0 Std" = rep[1,2], "b mean" = rep[2,1], "b Std" = rep[2,2], 
                   "k mean" = rep[3,1], "k Std" = rep[3,2], "m mean" = rep[4,1], "m Std" = rep[4,2], "CV Length Error" = rep[5,1], 
                   "CV Age Error" = CV_e, "alpha" = rep[6,1], "sigma" = rep[7,1], "beta" = NA, "shape" = NA, "scale" = NA,
                   "z" = NA, check.names = FALSE), digits = 4)
        } else if (input$likelihoods == "Exponential") {
          rv$type = "exp"
          parameters = list(r0 = r0, b = b, k = k, m = m, CV_Lt = CV_Lt, beta = beta, age_re = rep(mean(age, na.rm = TRUE), num))
          
          
          lower = c(-Inf, -Inf, -Inf, -Inf, exp(1) ^ (-10), exp(1) ^ (-100), rep(min(age), num))
          # upper bounds
          upper = c(Inf, Inf, Inf, Inf, exp(1) ^ (2), exp(1) ^ (100), rep(max(age), num) )
          obj = MakeADFun(data = data, parameters = parameters, random = "age_re", DLL = dll)
          opt = opt(obj, lower, upper)
          
          # return NULL if model didn't converge
          if (is.null(opt)) return(NULL)
          
          rep = summary(sdreport(obj))
          newRow = format(data.frame("Run name" = runname, "Starting r0" = r0, "Starting b" = b, "Starting k" = k, 
                       "Starting m" = m, "r0 mean" = rep[1,1], "r0 Std" = rep[1,2], "b mean" = rep[2,1], "b Std" = rep[2,2], 
                       "k mean" = rep[3,1], "k Std" = rep[3,2], "m mean" = rep[4,1], "m Std" = rep[4,2], "CV Length Error" = rep[5,1], 
                       "CV Age Error" = CV_e, "alpha" = NA, "sigma" = NA, "beta" = rep[6,1], "shape" = NA, "scale" = NA,
                       "z" = NA, check.names = FALSE), digits = 4)
        } else {
          rv$type = "gam"
          parameters = list(r0 = r0, b = b, k = k, m = m, CV_Lt = CV_Lt,
                            gam_shape = gam_shape, gam_scale = gam_scale, age_re = rep(mean(age, na.rm = TRUE), num))
          
          lower = c(-Inf, -Inf, -Inf, -Inf, exp(1) ^ (-10), 0, 0, rep(min(age), num))
          # upper bounds
          upper = c(Inf, Inf, Inf, Inf, exp(1) ^ (2), 100, 100, rep(max(age), num))
          obj = MakeADFun(data = data, parameters = parameters, random = "age_re", DLL = dll)
          opt = opt(obj, lower, upper)
          
          # return NULL if model didn't converge
          if (is.null(opt)) return(NULL)
          
          rep = summary(sdreport(obj))
          newRow = format(data.frame("Run name" = runname, "Starting r0" = r0, "Starting b" = b, "Starting k" = k, 
                   "Starting m" = m, "r0 mean" = rep[1,1], "r0 Std" = rep[1,2], "b mean" = rep[2,1], "b Std" = rep[2,2], 
                   "k mean" = rep[3,1], "k Std" = rep[3,2], "m mean" = rep[4,1], "m Std" = rep[4,2], "CV Length Error" = rep[5,1], 
                   "CV Age Error" = CV_e, "alpha" = NA, "sigma" = NA, "beta" = NA, "shape" = rep[6,1], "scale" = rep[7,1],
                   "z" = NA, check.names = FALSE), digits = 4)
        }
        rv$schnute_summaries = rbind(rv$schnute_summaries, newRow)
        rv$last_re_run = nrow(rv$schnute_summaries)
      } else { # input$model_re == "Von Bertlanffy"
        rv$last_re_run_type = "vb"
        dll = paste("vbre", input$likelihoods, sep = "_")
        linf = as.numeric(input$linf_re)
        kappa = as.numeric(input$kappa_re)
        t0 = as.numeric(input$t0_re)
        if (input$likelihoods == "Normal") {
          rv$type = "norm"
          parameters = list(linf = linf, kappa = kappa, t0 = t0, CV_Lt = CV_Lt, 
                            alpha = alpha, sigma_age = sigma_age, age_re = rep(mean(age, na.rm = TRUE), num))
          
          # lower bounds: linf, kappa, t0, CV_Lt, alpha, sigma_age, age_re
          lower = c(0.75 * max(len), 0.0001, -15, exp(1) ^ (-10), 0, exp(1) ^ (-100), rep(min(age), num))
          # upper bounds
          upper = c(3 * max(len), 1, 1, exp(1) ^ (2), 500, exp(1) ^ (100), rep(max(age), num))
          obj = MakeADFun(data = data, parameters = parameters, random = "age_re", DLL = dll)
          opt = opt(obj, lower, upper)
          
          # return NULL if model didn't converge
          if (is.null(opt)) return(NULL)
          
          rep = summary(sdreport(obj))
          newRow = format(data.frame("Run name" = runname, "Starting Linf" = linf, "Starting K" = kappa,
                   "Starting t0" = t0, "Linf mean" = rep[1,1], "Linf sd" = rep[1,2], "K mean" = rep[2,1], "K sd" = rep[2,2],
                   "t0 mean" = rep[3,1], "t0 sd" = rep[3,2], "CV Length Error" = rep[4,1], "CV Age Error" = CV_e, "alpha" = rep[5,1],
                   "sigma" = rep[6,1], "beta" = NA, "shape" = NA, "scale" = NA, "z" = NA, check.names = FALSE), digits = 4)
        } else if (input$likelihoods == "Exponential") {
          rv$type = "exp"
          parameters = list(linf = linf, kappa = kappa, t0 = t0, CV_Lt = CV_Lt, beta = beta, age_re = rep(mean(age, na.rm = TRUE), num))
          
          # lower bounds: linf, kappa, t0, CV_Lt, beta, age_re
          lower = c(0.75 * max(len), 0.0001, -15, exp(1) ^ (-10), exp(1) ^ (-100), rep(min(age), num))
          # upper bounds
          upper = c(3 * max(len), 1, 1, exp(1) ^ (2), exp(1) ^ (100), rep(max(age), num) )
          obj = MakeADFun(data = data, parameters = parameters, random = "age_re", DLL = dll)
          opt = opt(obj, lower, upper)
          
          # return NULL if model didn't converge
          if (is.null(opt)) return(NULL)
          
          rep = summary(sdreport(obj))
          newRow = format(data.frame("Run name" = runname, "Starting Linf" = linf, "Starting K" = kappa,
                   "Starting t0" = t0, "Linf mean" = rep[1,1], "Linf sd" = rep[1,2], "K mean" = rep[2,1], "K sd" = rep[2,2],
                   "t0 mean" = rep[3,1], "t0 sd" = rep[3,2], "CV Length Error" = rep[4,1], "CV Age Error" = CV_e, "alpha" = NA,
                   "sigma" = NA, "beta" = rep[5,1], "shape" = NA, "scale" = NA, "z" = NA, check.names = FALSE), digits = 4)
        } else { #gamma
          rv$type = "gam"
          parameters = list(linf = linf, kappa = kappa, t0 = t0, CV_Lt = CV_Lt, 
                            gam_shape = gam_shape, gam_scale = gam_scale, age_re = rep(mean(age, na.rm = TRUE), num))
          
          # lower bounds: linf, kappa, t0, CV_Lt, gam_shape, gam_scale, age_re
          lower = c(0.75 * max(len), 0.0001, -15, exp(1) ^ (-10), 0, 0, rep(min(age), num))
          # upper bounds
          upper = c(3 * max(len), 1, 1, exp(1) ^ (2), 100, 100, rep(max(age), num))
          obj = MakeADFun(data = data, parameters = parameters, random = "age_re", DLL = dll)
          opt = opt(obj, lower, upper)
          
          # return NULL if model didn't converge
          if (is.null(opt)) return(NULL)
          
          rep = summary(sdreport(obj))
          newRow = format(data.frame("Run name" = runname, "Starting Linf" = linf, "Starting K" = kappa,
                  "Starting t0" = t0, "Linf mean" = rep[1,1], "Linf sd" = rep[1,2], "K mean" = rep[2,1], "K sd" = rep[2,2],
                  "t0 mean" = rep[3,1], "t0 sd" = rep[3,2], "CV Length Error" = rep[4,1], "CV Age Error" = CV_e, "alpha" = NA,
                   "sigma" = NA, "beta" = NA, "shape" = rep[5,1], "scale" = rep[6,1], "z" = NA, check.names = FALSE), digits = 4)
        }
        rv$vb_summaries = rbind(rv$vb_summaries, newRow)
        rv$last_re_run = nrow(rv$vb_summaries)
      }
      
      return(rep)
    }
  })
  
  # nonlinear fit plot
  output$ModelNLS_Plot <- renderPlotly({
    if (is.null(rv$selected_data)) return()
    plot(rv, rep_nls(), "analyze", "nls")
  })
  
  # random effects plot
  output$ModelRE_Plot <- renderPlotly({
    if (is.null(rv$selected_data)) return()
    plot(rv, rep_re(), "analyze", "re")
  })
  
  # alerts the computed z value when user clicks get z value button
  # and updates the z value in the most recent RE model run
  observeEvent(input$get_z, {
    y = c(as.numeric(input$start_y), as.numeric(input$end_y))
    x = c(as.numeric(input$start_x), as.numeric(input$end_x))
    
    slope = tryCatch({
      lm(y ~ x)$coeff[[2]]
    },
    error = function(cond) {
      showModal(modalDialog(
        title = "Important message",
        "Check your input values.",
        easyClose = TRUE
      ))
      return(NULL)
    })
    
    if (!is.null(slope)) {
      slope = format(slope, digits = 4)
      rv[[paste0(rv$last_re_run_type, "_summaries")]][["z"]][rv$last_re_run] = slope
      showModal(modalDialog(
        title = "Important message",
        paste0("z = ", slope),
        easyClose = TRUE
      ))
      rv$z_plot_layer =  geom_line(data = data.frame(x = c(as.numeric(input$start_x), as.numeric(input$end_x)), 
                                                     y = c(as.numeric(input$start_y), as.numeric(input$end_y))),
                                   aes(x, y))
    }
  })
  
  observeEvent(input$get_start, {
    rv$get_end = FALSE
    rv$get_start = TRUE
  })
  
  observeEvent(input$get_end, {
    rv$get_start = FALSE
    rv$get_end = TRUE
  })
  
  # populate endpoints coordinates with histogram bar coordinates
  # when user clicks on the bars
  observe({
    coords = event_data("plotly_click", source = "hist_RE")
    if (rv$get_start) {
      updateTextInput(session, "start_x", value = coords$x)
      updateTextInput(session, "start_y", value = coords$y)   
    }
    if (rv$get_end) {
      updateTextInput(session, "end_x", value = coords$x)
      updateTextInput(session, "end_y", value = coords$y) 
    }
  })
  
  # age random effects histogram
  output$hist_RE <- renderPlotly({
    if (is.null(rv$selected_data) || is.null(rv$type)) return()
    num = length(rv$selected_data$Length)
    last_re_run_type = rv$last_re_run_type
    if (last_re_run_type == "vb" || last_re_run_type == "gompertz" || last_re_run_type == "logistic") {
      if (rv$type == "exp") {
        data = data.frame("Age" = rep_re()[6:(5+num), 1])
      } else {
        data = data.frame("Age" = rep_re()[7:(6+num), 1])
      }
    } else if (last_re_run_type == "linear") {
      if (rv$type == "exp") {
        data = data.frame("Age" = rep_re()[5:(4+num), 1])
      } else {
        data = data.frame("Age" = rep_re()[6:(5+num), 1])
      }
    } else { # last_re_run_type == "schnute"
      if (rv$type == "exp") {
        data = data.frame("Age" = rep_re()[7:(6+num), 1])
      } else {
        data = data.frame("Age" = rep_re()[8:(7+num), 1])
      }
    }

    
    bins = seq(min(data), max(data), length.out = input$bins + 1)
    
    p = ggplot(data, aes(Age)) + 
      geom_histogram(aes(fill = ..count..), breaks = bins) + 
      scale_fill_gradient("Count", low = "#071c07", high = "#35c435") +
      labs(title = "Age Distribution", x = "Age (yr)", y = "Count")
    
    p = p + rv$z_plot_layer
    
    ggplotly(p, tooltip = c("Age", "Count"), source = "hist_RE") %>% plotly::config(modeBarButtonsToRemove = 
      list("sendDataToCloud", "zoom2d", "pan2d", "select2d", "lasso2d", "zoomIn2d", "zoomOut2d", "autoScale2d",
           "resetScale2d", "hoverCompareCartesian", "hoverClosestCartesian"), displaylogo = FALSE, collaborate = FALSE)
  })
})