compile("vbre_Gamma.cpp")
compile("vbre_Exponential.cpp")
compile("vbre_Normal.cpp")
compile("vb_likelihood.cpp")
compile("vb_multiple_ages.cpp")
dyn.load(dynlib("vbre_Gamma"))
dyn.load(dynlib("vbre_Exponential"))
dyn.load(dynlib("vbre_Normal"))
dyn.load(dynlib("vb_likelihood"))
dyn.load(dynlib("vb_multiple_ages"))

# computes coefficient of variation for a vector
coeff_var <- function(v) {sd(v) / mean(v)}

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
plot <- function(rv, rep, mode) {
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
  if (mode == "analyze") {
    p = p + stat_function(fun = function(x) rep[1,1] * (1 - exp(-rep[2,1] * (x - rep[3,1]))),
                          aes(x = x, color = "Expected"), data = data.frame(x = c(0, age_max)))  
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
    summaries = data.frame("Run name" = character(), "Starting Linf" = double(), "Starting K" = double(),
                           "Starting t0" = double(), "Linf mean" = double(), "Linf sd" = double(), "K mean" = double(), "K sd" = double(),
                           "t0 mean" = double(), "t0 sd" = double(), "CV Length Error" = double(), "CV Age Error" = double(),
                           "alpha" = double(), "sigma" = double(), "beta" = double(), "shape" = double(), "scale" = double(), "z" = double(),
                           check.names = FALSE
    ),
    # filtered data table as estimation runs input
    selected_data = NULL,
    type = NULL,
    get_start = FALSE,
    get_end = FALSE,
    last_re_run = 0
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
      
      # swith to "Selected Data" tab
      updateTabsetPanel(session, "data", selected = "Selected Data") 
    }
  })
  
  # creates a scatter plot for selected data
  output$selected_data_plot <- renderPlotly({
    # do not plot when the filtered data table is not ready
    if (is.null(rv$selected_data)) return()
    plot(rv, NULL, "input")
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
  output$summaries <- renderDT(
    rv$summaries, extensions = 'Buttons',
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
  
  # run estimates for a nonlinear model and returns a report as a dataframe
  rep_nls <- eventReactive(input$run_nls, {
    if (!is.null(rv$selected_data)) {
      num_reads = ncol(rv$selected_data) - 4
      len = rv$selected_data$Length
      num = length(len)
      # age matrix: nrow = #reads, ncol = num
      age = t(rv$selected_data[c(-1, -2, -3, -4)])
      
      linf = as.numeric(input$linf_nls)
      kappa = as.numeric(input$kappa_nls)
      t0 = as.numeric(input$t0_nls)
      
      fit_data = strtoi(substring(input$fit_data, 1, 1))
      if (fit_data != 4) {
        # use R least square function for fit
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
          data = melt(rv$selected_data, id = c("Sample", "Area", "Sex", "Length"))
          age_use = data$value
          len_use = data$Length
        }
        
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
        newRow = format(data.frame("Run name" = input$runname_nls,
                                   "Starting Linf" = input$linf_nls,
                                   "Starting K" = input$kappa_nls,
                                   "Starting t0" = input$t0_nls,
                                   "Linf mean" = rep[1,1],
                                   "Linf sd" = rep[1,2],
                                   "K mean" = rep[2,1],
                                   "K sd" = rep[2,2],
                                   "t0 mean" = rep[3,1],
                                   "t0 sd" = rep[3,2],
                                   "CV Length Error" = NA,
                                   "CV Age Error" = NA,
                                   "alpha" = NA,
                                   "sigma" = NA,
                                   "beta" = NA,
                                   "shape" = NA,
                                   "scale" = NA,
                                   "z" = NA,
                                   check.names = FALSE), digits = 4)
      } else {
        #        if (fit_data == 4) {
        # fit to likelihood
        data = list(len = len, age = age)
        CV_Lt = as.numeric(input$CV_const)
        parameters = list(linf = linf, kappa = kappa, t0 = t0, CV_Lt = CV_Lt)
        obj = MakeADFun(data = data, parameters = parameters, DLL = "vb_likelihood")
        lower = c(0.75 * max(len), 0.0001, -15, exp(1) ^ (-10))
        upper = c(3 * max(len), 1, 1, exp(1) ^ (2))
        opt = opt(obj, lower, upper)
        
        # return NULL if model didn't converge
        if (is.null(opt)) return(NULL)
        
        rep = summary(sdreport(obj))
        newRow = format(data.frame("Run name" = input$runname_nls,
                                   "Starting Linf" = input$linf_nls,
                                   "Starting K" = input$kappa_nls,
                                   "Starting t0" = input$t0_nls,
                                   "Linf mean" = rep[1,1],
                                   "Linf sd" = rep[1,2],
                                   "K mean" = rep[2,1],
                                   "K sd" = rep[2,2],
                                   "t0 mean" = rep[3,1],
                                   "t0 sd" = rep[3,2],
                                   "CV Length Error" = rep[4,1],
                                   "CV Age Error" = NA,
                                   "alpha" = NA,
                                   "sigma" = NA,
                                   "beta" = NA,
                                   "shape" = NA,
                                   "scale" = NA,
                                   "z" = NA,
                                   check.names = FALSE), digits = 4)
        # } else { # fit_data == 5
        #   # fit to multiple reads (regression)
        #   data = list(len = len, age = age, num_reads = num_reads)
        #   parameters = list(linf = linf, kappa = kappa, t0 = t0)
        #   obj = MakeADFun(data = data, parameters = parameters, DLL = "vb_multiple_ages")
        #   opt = nlminb(obj$par, obj$fn, obj$gr)
        #   rv$rep_nls = summary(sdreport(obj))
        #   newRow = format(data.frame("Run name" = input$runname_nls,
        #                              "Starting Linf" = input$linf_nls,
        #                              "Starting K" = input$kappa_nls,
        #                              "Starting t0" = input$t0_nls,
        #                              "Linf mean" = rv$rep_nls[1,1],
        #                              "Linf sd" = rv$rep_nls[1,2],
        #                              "K mean" = rv$rep_nls[2,1],
        #                              "K sd" = rv$rep_nls[2,2],
        #                              "t0 mean" = rv$rep_nls[3,1],
        #                              "t0 sd" = rv$rep_nls[3,2],
        #                              "CV Length Error" = NA,
        #                              "CV Age Error" = NA,
        #                              "alpha" = NA,
        #                              "sigma" = NA,
        #                              "beta" = NA,
        #                              "shape" = NA,
        #                              "scale" = NA,
        #                              "z" = NA,
        #                              check.names = FALSE), digits = 4)
        # }
      }
      # updata the summary table with a new row from the run
      rv$summaries = rbind(rv$summaries, newRow)
      return(rep)
    }
  })
  
  # run estimates for a random effects model and returns a report as a dataframe
  rep_re <- eventReactive(input$run_re, {
    if (!is.null(rv$selected_data)) {
      num_reads = ncol(rv$selected_data) - 4
      len = rv$selected_data$Length
      num = length(len)
      # age matrix: nrow = #reads, ncol = num
      age = t(rv$selected_data[c(-1, -2, -3, -4)])
      
      dll = paste("vbre", input$likelihoods, sep = "_")
      
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
      linf = as.numeric(input$linf_re)
      kappa = as.numeric(input$kappa_re)
      t0 = as.numeric(input$t0_re)
      CV_Lt = as.numeric(input$CV_Lt)
      alpha = as.numeric(input$alpha)
      sigma_age = as.numeric(input$sigma_age)
      beta = as.numeric(input$beta)
      gam_shape = as.numeric(input$gam_shape)
      gam_scale = as.numeric(input$gam_scale)
      if (input$likelihoods == "Normal") {
        rv$type = "norm"
        parameters = list(linf = linf, kappa = kappa, t0 = t0, CV_Lt = CV_Lt, 
                          alpha = alpha, sigma_age = sigma_age, age_re = rep(mean(age), num))
        
        # lower bounds: linf, kappa, t0, CV_Lt, alpha, sigma_age, age_re
        lower = c(0.75 * max(len), 0.0001, -15, exp(1) ^ (-10), 0, exp(1) ^ (-100), rep(min(age), num))
        # upper bounds
        upper = c(3 * max(len), 1, 1, exp(1) ^ (2), 500, exp(1) ^ (100), rep(max(age), num))
        obj = MakeADFun(data = data, parameters = parameters, random = "age_re", DLL = dll)
        opt = opt(obj, lower, upper)
        
        # return NULL if model didn't converge
        if (is.null(opt)) return(NULL)
        
        rep = summary(sdreport(obj))
        newRow = format(data.frame("Run name" = input$runname_re,
                                   "Starting Linf" = input$linf_re,
                                   "Starting K" = input$kappa_re,
                                   "Starting t0" = input$t0_re,
                                   "Linf mean" = rep[1,1],
                                   "Linf sd" = rep[1,2],
                                   "K mean" = rep[2,1],
                                   "K sd" = rep[2,2],
                                   "t0 mean" = rep[3,1],
                                   "t0 sd" = rep[3,2],
                                   "CV Length Error" = rep[4,1],
                                   "CV Age Error" = CV_e,
                                   "alpha" = rep[5,1],
                                   "sigma" = rep[6,1],
                                   "beta" = NA,
                                   "shape" = NA,
                                   "scale" = NA,
                                   "z" = NA,
                                   check.names = FALSE), digits = 4)
      } else if (input$likelihoods == "Exponential") {
        rv$type = "exp"
        parameters = list(linf = linf, kappa = kappa, t0 = t0, CV_Lt = CV_Lt, beta = beta, age_re = rep(mean(age), num))
        
        # lower bounds: linf, kappa, t0, CV_Lt, beta, age_re
        lower = c(0.75 * max(len), 0.0001, -15, exp(1) ^ (-10), exp(1) ^ (-100), rep(min(age), num))
        # upper bounds
        upper = c(3 * max(len), 1, 1, exp(1) ^ (2), exp(1) ^ (100), rep(max(age), num) )
        obj = MakeADFun(data = data, parameters = parameters, random = "age_re", DLL = dll)
        opt = opt(obj, lower, upper)
        
        # return NULL if model didn't converge
        if (is.null(opt)) return(NULL)
        
        rep = summary(sdreport(obj))
        newRow = format(data.frame("Run name" = input$runname_re,
                                   "Starting Linf" = input$linf_re,
                                   "Starting K" = input$kappa_re,
                                   "Starting t0" = input$t0_re,
                                   "Linf mean" = rep[1,1],
                                   "Linf sd" = rep[1,2],
                                   "K mean" = rep[2,1],
                                   "K sd" = rep[2,2],
                                   "t0 mean" = rep[3,1],
                                   "t0 sd" = rep[3,2],
                                   "CV Length Error" = rep[4,1],
                                   "CV Age Error" = CV_e,
                                   "alpha" = NA,
                                   "sigma" = NA,
                                   "beta" = rep[5,1],
                                   "shape" = NA,
                                   "scale" = NA,
                                   "z" = NA,
                                   check.names = FALSE), digits = 4)
      } else { #gamma
        rv$type = "gam"
        parameters = list(linf = linf, kappa = kappa, t0 = t0, CV_Lt = CV_Lt, 
                          gam_shape = gam_shape, gam_scale = gam_scale, age_re = rep(mean(age), num))
        
        # lower bounds: linf, kappa, t0, CV_Lt, gam_shape, gam_scale, age_re
        lower = c(0.75 * max(len), 0.0001, -15, exp(1) ^ (-10), 0, 0, rep(min(age), num))
        # upper bounds
        upper = c(3 * max(len), 1, 1, exp(1) ^ (2), 100, 100, rep(max(age), num))
        obj = MakeADFun(data = data, parameters = parameters, random = "age_re", DLL = dll)
        opt = opt(obj, lower, upper)
        
        # return NULL if model didn't converge
        if (is.null(opt)) return(NULL)
        
        rep = summary(sdreport(obj))
        newRow = format(data.frame("Run name" = input$runname_re,
                                   "Starting Linf" = input$linf_re,
                                   "Starting K" = input$kappa_re,
                                   "Starting t0" = input$t0_re,
                                   "Linf mean" = rep[1,1],
                                   "Linf sd" = rep[1,2],
                                   "K mean" = rep[2,1],
                                   "K sd" = rep[2,2],
                                   "t0 mean" = rep[3,1],
                                   "t0 sd" = rep[3,2],
                                   "CV Length Error" = rep[4,1],
                                   "CV Age Error" = CV_e,
                                   "alpha" = NA,
                                   "sigma" = NA,
                                   "beta" = NA,
                                   "shape" = rep[5,1],
                                   "scale" = rep[6,1],
                                   "z" = NA,
                                   check.names = FALSE), digits = 4)
      }
      # updata the summary table with a new row from the run
      rv$summaries = rbind(rv$summaries, newRow)
      rv$last_re_run = nrow(rv$summaries)
      return(rep)
    }
  })
  
  # nonlinear fit plot
  output$ModelNLS_Plot <- renderPlotly({
    if (is.null(rv$selected_data)) return()
    plot(rv, rep_nls(), "analyze")
  })
  
  # random effects plot
  output$ModelRE_Plot <- renderPlotly({
    if (is.null(rv$selected_data)) return()
    plot(rv, rep_re(), "analyze")
  })
  
  # prints the computed z value when user clicks get z value button
  # and updates the z value in the most recent RE model run
  output$z_value <- renderText({
    if (input$get_z == 0) return()
    isolate({
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
        rv$summaries[["z"]][rv$last_re_run] = slope
        return(paste0("z = ", slope))
      }
    })
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
    if (rv$type == "exp") {
      data = data.frame("Age" = rep_re()[7:(5+num), 1])
    } else {
      data = data.frame("Age" = rep_re()[7:(6+num), 1])
    }
    
    bins = seq(min(data), max(data), length.out = input$bins + 1)
    
    p = ggplot(data, aes(Age)) + 
      geom_histogram(aes(fill = ..count..), breaks = bins) + 
      scale_fill_gradient("Count", low = "#071c07", high = "#35c435") +
      labs(title = "Age Distribution", x = "Age (yr)", y = "Count")
    ggplotly(p, tooltip = c("Age", "Count"), source = "hist_RE") %>% plotly::config(modeBarButtonsToRemove = 
      list("sendDataToCloud", "zoom2d", "pan2d", "select2d", "lasso2d", "zoomIn2d", "zoomOut2d", "autoScale2d",
           "resetScale2d", "hoverCompareCartesian", "hoverClosestCartesian"), displaylogo = FALSE, collaborate = FALSE)
  })
})