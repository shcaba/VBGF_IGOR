# shiny server functions

# to show mathjax in plot
addResourcePath("tmp", tempdir())
f = tempfile(fileext = ".html")

shinyServer(function(input, output, session) {

  hide(id = "loading-content", anim = TRUE, animType = "fade")    
  show("app-content")
  
  # set help content
  observe({
    indTab = NULL
    nextPage = NULL
    
      # Index help data into right tab
      if (input$menu == "tab1") {
        indTab = c(1, 2)
        nextPage = nextPage = c("tab2", "tab2-1")
      }
      if (input$menu == "tab2" && input$tab2_subpanel == "tab2-1") {
        if (input$oto_age_menu == "tab2-1-1") {
          indTab = 3
          nextPage = c("tab2", "tab2-1", "tab2-1-3")
        }
        if (input$oto_age_menu == "tab2-1-3") {
          indTab = 4
          nextPage = c("tab2", "tab2-2", "tab2-1-1") #tab2-1-1 reset very necessary!!!
        }
      }
      if (input$menu == "tab2" && input$tab2_subpanel == "tab2-2") {
        indTab = c(5, 6)
        nextPage = c("tab3", "tab3-1")
      }
      if (input$menu == "tab3") {
        if(input$data == "tab3-1") {
          indTab = 7
          nextPage = c("tab3", "tab3-2")
        } else {
          indTab = 8
          nextPage = c("tab4", "tab4-1")
        }
      }
      if (input$menu == "tab4" && input$fit_choice == "tab4-1") {
        indTab = 9
        nextPage = c("tab4", "tab4-2")
      }
      if (input$menu == "tab4" && input$fit_choice == "tab4-2") {
        if (input$re_results == "tab4-2-2") {
          indTab = 11
          nextPage = c("tab5", "tab4-2-1") #tab4-2-1 reset very necessary!!!
        } else {
          indTab = 10
          nextPage = c("tab4", "tab4-2", "tab4-2-2")
        }
      } 
      if (input$menu == "tab5") {
        indTab = 12
        nextPage = "tab1"
      }
      # set help content
      if (!is.null(indTab) && !is.null(nextPage)) {
        session$sendCustomMessage(type = 'setHelpContent', message = list(steps = toJSON(helpData[indTab, ]),
              nextPage = toJSON(nextPage)))  
      }
  })
  
  observeEvent(input$startHelp,{
    if(input$startHelp == 0) return()
    # Auto start help if we press nextpage
    session$sendCustomMessage(type = 'startHelp', message = list(""))
  })
  
  rv <- reactiveValues(
    # original oto age data uploaded by user
    oto_age_data = NULL,
    # same as oto_age_data but contains a "Selected" T/F column to indicate if data point is selected
    oto_age_data_filter_info = NULL,
    # original age length data uploaded by user
    age_length_data = NULL,
    # summaries for oto age analysis
    oto_age_summary = NULL,
    brkpt_summary = NULL,
    slope_summary = NULL,
    intercept_summary = NULL,
    
    oto_age_model = NULL,
    oto_age_model_name = NULL,
    jitter_plot = NULL,
    predicted_oto_age_data = NULL, # data table shown in (Run a new oto-age model) Fitted Values tab
    
    existing_oto_age_model = NULL, # oto-age model uploaded by user
    oto_data = NULL, # table that contains otowt, used with existing_oto_age_model to get ages
    predicted_oto_age_data_from_exisiting = NULL,
    
    # summary table for growth curve estimation runs
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
    # filtered age length data table as estimation runs input
    selected_data = NULL,
    rep_nls = NULL, # standard growth model
    rep_re = NULL, # random effects growth model
    type = NULL, # only for random effects model: normal, exponential, or gamma
    get_start = FALSE,
    get_end = FALSE,
    last_run_type = NULL, # vb, gompertz, linear, etc.
    last_re_run_type = NULL, # same as last_run_type
    last_re_run = 0,
    z_plot_layer = NULL
  )
  
  # download a sample age length data table
  output$download_age_length_data <- downloadHandler(
    filename = "sample_age_length_data.csv",
    content <- function(file) {
      file.copy("sample_downloads/sample_age_length_data.csv", file)
    }
  )
  
  # download a sample otowt age data table
  output$download_oto_age_data <- downloadHandler(
    filename = "sample_oto_age_data.csv",
    content <- function(file) {
      file.copy("sample_downloads/sample_oto_age_data.csv", file)
    }
  )
  
  # returns TRUE if uploaded oto age file has the columns required, FALSE otherwise
  # column order does not matter, and will not reorder
  check_oto_age_file <- function(col_names_upload, col_names) {
    if (length(col_names_upload) != length(col_names)) {
      return(FALSE)
    }
    for (i in 1 : length(col_names_upload)) {
      if (!col_names_upload[i] %in% col_names) {
        return(FALSE)
      } 
    }
    return(TRUE)
  }
  
  # check if data has the required columns, order does not matter, but will reorder if data format is correct
  # returns data with column order "Sample", "Species", "Area", "Sex", "Length", "Read1", "Read2"...
  # returns NULL if uploaded age length file does not have the columns required
  check_age_len_file <- function(data) {
    col_names_upload = colnames(data)
    col_names = c("Sample", "Species", "Area", "Sex", "Length")
    
    if (length(col_names_upload) <= length(col_names)) {
      return(NULL)
    }
    for (i in 1 : length(col_names)) {
      if (!col_names[i] %in% col_names_upload) {
        return(NULL)
      }
    }
    num_reads = length(col_names_upload) - length(col_names)
    # reads columns need to be labeled as Read1, Read2, etc.
    for (i in 1 : num_reads) {
      col = paste0("Read", i)
      if (!col %in% col_names_upload) {
        return(NULL)
      }
    }
    # reorder the columns
    reads = paste0("Read", seq_len(num_reads))
    cols_to_use = c(col_names, reads)
    return(data[cols_to_use])
  }
  
  # original age_length data table uploaded by user
  observe({
    req(input$age_length_file)
    data = try(read.csv(input$age_length_file$datapath), silent = TRUE)
    if (class(data)[[1]] == "try-error") {
      invalid_upload()
      shinyjs::reset("age_length_file")
    } else {
      data = check_age_len_file(data)
      if (is.null(data)) {
        invalid_upload()
        shinyjs::reset("age_length_file")
      } else {
        rv$age_length_data = data
      }
    }
  })
  
  # original otowt age data table uploaded by user
  observe({
    req(input$oto_age_file)
    data = try(read.csv(input$oto_age_file$datapath), silent = TRUE)
    if (class(data)[[1]] == "try-error") {
      invalid_upload()
      shinyjs::reset("oto_age_file")
    } else {
      col_names_upload = colnames(data)
      col_names = c("Sample", "Species", "Area", "Sex", "Length", "Age", "OtoWt")
      if (!check_oto_age_file(col_names_upload, col_names)) {
        invalid_upload()
        shinyjs::reset("oto_age_file")
      } else {
        rv$oto_age_data = data
      }
    }
  })
  
  # show the Age vs. Otolith Weight Data tab
  observeEvent(input$go_oto_age, updateTabsetPanel(session, "menu", selected = "tab2"))
  # show the Length vs. Age Data tab
  observeEvent(input$go_age_length, updateTabsetPanel(session, "menu", selected = "tab3"))
  
  # renders the original otowt age data table uploaded by user
  output$oto_age_data <- renderDT(
    rv$oto_age_data, filter = 'top', rownames = FALSE
  )
  
  # renders the original age length data table uploaded by user
  output$input_data <- renderDT(
    rv$age_length_data, rownames = FALSE
  )
  
  # computes an oto age model based on user choice
  observeEvent(input$run_oto_age, {
    if (is.null(rv$oto_age_data)) {
      showModal(modalDialog(
        title = "Important message",
        "Please upload a file first.",
        easyClose = TRUE
      ))
    } else {
      disable("oto_age_controls")
      # show Summary & Plot tab
      updateTabsetPanel(session, "oto_age_menu", selected = "tab2-1-2")
      
      rv$oto_age_summary = NULL
      rv$brkpt_summary = NULL
      rv$slope_summary = NULL
      rv$intercept_summary = NULL
      
      oto_age_data_filtered = rv$oto_age_data[input$oto_age_data_rows_all, ]
      Selected = rep("Not Selected", nrow(rv$oto_age_data))
      Selected[input$oto_age_data_rows_all] = "Selected"
      rv$oto_age_data_filter_info = cbind(rv$oto_age_data, "Selected" = Selected)
      
      plot_jitter = FALSE
      add_lengths = FALSE
      
      rv$oto_age_model_name = input$oto_age_model_type
      
      if (rv$oto_age_model_name == "Linear") {
        model = try(lm(Age ~ OtoWt, data = oto_age_data_filtered), silent = TRUE)
        if (class(model)[[1]] == "try-error") {
          convg_err()
          enable("oto_age_controls")
          return()
        }
        rv$oto_age_summary = summary(model)[["coefficients"]]
      } else if (rv$oto_age_model_name == "Quadratic") {
        model = try(lm(Age ~ poly(OtoWt, 2, raw = TRUE), data = na.omit(oto_age_data_filtered)), 
                    silent = TRUE)
        if (class(model)[[1]] == "try-error") {
          convg_err()
          enable("oto_age_controls")
          return()
        }
        rv$oto_age_summary = summary(model)[["coefficients"]]
      } else {
        if (rv$oto_age_model_name == "1-Breakpoint Piecewise") {
          names = "Breakpoint 1"
          if (input$use_quantile) {
            q1 = as.numeric(input$wt_bkpt1)
            if (q1 < 0 || q1 > 1) {
              invalid_quantile()
              enable("oto_age_controls")
              return()
            }
            wt_brkpts = quantile(oto_age_data_filtered$OtoWt, q1, na.rm = TRUE)
          } else {
            wt_brkpts = c(as.numeric(input$wt_bkpt1))
          }
          if (input$add_length) {
            lt_brkpts = c(as.numeric(input$lt_bkpt1))
            add_lengths = TRUE
          }
        } else { # rv$oto_age_model_name == "2-Breakpoints Piecewise"
          names = c("Breakpoint 1", "Breakpoint 2")
          if (input$use_quantile) {
           q1 = as.numeric(input$wt_bkpt1_2bp)
           q2 = as.numeric(input$wt_bkpt2_2bp)
           if (q1 < 0 || q1 > 1 || q2 < 0 || q2 > 1) {
             invalid_quantile()
             enable("oto_age_controls")
             return()
           }
           wt_brkpts = c(quantile(oto_age_data_filtered$OtoWt, q1, na.rm = TRUE), 
                         quantile(oto_age_data_filtered$OtoWt, q2, na.rm = TRUE))
          } else {
            wt_brkpts = c(as.numeric(input$wt_bkpt1_2bp), as.numeric(input$wt_bkpt2_2bp))  
          }
          if (input$add_length) {
            lt_brkpts = c(as.numeric(input$lt_bkpt1_2bp), as.numeric(input$lt_bkpt2_2bp))
            add_lengths = TRUE
          }
        }
        
        if (input$jitter) {
          if (rv$oto_age_model_name == "1-Breakpoint Piecewise") {
            plot_jitter = TRUE
          }
          jitter = as.numeric(input$num_jitter)
        } else {
          jitter = 0
        }
        model = oto_age_model(oto_age_data_filtered, wt_brkpts = wt_brkpts, lt_brkpts = lt_brkpts, 
                              jitter = jitter, add_lengths = add_lengths)
        
        if (is.null(model)) {
          convg_err()
          enable("oto_age_controls")
          return()
        }
        
        psi = summary(model)[["psi"]]
        rownames(psi) = names
        rv$brkpt_summary = psi
        rv$slope_summary = slope(model)
        rv$intercept_summary = intercept(model)
      }
      rv$oto_age_model = model
      enable("download_oto_age_model")
      pred_ages = predict(rv$oto_age_model, rv$oto_age_data)
      rv$predicted_oto_age_data = cbind(rv$oto_age_data, "Predicted Ages" = pred_ages)
  
      if (plot_jitter) {
        rv$jitter_plot = jitter_plot(as.numeric(input$num_jitter), rv$oto_age_model)
      } else {
        rv$jitter_plot = NULL
      }
      enable("oto_age_controls")
    }
  })

  # download a sample otowt age data table
  output$download_oto_age_model <- downloadHandler(
    filename = "oto_age_model.rds",
    content <- function(file) {
      saveRDS(rv$oto_age_model, file)
    }
  )
  
  # disable model download button when the model is null (nothing has been computed)
  observe({
    if (is.null(rv$oto_age_model)) {
      disable("download_oto_age_model")
    }
  })
  
  # prints oto age analysis summaries
  output$oto_age_summary <- renderTable(rv$oto_age_summary, rownames = TRUE, digits = 4)
  output$brkpt_summary <- renderTable(rv$brkpt_summary, rownames = TRUE, digits = 4)
  output$slope_summary <- renderTable(rv$slope_summary, rownames = TRUE, digits = 4)
  output$intercept_summary <- renderTable(rv$intercept_summary, rownames = TRUE, digits = 4)
  
  # renders the oto age plot and the model fit
  output$oto_age_plot <- renderUI({
    if (is.null(rv$oto_age_data) || is.null(rv$oto_age_model_name)) return()
    
    helpText("hello")
    if (rv$oto_age_model_name == "Quadratic") {
      quad = y ~ poly(x, 2, raw = TRUE)
      p = ggplot(data = rv$oto_age_data_filter_info) +
        geom_point(aes(x = OtoWt, y = Age, colour = Selected)) +
        geom_smooth(method = "lm", se = FALSE, formula = quad, color = "#0c3f10", aes(x = OtoWt, y = Age))
    } else {
      # linear/ 1 breakpoint / 2 breakpoints
      Expected = predict(rv$oto_age_model, rv$oto_age_data) # for all data points (unfiltered)
      data = cbind(rv$oto_age_data_filter_info, Expected = Expected)
      p = ggplot(data = data) +
        geom_point(aes(x = OtoWt, y = Age, colour = Selected)) +
        geom_line(aes(x = OtoWt, y = Expected), color = "blue")
    }
    p = p + labs(title = "Otolith Weight and Age Fit", x = "Otolith Weight (mg)", y = "Age (yr)") +
      scale_colour_discrete(name = "")
    return(ggplotly(p, tooltip = c("x", "y")) %>%
             # R square and AIC value as annotations
          layout(annotations = list(text = TeX(sprintf("R^2 = %.3f \\\\ \\text{AIC} = %.3f", AIC(rv$oto_age_model), summary(rv$oto_age_model)$r.squared)),
                                    xref = "paper",  # normalize coordinates
                                    yref = "paper",
                                    x = 0.95,  # appears in bottom right corner
                                    y = 0.05,
                                    showarrow = FALSE),
                 font = list(family = "'Times New Roman', 'arial', 'Raleway','verdana'")) %>%
          plotly::config(modeBarButtonsToRemove = list("sendDataToCloud", "zoom2d", "pan2d", "select2d", "lasso2d", "zoomIn2d", "zoomOut2d", "autoScale2d",
                "resetScale2d", "hoverCompareCartesian", "hoverClosestCartesian"),
                displaylogo = FALSE, collaborate = FALSE, mathjax = "cdn")) %>%
          saveWidget(f)  # to show mathjax
    # to show mathjax
    tags$iframe(
      src = file.path("tmp", basename(f)),
      width = "100%", 
      height = "400",
      scrolling = "no", 
      seamless = "seamless", 
      frameBorder = "0"
    )
  })
  
  # jitter plots
  output$jitter_intercept_plot <- renderPlotly({
    if (is.null(rv$jitter_plot)) return()
    return(rv$jitter_plot[[1]])
  })
  
  output$jitter_slope1_plot <- renderPlotly({
    if (is.null(rv$jitter_plot)) return()
    return(rv$jitter_plot[[2]])
  })
  
  output$jitter_slope2_plot <- renderPlotly({
    if (is.null(rv$jitter_plot)) return()
    return(rv$jitter_plot[[3]])
  })
  
  output$jitter_breakpoint_plot <- renderPlotly({
    if (is.null(rv$jitter_plot)) return()
    return(rv$jitter_plot[[4]])
  })
  
  output$jitter_AIC_plot <- renderPlotly({
    if (is.null(rv$jitter_plot)) return()
    return(rv$jitter_plot[[5]])
  })

  output$jitter_comp_plot <- renderPlotly({
    if (is.null(rv$jitter_plot)) return()
    return(rv$jitter_plot[[6]])
  })
  
  # renders the predicted age data table from oto age fitting(downloadable)
  output$predicted_age_preview <- renderDT(
    rv$predicted_oto_age_data,
    server = FALSE, # make sure it downloads all data, not just the ones in view
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
  
  # create an age length data table, where age is predicted from otowt,
  # and change tab focus to Length vs. Age Data tab
  observeEvent(input$set_age_length_data, {
    if (is.null(rv$oto_age_model)) {
      showModal(modalDialog(
        title = "Important message",
        "You have not run an Age vs. Otolith Weight analysis.",
        easyClose = TRUE
      ))
    } else {
      data = rv$predicted_oto_age_data
      colnames(data)[which(names(data) == "Predicted Ages")] = "Read1"
      rv$age_length_data = data.frame("Sample" = data$Sample, "Species" = data$Species,
                               "Area" = data$Area, "Sex" = data$Sex,
                               "Length" = data$Length, "Read1" = data$Read1)
      # show Length vs. Age Data tab
      updateTabsetPanel(session, "menu", selected = "tab3")
    }
  })
  
  # oto-age model uploaded by user
  observe({
    req(input$model_oto_file)
    model = try(readRDS(input$model_oto_file$datapath), silent = TRUE)
    if (!("lm" %in% class(model))) {
      invalid_upload()
      shinyjs::reset("model_oto_file")
    } else {
      # qualified model object has class "lm" or c("lm", "segmented")
      rv$existing_oto_age_model = model
    }
  })
  
  # oto file upload along side with existing_oto_age_model
  observe({
    req(input$oto_file)
    data = try(read.csv(input$oto_file$datapath), silent = TRUE)
    if (class(data)[[1]] == "try-error") {
      invalid_upload()
      shinyjs::reset("oto_file")
    } else {
      col_names_upload = colnames(data)
      col_names = c("Sample", "Species", "Area", "Sex", "Length", "OtoWt")
      if (!check_oto_age_file(col_names_upload, col_names)) {
        invalid_upload()
        shinyjs::reset("oto_file")
      } else {
        rv$oto_data = data
      }
    }
  })
  
  # download a sample otowt data table for fitting ages with existing model
  output$download_oto_data <- downloadHandler(
    filename = "sample_oto_data_for_existing_model.csv",
    content <- function(file) {
      file.copy("sample_downloads/sample_oto_data_for_existing_model.csv", file)
    }
  )
  
  # download a sample otowt age model
  output$download_oto_model <- downloadHandler(
    filename = "sample_oto_age_model.rds",
    content <- function(file) {
      file.copy("sample_downloads/sample_oto_age_model.rds", file)
    }
  )
  
  # compute the predicted age data table using existing oto age model
  observeEvent(input$run_existing_oto_model, {
    if (is.null(rv$existing_oto_age_model) || is.null(rv$oto_data)) {
      showModal(modalDialog(
        title = "Important message",
        "Please upload a file first.",
        easyClose = TRUE
      ))
      return()
    }
    pred_ages = predict(rv$existing_oto_age_model, rv$oto_data)
    rv$predicted_oto_age_data_from_exisiting = cbind(rv$oto_data, "Predicted Ages" = pred_ages)
  })
  
  # renders the predicted age data table using existing oto age model
  output$existing_model_predicted_data <- renderDT(
    rv$predicted_oto_age_data_from_exisiting,
    server = FALSE, # make sure it downloads all data, not just the ones in view
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
    rownames = FALSE,
    filter = 'top'
  )
  
  observeEvent(input$set_age_length_data_from_existing_model, {
    if (is.null(rv$predicted_oto_age_data_from_exisiting)) {
      showModal(modalDialog(
        title = "Important message",
        "You have not produced predicted ages with your existing model.",
        easyClose = TRUE
      ))
    } else {
      data = rv$predicted_oto_age_data_from_exisiting[input$existing_model_predicted_data_rows_all, ]
      colnames(data)[which(names(data) == "Predicted Ages")] = "Read1"
      rv$age_length_data = data.frame("Sample" = data$Sample, "Species" = data$Species,
                                      "Area" = data$Area, "Sex" = data$Sex,
                                      "Length" = data$Length, "Read1" = data$Read1)
      # show Length vs. Age Data tab
      updateTabsetPanel(session, "menu", selected = "tab3")
    }
  })
  
  # creates control panel to select reads, areas, and sex
  output$input_control <- renderUI({
    if (!is.null(rv$age_length_data)) {
      reads_choices = paste0("Read", seq_len(ncol(rv$age_length_data) - 5)) # only include the ReadX columns
      area_choices = unique(rv$age_length_data$Area)
      species_choices = unique(rv$age_length_data$Species)
      tagList(
        fluidRow(
          column(4, 
                 selectizeInput("reads_selected", "Reads to use", choices = reads_choices, selected = reads_choices, multiple = TRUE)
          ),
          column(4,
                 selectizeInput("areas_selected", "Areas to use", choices = area_choices, selected = area_choices, multiple = TRUE)
          ),
          column(4, 
                 selectizeInput("species_selected", "Species to use", choices = species_choices, selected = species_choices, multiple = TRUE)
          )
        ),
        fluidRow(
          column(3, selectInput("sex", "Sex", choices = c("Any", "F", "M"))),
          column(3, actionButton("preview_data", "Preview Selected Data", 
                                 icon = icon("wrench", lib = "glyphicon"), style = "margin-top:25px"))
        )
      ) 
    }
  })
  
  # show the Growth Curve Analysis tab when user clicks the analyze button
  observeEvent(input$go_analyze, updateTabsetPanel(session, "menu", selected = "tab4"))
  
  # renders the selected/filtered data table (downloadable)
  output$selected_data <- renderDT(
    rv$selected_data,
    server = FALSE, # make sure it downloads all data, not just the ones in view
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
      cols_to_use = c("Sample", "Species", "Area", "Sex", "Length", sort(input$reads_selected))
      if (input$sex == "Any") {
        rv$selected_data = subset(rv$age_length_data, Area %in% input$areas_selected &
                                   Species %in% input$species_selected)
      } else {
        rv$selected_data = subset(rv$age_length_data, Sex == input$sex & 
                                    Area %in% input$areas_selected & Species %in% input$species_selected)
      }
      rv$selected_data = rv$selected_data[cols_to_use]
      rv$last_run_type = NULL
      rv$last_re_run_type = NULL
      # swith to "Selected Data" tab
      updateTabsetPanel(session, "data", selected = "tab3-2") 
    }
  })
  
  # creates a scatter plot for selected data
  output$selected_data_plot <- renderPlotly({
    # do not plot when the filtered data table is not ready
    if (is.null(rv$selected_data)) return()
    plot_growth_curve(rv, "input", NULL)
  })
  
  # show the Growth Curve Summaries tab when user clicks the summaries button
  observeEvent(input$go_summaries, updateTabsetPanel(session, "menu", selected = "tab5"))
  
  # creates a select menu for the select read choice in nonlinear fit model 
  # based on the available reads in the filtered data table
  output$selected_read <- renderUI({
    if (!is.null(rv$selected_data) && input$fit_data == "1 - Selected Read") {
      reads_choices = names(rv$selected_data[c(-1, -2, -3, -4, -5)])
      tagList(
        helpText("Select a single read from the", em("filtered"), "data"),
        selectInput("read_selected", "Read to use", choices = reads_choices)
      )
    }
  })
  
  # creates a downloadable summary table for estimation runs
  output$vb_summaries <- renderDT(
    rv$vb_summaries, extensions = 'Buttons',
    server = FALSE, # make sure it downloads all data, not just the ones in view
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
    server = FALSE, # make sure it downloads all data, not just the ones in view
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
    server = FALSE, # make sure it downloads all data, not just the ones in view
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
    server = FALSE, # make sure it downloads all data, not just the ones in view
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
    server = FALSE, # make sure it downloads all data, not just the ones in view
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
  
  # run estimates for a nonlinear model
  observeEvent(input$run_nls, {
    if (!is.null(rv$selected_data)) {
      disable("age_len_std_control")
      num_reads = ncol(rv$selected_data) - 5
      len = rv$selected_data$Length
      num = length(len)
      # age matrix: nrow = #reads, ncol = num
      age = t(rv$selected_data[c(-1, -2, -3, -4, -5)])
      
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
        data = melt(rv$selected_data, id = c("Sample", "Species", "Area", "Sex", "Length"))
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
          if (is.null(model)) {
            enable("age_len_std_control")
            return()
          }
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
          if (is.null(opt)) {
            enable("age_len_std_control")
            return()
          }
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
          if (is.null(model)) {
            enable("age_len_std_control")
            return()
          }
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
          if (is.null(opt)) {
            enable("age_len_std_control")
            return()
          }
          
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
          if (is.null(model)) {
            enable("age_len_std_control")
            return()
          }
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
          if (is.null(opt)) {
            enable("age_len_std_control")
            return()
          }
          
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
          if (is.null(model)) {
            enable("age_len_std_control")
            return()
          }
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
          if (is.null(opt)) {
            enable("age_len_std_control")
            return()
          }
          
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
          if (is.null(model)) {
            enable("age_len_std_control")
            return()
          }
          
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
          if (is.null(opt)) {
            enable("age_len_std_control")
            return()
          }
          
          rep = summary(sdreport(obj))
          newRow = format(data.frame("Run name" = runname, "Starting Linf" = linf, "Starting K" = kappa,
                     "Starting t0" = t0, "Linf mean" = rep[1,1], "Linf sd" = rep[1,2], "K mean" = rep[2,1], "K sd" = rep[2,2],
                     "t0 mean" = rep[3,1], "t0 sd" = rep[3,2], "CV Length Error" = rep[4,1], "CV Age Error" = NA, "alpha" = NA, "sigma" = NA,
                     "beta" = NA, "shape" = NA, "scale" = NA, "z" = NA, check.names = FALSE), digits = 4)
        }
        rv$vb_summaries = rbind(rv$vb_summaries, newRow)
      }
      enable("age_len_std_control")
      rv$rep_nls = rep
    }
  })
  
  # run estimates for a random effects model
  observeEvent(input$run_re, {
    if (!is.null(rv$selected_data)) {
      disable("age_len_re_control")
      updateTextInput(session, "start_x", value = "")
      updateTextInput(session, "start_y", value = "")
      updateTextInput(session, "end_x", value = "")
      updateTextInput(session, "end_y", value = "") 
      rv$z_plot_layer = NULL
      rv$get_start = FALSE
      rv$get_end = FALSE
      
      num_reads = ncol(rv$selected_data) - 5
      len = rv$selected_data$Length
      num = length(len)
      # age matrix: nrow = #reads, ncol = num
      age = t(rv$selected_data[c(-1, -2, -3, -4, -5)])
      
      
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
          if (is.null(opt)) {
            enable("age_len_re_control")
            return()
          }
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
          if (is.null(opt)) {
            enable("age_len_re_control")
            return()
          }
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
          if (is.null(opt)) {
            enable("age_len_re_control")
            return()
          }
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
          if (is.null(opt)) {
            enable("age_len_re_control")
            return()
          }
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
          if (is.null(opt)) {
            enable("age_len_re_control")
            return()
          }
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
          if (is.null(opt)) {
            enable("age_len_re_control")
            return()
          }
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
          if (is.null(opt)) {
            enable("age_len_re_control")
            return()
          }
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
          if (is.null(opt)) {
            enable("age_len_re_control")
            return()
          }
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
          if (is.null(opt)) {
            enable("age_len_re_control")
            return()
          }
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
          if (is.null(opt)) {
            enable("age_len_re_control")
            return()
          }
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
          if (is.null(opt)) {
            enable("age_len_re_control")
            return()
          }
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
          if (is.null(opt)) {
            enable("age_len_re_control")
            return()
          }
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
          if (is.null(opt)) {
            enable("age_len_re_control")
            return()
          }
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
          if (is.null(opt)) {
            enable("age_len_re_control")
            return()
          }
          
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
          if (is.null(opt)) {
            enable("age_len_re_control")
            return()
          }
          
          rep = summary(sdreport(obj))
          newRow = format(data.frame("Run name" = runname, "Starting Linf" = linf, "Starting K" = kappa,
                  "Starting t0" = t0, "Linf mean" = rep[1,1], "Linf sd" = rep[1,2], "K mean" = rep[2,1], "K sd" = rep[2,2],
                  "t0 mean" = rep[3,1], "t0 sd" = rep[3,2], "CV Length Error" = rep[4,1], "CV Age Error" = CV_e, "alpha" = NA,
                   "sigma" = NA, "beta" = NA, "shape" = rep[5,1], "scale" = rep[6,1], "z" = NA, check.names = FALSE), digits = 4)
        }
        rv$vb_summaries = rbind(rv$vb_summaries, newRow)
        rv$last_re_run = nrow(rv$vb_summaries)
      }
      enable("age_len_re_control")
      rv$rep_re = rep
    }
  })
  
  # nonlinear fit plot
  output$ModelNLS_Plot <- renderPlotly({
    if (is.null(rv$selected_data) || is.null(rv$rep_nls)) return()
    plot_growth_curve(rv, "analyze", "nls")
  })
  
  # random effects plot
  output$ModelRE_Plot <- renderPlotly({
    if (is.null(rv$selected_data) || is.null(rv$rep_re)) return()
    plot_growth_curve(rv, "analyze", "re")
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
    if (is.null(rv$selected_data) || is.null(rv$type) || is.null(rv$rep_re)) return()
    num = length(rv$selected_data$Length)
    last_re_run_type = rv$last_re_run_type
    if (last_re_run_type == "vb" || last_re_run_type == "gompertz" || last_re_run_type == "logistic") {
      if (rv$type == "exp") {
        data = data.frame("Age" = rv$rep_re[6:(5+num), 1])
      } else {
        data = data.frame("Age" = rv$rep_re[7:(6+num), 1])
      }
    } else if (last_re_run_type == "linear") {
      if (rv$type == "exp") {
        data = data.frame("Age" = rv$rep_re[5:(4+num), 1])
      } else {
        data = data.frame("Age" = rv$rep_re[6:(5+num), 1])
      }
    } else { # last_re_run_type == "schnute"
      if (rv$type == "exp") {
        data = data.frame("Age" = rv$rep_re[7:(6+num), 1])
      } else {
        data = data.frame("Age" = rv$rep_re[8:(7+num), 1])
      }
    }

    
    bins = seq(min(data), max(data), length.out = input$bins + 1)
    
    p = ggplot(data, aes(Age)) + 
      geom_histogram(aes(fill = ..count..), breaks = bins) + 
      scale_fill_gradient("Count", low = "#071c07", high = "#35c435") +
      labs(title = "Age Distribution", x = "Age (yr)", y = "Count")
    
    p = p + rv$z_plot_layer
    
    ggplotly(p, tooltip = c("Age", "Count"), source = "hist_RE") %>%
      layout(font = list(family = "'Times New Roman', 'arial', 'Raleway','verdana'")) %>% 
      plotly::config(modeBarButtonsToRemove = 
      list("sendDataToCloud", "zoom2d", "pan2d", "select2d", "lasso2d", "zoomIn2d", "zoomOut2d", "autoScale2d",
           "resetScale2d", "hoverCompareCartesian", "hoverClosestCartesian"), displaylogo = FALSE, collaborate = FALSE)
  })
})