# libraries
library(shiny)
library(jsonlite)
library(shinycssloaders)
library(shinyjs)
library(segmented)
library(ggplot2)
library(RColorBrewer)
library(plotly)
library(DT)
library(reshape2)
library(TMB)
library(htmlwidgets)
library(htmltools)

source("oto_age_model_function.R")

compile("tmb_models/vbre_Gamma.cpp")
compile("tmb_models/vbre_Exponential.cpp")
compile("tmb_models/vbre_Normal.cpp")
compile("tmb_models/vb_likelihood.cpp")
compile("tmb_models/linear_likelihood.cpp")
compile("tmb_models/schnute_likelihood.cpp")
compile("tmb_models/logistic_likelihood.cpp")
compile("tmb_models/gompertz_likelihood.cpp")
compile("tmb_models/linear_Normal.cpp")
compile("tmb_models/linear_Exponential.cpp")
compile("tmb_models/linear_Gamma.cpp")
compile("tmb_models/gompertz_Normal.cpp")
compile("tmb_models/gompertz_Exponential.cpp")
compile("tmb_models/gompertz_Gamma.cpp")
compile("tmb_models/schnute_Normal.cpp")
compile("tmb_models/schnute_Exponential.cpp")
compile("tmb_models/schnute_Gamma.cpp")
compile("tmb_models/logistic_Normal.cpp")
compile("tmb_models/logistic_Exponential.cpp")
compile("tmb_models/logistic_Gamma.cpp")

dyn.load(dynlib("tmb_models/vbre_Gamma"))
dyn.load(dynlib("tmb_models/vbre_Exponential"))
dyn.load(dynlib("tmb_models/vbre_Normal"))
dyn.load(dynlib("tmb_models/vb_likelihood"))
dyn.load(dynlib("tmb_models/linear_likelihood"))
dyn.load(dynlib("tmb_models/schnute_likelihood"))
dyn.load(dynlib("tmb_models/logistic_likelihood"))
dyn.load(dynlib("tmb_models/gompertz_likelihood"))
dyn.load(dynlib("tmb_models/logistic_Normal"))
dyn.load(dynlib("tmb_models/logistic_Exponential"))
dyn.load(dynlib("tmb_models/logistic_Gamma"))
dyn.load(dynlib("tmb_models/schnute_Normal"))
dyn.load(dynlib("tmb_models/schnute_Exponential"))
dyn.load(dynlib("tmb_models/schnute_Gamma"))
dyn.load(dynlib("tmb_models/linear_Normal"))
dyn.load(dynlib("tmb_models/linear_Exponential"))
dyn.load(dynlib("tmb_models/linear_Gamma"))
dyn.load(dynlib("tmb_models/gompertz_Normal"))
dyn.load(dynlib("tmb_models/gompertz_Exponential"))
dyn.load(dynlib("tmb_models/gompertz_Gamma"))

helpData = read.csv("help.csv")

# computes coefficient of variation for a vector
coeff_var <- function(v) {sd(v, na.rm = TRUE) / mean(v, na.rm = TRUE)}

# report convergence failed message
convg_err <- function() {showModal(modalDialog(
  title = "Important message",
  "The model didn't converge. Check your starting values.",
  easyClose = TRUE
))}

# report invalid file upload
invalid_upload <- function() {
  showModal(
    modalDialog(
      title = "Important message",
      "The file you uploaded does not have the correct format. Refer to the sample file.",
      easyClose = TRUE
    )
  )
}

# report quantile value out of range
invalid_quantile <- function() {
  showModal(
    modalDialog(
      title = "Important message",
      "Enter a quantile value between 0 and 1, inclusive",
      easyClose = TRUE
    )
  )
}

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

# plot data points and the fitting curve for age and length
plot_growth_curve <- function(rv, mode, model, age_unit, len_unit) {
  age_max = max(rv$selected_data[c(-1, -2, -3, -4, -5)])
  # get current selected reads (the ones in the filtered data table)
  reads_choices = names(rv$selected_data[c(-1, -2, -3, -4, -5)])
  if (length(reads_choices) < 3) {
    # 2 greens that are not too pale
    palette = c("#41AB5D", "#238B45")
  } else {
    palette = brewer.pal(length(reads_choices), "Greens")
  }
  colors = c("black", palette)
  names(colors) = c("Expected", reads_choices)
  
  p = ggplot(data = rv$selected_data) +
     labs(title = "Age and Growth Fit", x = paste0("Age (", age_unit, ")"), y = paste0("Length (", len_unit, ")"))
    
  for (i in reads_choices) {
    p = p + geom_point(mapping = aes_string(Species = "Species", Area = "Area", Sex = "Sex", x = i, y= "Length", color = paste("'", i, "'", sep = "")))
  }
  if (mode == "analyze" ) {
    if (model == "nls") {
      model_name = rv$last_run_type
      rep = rv$rep_nls
    } else {
      model_name = rv$last_re_run_type
      rep = rv$rep_re
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
  return(ggplotly(p, tooltip = c("x", "y", "Species", "Area", "Sex"))  %>%
           layout(font = list(family = "'Times New Roman', 'arial', 'Raleway','verdana'")) %>%
           plotly::config(modeBarButtonsToRemove = 
                            list("sendDataToCloud", "zoom2d", "pan2d", "select2d", "lasso2d", "zoomIn2d", "zoomOut2d", "autoScale2d",
                                 "resetScale2d", "hoverCompareCartesian", "hoverClosestCartesian"), displaylogo = FALSE, collaborate = FALSE))
}
