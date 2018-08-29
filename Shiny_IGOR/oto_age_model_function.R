# Computes a piecewise model for the given data
# data: dataframe with OtoWt and Age columns at least
# jitter: number of times to jitter
# wt_brkpts: vector of initial guess for breakpoints by otolith weight
# returns a piecewise model of the given data, returns NULL on error
oto_age_model <- function(data, wt_brkpts = c(0, 0), jitter = 0) {
  # list of piecewise models of size jitter
  jitters = list()
  i = 0
  
  # computes the linear model later used for piecewise model
  linear_model = try(lm(Age ~ OtoWt, data = data), silent = TRUE)
  
  if (class(linear_model)[[1]] == "try-error") {
    return(NULL)
  }
  
  if (jitter > 0) {
    for (j in 1 : jitter) {
      jitters[[j]] = try(segmented(linear_model, seg.Z = ~OtoWt, psi = wt_brkpts, jt = TRUE), silent = TRUE)
      if (class(jitters[[j]])[1] == "try-error") {
        # discard missing values in finding minmum model AIC later
        jitters[[j]] = NA
      }
    }
    # choose the model with minimum AIC
    piecewise_model = try(jitters[[which.min(sapply(jitters, AIC))]], silent = TRUE)
  } else {
    repeat{
      piecewise_model = try(segmented(linear_model, seg.Z = ~OtoWt, psi = wt_brkpts), silent = TRUE)
      i = i + 1
      if (class(piecewise_model)[1] != "try-error" || i >= 10) break
    }
  }
  
  if (class(piecewise_model)[1] == "try-error") {
    piecewise_model = NULL
  } else if (jitter > 0) {
    piecewise_model[["jitter_all"]] = jitters
    piecewise_model[["no_jitter_model"]] = oto_age_model(data, wt_brkpts, lt_brkpts, jitter = 0, add_lengths)
  }
  
  return (piecewise_model)
}

# p: ggplot2 object
# returns a plotly object that has unnecessary buttons removed
clean_plot <- function(p) {
  l = list(orientation = "h", x = 0, y = 0) # horizontal legends on the bottom
  return(ggplotly(p, tooltip = c("x", "y", "yintercept")) %>%  
           layout(font = list(family = "'Times New Roman', 'arial', 'Raleway','verdana'"),
                  legend = l) %>%
    plotly::config(modeBarButtonsToRemove = list("sendDataToCloud", "zoom2d", "pan2d", "select2d", 
          "lasso2d", "zoomIn2d", "zoomOut2d", "autoScale2d", "resetScale2d", "hoverCompareCartesian", 
          "hoverClosestCartesian"), displaylogo = FALSE, collaborate = FALSE))
}

# returns a list of 6 plots: jittered intercept, slope1, slope2, breakpoint, AIC, 
# relative to unjittered values
jitter_plot <- function(jitter, oto_age_model) {
  no_jitter_model = oto_age_model$no_jitter_model
  if (is.null(no_jitter_model)) return()
  
  # pre: jitter > 0
  #      jitter_all contains a smallest AIC model
  #      no_jitter_model != NULL
  
  table = matrix(nrow = jitter, ncol = 5)
  # table looks like this
  # -------------------------------------------------------------
  #         intercept slope1 slope2 breakpoint AIC
  # jitter1
  # jitter2
  # ...
  
  x_vals = c(1 : jitter)
  jitter_all = oto_age_model$jitter_all
  for (i in 1 : length(jitter_all)) {
    model = jitter_all[[i]]
    if (!is.na(model)) {
      table[i, 1] = model$coefficients[[1]] # intercept
      table[i, 2] = slope(model)$OtoWt[1, 1] # slope1
      table[i, 3] = slope(model)$OtoWt[2, 1] # slope2
      table[i, 4] = model$psi[1, 2] # breakpoint
      table[i, 5] = AIC(model) # AIC
    }
  }
  
  # from no_jitter_model
  intercept = no_jitter_model$coefficients[[1]]
  slope1 = slope(no_jitter_model)$OtoWt[1, 1]
  slope2 = slope(no_jitter_model)$OtoWt[2, 1]
  breakpoint = no_jitter_model$psi[1, 2]
  AIC = AIC(no_jitter_model)
  
  p1 = ggplot(data = data.frame(x = x_vals, y = table[,1])) +
          geom_point(aes(x = x, y = y)) +
          labs(x = "jitter number", y = "intercept")
  min_intercept = min(table[,1], na.rm = TRUE)
  if (min_intercept >= 0.9 * intercept && max(table[,1], na.rm = TRUE) <= 1.1 * intercept) {
    p1 = p1 + ylim(0.9 * intercept, 1.1 * intercept)
  }
  p1 = p1 + geom_hline(aes(yintercept = min_intercept, color = "jitter minimum")) + 
          geom_hline(aes(yintercept = intercept, color = "base case")) +
          scale_color_manual(values = c("base case" = "green", "jitter minimum" = "red"), name = "")
  
  p2 = ggplot(data = data.frame(x = x_vals, y = table[,2])) +
          geom_point(aes(x = x, y = y)) +
          labs(x = "jitter number", y = "slope1")
  min_slope1 = min(table[,2], na.rm = TRUE)
  if (min_slope1 >= 0.9 * slope1 && max(table[,2], na.rm = TRUE) <= 1.1 * slope1) {
    p2 = p2 + ylim(0.9 * slope1, 1.1 * slope1)
  }
  p2 = p2 + geom_hline(aes(yintercept = min_slope1, color = "jitter minimum")) + 
          geom_hline(aes(yintercept = slope1, color = "base case")) +
          scale_color_manual(values = c("base case" = "green", "jitter minimum" = "red"), name = "")
  
  p3 = ggplot(data = data.frame(x = x_vals, y = table[,3])) +
          geom_point(aes(x = x, y = y)) +
          labs(x = "jitter number", y = "slope2")
  min_slope2 = min(table[,3], na.rm = TRUE)
  if (min_slope2 >= 0.9 * slope2 && max(table[,3], na.rm = TRUE) <= 1.1 * slope2) {
    p3 = p3 + ylim(0.9 * slope2, 1.1 * slope2)
  }
  p3 = p3 + geom_hline(aes(yintercept = min_slope2, color = "jitter minimum")) + 
          geom_hline(aes(yintercept = slope2, color = "base case")) +
          scale_color_manual(values = c("base case" = "green", "jitter minimum" = "red"), name = "")
  
  p4 = ggplot(data = data.frame(x = x_vals, y = table[,4])) +
          geom_point(aes(x = x, y = y)) +
          labs(x = "jitter number", y = "breakpoint")
  min_breakpoint = min(table[,4], na.rm = TRUE)
  if (min_breakpoint >= 0.9 * breakpoint && max(table[,4], na.rm = TRUE) <= 1.1 * breakpoint) {
    p4 = p4 + ylim(0.9 * breakpoint, 1.1 * breakpoint)
  }
  p4 = p4 + geom_hline(aes(yintercept = min_breakpoint, color = "jitter minimum")) + 
          geom_hline(aes(yintercept = breakpoint, color = "base case")) +
          scale_color_manual(values = c("base case" = "green", "jitter minimum" = "red"), name = "")
  
  p5 = ggplot(data = data.frame(x = x_vals, y = abs(table[,5] - AIC))) +
          geom_point(aes(x = x, y = y)) +
          labs(x = "jitter number", y = "absolute AIC difference")
  min_AIC_diff = min(abs(table[,5] - AIC), na.rm = TRUE)
  p5 = p5 + geom_hline(aes(yintercept = min_AIC_diff, color = "jitter minimum")) +
          geom_hline(aes(yintercept = 1.96, color = "significance level")) +
          scale_color_manual(values = c("jitter minimum" = "red", "significance level" = "blue"), name = "")
  
  p6 = ggplot(data = data.frame(x = no_jitter_model$fitted.values, y = oto_age_model$fitted.values)) +
         geom_point(aes(x = x, y = y)) +
         xlim(0, max(no_jitter_model$fitted.values) + 2) +
         ylim(0, max(oto_age_model$fitted.values) + 2) + 
         labs(x = "benchmark age values", y = "Alternative Parameter Age values") +
         geom_abline(intercept = 0, slope = 1, color = "blue")
  
  plots = list(clean_plot(p1), clean_plot(p2), clean_plot(p3), clean_plot(p4), clean_plot(p5), clean_plot(p6))
  
  return(plots)
}