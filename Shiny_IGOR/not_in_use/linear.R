rv$linear_summaries

# linear model
model = lm(len_use ~ age_use)
rep = summary(model)[["coefficients"]]
newRow = format(data.frame("Run name" = input$runname_nls,
                           "Intercept mean" = rep[1,1],
                           "Intercept Std" = rep[1,2],
                           "Slope mean" = rep[2,1],
                           "Slope Std" = rep[2,2],
                           "CV Length Error" = NA,
                           "CV Age Error" = NA,
                           "alpha" = NA,
                           "sigma" = NA,
                           "beta" = NA,
                           "shape" = NA,
                           "scale" = NA,
                           "z" = NA,
                           check.names = FALSE
                ), digits = 4
         )

data = list(len = len_use, age = age_use)
CV_Lt = as.numeric(input$CV_const)
parameters = list(intercept = input$intercept_nls, slope = input$slope_nls, CV_Lt = CV_Lt)
obj = MakeADFun(data = data, parameters = parameters, DLL = "linear_likelihood")
lower = c(-15, 0, exp(1) ^ (-10))
upper = c(3 * max(len_use), 1000, exp(1) ^ (2))
opt = opt(obj, lower, upper)

# return NULL if model didn't converge
if (is.null(opt)) return(NULL)

rep = summary(sdreport(obj))
newRow = format(data.frame("Run name" = input$runname_nls,
                           "Starting intercept" = input$intercept_nls,
                           "Starting slope" = input$slope_nls,
                           "Intercept mean" = rep[1,1],
                           "Intercept Std" = rep[1,2],
                           "Slope mean" = rep[2,1],
                           "Slope Std" = rep[2,2],
                           "CV Length Error" = rep[3,1],
                           "CV Age Error" = NA,
                           "alpha" = NA,
                           "sigma" = NA,
                           "beta" = NA,
                           "shape" = NA,
                           "scale" = NA,
                           "z" = NA,
                           check.names = FALSE), digits = 4)


rv$linear_summaries = rbind(rv$linear_summaries, newRow)