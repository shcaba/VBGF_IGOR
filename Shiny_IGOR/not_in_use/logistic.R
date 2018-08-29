rv$logistic_summaries
# Gompertz

model = nls(len ~ a / (1 + b * exp(-k * age)),
            data = list(age = age_use, len = len_use),
            start = list(a = input$log_a_nls, b = input$log_b_nls, k = input$log_k_nls),
            control = list(reltol=0.00000000001))

rep = summary(model)[["parameters"]]
newRow = format(data.frame("Run name" = input$runname_nls,
                           "Starting a" = input$log_a_nls,
                           "Starting b" = input$log_b_nls,
                           "Starting k" = input$log_k_nls,
                           "a mean" = rep[1,1],
                           "a Std" = rep[1,2],
                           "b mean" = rep[2,1],
                           "b Std" = rep[2,2],
                           "k mean" = rep[3,1],
                           "k Std" = rep[3,2],
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
parameters = list(a = input$log_a_nls, b = input$log_b_nls, k = input$log_k_nls, CV_Lt = CV_Lt)
obj = MakeADFun(data = data, parameters = parameters, DLL = "logistic_likelihood")
lower = c(-Inf, -Inf, -Inf, exp(1) ^ (-10))
upper = c(Inf, Inf, Inf, exp(1) ^ (2))
opt = opt(obj, lower, upper)

# return NULL if model didn't converge
if (is.null(opt)) return(NULL)

rep = summary(sdreport(obj))
newRow = format(data.frame("Run name" = input$runname_nls,
                           "Starting a" = input$log_a_nls,
                           "Starting b" = input$log_b_nls,
                           "Starting k" = input$log_k_nls,
                           "a mean" = rep[1,1],
                           "a Std" = rep[1,2],
                           "b mean" = rep[2,1],
                           "b Std" = rep[2,2],
                           "k mean" = rep[3,1],
                           "k Std" = rep[3,2],
                           "CV Length Error" = rep[4,1],
                           "CV Age Error" = NA,
                           "alpha" = NA,
                           "sigma" = NA,
                           "beta" = NA,
                           "shape" = NA,
                           "scale" = NA,
                           "z" = NA,
                           check.names = FALSE), digits = 4)



rv$logistic_summaries = rbind(rv$logistic_summaries, newRow)