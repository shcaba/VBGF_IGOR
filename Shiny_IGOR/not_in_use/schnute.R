rv$schnute_summaries
# Gompertz

model = nls(len ~ (r0 + b * exp(k * age)) ** m,
            data = list(age = age_use, len = len_use),
            start = list(r0 = input$sch_r0_nls, b = input$sch_b_nls, k = input$sch_k_nls, m = input$sch_m_nls),
            control = list(reltol=0.00000000001))

rep = summary(model)[["parameters"]]
newRow = format(data.frame("Run name" = input$runname_nls,
                           "Starting r0" = input$sch_r0_nls,
                           "Starting b" = input$sch_b_nls,
                           "Starting k" = input$sch_k_nls,
                           "Starting m" = input$sch_m_nls,
                           "r0 mean" = rep[1,1],
                           "r0 Std" = rep[1,2],
                           "b mean" = rep[2,1],
                           "b Std" = rep[2,2],
                           "k mean" = rep[3,1],
                           "k Std" = rep[3,2],
                           "m mean" = rep[4,1],
                           "m Std" = rep[4,2],
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
parameters = list(r0 = input$sch_r0_nls, b = input$sch_b_nls, k = input$sch_k_nls, m = input$sch_m_nls, CV_Lt = CV_Lt)
obj = MakeADFun(data = data, parameters = parameters, DLL = "schnute_likelihood")
lower = c(-Inf, -Inf, -Inf, -Inf, exp(1) ^ (-10))
upper = c(Inf, Inf, Inf, Inf, exp(1) ^ (2))
opt = opt(obj, lower, upper)

# return NULL if model didn't converge
if (is.null(opt)) return(NULL)

rep = summary(sdreport(obj))
newRow = format(data.frame("Run name" = input$runname_nls,
                           "Starting r0" = input$sch_r0_nls,
                           "Starting b" = input$sch_b_nls,
                           "Starting k" = input$sch_k_nls,
                           "Starting m" = input$sch_k_nls,
                           "r0 mean" = rep[1,1],
                           "r0 Std" = rep[1,2],
                           "b mean" = rep[2,1],
                           "b Std" = rep[2,2],
                           "k mean" = rep[3,1],
                           "k Std" = rep[3,2],
                           "m mean" = rep[4,1],
                           "m Std" = rep[4,2],
                           "CV Length Error" = rep[5,1],
                           "CV Age Error" = NA,
                           "alpha" = NA,
                           "sigma" = NA,
                           "beta" = NA,
                           "shape" = NA,
                           "scale" = NA,
                           "z" = NA,
                           check.names = FALSE), digits = 4)

rv$schnute_summaries = rbind(rv$schnute_summaries, newRow)