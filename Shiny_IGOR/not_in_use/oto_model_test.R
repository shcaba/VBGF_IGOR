library(ggplot2)

data = read.csv("sample_oto_age_data.csv")
result = oto_age_model(data = data, wt_brkpts = c(50), jitter = 5)
plots = jitter_plot(jitter = 5, oto_age_model = result)

plots[[1]]

downloaded_model = readRDS("~/Downloads/oto_age_model.rds")

model = try(lm(Age ~ poly(OtoWt, 2, raw = TRUE), data = na.omit(data)), 
            silent = TRUE)