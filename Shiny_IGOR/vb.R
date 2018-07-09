setwd("~/Rprojects/IGOR/model_wk2/")
library("TMB")
compile("vb_likelihood.cpp")
compile("vb_multiple_ages.cpp")

fish_re = read.csv("dat_re.csv")
num_reads = 3
len = fish_re$Length
num = length(len)
# age matrix: nrow = #reads, ncol = num
# age = t(matrix(c(fish_re$Read1, fish_re$Read2, fish_re$Read3), nrow = num, ncol = num_reads))
age = t(fish_re[c(-1, -2)])

age_max = max(age)

fit_data = 1

if (fit_data < 4) {
  # use R least square function for fit
  if (fit_data == 1) {
    age_use = age["Read1",]
  } else if (fit_data == 2) {
    # fit to average age
    # create an age vector that is the mean across all reads
    age_use = apply(age, 2, mean)
  } else {
    # fit to median age
    # create an age vector that is the median across all reads
    age_use = apply(age, 2, median)
  }
  model = nls(len ~ linf * (1 - exp(-kappa * (age - t0))),
              data = list(age = age_use, len = len),
              start = list(linf = max(len), kappa = 0.1, t0 = 0),
              control = list(reltol=0.00000000001)
  )
  summ = summary(model)
  linf = coef(model)[[1]]
  kappa = coef(model)[[2]]
  t0 = coef(model)[[3]]
 } else {
  if (fit_data == 4) {
    # fit to likelihood
    dyn.load(dynlib("vb_likelihood"))
    data = list(len = len, age = age)
    parameters = list(linf = 912, kappa = 0.04, t0 = 0, CV_Lt = 0.1)
    obj = MakeADFun(data = data, parameters = parameters, DLL = "vb_likelihood")
  } else { # fit_data == 5
    # fit to multiple reads (regression)
    dyn.load(dynlib("vb_multiple_ages"))
    data = list(len = len, age = age, num_reads = num_reads)
    parameters = list(linf = 518, kappa = 0.2, t0 = 0)
   obj = MakeADFun(data = data, parameters = parameters, DLL = "vb_multiple_ages")
  }
  opt = nlminb(obj$par, obj$fn, obj$gr)
  rep = summary(sdreport(obj))
  linf = rep[1,1]
  kappa = rep[2,1]
  t0 = rep[3,1]
}
plot(age[1,], len , xlab = "Age (yr)", ylab = "Length (mm)")
curve(linf * (1 - exp(-kappa * (x - t0))), from = 0 , to = age_max, col = 2, add=T)