library("TMB")
compile("vbre_Gamma.cpp")
compile("vbre_Exponential.cpp")
compile("vbre_Normal.cpp")

age_re_likelihood = "Exponential"
dll = paste("vbre", age_re_likelihood, sep = "_")
dyn.load(dynlib(dll))

fish_re = read.csv("sample.csv")
# num_reads = 3
# len = fish_re$Length
# num = length(len)
# #age = t(matrix(c(fish_re$Read1, fish_re$Read2, fish_re$Read3), nrow = num, ncol = num_reads))
# age = t(fish_re[c(-1, -2)])

num_reads = ncol(fish_re) - 4
len = fish_re$Length
num = length(len)
# age matrix: nrow = #reads, ncol = num
age = t(fish_re[c(-1, -2, -3, -4)])

coeff_var = function(v) {sd(v) / mean(v)}
cv = apply(age, 2, coeff_var)
CV_e = sqrt(sum(cv * cv) / num)

data = list(
  age = age,
  len = len,
  CV_e = CV_e,
  num_reads = num_reads
)

if (age_re_likelihood == "Normal") {
  parameters = list(linf = 912, kappa = 0.04, t0 = 0, CV_Lt = 0.1, alpha = 5, sigma_age = 1, age_re = rep(mean(age), num))
  
  # lower bounds: linf, kappa, t0, CV_Lt, alpha, sigma_age, age_re
  lower = c(0.75 * max(len), 0.0001, -15, exp(1) ^ (-10), 0, exp(1) ^ (-100), rep(min(age), num))
  # upper bounds
  upper = c(3 * max(len), 1, 1, exp(1) ^ (2), 500, exp(1) ^ (100), rep(max(age), num))
} else if (age_re_likelihood == "Exponential") {
  parameters = list(linf = 912, kappa = 0.04, t0 = 0, CV_Lt = 0.1, beta = 0.2, age_re = rep(mean(age), num))

  # lower bounds: linf, kappa, t0, CV_Lt, beta, age_re
  lower = c(0.75 * max(len), 0.0001, -15, exp(1) ^ (-10), exp(1) ^ (-100), rep(min(age), num))
  # upper bounds
  upper = c(3 * max(len), 1, 0, exp(1) ^ (2), exp(1) ^ (100), rep(max(age), num) )
} else { #gamma
  parameters = list(linf = 912, kappa = 0.04, t0 = 0, CV_Lt = 0.1, gam_shape = 5, gam_scale = 3, age_re = rep(mean(age), num))
  
  # lower bounds: linf, kappa, t0, CV_Lt, gam_shape, gam_scale, age_re
  lower = c(0.75 * max(len), 0.0001, -15, exp(1) ^ (-10), 0, 0, rep(min(age), num))
  # upper bounds
  upper = c(3 * max(len), 1, 0, exp(1) ^ (2), 100, 100, rep(max(age), num))
}

obj = MakeADFun(data = data, parameters = parameters, random = "age_re", DLL = dll)
opt = nlminb(obj$par, obj$fn, obj$gr, lower = lower, upper = upper)
rep = summary(sdreport(obj))

plot(age[1,], len, xlab = "Age (yr)", ylab = "Length (mm)", col = "green")
points(age[2,], len, xlab = "Age (yr)", ylab = "Length (mm)", col = "blue")
points(age[3,], len, xlab = "Age (yr)", ylab = "Length (mm)", col = "purple")
curve(rep[1,1] * (1 - exp(-rep[2,1] * (x - rep[3,1]))), from = 0 , to = max(age), col = 2, add = T)
