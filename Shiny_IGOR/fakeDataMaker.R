###### Simulate data
N = 100 # Number of observations
age_min = 0 # Minimum age in the population
age_max = 10 # Maximum age in the population

# Parameters
Obs_sd = 0.5 # Standard deviation
linf = 350 # Asymptoptic length/mean length at maximum age (mm)
kappa = .3 # Growth rate (yr^-1)
t0 = -1.5 # Estimated age at length 0 (yr)

model = "sch"

# Simulation
age1 = round(runif(N, age_min, age_max))
age2 = round(rnorm(N, age1, 0.27))
age3 = round(rnorm(N, age1, 0.3))
age4 = round(rnorm(N, age1, 0.22))
age5 = round(rnorm(N, age1, 0.18))
age6 = round(rnorm(N, age1, 0.16))
age7 = round(rnorm(N, age1, 0.2))

if (model == "gom") {
  len = round(rnorm(N, 800 * exp(-83 * exp(-0.5 * age1)), Obs_sd))
} else if (model == "sch") {
  len = round(rnorm(N, (10 + 5 * exp(0.5 * age1)) ** 0.5, Obs_sd))
  m = nls(len ~ (r0 + b * exp(k * age)) ** m,
      data = list(age = age1, len = len),
      start = list(r0 = 10, b = 5, k = 0.5, m = 0.5),
      control = list(reltol=0.00000000001))
  rep = summary(m)[["parameters"]]
} else {
  len = round(rnorm(N, linf * (1 - exp(-kappa * (age1 - t0))), Obs_sd))
}

sample = c(1:N)
area = sample( LETTERS[1:5], N, replace=TRUE, prob=c(0.17, 0.22, 0.25, 0.19, 0.17))
sex = sample(c('F','M'), N, replace = TRUE, prob = c(0.5, 0.5))

df = data.frame("Sample" = sample, "Area" = area, "Sex" = sex, "Length" = len,
                 "Read1" = age1, "Read2" = age2, "Read3" = age3, "Read4" = age4,
                 "Read5" = age5, "Read6" = age6, "Read7" = age7)

write.csv(df, paste0(model, "_fake.csv"), row.names = FALSE)