N = nrow(Hake.dat.otos)
sample = c(1:N)
sex = sample(c('F','M'), N, replace = TRUE, prob = c(0.5, 0.5))
area = sample( LETTERS[1:5], N, replace=TRUE, prob=c(0.17, 0.22, 0.25, 0.19, 0.17))

df = data.frame("Sample" = sample, "Species" = Hake.dat.otos$Species,
           "Area" = area, "Sex" = sex,
           "Length" = Hake.dat.otos$Lt, "Age" = Hake.dat.otos$Age, "OtoWt" = Hake.dat.otos$OtoWt)

write.csv(df, "sample_oto_age_data.csv", row.names = FALSE)