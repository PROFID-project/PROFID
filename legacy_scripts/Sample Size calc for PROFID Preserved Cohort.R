library(pmsampsize)

N <- 986
E <- 30 #number of SCD events
Mean_FU <- 56.57084 #mean follow-up time (time-to-SCD or time-to-death or time-to-censoring)
TT <- 55778.84 #person-years of follow-up time (time-to-SCD or time-to-death or time-to-censoring)

lnull <- (E*log(E/TT)) - E
max_R2 <- 1 - exp((2*lnull) / N)

Assumed_C_stat <- 0.6
D <- (5.50*(Assumed_C_stat - 0.5)) + (10.26*(Assumed_C_stat - 0.5)^3)
R2_D <- ((pi/8)*D^2) / (((pi^2)/6) + ((pi/8)*D^2))
R2_quingly <- ((-(pi^2)/6)*R2_D) / (((1-((pi^2)/6))*R2_D) - 1)
LR <- -E*log(1 - R2_quingly)
R2_cox <- 1 - exp(-(LR / N))


p <- 2:10
C_stat_SS <- rep(NA, length(p))
R2_SS <- rep(NA, length(p))
for (i in 1:length(p)) {
  
  C_stat_SS[i] <- pmsampsize(
    type = "s",
    rsquared = R2_cox,
    parameters = p[i],
    shrinkage = 0.9,
    rate = E / TT,
    timepoint = 12,
    meanfup = Mean_FU
  )$sample_size
  
  R2_SS[i] <- pmsampsize(
    type = "s",
    rsquared = 0.15*max_R2,
    parameters = p[i],
    shrinkage = 0.9,
    rate = E / TT,
    timepoint = 12,
    meanfup = Mean_FU
  )$sample_size
  
}
C_stat_SS
R2_SS






N <- 986
E <- 30 #number of SCD events
Mean_FU <- 56.57084 #mean follow-up time (time-to-SCD or time-to-death or time-to-censoring)
TT <- 55778.84 #person-years of follow-up time (time-to-SCD or time-to-death or time-to-censoring)

lnull <- (E*log(E/TT)) - E
max_R2 <- 1 - exp((2*lnull) / N)

Assumed_C_stat <- 0.75
D <- (5.50*(Assumed_C_stat - 0.5)) + (10.26*(Assumed_C_stat - 0.5)^3)
R2_D <- ((pi/8)*D^2) / (((pi^2)/6) + ((pi/8)*D^2))
R2_quingly <- ((-(pi^2)/6)*R2_D) / (((1-((pi^2)/6))*R2_D) - 1)
LR <- -E*log(1 - R2_quingly)
R2_cox <- 1 - exp(-(LR / N))


p <- 2:10
C_stat_SS <- rep(NA, length(p))
R2_SS <- rep(NA, length(p))
for (i in 1:length(p)) {
  
  C_stat_SS[i] <- pmsampsize(
    type = "s",
    rsquared = R2_cox,
    parameters = p[i],
    shrinkage = 0.9,
    rate = E / TT,
    timepoint = 12,
    meanfup = Mean_FU
  )$sample_size
  
  R2_SS[i] <- pmsampsize(
    type = "s",
    rsquared = 0.15*max_R2,
    parameters = p[i],
    shrinkage = 0.9,
    rate = E / TT,
    timepoint = 12,
    meanfup = Mean_FU
  )$sample_size
  
}
C_stat_SS
R2_SS
