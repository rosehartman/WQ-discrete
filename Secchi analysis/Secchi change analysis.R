require(dplyr)
require(sf)
require(ggplot2)
require(mgcv)
require(itsadug)
require(cenGAM)
require(brms)

Data_analysis<-readRDS("Secchi analysis/Discrete Secchi Analysis Data.Rds")

SC_gam1_NOAR <- bam(Secchi ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                          te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=Year_s) + 
                          s(Time_num_s, k=5), family=tobit1(right.threshold=470), data = Data_analysis, method="fREML", discrete=T, nthreads=4)
r1 <- start_value_rho(SC_gam1_NOAR, plot=TRUE)

#########Best Model####################

CC_gam1_AR <- bam(Secchi ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                        te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=Year_s) + 
                        s(Time_num_s, k=5), family=scat, rho=r1, AR.start=Start, data = Data_analysis, method="fREML", discrete=T, nthreads=4)

iterations <- 5e3
warmup <- iterations/4

SC_brm1 <- brm(Secchi ~ t2(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                      t2(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=Year_s) + 
                      s(Time_num_s, k=5), 
               family=student, 
               data = Data_analysis,
               prior=prior(normal(0,10), class="Intercept")+
                 prior(normal(0,5), class="b")+
                 prior(cauchy(0,5), class="sigma"),
               chains=1, cores=1,
               iter = iterations, warmup = warmup,
               backend = "cmdstanr", threads = threading(4))
