require(dplyr)
require(sf)
require(ggplot2)
require(mgcv)
require(itsadug)
require(cenGAM)
require(brms)
require(fitdistrplus)

Data_analysis<-readRDS("Secchi analysis/Discrete Secchi Analysis Data.Rds")


# Figure out censoring ----------------------------------------------------

secchi_cens<-Data_analysis%>%
  st_drop_geometry()%>%
  filter(is.na(Secchi_estimated) | !Secchi_estimated)%>%
  group_by(Source, Year)%>%
  summarise(max=max(Secchi), u99=quantile(Secchi, 0.99), N=n(), .groups="drop")

ggplot(secchi_cens, aes(x=Year))+
  geom_point(aes(y=max), color="red")+
  geom_point(aes(y=u99), color="blue")+
  facet_wrap(~Source)+
  ylab("Max (red) or 99% quantile (blue) Secchi Depth (cm)")+
  theme_bw()

# Test censoring approaches -----------------------------------------------

## Expectation maximum from baytrends R package

cens_data<-filter(Data_analysis, Source=="FMWT" & Year==2010)%>%
  mutate(censored=if_else(Secchi<200, "none", "right"))



hist(cens_data$Secchi)
cens_dist_data<-data.frame(left=cens_data$Secchi, right=if_else(censored=="right", NA_real_, cens_data$Secchi))
# Normal
fn <- fitdistcens(cens_dist_data, distr="norm")
plot(fn)
# Lognormal
fln <- fitdistcens(cens_dist_data, distr="lnorm")
plot(fln)
# Student's t
ft <- fitdistcens(cens_dist_data, distr=dstudent_t, start=list(mu=0, sigma=2, df=3), control=list(trace=1, report=1))
plot(ft)
# lognormal seems best

mu.lnorm <- fln$estimate[1]
sd.lnorm <- fln$estimate[2]

cens_data_fill<-cens_dist_data%>%
  mutate(right=if_else(is.na(right), Inf, right))

.ExpLNrCens <- function(l,u,mu,sigma,rCens=NA) {
  # computes expected value for right censored for log-normal
  df <- data.frame(l=l, u=u, mu=mu, row.names = NULL)
  
  if(is.na(rCens[1])) {
    df$rcens <- (is.finite(df$l) & is.infinite(df$u))
  } else {
    df$rcens <- rCens
  }
  
  if (any(df$rcens))  {
    df$zl[df$rcens] <- (log(df$l[df$rcens]) - df$mu[df$rcens])/sigma # 01Jun2018
    df$ec[df$rcens] <- exp(df$mu[df$rcens] + sigma**2/2)* (1-stats::pnorm(df$zl[df$rcens]-sigma)) / (1-stats::pnorm(df$zl[df$rcens]))
    df$l[df$rcens] <- df$ec[df$rcens]
    df$u[df$rcens] <- df$ec[df$rcens]
  }
  ExpLNrcens.return <- df[,c('l','u')]
  return(ExpLNrcens.return)
}

cens_data_fill$fill <- .ExpLNrCens(l=cens_data_fill$left, u=cens_data_fill$right, mu=mu.lnorm, sigma=sd.lnorm)[,1]

# This works, but it returns the same number for all censored observations 



# GAMs --------------------------------------------------------------------


SC_gam1_NOAR <- bam(Secchi ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                          te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=Year_s) + 
                          s(Time_num_s, k=5), family=tobit1(right.threshold=470), data = Data_analysis, method="fREML", discrete=T, nthreads=4)
r1 <- start_value_rho(SC_gam1_NOAR, plot=TRUE)


CC_gam1_AR <- bam(Secchi ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                        te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=Year_s) + 
                        s(Time_num_s, k=5), family=scat, rho=r1, AR.start=Start, data = Data_analysis, method="fREML", discrete=T, nthreads=4)



# brms --------------------------------------------------------------------

iterations <- 5e3
warmup <- iterations/4

SC_brm_cens <- brm(Secchi | cens(censored) ~ t2(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(10, 5)) + 
                 s(Time_num_s, k=5), 
               family=student, 
               data = Data_analysis,
               prior=prior(normal(0,10), class="Intercept")+
                 prior(normal(0,5), class="b")+
                 prior(cauchy(0,5), class="sigma"),
               chains=1, cores=1,
               iter = iterations, warmup = warmup,
               backend = "cmdstanr", threads = threading(4))

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
