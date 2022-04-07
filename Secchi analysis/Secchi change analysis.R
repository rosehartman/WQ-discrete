require(dplyr)
require(sf)
require(ggplot2)
require(mgcv)
require(itsadug)
require(cenGAM)
require(brms)
require(fitdistrplus)
require(ggridges)
require(patchwork)

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

cens_plot<-function(source){
  data_sum<-filter(secchi_cens, Source==source)
  
  p1<-ggplot(data_sum, aes(x=Year))+
    geom_point(aes(y=max), color="red")+
    geom_point(aes(y=u99), color="blue")+
    ggtitle(paste0(source, ": ", min(data_sum$Year), "-", max(data_sum$Year), 
                   " (", length(min(data_sum$Year):max(data_sum$Year)), " years)"))+
    ylab("Max (red) or 99% quantile (blue) Secchi Depth (cm)")+
    theme_bw()
  
  p2<-ggplot(filter(Data_analysis, Source==source), aes(x=Secchi, y=Year, group=Year))+
    geom_density_ridges()+
    theme_bw()
  p1|p2
}

# Specify censoring levels ------------------------------------------------
max_secchi<-max(Data_analysis$Secchi)*2

censored_years<-tibble(
  Source=c(rep("20mm", 27), rep("Baystudy", 41), rep("DJFMP", 45), rep("EMP", 46),
           rep("FMWT", 52), rep("SDO", 18), rep("SKT", 20), rep("SLS", 13), 
           rep("STN", 52), rep("Suisun", 42), rep("YBFMP", 20)),
  Year=c(1995:2021,
         1980:2020,
         1976:2020,
         1975:2020,
         1969:2020,
         2001:2018,
         2002:2021,
         2009:2021,
         1970:2021,
         1980:2021,
         1999:2018),
  Censor_level=c(rep(200, 27), 
                 rep(200, 40), 300,
                 rep(200, 36), rep(180, 9),
                 rep(200, 26), rep(max_secchi, 20),
                 rep(200, 44), rep(max_secchi, 8),
                 rep(max_secchi, 18),
                 rep(200, 20),
                 rep(200, 13),
                 rep(99, 11), rep(200, 33), rep(300, 8),
                 rep(max_secchi, 42),
                 rep(300, 17), rep(max_secchi, 3))
)


# Join censoring levels to data -------------------------------------------

Data_analysis_cens<-Data_analysis%>%
  left_join(censored_years, by=c("Source", "Year"))%>%
  mutate(Censored=if_else(Secchi>=Censor_level, "right", "none"),
         Secchi_cens=if_else(Censored=="right", Censor_level, Secchi))

# Test censoring approaches -----------------------------------------------

## Expectation maximum from baytrends R package

cens_data<-filter(Data_analysis, Source=="FMWT" & Year==2010)%>%
  mutate(censored=if_else(Secchi<200, "none", "right"))



hist(cens_data$Secchi)
cens_dist_data<-data.frame(left=Data_analysis_cens$Secchi_cens, right=if_else(Data_analysis_cens$Censored=="right", NA_real_, Data_analysis_cens$Secchi_cens))%>%
  filter(left>0)
# Normal
fn <- fitdistcens(cens_dist_data, distr="norm")
plot(fn)
# Lognormal
fln <- fitdistcens(cens_dist_data, distr="lnorm")
plot(fln)
# Student's t
ft <- fitdistcens(cens_dist_data, distr=dstudent_t, start=list(mu=0, sigma=2, df=3))
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


SC_gam1_NOAR <- bam(Secchi_cens ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                      te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=Year_s) + 
                      s(Time_num_s, k=5), family=tobit1(right.threshold=470), data = Data_analysis, method="fREML", discrete=T, nthreads=4)
r1 <- start_value_rho(SC_gam1_NOAR, plot=TRUE)


CC_gam1_AR <- bam(Secchi_cens ~ te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                    te(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=Year_s) + 
                    s(Time_num_s, k=5), family=scat, rho=r1, AR.start=Start, data = Data_analysis, method="fREML", discrete=T, nthreads=4)



# brms --------------------------------------------------------------------

iterations <- 5e3
warmup <- iterations/4

SC_brm_cens <- brm(Secchi_cens | cens(Censored) ~ t2(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("cs", "cc"), k=c(25, 13)) + 
                     t2(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("cs", "cc"), k=c(25, 13), by=Year_s) + 
                     s(Time_num_s, k=5), 
                   family=student, 
                   data = Data_analysis_cens,
                   prior=prior(normal(0,10), class="Intercept")+
                     prior(normal(0,5), class="b")+
                     prior(cauchy(0,5), class="sigma"),
                   chains=1, cores=1,
                   iter = iterations, warmup = warmup,
                   backend = "cmdstanr", threads = threading(15))

SC_brm1 <- brm(Secchi_cens ~ t2(Latitude_s, Longitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
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
