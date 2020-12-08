
require(dplyr)
require(tidyr)

Data<-readRDS("Temperature analysis/Discrete Temp Data.Rds")

n<-1e4

sin_sim<-function(years, cycles=10, sin_denom=3, exp_scale=5, exp_sub=7, lin_slope=300){
  x<-runif(years, -cycles*pi, cycles*pi)
  y<-sin(x)/sin_denom+exp((x)/exp_scale-exp_sub)+x/lin_slope
  return(tibble(x,y))
}

Time_sim<-function(n){
  n2<-n/2-0.5
  ((-1*(-n2:n2)^2)+n2^2)/n2^2
}

season<-cos(seq(-1*pi, pi, length.out=13))[1:12]
Year_means<-exp((1:51)/60)+1:51/500
#Time_means<-
Sim_data<-expand_grid(Month=1:12, Year=1970:2020, Location=1:50, Source=1:5)%>%
  mutate(Time_num_s=rnorm(n()),
         Year_num=as.integer(as.factor(Year)),
         Temp=rnorm(n(), Year_means[Year_num]+season[Month]*10))



# New method: Try simulating data (using the gratia simulate or predicted_sampled functions) from the global smoothing model to use for testing the model performance.
# Simulate with a balanced and imbalanced sampling design to examing that impact
# Also try simulating with all years vs. just years with similar average water temps to ensure model can differentiate real from spurious trends.

# First calculate year-to-year variability
# Start in 1975 when sampling became more regular across different regions
Data_CC4<-readRDS("Temperature analysis/Data_CC4.Rds")

require(brms)

model_var<-brm(Temperature ~ (1|Year_fac) + (1| Month) + (1|Station), data=Data_CC4,
               family=gaussian(),
               prior=prior(normal(0,5), class="Intercept")+
                 prior(cauchy(0,5), class="sigma")+
                 prior(cauchy(0,5), class="sd"),
               control=list(adapt_delta=0.85, max_treedepth=15),
               iter=8e3, warmup=2e3, cores=3, chains=3)


model_sum<-summary(model)
year_sd<-model_sum$p.coeff[!names(model_sum$p.coeff) %in% c("(Intercept)", "Year_fac1970", "Year_fac1972", "Year_fac1974")]

gamsim<-function(model, newdata, years=1970:2020, year_slope=0.02, year_sd=year_sd, nsim=10){
  mu<-predict(model, newdata=newdata, type="response", discrete=T, n.threads=4)
  scale <- model[["sig2"]]
  
  sim_fun<-function(){
    year_add<-tibble(Year=years, add=rnorm(n=Year, mean=(Year-min(Year))*year_slope, sd=year_sd))
    
    newdata2<-newdata%>%
      left_join(year_add, by="Year")%>%
      mutate(pred=mu+add)
    
    sims <- rnorm(n=newdata2$pred, mean=newdata2$pred, sd=scale)
    return(list(sims, year_add))
  }
  
  out<-replicate(nsim, sim_fun, simplify=F)
  return(out)
  
}