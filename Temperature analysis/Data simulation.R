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
Time_means<-
Sim_data<-expand_grid(Month=1:12, Year=1970:2020, Location=1:50, Source=1:5)%>%
  mutate(Time_num_s=rnorm(n()),
         Year_num=as.integer(as.factor(Year)),
         Temp=rnorm(n(), Year_means[Year_num]+season[Month]*10+))
