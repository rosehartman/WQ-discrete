model_validation<-function(model, response, tag_annotation="A"){
  require(mgcv)
  require(patchwork)
  require(ggplot2)
  require(tibble)
  require(itsadug)
  
  sum<-summary(model)
  resids <- resid_gam(model, incl_na=TRUE)
  file<-tempfile()
  png(file)
  quantiles<-qq.gam(model)
  dev.off()
  file.remove(file)
  qq_data<-as_tibble(qqplot(quantiles, resids, plot.it=FALSE))
  
  p_qq<-ggplot(qq_data, aes(x=x, y=y))+
    geom_abline(intercept=0, slope=1, color="firebrick3", size=1)+
    geom_point(size=0.5)+
    xlab("Theoretical quantiles")+
    ylab("Residuals")+
    theme_bw()
  
  p_hist<-ggplot(data.frame(Residuals=resids), aes(x=Residuals))+
    geom_histogram()+
    scale_y_continuous(expand=expansion(mult=c(0,0.1)))+
    theme_bw()
  
  linpred<-napredict(model$na.action, model$linear.predictors)
  
  p_pred_resid<-ggplot(tibble(Predictor=linpred, Residuals=resids), aes(x=Predictor, y=Residuals))+
    geom_point(size=0.5)+
    xlab("Linear predictor")+
    theme_bw()
  
  fitted<-predict(model, type = "response")
  
  p_fitted_resid<-ggplot(data=tibble(Fitted=fitted, Response=response), aes(x=Fitted, y=Response))+
    geom_point(size=0.5)+
    annotate("label", x=10, y=30, label=paste0("R^2 == ", round(sum$r.sq, 4)), parse=T)+
    xlab("Fitted values")+
    theme_bw()
  
  p_check<-(p_qq|p_pred_resid)/(p_hist|p_fitted_resid)+plot_annotation(tag_levels=tag_annotation)
  return(p_check)
}

# Function to rasterize all dates. Creates a 3D raster Latitude x Longitude x Date 
Rasterize_all <- function(data, var, region=Delta, out_crs=4326, n=100){
  require(stars)
  require(purrr)
  require(sf)
  var<-rlang::enquo(var)
  rlang::as_name(var)
  preds<-map(unique(data$Date), function(x) st_rasterize(data%>%
                                                           filter(Date==x)%>%
                                                           dplyr::select(!!var), 
                                                         template=st_as_stars(st_bbox(region), dx=diff(st_bbox(region)[c(1, 3)])/n, dy=diff(st_bbox(region)[c(2, 4)])/n, values = NA_real_))%>%
               st_warp(crs=out_crs))
  
  # Then bind all dates together into 1 raster
  out <- exec(c, !!!preds, along=list(Date=unique(data$Date)))
  return(out)
}

ST_variogram<-function(model, data, cores){
  require(spacetime)
  require(sp)
  require(gstat)
  require(ggplot2)
  require(dplyr)
  require(itsadug)
  require(sf)
  require(patchwork)
  
  norm_resids<-resid_gam(model, incl_na=TRUE)
  
  Data_vario<-model%>%
    mutate(Resid=norm_resids$norm_res)
  
  Data_coords<-Data_vario%>%
    st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)%>%
    st_transform(crs=26910)%>%
    st_coordinates()%>%
    as_tibble()%>%
    mutate(across(c(X,Y), ~(.x-mean(.x))/1000))
  
  Data_vario<-bind_cols(Data_vario%>%
                          select(Date, Resid), Data_coords)
  sp<-SpatialPoints(coords=data.frame(X=Data_vario$X, Y=Data_vario$Y))
  sp2<-STIDF(sp, time=Data_vario$Date, 
             data=data.frame(Residuals=Data_vario$Resid))
  vario<-variogramST(Residuals~1, data=sp2, tunit="weeks", cores=cores, tlags=(30/7)*0:10)
  
  vario_plot<-vario%>%
    mutate(monthlag=as.integer(as.factor(timelag))-0.5)
  
  p_time<-ggplot(vario_plot, aes(x=monthlag, y=gamma, color=spacelag, group=spacelag))+
    geom_line()+
    geom_point()+
    scale_color_viridis_c(name="Distance (km)")+
    scale_x_continuous(breaks=c(2,4,6,8,10))+
    xlab("Time difference (months)")+
    theme_bw()+
    theme(legend.justification = "left")
  
  p_space<-ggplot(vario_plot, aes(x=spacelag, y=gamma, color=monthlag, group=monthlag))+
    geom_line()+
    geom_point()+
    scale_color_viridis_c(breaks=c(2,4,6,8,10), name="Time difference\n(months)")+
    xlab("Distance (km)")+
    theme_bw()+
    theme(legend.justification = "left")
  
  p_variogram<-p_time/p_space+plot_annotation(tag_levels="A")
  
  return(p_variogram)
}