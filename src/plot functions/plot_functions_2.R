#########################   GRID MAP ##################################
cols = colora(2,15,show=F)
color_low =  cols[1]
color_high = cols[2]
gridMap <- function(initial_date,final_date,every,file_name,chosen_variable_name){
	
	filter_date_list = filter_date(data_agc_lomb,initial_date,final_date,every)
	data_from_to = filter_date_list[[1]]
	len_time = filter_date_list[[2]]
	chosen_variable = as.numeric(data_from_to[,chosen_variable_name])
	

	gradient_map <- ggplot() +
		geom_hline( yintercept = data_from_to$Latitude, color = "black", size = 0.1) +  # Linee orizzontali
		geom_vline(xintercept = data_from_to$Longitude, color = "black", size = 0.1) + # Linee verticali
		geom_tile(data = data_from_to, aes(x = Longitude, y = Latitude, fill = chosen_variable), 
				  colour = "grey50",width = 1, height = 1,alpha = 1) +
		scale_fill_gradient(low = color_low, high = color_high,na.value = "gray") +
		labs(title = chosen_variable_name)+
		guides(color = guide_legend(title = chosen_variable_name))
	
	return(animator(file_name,gradient_map,data_from_to,len_time))
	

}

#########################   STATION MAP ###############################
color_comuni_lombardia = "red"
color_empty = "white"
color_station = "yellow"
color_fill = "forestgreen"
color_station = "darkred"

stationPlot <- function(){
	# crea mappa lombardia
	mappa_migliorata <- ggplot() +
		# background_image(img)+
		# l'ordine è importante!
		geom_sf(data = altre_regioni, fill = color_empty ,color = color_fill, size = 1,alpha=0.4, show.legend = FALSE) +
		scale_fill_manual(values = c(color_station,color_empty),na.value = color_empty) +  # Define colors for inside/outside stations
		
		geom_sf(data = lombardia,aes(fill = station_inside), color = color_comuni_lombardia, size = 0.1,alpha=0.7, show.legend = FALSE) +
		scale_fill_manual(values = c(color_station, color_empty),na.value = color_fill) +  # Define colors for inside/outside stations
		
		coord_sf(xlim = range(sites$longitude) + pad, ylim = range(sites$latitude) + pad, expand = FALSE)+
		geom_point(data = sites, aes(x = longitude, y = latitude), size = 1, shape = 23, fill = color_station) +
		
		theme_bw()+
		theme(panel.grid = element_blank())+
		labs(title = "stations positions in Lombardy")
	
	
	print(mappa_migliorata)
	return(mappa_migliorata)
}
#########################   CIRCLES MAP  #################################


circlesPlot <- function(initial_date,final_date,every,file_name,chosen_var_name){
	
	filter_date_list = filter_date(df_agri,initial_date,final_date,every)
	data_from_to = filter_date_list[[1]]
	len_time = filter_date_list[[2]]
	chosen_var = as.numeric(data_from_to[,chosen_var_name])
	
	
	
	mappa_expanding <- stationPlot()+
		geom_point(data =data_from_to, aes(x=Longitude,y=Latitude),
				   size=lerp_pm10_radius(chosen_var),
				   color = lerp_pm10_color(chosen_var), alpha=0.6)+
		labs(title=paste0("measurments of ",chosen_var_name))+
		guides(color = guide_legend(title = chosen_var_name))

	return(animator(file_name,mappa_expanding,data_from_to,len_time))

}


#########################   WIND MAP  #################################
color_background_map = "forestgreen"

windPlot <- function(initial_date,final_date,every,file_name){
	
	filter_date_list = filter_date(data_agc_lomb,initial_date,final_date,every)
	data_from_to = filter_date_list[[1]]
	len_time = filter_date_list[[2]]

	wind_arrows <- data.frame(
		longitude = data_from_to$Longitude,
		latitude = data_from_to$Latitude,
		direction = cardinal_to_degree(data_from_to$WE_mode_wind_direction_10m),
		intensity = as.numeric(data_from_to$WE_wind_speed_10m_mean)
	)
	
	# Calcola le coordinate di fine delle frecce in base alla direzione 
	# l'intensità verra colorata invece di cambiare in lunghezza
	wind_arrows$end_longitude =( wind_arrows$longitude + sin(wind_arrows$direction)/10)
	wind_arrows$end_latitude = ( wind_arrows$latitude + cos(wind_arrows$direction)/10)
	# put a 0 in the NA
	wind_arrows[is.na(wind_arrows)] <- 0
	
	time_vect = c()
	for(i in 1:1053) {
		time_vect <- c(time_vect,seq(1,len_time))
	}
	wind_arrows$Time = time_vect
	
	
	mappa_wind <- ggplot(data = wind_arrows) +
		
		geom_sf(data = shp_map, fill = color_background_map , color = "black", size = 0.5,alpha=0.6)+
		coord_sf(xlim = range(na.omit(wind_arrows$longitude))+pad,
				 ylim = range(na.omit(wind_arrows$latitude))+pad, expand = FALSE)+
		
		geom_segment(data = wind_arrows,
					 aes(x = longitude, y = latitude, xend = end_longitude, yend = end_latitude,
					 	color = intensity),
					 arrow = arrow(type = "closed", length = unit(0.08, "inches"), ends = "last"),
					 lineend = "round", size = 0.3,alpha=0.9 )+
		theme_bw()+
		theme(panel.grid = element_blank())+
		labs(title="Wind map")
	
	return(animator(file_name,mappa_wind,wind_arrows,len_time))

}

#########################   XY PLOT  #####################################
cols = colora(6,970,show=F)
xyPlot <- function(initial_date,final_date,every,file_name,var1_name,var2_name,size_name,colors_factor_name){

	filter_date_list = filter_date(df_agri,initial_date,final_date,every)
	data_from_to = filter_date_list[[1]]
	len_time = filter_date_list[[2]]
	
	
	var1 = as.numeric(data_from_to[,var1_name])
	var2 = as.numeric(data_from_to[,var2_name])
	if(class(size_name)!="numeric"){
		size = as.numeric(data_from_to[,size_name])
	}else{
		size = size_name
	}
	colors_factor = as.factor(data_from_to[,colors_factor_name])
	
	if(class(size_name)!="numeric"){
		p <- ggplot(
			data_from_to, 
			aes(x = var1, y=var2, size = size, colour = colors_factor)) +
			geom_point(alpha = 1) +
			scale_color_manual(values = cols) +
			#scale_color_viridis_d() +
			scale_size(range = c(2, 12)) +
			labs(x = var1_name, y =var2_name)+
			guides(size = guide_legend(title = size_name), color = "none")+
			theme_bw()
	}else{
		p <- ggplot(
			data_from_to, 
			aes(x = var1, y=var2, size = size, colour = colors_factor)) +
			geom_point(alpha = 1,show.legend = FALSE) +
			scale_color_viridis_d() +
			scale_size(range = c(2, 12)) +
			labs(x = var1_name, y =var2_name)+
			theme_bw()
	}
	
	return(animator(file_name,p,data_from_to,len_time))
	
}
	



#########################   TREND PLOTS-stations  ########################
cols = colora(6,970,show=F)

trendStationYear <- function(chosen_station,initial_date,final_date,file_name,chosen_variable_name){
	df_st = df_stat[[stations[chosen_station]]]
	stations = unique(df_agri$IDStations)
	
	
	filter_date_list = filter_date(df_st,initial_date,final_date,"day")
	data_from_to = filter_date_list[[1]]
	len_time = filter_date_list[[2]]
	
	
	data_from_to$month_day = as.Date(paste0("1990-",substr(data_from_to$Time,6,10)))
	data_from_to$t = as.numeric(format(data_from_to$Time,"%j"))
	chosen_variable = as.numeric(data_from_to[,chosen_variable_name])
	
	# Crea il grafico ggplot
	time_trend <- ggplot(data_from_to,aes(x = month_day, 
										  y = chosen_variable,
										  group=year(Time), 
										  color = as.factor(year(Time)))) +
		
		geom_line() +
		labs(x = "Year", y = chosen_variable_name, title = paste0("Station ", stations[chosen_station], ", all years")) +
		ylim(range(na.omit(chosen_variable))) +
		scale_x_date(date_labels = "%b",date_breaks = "1 month")+
		scale_color_manual(values = cols) +
		theme(legend.position = "top")+
		guides(color = guide_legend(title = "Years"))+
		theme_bw()+
		theme(panel.grid = element_blank()) 
	
	print(time_trend)
	
	if(file_name!="None"){
		trend_animate <- time_trend + transition_reveal(data_from_to$t)+ ggtitle(data_from_to$t) + geom_point() 
		
		output_file <- paste0("./gifs/",file_name,".mp4")
		anim_save(output_file,trend_animate,
				  height = 1080, width = 1920,
				  duration = len_time,
				  fps = 10,
				  renderer = av_renderer(),
				  res = 200, type = "cairo")
	}
	
	return(time_trend)
	
}


#########################  TREND PLOTS-years #############################
cols = colora(6,970,show=F)

trendYearStation <- function(initial_date,final_date,chosen_stations,file_name,chosen_variable_name){
	
	data_year<-df_agri %>% subset(IDStations %in% stations[chosen_stations])
	
	filter_date_list = filter_date(data_year,initial_date,final_date,"day")
	data_from_to = filter_date_list[[1]]
	len_time = filter_date_list[[2]]
	

	chosen_variable = as.numeric(data_from_to[,chosen_variable_name])
	
	
	# Crea il grafico ggplot
	station_trend <- ggplot(data_from_to,aes(x = Time, 
												 y = chosen_variable,
												 group=IDStations, 
												 color = as.factor(IDStations))) +
		
		geom_line() +
		labs(x = "Stations", y = chosen_variable_name, title = paste0("Year: ",year(initial_date), ", all stations")) +
		ylim(range(na.omit(chosen_variable))) +
		scale_color_manual(values = cols) +
		theme(legend.position = "top")+
		scale_x_date(date_labels = "%b",date_breaks = "1 month")+
		theme_bw()+
		theme(panel.grid = element_blank()) +
		guides(color = guide_legend(title = "Stations"))
	
	print(station_trend)
	
	if(file_name!="None"){
		station_trend_animate <- station_trend + transition_reveal(data_from_to$Time)+ ggtitle(data_from_to$Time) + geom_point() 
		output_file <- paste0("./gifs/",file_name,".mp4")
		anim_save(output_file,station_trend_animate,
				  		  height = 1080, width = 1920, 
						  duration = (len_time%/%365),
				  		  fps = 10, 
				  		  renderer = av_renderer(),
				  		  res = 200, type = "cairo")
	
	}
	
	return(station_trend)
	
	
}





