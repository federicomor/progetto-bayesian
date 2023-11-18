library(maps)
library(ggplot2)
library(sf)
library(lubridate)
#devtools::install_github("dgrtwo/gganimate")
library(gganimate)
library(av)
library(ggpubr)
library(jpeg)

##### OLD
# load("../data/data_agc_lomb.Rdata")
##### NEW
load("../data/data_agc_lomb_part1.RData")
load("../data/data_agc_lomb_part2.RData")
AGC_Dataset = rbind(parte1, parte2)

head(AGC_Dataset)
data_agc_lomb=AGC_Dataset
data_agc_lomb$Time = as.Date(data_agc_lomb$Time)
df_agri$Time = as.Date(df_agri$Time)

sites <- data.frame(
	longitude = unique(df_agri$Longitude), 
	latitude = unique(df_agri$Latitude))

stations = unique(df_agri$IDStations)
cols = colora(6,970,show=F)
pad <- 0.2 * c(-1, 1)
DATE_FORMAT = "%Y-%m-%d"
########################

# carica file shp
details <- 3
details_altre_regioni <- 2

confini  <- c("Emilia-Romagna","Piemonte","Lombardia","Trentino-Alto Adige","Veneto")
regioni_italiane <- st_read(paste0("italia/gadm40_ITA_",details,".shp"))
regioni_italiane_2 <- st_read(paste0("italia/gadm40_ITA_",details_altre_regioni,".shp"))

lombardia <- regioni_italiane[regioni_italiane$NAME_1 == "Lombardia",]
altre_regioni <- regioni_italiane_2[regioni_italiane_2$NAME_1 %in% confini,]

shp_map <- st_read(paste0("italia/gadm40_ITA_",1,".shp"))
#######################################################################
#########################   UTILITIES #################################
#######################################################################
# Check if a station is inside each subregion
lombardia$station_inside <- sapply(lombardia$geometry, function(x) any(st_intersects(st_as_sf(sites, coords = c("longitude", "latitude")), x)))

altre_regioni$station_inside <- sapply(altre_regioni$geometry, function(x) any(st_intersects(st_as_sf(sites, coords = c("longitude", "latitude")), x)))


cardinal_to_degree <- function(cardinal) {
	cardinal_degrees <- c(N = 0, NNE = 22.5, NE = 45, ENE = 67.5, E = 90,
						  ESE = 112.5, SE = 135, SSE = 157.5,
						  S = 180, SSW = 202.5, SW = 225, WSW = 247.5,
						  W = 270, WNW = 292.5, NW = 315, NNW = 337.5)
	return(cardinal_degrees[cardinal])
}

# LERP FUNCTIONS
lerp_pm10_radius <- function(val){
	radius_0 = 1
	radius_f = 10
	min_max = range(val[!is.na(val)])
	min_val = min_max[1]
	max_val = min_max[2]
	
	# linear, you can use any function!
	return_val = radius_0 + (radius_f - radius_0) * ((val - min_val) / (max_val - min_val))
	return(return_val) 
}

lerp_pm10_color <- function(val) {
	color_0 <- col2rgb("#008F39")
	color_f <- col2rgb("#F80000")
	
	min_max <- range(val[!is.na(val)])
	min_val <- min_max[1]
	max_val <- min_max[2]
	
	return_vals <- vector("list", length = length(val))
	
	for (i in seq_along(val)) {
		if (!is.nan(val[i])) {
			return_val <- color_0 + (color_f - color_0) * ((val[i] - min_val) / (max_val - min_val))
			return_vals[[i]] <- return_val
		} else {
			return_vals[[i]] <- rep(255, 3)  # Set RGB values to 255 (white)
		}
	}
	
	# Convert RGB components to color strings
	colors <- lapply(return_vals, function(rgb_vec) {
		rgb_string <- rgb(rgb_vec[1], rgb_vec[2], rgb_vec[3], maxColorValue = 255)
		return(rgb_string)
	})
	
	return(colors) 
}


filter_date <- function(dataframe,initial_date,final_date,every){
	
	date_time = as.Date(seq.Date(from=as.Date(initial_date,DATE_FORMAT),
								 to=as.Date(final_date,DATE_FORMAT),by=every))
	filtered_dataset = dataframe[which(as.Date(dataframe$Time) %in% date_time),]

	return(list(filtered_dataset,length(date_time)))
	
}
#######################################################################
#########################   GRID MAP ##################################
#######################################################################
gridMap <- function(initial_date,final_date,every,file_name,color_low,color_high,chosen_variable_name){
	
	filter_date_list = filter_date(data_agc_lomb,initial_date,final_date,every)
	data_from_to = filter_date_list[[1]]
	len_time = filter_date_list[[2]]

	
	chosen_variable = data_from_to[,chosen_variable_name]
	
	# aggiungi grid for long and lat
	grid <- ggplot() +
		geom_hline( yintercept = data_from_to$Latitude, color = "black", size = 0.1) +  # Linee orizzontali
		geom_vline(xintercept = data_from_to$Longitude, color = "black", size = 0.1)  # Linee verticali
	
	
	# Aggiungi griglia colorata per altitudine
	
	gradient_map <- 
		grid+
		geom_tile(data = data_from_to, aes(x = Longitude, y = Latitude, fill = chosen_variable), 
				  colour = "grey50",width = 1, height = 1,alpha = 1) +
		scale_fill_gradient(low = color_low, high = color_high,na.value = "gray") 
	
	
	if (file_name == "None"){
		gradient_map <- gradient_map +facet_wrap(~Time)
		print(gradient_map)
		return(gradient_map)
	} else {
		gradient_map <- gradient_map +  
			ggtitle(data_from_to$Time) +
			transition_time(data_from_to$Time)+
			labs(title = paste0(every,": {frame_time}"))
		
		output_file <- paste0("./gifs/",file_name,".mp4")
		anim_save(output_file,gradient_map,duration = len_time,fps = 10, renderer = av_renderer())
	}

}
#######################################################################
#########################   STATION MAP ###############################
#######################################################################
stationPlot <- function(){
	# crea mappa lombardia
	mappa_migliorata <- ggplot() +
		# background_image(img)+
		# l'ordine è importante!
		geom_sf(data = altre_regioni, fill = "white" ,color = "lightgreen", size = 1) +
		#scale_fill_manual(values = c("gold", "white"),na.value = "white") +  # Define colors for inside/outside stations
		geom_sf(data = lombardia,aes(fill = station_inside), color = "red", size = 1) +
		scale_fill_manual(values = c("yellow", "white"),na.value = "lightblue") +  # Define colors for inside/outside stations
		coord_sf(xlim = range(sites$longitude) + pad, ylim = range(sites$latitude) + pad, expand = FALSE)+
		theme(legend.position = "none")+
		theme(panel.grid = element_blank())+
		theme_bw()
	
	# aggiungi stazioni
	mappa_migliorata <- mappa_migliorata + 
		geom_point(data = sites, aes(x = longitude, y = latitude), size = 1, shape = 23, fill = "darkred") 
	
	print(mappa_migliorata)
	return(mappa_migliorata)
}
#######################################################################
#########################   WIND MAP  #################################
#######################################################################

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
	wind_arrows$t = time_vect
	
	
	mappa_wind <- ggplot(data = wind_arrows) +
		geom_sf(data = shp_map, fill = "lightgreen", color = "black", size = 0.5)+
		coord_sf(xlim = range(na.omit(wind_arrows$longitude))+pad,
				 ylim = range(na.omit(wind_arrows$latitude))+pad, expand = FALSE)+
		
		geom_segment(data = wind_arrows,
					 aes(x = longitude, y = latitude, xend = end_longitude, yend = end_latitude,
					 	color = intensity),
					 arrow = arrow(type = "closed", length = unit(0.08, "inches"), ends = "last"),
					 lineend = "round", size = 0.3,alpha=0.9 )+
		theme_bw()
	
	if (file_name == "None"){
		mappa_wind <- mappa_wind +facet_wrap(~t)
		print(mappa_wind)
	} else {
		mappa_animata <- mappa_wind +  
			ggtitle(wind_arrows$t) +
			transition_time(wind_arrows$t)+
			labs(title = paste0(every,": {frame_time}"))
		
		output_file <- paste0("./gifs/",file_name,".mp4")
		anim_save(output_file,mappa_animata,duration = len_time,fps = 10, renderer = av_renderer(), width = 1280, height = 720)
	}
}
##########################################################################
#########################   XY PLOT  #####################################
##########################################################################
xyPlot <- function(initial_date,final_date,every,file_name,var1_name,var2_name,size_name,colors_factor_name){

	filter_date_list = filter_date(df_agri,initial_date,final_date,every)
	data_from_to = filter_date_list[[1]]
	len_time = filter_date_list[[2]]
	
	#data_from_to$Time = as.Date(data_from_to$Time )
	
	var1 = data_from_to[,var1_name]
	var2 = data_from_to[,var2_name]
	if(class(size_name)!="numeric"){
		size = data_from_to[,size_name]	
	}else{
		size = size_name
	}
	colors_factor = data_from_to[,colors_factor_name]
	
	p <- ggplot(
		data_from_to, 
		aes(x = var1, y=var2, size = size, colour = colors_factor)) +
		
		geom_point(alpha = 0.7,show.legend = FALSE) +
		scale_color_viridis_d() +
		scale_size(range = c(2, 12)) +
		labs(x = var1_name, y =var2_name)+
		theme_bw()
	if(class(size_name)!="numeric"){
		p+guides(size = guide_legend(title = size_name), color = "none")
	}
	
	
	if (file_name == "None"){
		p <- p +facet_wrap(~Time)
		print(p)
		return(p)
	} else {
		mappa_animata <- p+ transition_time(data_from_to$Time) +
			ggtitle(data_from_to$Time) +
			labs(paste0(every,": {frame_time}"))+
			shadow_wake(wake_length = 0.1, alpha = FALSE)
		
		output_file <- paste0("./gifs/",file_name,".mp4")
		anim_save(output_file,mappa_animata,duration = len_time,fps = 10, renderer = av_renderer())
	}
}

##########################################################################
#########################   TREND PLOTS-stations  ########################
##########################################################################

trendStationYear <- function(chosen_station,initial_date,final_date,file_name,chosen_variable_name){
	st = stations[chosen_station]
	df_st = df_stat[[st]]
	stations = unique(df_agri$IDStations)
	cols = colora(6,970,show=F)
	
	filter_date_list = filter_date(df_st,initial_date,final_date,every)
	data_from_to = filter_date_list[[1]]
	len_time = filter_date_list[[2]]


	data_from_to$month_day = as.Date(paste0("1990-",substr(data_from_to$Time,6,10)))
	data_from_to$t = format(data_from_to$Time,"%j")
	chosen_variable = data_from_to[,chosen_variable_name]
	
	# Crea il grafico ggplot
	time_trend <- ggplot(data_from_to,aes(x = month_day, 
										  y = chosen_variable,
										  group=year(Time), 
										  color = as.factor(year(Time)))) +
		
		geom_line() +
		
		labs(x = "Year", y = chosen_variable_name, title = paste0("Station ", st, ", all years")) +
		ylim(range(na.omit(chosen_variable))) +
		scale_x_date(date_labels = "%b",date_breaks = "1 month")+
		
		scale_color_manual(values = cols) +
		theme(legend.position = "top")+
		guides(color = guide_legend(title = "Years"))+
		theme_bw()+
		theme(panel.grid = element_blank()) 
	
	print(time_trend)
	
	if(file_name!="None"){
	#	trend_animate <- time_trend + transition_reveal(as.numeric(t))+ ggtitle(as.numeric(t)) + geom_point() 
		trend_animate <- time_trend + transition_reveal(data_from_to$Time)+ ggtitle(data_from_to$Time) + geom_point()
		output_file <- paste0("./gifs/",file_name,".mp4")
		anim_save(output_file,trend_animate,duration = len_time,fps = 10, renderer = av_renderer())
	}
	
	return(time_trend)
	
}

##########################################################################
#########################  TREND PLOTS-years #############################
##########################################################################
trendYearStation <- function(initial_date,final_date,chosen_stations,file_name,chosen_variable_name){
	
	st = stations[chosen_stations]
	data_year<-df_agri %>% subset(IDStations %in% st)
	
	filter_date_list = filter_date(data_year,initial_date,final_date,every)
	data_from_to = filter_date_list[[1]]
	len_time = filter_date_list[[2]]
	

	chosen_variable = data_from_to[,chosen_variable_name]
	
	
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
		anim_save(output_file,station_trend_animate,duration = len_time,fps = 10, renderer = av_renderer())
	}
	
	return(station_trend)
	
	
}

##########################################################################
#########################   CIRCLES MAP  #################################
##########################################################################







