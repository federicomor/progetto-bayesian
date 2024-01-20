file_type=".mp4"
folder = "plot functions/gifs/"
#########################   UTILITIES #################################

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
lerp_pm10_radius <- function(val,rmin=1,rmax=10){
	radius_0 = rmin
	radius_f = rmax
	min_max = range(val[!is.na(val)])
	min_val = min_max[1]
	max_val = min_max[2]
	
	# linear, you can use any function!
	return_val = radius_0 + (radius_f - radius_0) * ((val - min_val) / (max_val - min_val))
	return(return_val) 
}

lerp_pm10_color <- function(val,colmin="yellow",colmax="#F80000") {
	color_0 <- col2rgb(colmin)
	color_f <- col2rgb(colmax)
	
	min_max <- range(val[!is.na(val)])
	min_val <- min_max[1]
	max_val <- min_max[2]
	
	return_vals <- vector("list", length = length(val))
	
	for (i in seq_along(val)) {
		if (!is.na(val[i]) & !is.nan(val[i])) {
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
	if(every == ""){
		date_time = as.Date(c(initial_date,final_date))
	}else{
		date_time = as.Date(seq.Date(from=as.Date(initial_date,DATE_FORMAT),
									 to=as.Date(final_date,DATE_FORMAT),by=every))
	}
	
	filtered_dataset = dataframe[which(as.Date(dataframe$Time) %in% date_time),]
	
	return(list(filtered_dataset,length(date_time)))
	
}


animator <- function(file_name,plotgg, df,len_time,w,h,shadow=FALSE){
	if (file_name == "None"){
		plotgg <- plotgg +facet_wrap(~Time)
		#print(plotgg)
	} else {
		mappa_animata <- plotgg +  
			ggtitle(df$Time) +
			transition_time(df$Time)+
			labs(title = paste0(every,": {frame_time}"))
		if(shadow){
			mappa_animata <- mappa_animata+shadow_wake(wake_length = 0.1, alpha = FALSE)
		}
		output_file <- paste0(folder,file_name,file_type)
		
		anim_save(output_file,mappa_animata,
				  
				  duration = len_time,
				  fps = 10, 
				  width = w, height = h, 
				  renderer = av_renderer(),
				  res = 200, type = "cairo")
	}
	return(plotgg)
}

trend_animator <- function(file_name,plotgg, time,len_time){
	
	if(file_name!="None"){
		plot_trend_animate <- plotgg + transition_reveal(time)+ ggtitle(time) + geom_point() 
		output_file <- paste0(folder,file_name,file_type)
		anim_save(output_file,plot_trend_animate,
				  height = 1080, width = 1920, 
				  duration = len_time,
				  fps = 10, 
				  renderer = av_renderer(),
				  res = 200, type = "cairo")
		
	}
	#print(plotgg)
	return(plotgg)
	
	
}
