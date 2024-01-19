####################################
## Include Me in your scripts with
## "source("include.R")
####################################

# morally always needed
library(ggplot2)
library(grDevices)
library(png)

Mode <- function(x,verbose=0) {
	ux <- unique(x)
	if(verbose){
		print(table(x))
	}
	ux[which.max(tabulate(match(x, ux)))]
}

########################
# useful aliases
########################
size = dim
end = len = lenght = length
assert = stopifnot
ln = log

extrema = function(vec){
	return(c(min(vec),max(vec)))
}

count_na = function(vec){
	na_vals = sum(as.numeric(is.na(vec)))
	return(na_vals)
}
na_count = quanti_na = count_na

########################
# loading procedure, with feedback
########################
library(crayon)
library(hash)

h="❤"
# load("../data/df_agri.Rdata")
# df_agri_full = df_agri
# rm(df_agri)
cat(crayon::strikethrough(
	crayon::cyan(h,"Loaded agrimonia dataset. Available as"),
	crayon::red("df_agri_full.\n")))
# cat(crayon::italic("Consider to remove it with rm(df_agri_full).\n"))
cat(crayon::italic("Removed it as useless now.\n\n"))


load("../data/Cleaned_data.Rdata")
cat(crayon::cyan(h,"Loaded cleaned dataset. Available as"),crayon::red("df_agri.\n"))

indeces_2018 = df_agri$Time>=as.Date("2018-01-01") & df_agri$Time<=as.Date("2018-12-31")
assert(sum(indeces_2018==TRUE) == 365*length(unique(df_agri$IDStations)))
df_2018 = df_agri[indeces_2018,]
cat(crayon::cyan(h,"Loaded 2018 dataset. Available as"),crayon::red("df_2018.\n"))
rm(indeces_2018)

df_agri$Time <- as.Date(df_agri$Time, "%Y-%m-%d")
df_2018$Time <- as.Date(df_2018$Time, "%Y-%m-%d")

cat(crayon::italic("Converted column Time (of all dfs) to Date variable type, format: year-month-day.\n"))
DATE_FORMAT = "%Y-%m-%d"
cat(crayon::italic("Created DATE_FORMAT = \"%Y-%m-%d\" variable for date comparisons.\n"))
cat(crayon::italic("Usage example: df_agri$Time >= as.Date(\"2017-01-01\",DATE_FORMAT).\n"))
cat(crayon::italic("Actually also this works: df_agri$Time >= as.Date(\"2017-01-01\").\n\n"))

load("../data/df_weekly.Rdata") # here it is not log-trasformed, we asser it
if(!(min(df_weekly$AQ_pm10)>2 && max(df_weekly$AQ_pm10)>78)){
	warning("Dataset seems already log-transformed when it shouldn't be.")
}
# extrema(df_weekly$AQ_pm10)
# [1] 2.428571 79.857143

# df_weekly_no_log_transf = df_weekly # so here we make a copy of it, the original one
df_weekly$AQ_pm10 = log(df_weekly$AQ_pm10) # here we overwrite df_weekly to be log-transformed
# extrema(df_weekly$AQ_pm10)
# [1] 0.8873032 4.3802393

cat(crayon::cyan(h,"Loaded weekly divided dataset(s). Available as"),
	crayon::red("df_weekly.\n"))

df_weekly$IDStations[which(df_weekly$IDStations=="STA-CH0011A")] = "STA.CH0011A"
df_weekly$IDStations[which(df_weekly$IDStations=="STA-CH0033A")] = "STA.CH0033A"
df_weekly$IDStations[which(df_weekly$IDStations=="STA-CH0043A")] = "STA.CH0043A"
cat(crayon::italic("Uniformed names of stations (some were STA-ecc and some STA.ecc; now are all STA.ecc).\nOnly of df_weekly and the following df_weekly_scaled_centered (ie df_wsc).\n"))
cat(crayon::italic("This name change was also needed for the graph cluster plot function.\n\n"))




df_weekly_scaled_centered = df_weekly
# colnames(df_weekly_scaled_centered)
numerical_covariates_to_scale = c(3,4,6,8:10,12:20,22:37)
df_weekly_scaled_centered[,numerical_covariates_to_scale] = scale(df_weekly_scaled_centered[,numerical_covariates_to_scale],center=TRUE,scale=TRUE)
df_weekly_scaled_centered$AQ_pm10 = scale(df_weekly_scaled_centered$AQ_pm10,center=TRUE,scale=FALSE)
cat(crayon::cyan(h,"Created scaled df_weekly dataset. Available as"),
	crayon::red("df_weekly_scaled_centered"),
	crayon::cyan("(or"),
	crayon::red("df_wsc).\n")) # add centering on the PM10
cat(crayon::italic(
	"Scaled variables were c(3,4,6,8:10,12:20,22:37)
Untouched variables were
   (col 1) X, 
   (col 2) IDStations,
   (col 5) Time,
  (col 11) WE_mode_wind_direction_10m,
  (col 21) WE_mode_wind_direction_100,
  (col 38) day,
  (col 39) week
------------------------------------------------------
  (col 7) AQ_pm10 has been centered, not scaled
(col 3&4) Latitude and Longitude have also been scaled
          (fits were better in this way)\n\n"))

# par(mfrow=c(2,2))
# par(mar=c(9,2,1,1))
# boxplot(df_weekly_scaled_centered[,c(3,4,6,8:10,12:20,22:37)],title="scaled variables",las=2,cex=0.4)
# par(mar=c(3,3,2,2))
# plot(df_weekly_scaled_centered$Latitude,df_weekly_scaled_centered$Longitude,main="std coordinates (look at axes)")
# boxplot(df_weekly_scaled_centered$AQ_pm10,main="AQ_pm10 values")
# boxplot(df_weekly_scaled_centered[,c(7,38,39)],title="remained variables",las=2,cex=0.4)

df_wsc = df_weekly_scaled_centered # short alias


########################
# df_stat split
########################
create_df_stat = function(df_data){
	df_stat_ = hash()
	if("IDStations" %in% colnames(df_data)){
		stations_ = unique(df_data$IDStations)
		for (st in stations_){
			df_stat_[[st]] = df_data[which(df_data$IDStations == st),]
		}
		rm(st)
		rm(stations_)
		return(df_stat_)
	}
	else{
		cat(crayon::red("Error: missing IDStations field in df.\n"))
		return(NA)
	}
}

cat(crayon::cyan(h,"Created stations split function Available as"),crayon::red("create_df_stat(df).\n"))
cat(crayon::italic("Use it as my_df_stat = create_df_stat(df_2018).\nThen for example my_df_stat[[\"1264\"]] retrieves the dataset for station 1264.\n\n"))

########################
# function to get colors for plotting
########################
library(RColorBrewer)
fun_colori = function(len=2, seed=33, show=1, seed_div = "Set3"){
	hcols_ = hcl.pals()
	if(seed=="rand"){
		seed = round(runif(1,0,115))
		col.ramp_ = hcl.colors(len,palette=hcols_[seed%%115+1])
	}
	if(seed=="div"){ # create a divergent palette
		col.ramp_ = brewer.pal(len, seed_div)
		# possible seed_div choices (and max len supported)
		# Set3	    12
		# Paired    12
		# Pastel1   9
		# Set1	    9
		# Accent    8
		# Dark2     8
		# Pastel2   8
		# Set2	    8
	}
	else{
		col.ramp_ = hcl.colors(len,palette=hcols_[seed%%115+1])
	}
	if(show==1){
		dati_ <- matrix(1:100, ncol = 1)
		par(mar=rep(2,4))
		image(dati_, col = col.ramp_, axes = FALSE)
		title(main=paste("palette",seed,"of",len,"colors"))
		# title(main=paste("palette",seed))
	}
	return(col.ramp_)
	
}
colori_fun = colorami = colora = fun_colori # aliases
# usage: cols = colora(7) for getting a palette with 7 colors
# usage: cols = colora(7,456) for getting a palette of seed 456 with 7 colors

cat(crayon::cyan(h,"Created function to get color palettes. Available as"),crayon::red("colora(len, seed, show).\n"))
cat(crayon::italic("Try for example colora(10,56,1).\n"))


########################
# function to show the hash content
########################
# what_is <- new.env(hash = TRUE, parent = emptyenv(), size = NA)
what_is = hash()
# usage: explain("AQ_pm10") show AQ_pm10 explanation
# usage: explain("pm10") yes works also with partial names
# usage: explain("pm") here explains both pm 10 and 2.5

spiegami = function(pattern) {
	matching_vars <- names(what_is)[grep(pattern, names(what_is))]
	# chat gpt help here :/	
	if (length(matching_vars) == 0) {
		cat("no matching variable name found\n")
	} else {
		for (var in matching_vars) {
			cat(crayon::red(var), ":", what_is[[var]], "\n")
		}
	}
}
spiega = tellme = tell_me = tell = explain = spiegami # aliases
# spiegami = function(colonna){
# 	cat(what_is[[colonna]])
# }


what_is[["AQ"]] = "Air quality - Air pollutants concentrations"
what_is[["WE"]] = "Wheather - Estimates of atmospheric variables"
what_is[["EM"]] = "Emission - Emission inventories"
what_is[["LI"]] = "Livestock - Livestock inventories"
what_is[["LA"]] = "Land cover -  Estimates of land cover variables"

what_is[["IDStations"]] = "Unique ID (uid) of station present in sensor network [For detail see Station_Registry.csv file]. \n[ - ]"
what_is[["Latitude"]] = "Latitude of station \n[ Degrees ]"
what_is[["Longitude"]] = "Longitude of station \n[ Degrees ]"
what_is[["Time"]] = "Date \n[ yyyy-mm-dd ]"
what_is[["Altitude"]] = "Elevations of the stations \n[ m ]"
what_is[["AQ_pm10"]] = "Concentration of particles with an aerodynamic diameter of less than 10 micrometers (µm) \n[ μg m^(-3) ]"
what_is[["AQ_pm25"]] = "Concentration of particles with an aerodynamic diameter of less than 2.5 micrometers (µm) \n[ μg m^(-3) ]"
what_is[["AQ_co"]] = "Concentration of carbon monoxide \n[ μg m^(-3) ]"
what_is[["AQ_nh3"]] = "Concentration of ammonia \n[ μg m^(-3) ]"
what_is[["AQ_nox"]] = "Concentration of several oxides of nitrogen most of which are produced in combustion and are considered to be atmospheric pollutants \n[ μg m^(-3) ]"
what_is[["AQ_no2"]] = "Concentration of nitrogen dioxide \n[ μg m^(-3) ]"
what_is[["AQ_so2"]] = "Concentration of sulfur dioxide \n[ μg m^(-3) ]"
what_is[["WE_temp_2m"]] = "Temperature of air at 2 m above the surface of land or sea or inland water \n[ Celsius ]"
what_is[["WE_wind_speed_10m_mean"]] = "Mean horizontal wind speed at a height of 10 m above the surface of the Earth \n[ m/s ]"
what_is[["WE_wind_speed_10m_max"]] = "Maximum horizontal wind speed at a height of 10 m above the surface of the Earth \n[ m/s ]"
what_is[["WE_mode_wind_direction_10m"]] = "Direction of horizontal wind at a height of 10 m above the surface of the Earth \n[ categorical ]"
what_is[["WE_tot_precipitation"]] = "The accumulated liquid and frozen water (comprising rain and snow) that falls to the Earth’s surface \n[ m ]"
what_is[["WE_precipication_t"]] = "The type of precipitation at the surface at the specified time. \n[ categorical ]"
what_is[["WE_surface_pressure"]] = "The pressure (force per unit area) of the atmosphere at the surface of land and sea and inland water \n[ Pa ]"
what_is[["WE_solar_radiation"]] = "Amount of solar radiation that reaches a horizontal plane at the surface of the Earth (both direct and diffuse) minus the amount reflected by the Earth’s surface \n[ J m^(-2) ]"
what_is[["WE_wind_speed_100m_mean"]] = "Mean horizontal wind speed at a height of 100 m above the surface of the Earth \n[ m/s ]"
what_is[["WE_wind_speed_100m_max"]] = "Max horizontal wind speed at a height of 100 m above the surface of the Earth \n[ m/s ]"
what_is[["WE_mode_wind_direction_100m"]] = "Direction of horizontal wind at a height of 100 m above the surface of the Earth \n[ categorical ]"
what_is[["WE_blh_layer_max"]] = "The maximum depth of air next to the Earth’s surface which is most affected by the resistance to the transfer of momentum or heat or moisture across the surface \n[ m ]"
what_is[["WE_blh_layer_min"]] = "The minimum depth of air next to the Earth’s surface which is most affected by the resistance to the transfer of momentum or heat or moisture across the surface \n[ m ]"
what_is[["WE_rh_min"]] = "Minimum of amount of water vapour present in air expressed as a percentage of the amount needed for saturation at the same temperature \n[ % ]"
what_is[["WE_rh_mean"]] = "Mean of amount of water vapour present in air expressed as a percentage of the amount needed for saturation at the same temperature \n[ % ]"
what_is[["WE_rh_max"]] = "Maximum of amount of water vapour present in air expressed as a percentage of the amount needed for saturation at the same temperature \n[ % ]"
what_is[["EM_nh3_livestock_mm"]] = "Emissions of NH3 originating from the livestock sector within the manure management \n[ mg m^(-2) ]"
what_is[["EM_nh3_agr_soils"]] = "Emissions of NH3 originating from the agriculture soils \n[ mg m^(-2) ]"
what_is[["EM_nh3_agr_waste_burn"]] = "Emissions of NH3 originating from the burning of agriculture waste \n[ mg m^(-2) ]"
what_is[["EM_nh3_sum"]] = "Total emissions of NH3 across all sectors (anthropogenic) \n[ mg m^(-2) ]"
what_is[["EM_nox_traffic"]] = "Emissions of NOx from the on road transportation \n[ mg m^(-2) ]"
what_is[["EM_nox_sum"]] = "Emissions of NOx across all sectors (anthroprogenic) \n[ mg m^(-2) ]"
what_is[["EM_so2_sum"]] = "Total emissions of SO2 across all sectors (anthroprogenic) \n[ mg m^(-2) ]"
what_is[["LI_pigs"]] = "Density of pigs \n[ Number km^(-2) ]"
what_is[["LI_bovine"]] = "Density of bovines \n[ Number km^(-2) ]"
what_is[["LI_pigs_v2"]] = "Density of pigs \n[ Number km^(-2) ]"
what_is[["LI_bovine_v2"]] = "Density of bovines \n[ Number km^(-2) ]"
what_is[["LA_hvi"]] = "One-half of the total green leaf area per unit horizontal ground surface area for high vegetation type \n[ unitless (m^2/m^2) ]"
what_is[["LA_lvi"]] = "One-half of the total green leaf area per unit horizontal ground surface area for low vegetation type. \n[ unitless (m^2/m^2) ]"
what_is[["LA_land_use"]] = "Land use across 44 sectors \n[ categorical ]"
what_is[["LA_soil_use"]] = "Lombardy soil use across 21 sectors \n[ categorical ]"

cat(crayon::cyan(h,"Created utility to explain covariates. Available as"),crayon::red("spiega(string).\n"))
cat(crayon::italic("Try for example spiega(\"wind\").\n\n"))

# rm(h)

