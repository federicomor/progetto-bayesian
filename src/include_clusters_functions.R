cat(crayon::cyan("\nGenerated these functions:\n"))

################################################################################
#### mode correction
mode_correct_clusters = function(cl_old, cl_cur,verbose=0){ # and returns cl_new the updated cl_cur
	if(is.null(cl_old)){
		warning("The cl_old vector is NULL (probably it's just iteration 1 initialization).\nReturning original clustering.\n")
		return(cl_cur) # nothing to correct
	}
	
	cl_new = cl_cur
	new_labels = c()
	tied_labels = c() # solve them later
	already_taken_labels = c() # solve them later
	
	for (i in 1:max(cl_cur)){
		indici_label_i = which(cl_cur == i)
		freq = table(cl_old[indici_label_i])
		better_labels = as.numeric(names(freq[freq == max(freq)]))
		num_mode = length(better_labels)
		
		
		if(num_mode>=2){ # per risolvere i pareggi
			# warning("Tie case happened.\n")
			tied_labels = c(tied_labels,i)
			
		} else { # una sola moda, ma c'è da vedere se non è già presa
			if( !(better_labels %in% new_labels) ){
				better_label = better_labels
				new_labels = c(new_labels,better_label)
				cl_new[indici_label_i] = better_label
			} else {
				already_taken_labels = c(already_taken_labels,i)
			}
		}
		# cat(new_labels,"\n")
	}
	max_it = 100
	it = 1
	if(verbose==1){
		cat("tied_labels =",tied_labels,"\n")
		cat("already_taken_labels =",already_taken_labels,"\n")
	}
	
	while(length(tied_labels)!=0 && it<max_it){
		for(i in tied_labels){
			indici_label_i = which(cl_cur == i)
			freq = table(cl_old[indici_label_i])
			better_labels = as.numeric(names(freq[freq == max(freq)]))
			
			better_labels = setdiff(better_labels,new_labels)
			if(length(better_labels)==1){
				# risolvi il pareggio ed eilimina label i da tied_labels
				better_label = better_labels
				new_labels = c(new_labels,better_label)
				cl_new[indici_label_i] = better_label
				tied_labels = setdiff(tied_labels,i)
			}
			else{
				# prendi il primo label, tanto comunque non c'erano in new_labels
				# perché abbiam fatto prima il set diff, quindi sono label liberi
				better_label = better_labels[1]
				new_labels = c(new_labels,better_label)
				cl_new[indici_label_i] = better_label
				tied_labels = setdiff(tied_labels,i)
			}
		}
		it = it+1
	}
	for(i in already_taken_labels){
		indici_label_i = which(cl_cur == i)
		# cat(indici_label_i,"\n")
		for(k in 1:length(unique(cl_cur))){
			# cat(k,"\n")
			if( !(k %in% new_labels) ){
				new_labels = c(new_labels,k)
				cl_new[indici_label_i] = k
				break
			}
		}
	}
	if(verbose==1){
		cat("cl_old =",cl_old,"\n")
		cat("cl_cur =",cl_cur,"\n")
		cat("cl_new =",cl_new,"\n")
	}
	
	if(it==max_it){
		warning("Something went wrong. Returning original clustering.\n")
		return(cl_cur)
	}
	if (length(unique(cl_new))!=length(unique(cl_new)) ||
		length(levels(factor(cl_new)))!=length(levels(factor(cl_cur))) ){
		warning("Something went wrong. Returning original clustering.\n")
		return(cl_cur)
	}

	if(!all(cl_new==cl_cur)){
		cat(crayon::italic("Some change was made!"))
	}
	return(cl_new)
	
}
cat(crayon::red("- mode_correct_clusters(cl_old, cl_cur)\n"))

cat(crayon::blue("Test sets\n"))
### test sets

cat(crayon::blue("- Ettore test\n"))
cl_old = c(1,2,1,1,2,2,3,3)
cl_cur = c(1,1,2,2,1,1,3,3)
check_ = c(2,2,1,1,2,2,3,3)
mode_correct_clusters(cl_old,cl_cur,verbose=1)
# mode_correct_clusters(cl_old,cl_cur) == check_

cat(crayon::blue("- No change case\n"))
cl_old = c(1,2,1,1,2,2,3,3)
cl_cur = c(1,2,1,1,2,2,3,3)
mode_correct_clusters(cl_old,cl_cur,verbose=1)


cat(crayon::blue("- Tie case\n"))
# tie case
#   units  1 2 3 4 5 6 7 8 9 10
cl_old = c(1,2,1,1,2,1,3,3)
cl_cur = c(1,1,2,2,1,1,3,3)
# ora (pareggio) vedere sia il cluster 1 è entrato in 2 ma anche viceversa
check1 = c(2,2,1,1,2,2,3,3)
check2 = c(1,1,2,2,1,1,3,3)
mode_correct_clusters(cl_old,cl_cur,verbose=1)
# mode_correct_clusters(cl_old,cl_cur) == check1
# mode_correct_clusters(cl_old,cl_cur) == check2

cat(crayon::blue("- Double tie case\n"))
# tie case double
cl_old = c(1,2,1,1,2,1,3,3,4,4,5,5)
cl_cur = c(1,1,2,2,1,1,3,3,4,5,4,5)
# ora (pareggio) vedere sia il cluster 1 è entrato in 2 ma anche viceversa
# e idem per cluster 4 e 5. Vediamo cosa ritorna
mode_correct_clusters(cl_old,cl_cur,verbose=1)

cat(crayon::blue("- collapsed clusters test\n"))
# case of collapsed clusters
cl_old = c(1,1,1,2,2,2,1,1,3,4,4,4,4)
cl_cur = c(1,1,1,1,1,1,1,1,1,2,2,2,2) 
mode_correct_clusters(cl_old,cl_cur,verbose=1) # perfect!
# it stores the past value, 4, skipping the others; at least, that was my desired beahaviour

cat(crayon::blue("- orphan value test\n"))
# case of orphan value which was in a bigger cluster
cl_old = c(1,1,1,1,1,2,2,2,2,2,2,2)
cl_cur = c(2,2,2,2,2,1,1,1,1,1,1,3)
mode_correct_clusters(cl_old,cl_cur,verbose=1)

cat(crayon::blue("- complex case of cl_old\n"))
cl_old = c(2,2,2,2,2,2,2,2,2,2,2,2)
cl_cur = c(2,2,2,2,2,1,1,1,1,1,1,3)
mode_correct_clusters(cl_old,cl_cur,verbose=1)



################################################################################
### graph plot

library(igraph)
# g <- sample_gnp(10, 3/10)
# g_mst <- mst(g)
# 
# print(g)
# plot(g)
# plot(g_mst,edge.arrow.size = 0.5)
#
# plot(make_graph(c(1, 2, 2, 3, 3, 4, 5, 6), directed = FALSE))
# 
# #            lati:  A->B      B->C      C->D     D->As
# gg = make_graph(c("A", "B", "B", "C", "C", "D","D","A"), directed = FALSE)
# plot(gg)
# plot(mst(gg))
# 
# get.edgelist(gg)
# get.edgelist(mst(gg))


# si chiama my_sites sennò dava conflitto con sites dei plot :/
my_sites = data.frame(
	id = unique(df_weekly$IDStations),
	x = unique(df_weekly$Longitude),
	y = unique(df_weekly$Latitude)
)

stations = unique(df_weekly$IDStations)

graph_from_string <- function(x) {
	e <- str2expression(strsplit(x, ",")[[1]])
	do.call(igraph:::graph_from_literal_i, list(e))
}

stations_distance = function(st_1,st_2){
	coords_1 = my_sites[which(my_sites$id == st_1),c("x","y")]
	coords_2 = my_sites[which(my_sites$id == st_2),c("x","y")]
	dist2 = (coords_1$x - coords_2$x)^2 + (coords_1$y - coords_2$y)^2
	return(dist2)
}

cat(crayon::red("- assemble_edges(clusters)\n"))
assemble_edges = function(clusters,need_to_debug=0){
	edges_temp_list = NULL
	edges = data.frame(
		x = c(),
		y = c(),
		xend = c(),
		yend = c(),
		cluster = c()
	)
	for(cl in 1:max(clusters)){
		stations_here = stations[which(clusters==cl)]
		
		# graph_from_literal( A:B:C:D -- A:B:C:D )
		# crea un grafo con tutte le connessioni
		x = paste(paste(stations_here, collapse = ":"),"--",paste(stations_here, collapse = ":"))
		gg=graph_from_string(x)
		# plot(gg)
		# plot(mst(gg))
		
		# define pesi
		elist = get.edgelist(gg)
		pesi=c()
		if(need_to_debug==1){
			cat("dealing with cluster",cl,"\n")
		}
		if(size(elist)[1]>0){
			for(i in 1:size(elist)[1]){
				pesi = c(pesi,stations_distance(elist[i,1],elist[i,2]))
			}
			elist = get.edgelist(mst(gg,weights = pesi))
		} else{
			elist = get.edgelist(mst(gg))
		}
		
		edges_temp = data.frame(
			x = c(),
			y = c(),
			xend = c(),
			yend = c(),
			cluster = c()
		)
		if(size(elist)[1]>0){
			for (i in 1:size(elist)[1]){
				# c'è un problema con STA- nei nomi delle stazioni, che vengono presi come STA lato ecc
				edges_temp[i,c("x","y")] = my_sites[which(my_sites$id==elist[i,1]),c("x","y")]
				edges_temp[i,c("xend","yend")] = my_sites[which(my_sites$id==elist[i,2]),c("x","y")]
				edges_temp[i,"cluster"] = cl
			}
		} else {
			edges_temp[1,c("x","y")] = my_sites[which(my_sites$id==stations_here),c("x","y")]
			edges_temp[1,c("xend","yend")] = my_sites[which(my_sites$id==stations_here),c("x","y")]
			edges_temp[1,"cluster"] = cl
		}
		edges_temp_list[[cl]] = edges_temp
		edges = rbind(edges,edges_temp)
		# cat("current size",size(edges),"\n")
	}
	edges$cluster = factor(edges$cluster)
	return(edges)
}

cat(crayon::red("- assemble_edges_list(clusters)\n"))
assemble_edges_list = function(clusters,need_to_debug=0){
	edges_temp_list = NULL
	edges = data.frame(
		x = c(),
		y = c(),
		xend = c(),
		yend = c(),
		cluster = c()
	)
	for(cl in 1:max(clusters)){
		stations_here = stations[which(clusters==cl)]
		
		# graph_from_literal( A:B:C:D -- A:B:C:D )
		# crea un grafo con tutte le connessioni
		x = paste(paste(stations_here, collapse = ":"),"--",paste(stations_here, collapse = ":"))
		gg=graph_from_string(x)
		# plot(gg)
		# plot(mst(gg))
		
		# define pesi
		elist = get.edgelist(gg)
		pesi=c()
		if(need_to_debug==1){
			cat("dealing with cluster",cl,"\n")
		}
		if(size(elist)[1]>0){
			for(i in 1:size(elist)[1]){
				pesi = c(pesi,stations_distance(elist[i,1],elist[i,2]))
			}
			elist = get.edgelist(mst(gg,weights = pesi))
		} else{
			elist = get.edgelist(mst(gg))
		}
		
		edges_temp = data.frame(
			x = c(),
			y = c(),
			xend = c(),
			yend = c(),
			cluster = c()
		)
		if(size(elist)[1]>0){
			for (i in 1:size(elist)[1]){
				# c'è un problema con STA- nei nomi delle stazioni, che vengono presi come STA lato ecc
				edges_temp[i,c("x","y")] = my_sites[which(my_sites$id==elist[i,1]),c("x","y")]
				edges_temp[i,c("xend","yend")] = my_sites[which(my_sites$id==elist[i,2]),c("x","y")]
				edges_temp[i,"cluster"] = cl
			}
		} else {
			edges_temp[1,c("x","y")] = my_sites[which(my_sites$id==stations_here),c("x","y")]
			edges_temp[1,c("xend","yend")] = my_sites[which(my_sites$id==stations_here),c("x","y")]
			edges_temp[1,"cluster"] = cl
		}
		edges_temp_list[[cl]] = edges_temp
		edges = rbind(edges,edges_temp)
		# cat("current size",size(edges),"\n")
	}
	edges$cluster = factor(edges$cluster)
	return(edges_temp_list)
}


# test zone
# time = 4
# df_cluster_cut = df_cluster[df_cluster$Time==time,]
# clusters = df_cluster_cut$clusters
# edges = assemble_edges(clusters)


cat(crayon::bold("plotter library required!\n"))


cat(crayon::red("- get_graph_plot(df_cluster_cut)\n"))
get_graph_plot = function(df_cluster_cut){ # already mode_corrected
	clusters_now = df_cluster_cut$clusters
	edges_list = assemble_edges_list(clusters_now)
	
	p  <-  DefaultPlot()+
		geom_point(data = df_cluster_cut, aes(x = Longitude, y = Latitude,
											  # color = cols[factor(clusters)]), size = 2)+
											  color = cols[clusters_now]), size = 2)+
		labs(title = paste("Cluster map - time",time))
	
	q = p
	for(cl in 1:len(edges_list)){
		edges_to_plot = edges_list[[cl]]
		q = q + geom_segment(aes(x = x, y = y, xend = xend, yend = yend,
								 color = cols[as.numeric(cluster)]),
							 # color = paste0("cl",cl)),
							 linewidth=1.2,
							 data = edges_to_plot,show.legend = TRUE)+
			guides(color = guide_legend(title = "Clusters"))
	}
	
	q = q +
		scale_colour_identity(guide="legend",labels=paste0("cl",1:len(edges_list)),
							  breaks=cols[1:len(edges_list)])
}


cat(crayon::red("- get_hist_fill_plot(df_cluster_cut)\n"))
get_hist_fill_plot = function(df_cluster_cut,verbose=0){
	clusters_now = df_cluster_cut$clusters
	n_clusters = max(clusters_now)
	ycurrent = y[,paste0("w",time)]
	
	if(verbose==1){
	cat(crayon::red("Time",time,"\n"))
		for (cl in 1:n_clusters){
			cat("Cluster",cl,"- size",length(ycurrent[which(clusters_now==cl)]),
				"- mean",mean(ycurrent[which(clusters_now==cl)]),"\n")
			
		}
	}
	
	clust_vals = clusters_now[1:105]
	df_temp = data.frame(clusters=clust_vals,ycurrent=ycurrent)
	
	pad = 2	
	p = ggplot(df_temp, aes(ycurrent,
							fill = cols[clust_vals]
							# color = cols[clust_vals]
	))+
		
		geom_histogram(alpha=0.3,
					   # fill="white",
					   position="identity")+ # to have histograms
		# geom_density(alpha = 0.3)+ # to have the kernel/density estimation
		
		ggtitle(paste("Time",time))+
		# labs(title = paste("Cluster map - time",time))+
		guides(fill = guide_legend(title = "Clusters"))+
		
		# theme_classic()
		theme_bw()+
		# xlim(extrema(ycurrent)+c(-pad,pad))+
		
		scale_fill_identity(guide="legend",labels=paste0("cl",1:max(clust_vals)),
							breaks=cols[1:max(clust_vals)])+
		# scale_color_identity(guide="legend",labels=paste0("cl",1:max(clust_vals)),
		# breaks=cols[1:max(clust_vals)])
		xlab("log(PM10) values")+
		xlim(c(0.5,5))
	
}


cat(crayon::red("- get_hist_color_plot(df_cluster_cut)\n"))
get_hist_color_plot = function(df_cluster_cut,verbose=0){
	clusters_now = df_cluster_cut$clusters
	n_clusters = max(clusters_now)
	ycurrent = y[,paste0("w",time)]
	
	if(verbose==1){
	cat(crayon::red("Time",time,"\n"))
		for (cl in 1:n_clusters){
			cat("Cluster",cl,"- size",length(ycurrent[which(clusters_now==cl)]),
				"- mean",mean(ycurrent[which(clusters_now==cl)]),"\n")
			
		}
	}
	
	clust_vals = clusters_now[1:105]
	df_temp = data.frame(clusters=clust_vals,ycurrent=ycurrent)
	
	pad = 2	
	p = ggplot(df_temp, aes(ycurrent,
							# fill = cols[clust_vals]
							color = cols[clust_vals]
	))+
		
		geom_histogram(alpha=0.5,
					   fill="white",
					   position="identity")+ # to have histograms
		# geom_density(alpha = 0.3)+ # to have the kernel/density estimation
		
		ggtitle(paste("Time",time))+
		# labs(title = paste("Cluster map - time",time))+
		guides(color = guide_legend(title = "Clusters"))+
		
		# theme_classic()
		theme_bw()+
		# xlim(extrema(ycurrent)+c(-pad,pad))+
		
		# scale_fill_identity(guide="legend",labels=paste0("cl",1:max(clust_vals)),
		# breaks=cols[1:max(clust_vals)])
		scale_color_identity(guide="legend",labels=paste0("cl",1:max(clust_vals)),
							 breaks=cols[1:max(clust_vals)])+
		xlab("log(PM10) values")+
		xlim(c(0.5,5))
}


cat(crayon::red("- get_hist_continuos_plot(df_cluster_cut)\n"))
get_hist_continuos_plot = function(df_cluster_cut,verbose=1){
	clusters_now = df_cluster_cut$clusters
	n_clusters = max(clusters_now)
	ycurrent = y[,paste0("w",time)]
	cat(crayon::red("Time",time,"\n"))
	
	if(verbose==1){
		for (cl in 1:n_clusters){
			cat("Cluster",cl,"- size",length(ycurrent[which(clusters_now==cl)]),
				"- mean",mean(ycurrent[which(clusters_now==cl)]),"\n")
			
		}
	}
	
	clust_vals = clusters_now[1:105]
	df_temp = data.frame(clust_vals=clust_vals,ycurrent=ycurrent)
	
	pad = 2	
	# p = ggplot(df_temp, aes(ycurrent, fill = factor(clust_vals))) +
	p = ggplot(df_temp, aes(ycurrent, fill = cols[clust_vals])) +
		
		# scale_fill_manual(values = colora(n_clusters,77,0),name="Cluster")+
		scale_fill_identity(guide="legend",labels=paste0("cl",1:max(clust_vals)),
							breaks=cols[1:max(clust_vals)])+
		guides(fill = guide_legend(title = "Clusters"))+
		
		geom_density(alpha = 0.3)+
		ggtitle(paste("Time",time))+
		theme_bw()+
		xlab("log(PM10) values")+
		ylab("")+
		xlim(c(0.5,5))
}


cat(crayon::red("- plot_graph_and_hist(df_cluster_cut)\n"))
plot_graph_and_hist = function(df_cluster_cut){
	
# GRAPH #######################
q_graph = get_graph_plot(df_cluster_cut)

# HIST #######################
# by hand as we have to remove the legend here, while the function produces it
n_clusters = max(clusters_now)
ycurrent = y[,paste0("w",time)]
clust_vals = clusters_now[1:105]
df_temp = data.frame(clusters=clust_vals,ycurrent=ycurrent)
p = ggplot(df_temp, aes(ycurrent,
						color = cols[clust_vals]
))+
	geom_histogram(alpha=0.5,
				   fill="white",
				   position="identity")+ 
	ggtitle(paste("Time",time))+
	# guides(color = guide_legend(title = "Clusters"))+
	theme_bw()+
	theme(legend.position = "none")+
	scale_color_identity(guide="legend",labels=paste0("cl",1:max(clust_vals)),
						 breaks=cols[1:max(clust_vals)])+
	xlab("log(PM10) values")+
	xlim(c(0,5))
# xlim(extrema(df_weekly$AQ_pm10))
p_hist = p

grid.arrange(q_graph, p_hist, ncol=2,widths=c(1.5,1))
}
