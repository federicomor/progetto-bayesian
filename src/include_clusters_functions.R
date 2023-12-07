cat(crayon::cyan("\nGenerated these functions:\n"))

################################################################################
#### mode correction
mode_correct_clusters = function(cl_old, cl_cur){ # and returns cl_new the updated cl_cur
	if(is.null(cl_old)){
		warning("The cluster_old vector is NULL.")
		return(cl_cur) # nothing to correct
	}
	
	cl_new = cl_cur
	new_labels = c()
	tied_labels = c() # solve them later
	
	for (i in 1:max(cl_cur)){
		indici_label_i = which(cl_cur == i)
		freq = table(cl_old[indici_label_i])
		better_labels = as.numeric(names(freq[freq == max(freq)]))
		num_mode = length(better_labels)
		
		if(num_mode>=2){ # per risolvere i pareggi
			# warning("Tie case happened.\n")
			
			tied_labels = c(tied_labels,i)
		} else{
			better_label = better_labels
			new_labels = c(new_labels,better_label)
			cl_new[indici_label_i] = better_label
		}
		# cat(new_labels,"\n")
	}
	max_it = 100
	it = 1
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
	if(it==max_it){
		warning("Tie cases not solved :'/")
		return(cl_cur)
	} else {
		warning("Tie cases solved :)")
		return(cl_new)
	}
	
}
cat(crayon::red("- mode_correct_clusters(cl_old, cl_cur)\n"))

### test sets
# cl_old = c(1,2,1,1,2,2,3,3)
# cl_cur = c(1,1,2,2,1,1,3,3)
# check_ = c(2,2,1,1,2,2,3,3)
# mode_correct_clusters(cl_old,cl_cur)
# mode_correct_clusters(cl_old,cl_cur) == check_
# 
# # tie case
# #   units  1 2 3 4 5 6 7 8 9 10
# cl_old = c(1,2,1,1,2,1,3,3)
# cl_cur = c(1,1,2,2,1,1,3,3)
# # ora (pareggio) vedere sia il cluster 1 è entrato in 2 ma anche viceversa
# check1 = c(2,2,1,1,2,2,3,3)
# check2 = c(1,1,2,2,1,1,3,3)
# mode_correct_clusters(cl_old,cl_cur)
# mode_correct_clusters(cl_old,cl_cur) == check1
# mode_correct_clusters(cl_old,cl_cur) == check2
# 
# # tie case double
# cl_old = c(1,2,1,1,2,1,3,3,4,4,5,5)
# cl_cur = c(1,1,2,2,1,1,3,3,4,5,4,5)
# # ora (pareggio) vedere sia il cluster 1 è entrato in 2 ma anche viceversa
# # e idem per cluster 4 e 5. Vediamo cosa ritorna
# mode_correct_clusters(cl_old,cl_cur)
# 
# cl_old = c(1,1,1,2,2,2,1,1,3,4,4,4,4) 
# cl_cur = c(1,1,1,1,1,1,1,1,1,2,2,2,2) # case of collapsed clusters
# mode_correct_clusters(cl_old,cl_cur) # perfect!
# # it stores the past value, 4, skipping the others; at least, that was my desired beahaviour



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
		edges = rbind(edges,edges_temp)
		# cat("current size",size(edges),"\n")
	}
	edges$cluster = factor(edges$cluster)
	return(edges)
}

# test zone
# time = 4
# df_cluster_cut = df_cluster[df_cluster$Time==time,]
# clusters = df_cluster_cut$clusters
# edges = assemble_edges(clusters)


