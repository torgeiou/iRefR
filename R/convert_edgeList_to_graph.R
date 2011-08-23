####
# Generation of graphical objects from iRefIndex/edgeList tables:
####
convert_edgeList_to_graph = function(edgeList, directionality="undirected", graphical_package="igraph") {
	if (graphical_package == "graph") {
		# 1. Remove polymers and other self-loop edges:
		position_self_loop_edges = which(edgeList[,1]==edgeList[,2])
		nodes_in_loops = unique(edgeList[position_self_loop_edges,1])
		if (length(nodes_in_loops)>0) {
			print("Polymers and loop edges will be removed. If you want to work with a multigraph (loops included), use igraph as graphical_package.")
			edgeList = edgeList[-position_self_loop_edges,]	# get edgeList without loops
			#nodes_in_loops = setdiff(nodes_in_loops, unique(c(edgeList[,1], edgeList[,2])))	# only new nodes (not already in graph) left for this list
		}

		# 2. Generate graph:
		edges = edgeList[,1:2]
		weights = as.integer(edgeList[,3])
		the_graph = ftM2graphNEL(edges, W=weights, edgemode=directionality)
		#if (length(nodes_in_loops)>0) {
		#	for (i in nodes_in_loops) {
		#		the_graph = addNode(i, the_graph)
		#	}
		#}
	} else {
		edgel = as.data.frame(edgeList)
		colnames(edgel) = c("V1", "V2", "weight")
		if (directionality=="directed") {
			direction_flag = "TRUE"
		} else {
			direction_flag = "FALSE"
		}
		the_graph = graph.data.frame(edgel, directed=direction_flag)
	}

	output = the_graph
}
