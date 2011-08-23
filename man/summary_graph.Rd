    % File /man/summary_graph.Rd
     \name{summary_graph}
     \alias{summary_graph}
     \title{Get Summary Information for a Graph}
     \description{
       Get statistical information about the number of nodes, edges, degree distribution, nodes per connected component and graphs of connected components in an R Graph.

Here "iRefR" uses the "graph", "igraph" and "RBGL" libraries.
     }
     \usage{
     summary_graph(graph_object, graphical_package)
     }
     \arguments{
       \item{graph_object}{R graph to summarize.}
       \item{graphical_package}{Either "graph" or "igraph", depending on the R package the user wants to use. Default="igraph".}
     }
     \value{
     \item{output}{R list containing: Number of nodes and edges, degree distribution, nodes per connected component and graphs per connected component. A plot of the graph is also produced if the "igraph" option is chosen.}
     }
     \author{Antonio Mora <a.m.ortiz@biotek.uio.no>}
     \examples{
     ## get tables
     irefindex_curr_ecoli = get_irefindex("562", "current", tempdir())
     INTACT_curr_ecoli = select_database("intact", irefindex_curr_ecoli, "this_database")
     edgeList_ecoli = convert_MITAB_to_edgeList(INTACT_curr_ecoli, "default", "spoke")
     graph_ecoli = convert_edgeList_to_graph(edgeList_ecoli, graphical_package="igraph")

     ## execute function
     if(interactive()) {
     	summary = summary_graph(graph_ecoli, "igraph")
     	summary$nodes_and_edges
     	table(summary$degree_distribution)
     	x11(); hist(summary$degree_distribution)

     	table(summary$nodes_per_connected_component)
     	x11(); barplot(summary$nodes_per_connected_component, xlab="Size Components", ylab="Number Components")
     	title(main="Distribution of Connected Components")

     	summary$graph_per_module[[1]]
     }
     }
