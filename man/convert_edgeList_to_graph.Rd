    % File /man/convert_edgeList_to_graph.Rd
     \name{convert_edgeList_to_graph}
     \alias{convert_edgeList_to_graph}
     \title{Convert iRefIndex/edgeList table to Graph}
     \description{
       Create R graphical objects of an iRefIndex/edgeList interaction table.
     }
     \usage{
     convert_edgeList_to_graph(edgeList, directionality, graphical_package)
     }
     \arguments{
       \item{edgeList}{iRefIndex/edgeList R table.}
       \item{directionality}{Either "directed" for directed graphs or "undirected" for undirected graphs. The chosen value must be the same used to generate the edgeList table. Default="undirected".}
       \item{graphical_package}{Either "graph" or "igraph", depending on the R package the user wants to use. If you use "graph", no loop edges will be allowed. Default="igraph".}
     }
     \value{
     \item{output}{R graphical object corresponding to the original edgeList.}
     }
     \author{Antonio Mora <a.m.ortiz@biotek.uio.no>}
     \examples{
     ## get tables
     irefindex_curr_ecoli = get_irefindex("562", "current", tempdir())
     all_INTACT = select_database("intact", irefindex_curr_ecoli)
     binary_INTACT = select_interaction_type("binary", all_INTACT)
     complex_INTACT = select_interaction_type("complex", all_INTACT)

     edgeList_binary_INTACT = convert_MITAB_to_edgeList(binary_INTACT)
     edgeList_complex_INTACT_s = convert_MITAB_to_edgeList(complex_INTACT, "default", "spoke")
     edgeList_complex_INTACT_m = convert_MITAB_to_edgeList(complex_INTACT, "default", "matrix")
     edgeList_all_INTACT = convert_MITAB_to_edgeList(all_INTACT, "default", "spoke")
     edgeList_binary_INTACT_dir = convert_MITAB_to_edgeList(binary_INTACT, "default", "bipartite", "yes", "directed")
     edgeList_complex_INTACT_sdir = convert_MITAB_to_edgeList(complex_INTACT, "default", "spoke", "yes", "directed")

     ## execute function
     graph_binary_INTACT = convert_edgeList_to_graph(edgeList_binary_INTACT, "undirected", "igraph")
     graph_complex_INTACT_s = convert_edgeList_to_graph(edgeList_complex_INTACT_s)
     graph_complex_INTACT_m = convert_edgeList_to_graph(edgeList_complex_INTACT_m)

     \dontrun{
     graph_all_INTACT = convert_edgeList_to_graph(edgeList_all_INTACT)

     graph_binary_INTACT_dir = convert_edgeList_to_graph(edgeList_binary_INTACT_dir, "directed")
     graph_complex_INTACT_sdir = convert_edgeList_to_graph(edgeList_complex_INTACT_sdir, "directed")
     # For the binary case, edge number is smaller than the undirected case because only interactions
     # with bait information are used. For the complex spoke case, the edge number is bigger than the
     # undirected case because edges with two baits give raise to two directed edges and only one
     # undirected.

     ## Additional info if you are using "graph":
     library(Rgraphviz)
     edgeList_complex = convert_MITAB_to_edgeList(complex_INTACT, node_names="uids")
     graph_complex = convert_edgeList_to_graph(edgeList_complex, graphical_package="graph")
     plot(graph_complex, "neato")

     nodes(graph_complex)
     edges(graph_complex)
     degree(graph_complex)
     numNoEdges(graph_complex)
     mostEdges(graph_complex)
     adj(graph_complex,"uniprotkb:P28043")
     acc(graph_complex, "uniprotkb:P28043")
     sg <- subGraph(nodes(graph_complex)[1:20], graph_complex)
     plot(sg, "neato")

     ## Additional info if you are using "igraph":
     V(graph_binary_INTACT)$name
     V(graph_binary_INTACT)[igraph::degree(graph_binary_INTACT) == 1]$name
     par(mfrow=c(2,3))
     plot(graph_binary_INTACT, layout=layout.fruchterman.reingold, vertex.size=3, vertex.color="green", frame=TRUE, main="Binary", vertex.label=NA)
     plot(graph_complex_INTACT_s, layout=layout.fruchterman.reingold, vertex.size=3, vertex.color="green", frame=TRUE, main="Complex Spoke", vertex.label=NA)
     plot(graph_complex_INTACT_m, layout=layout.fruchterman.reingold, vertex.size=3, vertex.color="green", frame=TRUE, main="Complex Matrix", vertex.label=NA)
     plot(graph_all_INTACT, layout=layout.fruchterman.reingold, vertex.size=3, vertex.color="green", frame=TRUE, main="Binary and Complex Spoke", vertex.label=NA)
     plot(graph_binary_INTACT_dir, layout=layout.fruchterman.reingold, vertex.size=3, edge.arrow.size=0.3, vertex.color="green", frame=TRUE, main="Binary Directed", vertex.label=NA)
     plot(graph_complex_INTACT_sdir, layout=layout.fruchterman.reingold, vertex.size=3, edge.arrow.size=0.3, vertex.color="green", frame=TRUE, main="Complex Spoke Directed", vertex.label=NA)
     x11(); par(mfrow=c(1,2))
     plot(degree.distribution(graph_binary_INTACT), main="Degree Distribution (Binary, INTACT, E.coli)")
     plot(degree.distribution(graph_binary_INTACT),log="xy", main="Log-Log Degree Distribution (Binary, INTACT, E.coli)")
     }
     }
