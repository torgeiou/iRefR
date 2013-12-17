    % File /man/convert_graph_to_edgeList.Rd
     \name{convert_graph_to_edgeList}
     \alias{convert_graph_to_edgeList}
     \title{Convert Graph to iRefIndex/edgeList table}
     \description{
       Generate an iRefIndex/edgeList interaction table from an R graphical object.
     }
     \usage{
     convert_graph_to_edgeList(graph, directionality, graphical_package)
     }
     \arguments{
       \item{graph}{R graphical object.}
       \item{directionality}{Either "directed" for directed graphs or "undirected" for undirected graphs. Default="undirected".}
       \item{graphical_package}{Either "graph" or "igraph", depending on the R package the user used to generate the graph. If "graph" was used, no loop edges were allowed and can not be reconstructed. Default="igraph".}
     }
     \value{
     \item{output}{iRefIndex/edgeList R table corresponding to the original graph.}
     }
     \author{Antonio Mora <antoniocmora@gmail.com>}
     \examples{
     ## get tables
     irefindex_curr_ecoli = get_irefindex("562", "13.0", tempdir())
     all_INTACT = select_database("intact", irefindex_curr_ecoli)
     binary_INTACT = select_interaction_type("binary", all_INTACT)
     complex_INTACT = select_interaction_type("complex", all_INTACT)

     ## generate graphs
     edgeList_binary_INTACT = convert_MITAB_to_edgeList(binary_INTACT)
     edgeList_binary_INTACT_dir = convert_MITAB_to_edgeList(binary_INTACT, "default",
 "bipartite", "yes", "directed")
     edgeList_complex_INTACT_s = convert_MITAB_to_edgeList(complex_INTACT, "default",
 "spoke")
     edgeList_complex_INTACT_sdir = convert_MITAB_to_edgeList(complex_INTACT, "default",
 "spoke", "yes", "directed")

     graph_binary_INTACT = convert_edgeList_to_graph(edgeList_binary_INTACT,
 "undirected")
     graph_binary_INTACT_dir = convert_edgeList_to_graph(edgeList_binary_INTACT_dir,
 "directed")
     graph_complex_INTACT_s = convert_edgeList_to_graph(edgeList_complex_INTACT_s)
     graph_complex_INTACT_sdir = convert_edgeList_to_graph(edgeList_complex_INTACT_sdir,
 "directed")

     ## execute function
     reconstructed_edgeList_binary_INTACT = convert_graph_to_edgeList(
 graph_binary_INTACT)
     setequal(dim(edgeList_binary_INTACT), dim(reconstructed_edgeList_binary_INTACT))

     reconstructed_edgeList_binary_INTACT_dir = convert_graph_to_edgeList(
 graph_binary_INTACT_dir, "directed")
     setequal(dim(edgeList_binary_INTACT_dir), dim(
 reconstructed_edgeList_binary_INTACT_dir))

     reconstructed_edgeList_complex_INTACT_s = convert_graph_to_edgeList(
 graph_complex_INTACT_s)
     setequal(dim(edgeList_complex_INTACT_s), dim(
 reconstructed_edgeList_complex_INTACT_s))

     reconstructed_edgeList_complex_INTACT_sdir = convert_graph_to_edgeList(
 graph_complex_INTACT_sdir, "directed")
     setequal(dim(edgeList_complex_INTACT_sdir), dim(
 reconstructed_edgeList_complex_INTACT_sdir))
     # NOTE: All four tests give "TRUE" when using "igraph" as graphical package, while
     # give "FALSE" when using "graph", due to the lack of polymer information in the
     # "graph" option.
     }
