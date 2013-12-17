    % File /man/convert_MITAB_to_edgeList.Rd
     \name{convert_MITAB_to_edgeList}
     \alias{convert_MITAB_to_edgeList}
     \title{Convert MITAB format to edgeList format}
     \description{
       Convert a table from MITAB format to edgeList format.

On iROGs and icROGs:
"iRefIndex" guarantees the non-redundancy of protein information by assigning a different protein identifier (called Redundant Object Group, ROG) to every different protein sequence. This is called the "non-canonical" representation. At the same time, groups of non-redundant proteins might be different isoforms of a given protein, and, in this case, the identifier of one protein of the group (called cROG or canonical ROG) is chosen to represent the entire group of similar proteins. This is called the "canonical" representation of proteins. The iROG and icROG IDs correspond to integer representations of the ROG and cROG, respectively. Most of the "iRefR" functions work with iROGs and icROGs.

"iRefIndex" also stores the Unique Identifier, UID, which is a protein identifier that represents a protein with its UniProt identifier wherever possible. If this is unknown, uses its RefSeq identifier and, if none of these exists, it uses the ROG ID.
     }
     \usage{
     convert_MITAB_to_edgeList(MITAB_table, edge_weight_col, complex_rep, canonical_rep,
 directionality, node_names, multi_edge)
     }
     \arguments{
       \item{MITAB_table}{iRefIndex/MITAB R table.}
       \item{edge_weight_col}{The number of the column of the MITAB_table that includes some quantitative score of an interaction that we would like to use as edge weight in a graphical representation of the table. In case the user is not interested in weighting edges, the "default" option sets the value of "1" for all edge weights. Default="default".}
       \item{complex_rep}{Three possible ways of representing a protein complex: "bipartite" is the format followed by iRefIndex, consisting in pairwise interactions between a node that represents the complex and each of the protein members of the complex. "matrix" generates the matrix model of the complex, i.e., all possible pairwise interactions among the members of the complex. "spoke" generates the spoke model of the complex, i.e., pairwise interactions between one single chosen protein of the complex and each of the other members of the complex. Default="spoke".}
       \item{canonical_rep}{Either "yes" for working with the canonical representation (icROGs) or "no" for working with the non-canonical representation (iROGs). Default="yes".}
       \item{directionality}{Either "directed" if we want to produce an ordered list of directed edges (where edge A-B != edge B-A), or "undirected" in the opposite case (edge A-B = edge B-A). "Directed" lists are only possible for binary interactions and spoke models of complex interactions. In case of binary interactions, directed edges point from the bait to a prey or to another bait. In case of spoke models of complexes, the bait is chosen as the center of the spoke model; if there is no bait, the first protein of the list is arbitrarily chosen as a center. Default="undirected".}
       \item{node_names}{Either "uids" if we want to identify the proteins according to the iRefIndex's Unique Identifier, UID, or "rogs" if we want to use the iRefIndex ROG identifiers. Default="rogs".}
       \item{multi_edge}{Either "yes" if we want to observe all edges with different RIG IDs between two nodes (for example, edges between two nodes belonging to different complexes or to a binary and a complex interaction), or "no" if we want to combine such edges and work with maximum one edge between every pair of nodes. Default="no".}
     }
     \value{
     \item{edgeList}{iRefIndex/edgeList R table with the previously specified parameters.}
     }
     \author{Antonio Mora <antoniocmora@gmail.com>}
     \examples{
     ## get tables
     irefindex_curr_ecoli = get_irefindex("562", "13.0", tempdir())
     all_INTACT = select_database("intact", irefindex_curr_ecoli, "this_database")
     binary_INTACT = select_interaction_type("binary", all_INTACT)
     complex_INTACT = select_interaction_type("complex", all_INTACT)

     ## execute function
     # NOTE: The matrix model is more time-consuming.
     edgeList_binary_INTACT = convert_MITAB_to_edgeList(binary_INTACT)
     edgeList_complex_INTACT_s = convert_MITAB_to_edgeList(complex_INTACT, "default",
 "spoke")

     \dontrun{
     edgeList_complex_INTACT_m = convert_MITAB_to_edgeList(complex_INTACT, "default",
 "matrix")
     edgeList_all_INTACT = convert_MITAB_to_edgeList(all_INTACT, "default", "spoke")

     edgeList_binary_INTACT_dir = convert_MITAB_to_edgeList(binary_INTACT, "default",
 "bipartite", "yes", "directed")
     edgeList_complex_INTACT_sdir = convert_MITAB_to_edgeList(complex_INTACT, "default",
 "spoke", "yes", "directed")
     edgeList_all_INTACT_dir = convert_MITAB_to_edgeList(all_INTACT, "default", "spoke",
 "yes", "directed")

     # Multi-edge interaction (C|X case) using UIDs:
     edgeList_irefindex = convert_MITAB_to_edgeList(irefindex_curr_ecoli, node_names=
 "uids", multi_edge="no")
     chosen_int = 1		# Choose first multiple edge in the data frame
     rigids_multi_edge = edgeList_irefindex[grep("C\\|X", edgeList_irefindex[,4])[[
 chosen_int]], 3]
     rigids_multi_edge = strsplit(rigids_multi_edge, "|", fixed="TRUE")[[1]]

     #edgeList_irefindex = convert_MITAB_to_edgeList(irefindex_curr_ecoli, node_names=
     #"uids", multi_edge="yes")
     index_interesting_multi_edge = unlist(lapply(rigids_multi_edge, grep,
 edgeList_irefindex[,3]))
     edgeList_interesting_complex = edgeList_irefindex[index_interesting_multi_edge,]
     graph_irefindex = convert_edgeList_to_graph(edgeList_interesting_complex)
     list_edge_colors = NULL
     list_edge_colors[which(edgeList_interesting_complex[,3]=="625881")] = "blue"
     list_edge_colors[which(edgeList_interesting_complex[,3]=="634347")] = "red"
     list_edge_colors[which(edgeList_interesting_complex[,3]=="922617")] = "green"
     list_edge_colors[which(edgeList_interesting_complex[,3]=="982855")] = "purple"
     graph_irefindex <- set.edge.attribute(graph_irefindex, "color", value=
 list_edge_colors)

     plot(graph_irefindex, layout=layout.kamada.kawai, vertex.label=
 V(graph_irefindex)$name, vertex.size=7, main="multi_edge interactions")
     legend("bottomright", c("complex1", "binary", "complex2", "complex3"),
 fill=c("blue", "red", "green", "purple"))
     }
     }
