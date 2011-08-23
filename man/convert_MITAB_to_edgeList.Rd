    % File /man/convert_MITAB_to_edgeList.Rd
     \name{convert_MITAB_to_edgeList}
     \alias{convert_MITAB_to_edgeList}
     \title{Convert MITAB format to edgeList format}
     \description{
       Convert a table from MITAB format to edgeList format.

On iROGs and icROGs:
"iRefIndex" guarantees the non-redundancy of protein information by assigning a different protein identifier (called Redundant Object Group, ROG) to every different protein sequence. This is called the "non-canonical" representation. At the same time, groups of non-redundant proteins might be different isoforms of a given protein, and, in this case, the identifier of one protein of the group (called cROG or canonical ROG) is chosen to represent the entire group of similar proteins. This is called the "canonical" representation of proteins. The iROG and icROG IDs correspond to integer representations of the ROG and cROG, respectively. Most of the "iRefR" functions work with iROGs and icROGs.
     }
     \usage{
     convert_MITAB_to_edgeList(MITAB_table, edge_weight_col, complex_rep, canonical_rep, directionality)
     }
     \arguments{
       \item{MITAB_table}{iRefIndex/MITAB R table.}
       \item{edge_weight_col}{The number of the column of the MITAB_table that includes some quantitative score of an interaction that we would like to use as edge weight in a graphical representation of the table. In case the user is not interested in weighting edges, the "default" option sets the value of "1" for all edge weights. Default="default".}
       \item{complex_rep}{Three possible ways of representing a protein complex: "bipartite" is the format followed by iRefIndex, consisting in pairwise interactions between a node that represents the complex and each of the protein members of the complex. "matrix" generates the matrix model of the complex, i.e., all possible pairwise interactions among the members of the complex. "spoke" generates the spoke model of the complex, i.e., pairwise interactions between one single chosen protein of the complex and each of the other members of the complex. Default="spoke".}
       \item{canonical_rep}{Either "yes" for working with the canonical representation (icROGs) or "no" for working with the non-canonical representation (iROGs). Default="yes".}
       \item{directionality}{Either "directed" if we want to produce an ordered list of directed edges (where edge A-B != edge B-A), or "undirected" in the opposite case (edge A-B = edge B-A). "Directed" lists are only possible for binary interactions and spoke models of complex interactions. In case of binary interactions, directed edges point from the bait to a prey or to another bait. In case of spoke models of complexes, the bait is chosen as the center of the spoke model; if there is no bait, the first protein of the list is arbitrarily chosen as a center. Default="undirected".}
     }
     \value{
     \item{edgeList}{iRefIndex/edgeList R table with the previously specified parameters.}
     }
     \author{Antonio Mora <a.m.ortiz@biotek.uio.no>}
     \examples{
     ## get tables
     irefindex_curr_ecoli = get_irefindex("562", "current", tempdir())
     all_INTACT = select_database("intact", irefindex_curr_ecoli, "this_database")
     binary_INTACT = select_interaction_type("binary", all_INTACT)
     complex_INTACT = select_interaction_type("complex", all_INTACT)

     ## execute function
     # NOTE: The matrix model is more time-consuming. Generating an edgeList for the whole human iRefIndex using a matrix model can take less than 2 minutes.
     edgeList_binary_INTACT = convert_MITAB_to_edgeList(binary_INTACT)
     edgeList_complex_INTACT_s = convert_MITAB_to_edgeList(complex_INTACT, "default", "spoke")
     edgeList_complex_INTACT_m = convert_MITAB_to_edgeList(complex_INTACT, "default", "matrix")
     edgeList_all_INTACT = convert_MITAB_to_edgeList(all_INTACT, "default", "spoke")

     edgeList_binary_INTACT_dir = convert_MITAB_to_edgeList(binary_INTACT, "default", "bipartite", "yes", "directed")
     edgeList_complex_INTACT_sdir = convert_MITAB_to_edgeList(complex_INTACT, "default", "spoke", "yes", "directed")
     edgeList_all_INTACT_dir = convert_MITAB_to_edgeList(all_INTACT, "default", "spoke", "yes", "directed")
     }
