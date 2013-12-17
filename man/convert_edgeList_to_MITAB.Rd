    % File /man/convert_edgeList_to_MITAB.Rd
     \name{convert_edgeList_to_MITAB}
     \alias{convert_edgeList_to_MITAB}
     \title{Convert edgeList format to MITAB format}
     \description{
       Convert a table from edgeList format to MITAB format.

On iROGs and icROGs:
"iRefIndex" guarantees the non-redundancy of protein information by assigning a different protein identifier (called Redundant Object Group, ROG) to every different protein sequence. This is called the "non-canonical" representation. At the same time, groups of non-redundant proteins might be different isoforms of a given protein, and, in this case, the identifier of one protein of the group (called cROG or canonical ROG) is chosen to represent the entire group of similar proteins. This is called the "canonical" representation of proteins. The iROG and icROG IDs correspond to integer representations of the ROG and cROG, respectively. Most of the "iRefR" functions work with iROGs and icROGs.
     }
     \usage{
     convert_edgeList_to_MITAB(edgeList, MITAB_table, canonical_rep, int_type)
     }
     \arguments{
       \item{edgeList}{iRefIndex/edgeList R table.}
       \item{MITAB_table}{Full iRefIndex/MITAB R table for a given organism.}
       \item{canonical_rep}{Either "yes" for working with the canonical representation (icROGs) or "no" for working with the non-canonical representation (iROGs). Default="yes".}
       \item{int_type}{Interaction Type: "binary" for binary interactions, "complex" for complexes, "polymer" for polymers and "all" for the three subsets. Default="all".}
     }
     \value{
     \item{MITAB_table}{iRefIndex/MITAB R table corresponding to the original edgeList.}
     }
     \author{Antonio Mora <antoniocmora@gmail.com>}
     \examples{
     ## get tables
     irefindex_curr_ecoli = get_irefindex("562", "13.0", tempdir())
     all_INTACT = select_database("intact", irefindex_curr_ecoli, "this_database")
     binary_INTACT = select_interaction_type("binary", all_INTACT)
     complex_INTACT = select_interaction_type("complex", all_INTACT)
     edgeList_binary_INTACT = convert_MITAB_to_edgeList(binary_INTACT)
     edgeList_complex_INTACT_s = convert_MITAB_to_edgeList(complex_INTACT, "default",
 "spoke")
     edgeList_all_INTACT = convert_MITAB_to_edgeList(all_INTACT, "default", "spoke")

     ## execute function
     new_binary = convert_edgeList_to_MITAB(edgeList_binary_INTACT, irefindex_curr_ecoli,
 "yes", "binary")
     reconstructed_binary_INTACT = select_database("intact", new_binary, "this_database")
     setequal(dim(binary_INTACT), dim(reconstructed_binary_INTACT))

     new_complex_s = convert_edgeList_to_MITAB(edgeList_complex_INTACT_s,
 irefindex_curr_ecoli, "yes", "complex")
     reconstructed_complex_INTACT_s = select_database("intact", new_complex_s,
 "this_database")
     setequal(dim(complex_INTACT), dim(reconstructed_complex_INTACT_s))

     new_all = convert_edgeList_to_MITAB(edgeList_all_INTACT, irefindex_curr_ecoli)
     reconstructed_all_INTACT = select_database("intact", new_all, "this_database")
     setequal(dim(all_INTACT), dim(reconstructed_all_INTACT))

     # Note that the dimensions would not be equal for the directed case at previous
     # examples because only records with baits can be reconstructed. Directed edgeLists
     # only include information with baits, so that all other cases are ignored during
     # edgeList generation and there is no way to go back to the original MITAB.
     }
