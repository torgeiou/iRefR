    % File /man/merge_complexes_lists.Rd
     \name{merge_complexes_lists}
     \alias{merge_complexes_lists}
     \title{Merge tables in complexList Format}
     \description{
       Merge a list of tables in complexList format, written in a prioritized order, and produce a single table where redundancies where removed according to the specified priorities.

The "complexList" format is one of the three table formats used by "iRefR", together with the MITAB format and the edgeList. The "complexList" format is a table with two columns, where the first column corresponds to an identifier and the second column to a group of proteins, commonly a complex, written in a comma-separated string.
     }
     \usage{
     merge_complexes_lists(list_of_complexLists)
     }
     \arguments{
       \item{list_of_complexLists}{List of tables in the "complexList" format. All of these tables must share the same representation of proteins, either canonical or non-canonical. Mixes would produce meaningless results. Repeated complexes inside individual tables would be removed. At the same time, when some group of proteins is repeated between two or more tables, it will be removed from all tables with less priority and will only stay in the table with the highest priority. The highest priority would be assigned to the first written table, while the smallest priority would go to the last one.}
     }
     \value{
     \item{result}{Table in the "complexList" format, with the prioritized merging of a set of tables in the "complexList" format.}
     }
     \author{Antonio Mora <a.m.ortiz@biotek.uio.no>}
     \examples{
     ## get tables
     irefindex_curr_ecoli = get_irefindex("562", "current", tempdir())
     iRef_binary = select_interaction_type("binary", irefindex_curr_ecoli)
     iRef_complex = select_interaction_type("complex", irefindex_curr_ecoli)
     complex_complexList = convert_MITAB_to_complexList(iRef_complex, "yes", "no")
     binary_complexList = convert_MITAB_to_complexList(iRef_binary, "yes", "yes")

     ## execute function
     consolidated_complexList = merge_complexes_lists(list(complex_complexList, binary_complexList))
     number_complexes_repeated_inside_lists_or_between_lists = dim(complex_complexList)[1] + dim(binary_complexList)[1] - dim(consolidated_complexList)[1]
     }
