    % File /man/convert_MITAB_to_complexList.Rd
     \name{convert_MITAB_to_complexList}
     \alias{convert_MITAB_to_complexList}
     \title{Convert MITAB format to complexList format}
     \description{
       Convert a table from MITAB format to complexList format. To do that, convert all interactions listed as complexes in the MITAB table and, optionally, check if some groups of binary interactions can be interpreted as complexes (possibly misrepresented or binary-represented), according to an "iRefR" algorithm that searches for groups of proteins coming from co-precipitation experiments, sharing the same PMID, source database and bait (when available).

This function might generate repeated complexes; f.ex. when the same complex is generated from two different sources. If the user wants a list of unique complexes, the "unique" function must be applied after generating the complexList.

The "complexList" format is one of the three table formats used by "iRefR", together with the MITAB format and the edgeList. The "complexList" format is a table with two columns, where the first column corresponds to an identifier and the second column to a group of proteins, commonly a complex, written as a comma-separated string. This representation has less information than MITAB but simplifies the visualization of groups of proteins, such as complexes.

On iROGs and icROGs:
"iRefIndex" guarantees the non-redundancy of protein information by assigning a different protein identifier (called Redundant Object Group, ROG) to every different protein sequence. This is called the "non-canonical" representation. At the same time, groups of non-redundant proteins might be different isoforms of a given protein, and, in this case, the identifier of one protein of the group (called cROG or canonical ROG) is chosen to represent the entire group of similar proteins. This is called the "canonical" representation of proteins. The iROG and icROG IDs correspond to integer representations of the ROG and cROG, respectively. Most of the "iRefR" functions work with iROGs and icROGs.
     }
     \usage{
     convert_MITAB_to_complexList(MITAB_table, canonical_rep,
 include_generated_complexes, list_methods, bait_use, node_names)
     }
     \arguments{
       \item{MITAB_table}{iRefIndex/MITAB R table.}
       \item{canonical_rep}{Either "yes" for working with the canonical representation (icROGs) or "no" for working with the non-canonical representation (iROGs). Default="yes".}
       \item{include_generated_complexes}{Parameter to include a set of complexes generated from binary-represented interactions that might be spoke-represented or misrepresented complexes, according to an "iRefR" algorithm that searches for groups of proteins coming from co-precipitation experiments, sharing the same PMID, source database and bait (when available). The user can choose between ignoring those records for which bait-prey information is not available or generating complexes from proteins with no bait information but 3 or more neighbors. To be used when the user can accept hypothetical complexes in the analysis. Default="no".}
       \item{list_methods}{Vector with a list of experimental methods that the user assumes should offer co-complexing information (when regenerating complexes). Methods should be described according to the MI-ontology format (see: http://www.ebi.ac.uk/ontology-lookup/browse.do?ontName=MI under "interaction detection method"). The user can either introduce a vector with the MI terms corresponding to the methods of interest or use the "default" parameter. The methods listed under "default" are: MI:0004, MI:0006, MI:0007, MI:0019, MI:0027, MI:0028, MI:0029, MI:0059, MI:0061, MI:0071, MI:0096, MI:0114, MI:0226, MI:0227, MI:0401, MI:0437, MI:0676, MI:0858 and MI:0963. The description of these methods can be found at the above-mentioned web address. Some other methods the user might be interested to add might be: MI:0067, MI:0038, MI:0104, MI:0826, MI:0069, MI:0894, MI:0893, MI:0938, MI:0099, MI:0888, MI:0966, MI:0825, MI:0824, MI:0025, MI:0928, MI:1022, MI:0982, MI:0807, MI:0276, MI:0404, MI:0808, MI:0412 and MI:0413. Default="default".}
       \item{bait_use}{Decision of either using only records with bait-prey information ("only_baits") or including records with no bait-prey information ("include_non_baits"), when regenerating complexes from binaries that share the same selected interaction method, PMID and source database. If the first case is chosen, all complexes will be generated as a spoke model around the bait protein. If the second case is chosen, complexes will be additionally generated as spoke models around every protein with 3 or more neighbors, for groups of proteins without bait-prey information. Default="only_baits".}
       \item{node_names}{Decision of either using iRefIndex ROG Identifiers or iRefIndex UIDs. Default="rogs".}
     }
     \value{
     \item{complexList}{iRefIndex/complexList R table corresponding to the original MITAB table.}
     }
     \author{Antonio Mora <antoniocmora@gmail.com>}
     \examples{
     \dontrun{
     ## get tables
     irefindex_curr_ecoli = get_irefindex("562", "13.0", tempdir())
     iRef_binary = select_interaction_type("binary", irefindex_curr_ecoli)

     ## execute function
     known_complexes_ecoli_complexList = convert_MITAB_to_complexList(
 irefindex_curr_ecoli, canonical_rep="yes", include_generated_complexes="no")

     ## more examples
     generated_complexes_ecoli_complexList = convert_MITAB_to_complexList(iRef_binary,
 canonical_rep="yes", include_generated_complexes="yes", list_methods="default")
     all_complexes_ecoli_complexList = convert_MITAB_to_complexList(irefindex_curr_ecoli,
 canonical_rep="yes", include_generated_complexes="yes", list_methods="default")
     number_unique_complexes_all = length(unique(all_complexes_ecoli_complexList[,2]))

     known_complexes_ecoli_complexList_extended = convert_MITAB_to_complexList(
 irefindex_curr_ecoli, canonical_rep="yes", include_generated_complexes="no",
 list_methods="default", bait_use="include_non_baits")
     generated_complexes_ecoli_complexList_extended = convert_MITAB_to_complexList(
 iRef_binary, canonical_rep="yes", include_generated_complexes="yes", list_methods=
 "default", bait_use="include_non_baits")
     all_complexes_ecoli_complexList_extended = convert_MITAB_to_complexList(
 irefindex_curr_ecoli, canonical_rep="yes", include_generated_complexes="yes",
 list_methods="default", bait_use="include_non_baits")
     number_unique_complexes_all_extended = length(unique(
 all_complexes_ecoli_complexList_extended[,2]))
     }
     }
