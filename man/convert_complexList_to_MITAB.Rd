    % File /man/convert_complexList_to_MITAB.Rd
     \name{convert_complexList_to_MITAB}
     \alias{convert_complexList_to_MITAB}
     \title{Convert complexList format to MITAB format}
     \description{
       Convert a table from complexList format to MITAB format.

The "complexList" format is one of the three table formats used by "iRefR", together with the MITAB format and the edgeList. The "complexList" format is a table with two columns, where the first column corresponds to an identifier and the second column to a group of proteins, commonly a complex, written as a comma-separated string. This representation has less information than MITAB but simplifies the visualization of groups of proteins, such as complexes.

On iROGs and icROGs:
"iRefIndex" guarantees the non-redundancy of protein information by assigning a different protein identifier (called Redundant Object Group, ROG) to every different protein sequence. This is called the "non-canonical" representation. At the same time, groups of non-redundant proteins might be different isoforms of a given protein, and, in this case, the identifier of one protein of the group (called cROG or canonical ROG) is chosen to represent the entire group of similar proteins. This is called the "canonical" representation of proteins. The iROG and icROG IDs correspond to integer representations of the ROG and cROG, respectively. Most of the "iRefR" functions work with iROGs and icROGs.
     }
     \usage{
     convert_complexList_to_MITAB(complexList, MITAB_table, canonical_rep,
 include_generated_complexes)
     }
     \arguments{
       \item{complexList}{iRefIndex/complexList R table.}
       \item{MITAB_table}{Full iRefIndex/MITAB R table for a given organism.}
       \item{canonical_rep}{Either "yes" for working with the canonical representation (icROGs) or "no" for working with the non-canonical representation (iROGs). Default="yes".}
       \item{include_generated_complexes}{Parameter to include the set of complexes generated from binary-represented interactions that might be misrepresented complexes, according to an "iRefR" algorithm that searches for groups of proteins coming from co-precipitation experiments, sharing the same PMID, source database and bait (when baits are available). Such set of complexes might be in the complexList only if the user generated them using "convert_MITAB_to_complexList", and can be easily identified by their complex ID constructed from their source database, PMID and experimental method names. If it is set to "yes", the user will get the original binary interactions these complexes were generated from. Default="no".}
     }
     \value{
     \item{MITAB_table}{iRefIndex/MITAB R table corresponding to the original complexList table.}
     }
     \author{Antonio Mora <antoniocmora@gmail.com>}
     \examples{
     \dontrun{
     ## get tables
     irefindex_curr_ecoli = get_irefindex("562", "current", tempdir())
     iRef_binary = select_interaction_type("binary", irefindex_curr_ecoli)
     iRef_complex = select_interaction_type("complex", irefindex_curr_ecoli)
     known_complexes_ecoli_complexList = convert_MITAB_to_complexList(iRef_complex,
 canonical_rep="yes", include_generated_complexes="no", list_methods="default")

     ## execute function
     reconstructed_known_ecoli_MITAB = convert_complexList_to_MITAB(
 known_complexes_ecoli_complexList, irefindex_curr_ecoli, "yes", "no")
     setequal(dim(iRef_complex), dim(reconstructed_known_ecoli_MITAB))

     ## more examples
     generated_complexes_ecoli_complexList = convert_MITAB_to_complexList(iRef_binary,
 canonical_rep="yes", include_generated_complexes="yes", list_methods="default")
     reconstructed_generated_ecoli_MITAB = convert_complexList_to_MITAB(
 generated_complexes_ecoli_complexList, irefindex_curr_ecoli, "yes", "yes")
     reconstructed_complexList = convert_MITAB_to_complexList(
 reconstructed_generated_ecoli_MITAB, canonical_rep="yes", include_generated_complexes=
 "yes")
     setequal(dim(generated_complexes_ecoli_complexList), dim(reconstructed_complexList))

     all_complexes_ecoli_complexList = convert_MITAB_to_complexList(irefindex_curr_ecoli,
 canonical_rep="yes", include_generated_complexes="yes", list_methods="default")
     reconstructed_all_ecoli_MITAB = convert_complexList_to_MITAB(
 all_complexes_ecoli_complexList, irefindex_curr_ecoli, "yes", "yes")
     reconstructed_all_complexList = convert_MITAB_to_complexList(
 reconstructed_all_ecoli_MITAB, canonical_rep="yes", include_generated_complexes="yes")
     setequal(dim(all_complexes_ecoli_complexList), dim(reconstructed_all_complexList))
     }
     }
