    % File /man/summary_protein.Rd
     \name{summary_protein}
     \alias{summary_protein}
     \title{Get Summary Information for a given Protein}
     \description{
       Get information on interactors and complexes, for a given protein.
     }
     \usage{
     summary_protein(id_type, id_value, MITAB_table)
     }
     \arguments{
       \item{id_type}{Protein Identifier type ("irogid" for iRefIndex's protein iROG ID; "icrogid" for the canonical version of the iRefIndex's protein iROG ID).}
       \item{id_value}{Value of the iROG or icROG to summarize.}
       \item{MITAB_table}{MITAB/iRefIndex R table.}
     }
     \value{
     \item{output}{R list containing: Number of interactions per interaction type, interactors in binary interactions, interactors in complex interactions and interactors in polymer interactions.}
     }
     \author{Antonio Mora <a.m.ortiz@biotek.uio.no>}
     \examples{
     \dontrun{
     ## get tables
     irefindex_80_mouse = get_irefindex("10090", "8.0", tempdir())

     ## execute function
     output_1 = summary_protein("irogid", "4374882", irefindex_80_mouse)
     output_2 = summary_protein("icrogid", "4374882", irefindex_80_mouse)
     output_3 = summary_protein("icrogid", "2892004", irefindex_80_mouse)
     output_4 = summary_protein("icrogid", "1365685", irefindex_80_mouse)
     }
     }
