    % File /man/create_id_conversion_table.Rd
     \name{create_id_conversion_table}
     \alias{create_id_conversion_table}
     \title{Create a Protein ID Conversion Table}
     \description{
       Use iRefIndex to generate a lookup table to convert protein IDs between iRefIndex's iROG IDs, icROG IDs, and other protein identifiers such as Gene IDs, RefSeq IDs, PDBs and UniProt IDs.

On iROGs and icROGs:
"iRefIndex" guarantees the non-redundancy of protein information by assigning a different protein identifier (called Redundant Object Group, ROG) to every different protein sequence. This is called the "non-canonical" representation. At the same time, groups of non-redundant proteins might be different isoforms of a given protein, and, in this case, the identifier of one protein of the group (called cROG or canonical ROG) is chosen to represent the entire group of similar proteins. This is called the "canonical" representation of proteins. The iROG and icROG IDs correspond to integer representations of the ROG and cROG, respectively. Most of the "iRefR" functions work with iROGs and icROGs.
     }
     \usage{
     create_id_conversion_table(MITAB_table, data_folder, output_filename, IDs_to_include)
     }
     \arguments{
       \item{MITAB_table}{iRefIndex/MITAB R table.}
       \item{data_folder}{Folder to save the text file: type "data" to save it in the "iRefR/data" directory, "home" to save it in the "R.home()" directory, or any other destination folder. Default = getwd().}
       \item{output_filename}{File name for the conversion table files. Extensions ".txt" and ".RData" would be automatically added. Default="id_conversion_table"}
       \item{IDs_to_include}{Vector with the protein identifiers you want to convert from/to. iRefIndex iROGs and canonical iROGs are always included in the table. In addition, currently supports: UniProt ("uniprotkb"), RefSeq ("refseq"), GeneID ("entrezgene/locuslink"), Protein Data Bank ("PDB"), c("uniprotkb", "refseq", "entrezgene/locuslink") ("default") and all of them ("all"). Default = "default".}
     }
     \value{
     \item{output}{Lookup table including iROG ID, icROG ID, and the other specified protein IDs.}
     }
     \author{Antonio Mora <antoniocmora@gmail.com>}
     \examples{
     ## get table
     irefindex_curr_ecoli = get_irefindex("562", "13.0", tempdir())

     ## execute function
     id_conversion_table_ecoli = create_id_conversion_table(irefindex_curr_ecoli,
 tempdir(), "id_conversion_table_562_a")
     }
