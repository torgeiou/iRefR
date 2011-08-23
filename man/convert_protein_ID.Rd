    % File /man/convert_protein_ID.Rd
     \name{convert_protein_ID}
     \alias{convert_protein_ID}
     \title{Convert Protein Identifier}
     \description{
       Convert protein identifiers from one of the following formats: iROG ID, icROG ID, gene ID, RefSeq and UniProt, to any other protein identifier in the list.

On iROGs and icROGs:
"iRefIndex" guarantees the non-redundancy of protein information by assigning a different protein identifier (called Redundant Object Group, ROG) to every different protein sequence. This is called the "non-canonical" representation. At the same time, groups of non-redundant proteins might be different isoforms of a given protein, and, in this case, the identifier of one protein of the group (called cROG or canonical ROG) is chosen to represent the entire group of similar proteins. This is called the "canonical" representation of proteins. The iROG and icROG IDs correspond to integer representations of the ROG and cROG, respectively. Most of the "iRefR" functions work with iROGs and icROGs.
     }
     \usage{
     convert_protein_ID(input_id_type, input_id_value, output_id_type, id_conversion_table)
     }
     \arguments{
       \item{input_id_type}{Input Type of Protein Identifier: iRefIndex non-canonical ("irogid"), iRefIndex canonical ("icrogid"), UniProt ("uniprotkb"), RefSeq ("refseq"), and GeneID ("entrezgene/locuslink").}
       \item{input_id_value}{Input Protein Identifier (value).}
       \item{output_id_type}{Output Type of Protein Identifier: iRefIndex non-canonical ("irogid"), iRefIndex canonical ("icrogid"), UniProt ("uniprotkb"), RefSeq ("refseq"), and GeneID ("entrezgene/locuslink").}
       \item{id_conversion_table}{Variable storing the ID Conversion Table (generated using "create_id_conversion_table()").}
     }
     \value{
     \item{output_id_value}{Output Protein Identifier (value).}
     }
     \author{Antonio Mora <a.m.ortiz@biotek.uio.no>}
     \examples{
     ## get tables
     irefindex_80_ecoli = get_irefindex("562", "8.0", tempdir())
     id_conversion_table_ecoli = create_id_conversion_table(irefindex_80_ecoli, tempdir(), "id_conversion_table_562")

     ## converting from geneID to uniprotID
     uniprot_id_value = convert_protein_ID("entrezgene/locuslink", "2986728", "uniprotkb", id_conversion_table_ecoli)

     ## converting a list of iROG IDs (from version 8.0) to RefSeq IDs
     irog_list = c("9917588", "9917576", "9917577", "369724", "593426", "3033592")
     refseq_list = list()
     for (i in irog_list) {
     	refseq_list[[i]] = convert_protein_ID("irogid", i, "refseq", id_conversion_table_ecoli)
     }
     }
