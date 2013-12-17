    % File /man/select_confidence.Rd
     \name{select_confidence}
     \alias{select_confidence}
     \title{Select the subset of iRefIndex corresponding to a set of confidence scores}
     \description{
       Select the subset of records of iRefIndex corresponding to a given set of confidence scores.

The lpr score (lowest pmid re-use) is the lowest number of distinct interactions (irigids) that any PMID (supporting the interaction in this row) is used to support. A value of one indicates that at least one of the PMIDs supporting this interaction has never been used to support any other interaction. This likely indicates that only one interaction was described by that reference and that the present interaction is not derived from high throughput methods.

The hpr score (highest pmid re-use) is the highest number of interactions that any PMID (supporting the interaction in this row) is used to support. A high value (e.g. greater than 50) indicates that one PMID describes at least 50 other interactions and it is more likely that high-throughput methods were used.

The np score (number pmids) is the total number of unique PMIDs used to support the interaction described in this row. 
     }
     \usage{
     select_confidence(confidence_type, confidence_value, MITAB_table)
     }
     \arguments{
       \item{confidence_type}{Name of the confidence score for an interaction record. Options: lowest pmid re-use ("lpr"), highest pmid re-use ("hpr"), number pmids ("np"). Default: "lpr".}
       \item{confidence_value}{Vector with the values of the score whose interaction records we want to retrieve.}
       \item{MITAB_table}{iRefIndex/MITAB R table.}
     }
     \value{
     \item{MITAB_output}{iRefIndex/MITAB R table with the previously specified features.}
     }
     \author{Antonio Mora <antoniocmora@gmail.com>}
     \examples{
     ## get tables
     irefindex_curr_ecoli = get_irefindex("562", "13.0", tempdir())
     iRef_complex = select_interaction_type("complex", irefindex_curr_ecoli)

     ## execute function
     high_confidence_complexes = select_confidence("lpr", c(1:3, 10:15), iRef_complex)
     }
