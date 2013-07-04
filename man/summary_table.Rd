    % File /man/summary_table.Rd
     \name{summary_table}
     \alias{summary_table}
     \title{Get Summary Information for a MITAB Table}
     \description{
       Get statistical information about some of the columns of a given iRefIndex/MITAB table.
     }
     \usage{
     summary_table(MITAB_table)
     }
     \arguments{
       \item{MITAB_table}{MITAB/iRefIndex R table to summarize.}
     }
     \value{
     \item{output}{R list containing: Number of non-canonical interactions, number of canonical interactions, number of non-canonical proteins, number of canonical proteins, number of publications, number of experimental methods, distribution of source databases, distribution of interaction types, distribution of number of participants.}
     }
     \author{Antonio Mora <antoniocmora@gmail.com>}
     \examples{
     ## get tables
     irefindex_80_ecoli = get_irefindex("562", "8.0", tempdir())

     ## execute function
     general_ecoli_statistics = summary_table(irefindex_80_ecoli)
     intact_statistics = summary_table(select_database("intact", irefindex_80_ecoli,
 "this_database"))
     complexes_statistics = summary_table(select_interaction_type("complex",
 irefindex_80_ecoli))
     low_thruput_statistics = summary_table(select_confidence("lpr", c(1,2),
 irefindex_80_ecoli))
     }
