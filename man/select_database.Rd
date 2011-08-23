    % File /man/select_database.Rd
     \name{select_database}
     \alias{select_database}
     \title{Select iRefIndex Records belonging to given Primary Interaction Databases}
     \description{
       Select the subset of records in iRefIndex belonging to one or a group of Primary Interaction Databases. In addition, select all records except for that database or those databases.
     }
     \usage{
     select_database(database, MITAB_table, flag)
     }
     \arguments{
       \item{database}{Name of the protein interaction DB or vector of databases to extract from iRefIndex (options: "intact", "bind", "grid", "mint", "ophid", "MPPI", "HPRD", "dip", "CORUM", "MPACT").}
       \item{MITAB_table}{iRefIndex/MITAB R table.}
       \item{flag}{Either "this_database" to get all "database" records, or "except_this_database" to get all records except for the ones belonging to "database". Default="this_database"}
     }
     \value{
     \item{output}{iRefIndex/MITAB R table with the previously specified features.}
     }
     \author{Antonio Mora <a.m.ortiz@biotek.uio.no>}
     \examples{
     ## get tables
     irefindex_curr_ecoli = get_irefindex("562", "current", tempdir())
     iRef_binary = select_interaction_type("binary", irefindex_curr_ecoli)
     iRef_complex = select_interaction_type("complex", irefindex_curr_ecoli)

     ## execute function
     binary_INTACT = select_database("intact", iRef_binary, "this_database")
     binary_non_INTACT = select_database("intact", iRef_binary, "except_this_database")
     complex_INTACT_mint = select_database(c("intact", "mint"), iRef_complex, "this_database")
     complex_non_INTACT_mint = select_database(c("intact", "mint"), iRef_complex, "except_this_database")
     }