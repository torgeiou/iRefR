    % File /man/select_protein.Rd
     \name{select_protein}
     \alias{select_protein}
     \title{Get all iRefIndex records for a given Protein}
     \description{
       Get all iRefIndex records for a given protein.
     }
     \usage{
     select_protein(id_type, id_value, MITAB_table, complex_info)
     }
     \arguments{
       \item{id_type}{Protein ID type ("irog" for iRefIndex's protein iROG ID; "icrog" for the canonical version of the iRefIndex's protein iROG ID).}
       \item{id_value}{A value or a vector of values of the iROG or icROG to search for in a given iRefIndex/MITAB table.}
       \item{MITAB_table}{iRefIndex/MITAB R table.}
       \item{complex_info}{Flag that gives the user the option to retrieve either all the records of the complex to which such protein belongs to ("full_complex"), or only the record of the complex where the given protein is present ("not_full_complex"). Default="full_complex".}
     }
     \value{
     \item{output}{iRefIndex/MITAB R table.}
     }
     \author{Antonio Mora <a.m.ortiz@biotek.uio.no>}
     \examples{
     ## get tables
     irefindex_80_mouse = get_irefindex("10090", "8.0", tempdir())

     ## single queries
     output_1 = select_protein("irogid", "4374882", irefindex_80_mouse, "not_full_complex")
     output_2 = select_protein("icrogid", "4374882", irefindex_80_mouse, "not_full_complex")
     output_3 = select_protein("icrogid", "2892004", irefindex_80_mouse, "not_full_complex")
     output_4 = select_protein("icrogid", "2892004", irefindex_80_mouse, "full_complex")
     #output_5 = select_protein("icrogid", "1365685", irefindex_80_mouse, "not_full_complex")
     #output_6 = select_protein("icrogid", "1365685", irefindex_80_mouse, "full_complex")

     ## getting list of tables
     irog_list = c("1828087", "2892004", "1365685")
     table_list = list()
     for (i in irog_list) {
     	table_list[[i]] = select_protein("irogid", i, irefindex_80_mouse, "not_full_complex")
     }

     ## getting single table
     irog_list = c("1828087", "2892004", "1365685")
     table_single = select_protein("irogid", irog_list, irefindex_80_mouse, "not_full_complex")
     }
