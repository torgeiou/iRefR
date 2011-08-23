    % File /man/select_interaction_type.Rd
     \name{select_interaction_type}
     \alias{select_interaction_type}
     \title{Select Binary||Polymer||Complex in iRefIndex}
     \description{
       Select the subset of records in iRefIndex belonging to a certain interaction type (binary, polymer, complex).
     }
     \usage{
     select_interaction_type(int_type, MITAB_table)
     }
     \arguments{
       \item{int_type}{Interaction type ("binary", "complex" or "polymer").}
       \item{MITAB_table}{iRefIndex/MITAB R table.}
     }
     \value{
     \item{output}{R Table containing the requested subset of iRefIndex.}
     }
     \author{Antonio Mora <a.m.ortiz@biotek.uio.no>}
     \examples{
     ## get tables
     irefindex_curr_ecoli = get_irefindex("562", "current", tempdir())

     ## execute function
     iRef_binary = select_interaction_type("binary", irefindex_curr_ecoli)
     iRef_complex = select_interaction_type("complex", irefindex_curr_ecoli)
     iRef_polymer = select_interaction_type("polymer", irefindex_curr_ecoli)
     
     ## get results
     #pie_int_type = pie(c(dim(iRef_binary)[1], dim(iRef_polymer)[1], dim(iRef_complex)[1]),labels=c("Binaries","Polymers","Complexes"), main="iRefIndex Records per Interaction Type")
     #edit(iRef_complex[1:10,])
     }
