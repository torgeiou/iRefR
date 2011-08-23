    % File /man/get_irefindex.Rd
     \name{get_irefindex}
     \alias{get_irefindex}
     \title{Get iRefIndex}
     \description{
       Get a given distribution of the iRefIndex database from the ftp site.

"iRefIndex" is one of the main efforts to consolidate the information coming from different protein interaction databases in one single repository avoiding redundant information. It has been defined as "an index of protein interactions available in a number of primary interaction databases including BIND, BioGRID, CORUM, DIP, HPRD, IntAct, MINT, Mpact, MPPI and OPHID... This index allows the user to search for a protein and retrieve a non-redundant list of interactors for that protein".
     }
     \usage{
     get_irefindex(tax_id, iref_version, data_folder)
     }
     \arguments{
       \item{tax_id}{Taxon ID of the organism ("10090" for mouse, "10116" for rat, "4932" for S.Cerevisiae, "562" for E.Coli, "6239" for C.Elegans, "7227" for D.Melanogaster, "9606" for H.Sapiens, and "All" for all organisms). Default="All".}
       \item{iref_version}{iRefIndex distribution ("7.0", "8.0" or "current"). Default="current".}
       \item{data_folder}{Folder to save the iRefIndex flat file: type "data" to save it in the "iRefR/data" directory, "home" to save it in the "R.home()" directory, or any other destination folder. Default=getwd().}
     }
     \value{
     \item{irefindex_tab}{R Table (data frame) containing the desired iRefIndex distribution.}
     }
     \author{Antonio Mora <a.m.ortiz@biotek.uio.no>}
     \examples{
     ## execute function
     irefindex_80_ecoli = get_irefindex(tax_id="562", iref_version="8.0", data_folder=tempdir())
     }
