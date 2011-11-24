####
# Get the desired distribution of iRefIndex from the ftp site:
####
get_irefindex = function(tax_id="All", iref_version="current", data_folder=getwd()) {

	# 1. Get data folder:
	if (data_folder == "data") {
		datafolder = system.file("data", package = "iRefR")
	} else if (data_folder == "home") {
		datafolder = R.home()
	} else {
		datafolder = data_folder
	}

	# 2. Get release dates and full URLs:
	if (iref_version == "current") {
		iref_version = "9.0"
		release_date = "10182011"
		url = paste("ftp://ftp.no.embnet.org/irefindex/data/current/psimi_tab/MITAB2.6/", tax_id, ".mitab.", release_date, ".txt.zip", sep="")
	} else {
		if (iref_version == "7.0") {
			release_date = "11042010"
		}
		if (iref_version == "8.0") {
			release_date = "01192011"
		}
		if (iref_version == "9.0") {
			release_date = "10182011"
		}
		url = paste("ftp://ftp.no.embnet.org/irefindex/data/archive/release_", iref_version, "/psimi_tab/MITAB2.6/", tax_id, ".mitab.", release_date, ".txt.zip", sep="")
	}

	# 3. Check if file already exists. Otherwise, download and save:
	file_location = paste(datafolder, "/", tax_id,".mitab.",release_date,".txt",sep="")
	if (file.exists(file_location) == TRUE) {
		cat("Reading available iRefIndex file...\n")
		irefindex_tab = unique(read.table(file_location, header=TRUE, comment.char="", sep='\t', quote=""))
	} else {
		cat("Downloading iRefIndex file...\n")
		zipfile = paste(datafolder, "/", tax_id, ".mitab.", release_date, ".txt.zip", sep="")
		download.file(url, destfile=zipfile)
		unzip(zipfile, exdir=datafolder)
		file.remove(zipfile)
		cat("Reading downloaded file...\n")
		txtfile = paste(datafolder, "/", tax_id,".mitab.",release_date,".txt", sep="")
		irefindex_tab = unique(read.table(txtfile, header=TRUE, comment.char="", sep='\t', quote=""))
		cat("File has been saved as:\n")
		cat(paste(txtfile, "\n"))
		save(file = paste(datafolder, "/", tax_id,".mitab.",release_date,".RData",sep=""), list = "irefindex_tab")
	}

	irefindex_tab
}
