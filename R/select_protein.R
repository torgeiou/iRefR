####
# Select all iRefIndex records for a given protein:
####
select_protein = function(id_type, id_value, MITAB_table, complex_info) {
	MITAB_output = NULL
	for (i in id_value) {
		if (id_type == "irogid") {
			all_interactions_A = MITAB_table[which(MITAB_table$irogida == i),]
			all_interactions_B = MITAB_table[which(MITAB_table$irogidb == i),]
		}
		if (id_type == "icrogid") {
			all_interactions_A = MITAB_table[which(MITAB_table$icrogida == i),]
			all_interactions_B = MITAB_table[which(MITAB_table$icrogidb == i),]
		}
		all_interactions = unique(rbind(all_interactions_A, all_interactions_B))
		MITAB_output = rbind(MITAB_output, all_interactions)
	}

	if (complex_info == "not_full_complex") {
		MITAB_output = MITAB_output
	}
	if (complex_info == "full_complex") {
		tmp = select_interaction_type("complex", MITAB_output)
		irigid_list = unique(tmp$irigid)
		full_complexes = NULL
		for (i in irigid_list) {
			this_irigid = MITAB_table[which(MITAB_table$irigid == i),]
			full_complexes = rbind(full_complexes, this_irigid)
		}
		MITAB_output = unique(rbind(MITAB_output, full_complexes))
	}

	output = MITAB_output
}
