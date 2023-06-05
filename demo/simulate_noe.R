# names of PDB files in ensemble
pdb_names <- "ensemble.pdb"

# exclude groups with multiple atoms
exclude_multiatom_groups <- FALSE

# distance cutoff in Angstroms
dist_cutoff <- 5

if (! "ensemble_coord" %in% ls()) {
	ensemble_coord <- ke::read_ensemble(pdb_names, proton_only=FALSE)[,,,drop=FALSE]
}

# groupings
grouping_list <- list(
	g1=list(seq_len(dim(ensemble_coord)[3])),
	g2=as.list(seq_len(dim(ensemble_coord)[3]))
)

proton_groups <- ke::coord_proton_groups(ensemble_coord[,,1])

if (exclude_multiatom_groups) {

	excluded_groups <- unique(proton_groups[duplicated(proton_groups)])
	
	cat("Excluding protons bound to these heavy atoms:", sep="\n")
	print(excluded_groups)
	
	proton_groups <- proton_groups[!proton_groups %in% excluded_groups]
}

# create list with each element giving the protons in the group
group_list <- as.list(tapply(names(proton_groups), unname(proton_groups), c))
# reorder group list to match input order
group_list <- group_list[unique(proton_groups)]

# get mean coordinates of the atoms in each group
group_mean_coord <- array(NA_real_, dim=c(3, length(group_list), dim(ensemble_coord)[3]), dimnames=list(NULL, names(group_list), dimnames(ensemble_coord)[[3]]))
for (i in seq_along(group_list)) {
	group_mean_coord[,i,] <- apply(ensemble_coord[,group_list[[i]],,drop=FALSE], c(1,3), mean)
}

# calculate a distance matrix for every member of the ensemble (stored in a 3D array)
dist_array <- apply(group_mean_coord, 3, function(x) as.matrix(dist(t(x))))
dim(dist_array) <- c(dim(group_mean_coord)[2], dim(group_mean_coord)[2], dim(group_mean_coord)[3])
dimnames(dist_array) <- list(dimnames(group_mean_coord)[[2]], dimnames(group_mean_coord)[[2]], NULL)

# average over ensemble members to get average distances in the ensemble
mean_dist_mat <- rowMeans(dist_array, dims=2)

# determine set of atoms that are within the distance cutoff
eval_mat <- mean_dist_mat < dist_cutoff & upper.tri(mean_dist_mat)
eval_idx <- which(eval_mat, arr.ind=TRUE)
group_pairs <- matrix(names(group_list)[eval_idx], ncol=2)

# calculate array of internuclear vectors
group_r_array <- ke::coord_array_to_r_array(aperm(group_mean_coord, c(2,1,3)), group_pairs)

# calculate array if dipole-dipole interaction tensor vectors
group_d_array <- ke::r_array_to_d_array(group_r_array)

# calculat matrix of radii
group_r_mat <- sqrt(rowSums(group_r_array^2, dims=2))
# alternate approach directly from d_array
group_r_mat_alt <- rowSums(group_d_array^2, dims=2)^(-1/6)

# calculate r to the minus 3 and 6 averages
group_rm3_mean <- rowMeans(group_r_mat^-3)^(-1/3)
group_rm6_mean <- rowMeans(group_r_mat^-6)^(-1/6)

# calculate matrix of radii ranges
group_r_range <- t(apply(group_r_mat, 1, range))
colnames(group_r_range) <- c("rmin", "rmax")

# calculate matrix of radii standard deviations
group_r_sd <- apply(group_r_mat, 1, sd)

# calculate order parameter using normalized array of internuclear vectors
group_rnorm_array <- group_r_array/as.vector(group_r_mat)
group_dnorm_array <- ke::r_array_to_d_array(group_rnorm_array)
group_s2 <- rowSums(colMeans(aperm(group_dnorm_array, c(2,1,3)))^2)
# alternate approach directly from d_array
group_dnorm_array_alt <- group_d_array/as.vector(sqrt(rowSums(group_d_array^2, dims=2)))
group_s2_alt <- rowSums(colMeans(aperm(group_dnorm_array_alt, c(2,1,3)))^2)

# determine which pairs represent interactions between singletons (single atom per group)
singleton_idx <- sapply(group_list, length) == 1
singleton_pair_idx <- singleton_idx[group_pairs[,1]] & singleton_idx[group_pairs[,2]]
singleton_atom_pairs <- matrix(unlist(group_list[group_pairs[singleton_pair_idx,]]), ncol=2)
colnames(singleton_atom_pairs) <- paste0("atom", 1:2)
singleton_idx_pairs <- matrix(match(singleton_atom_pairs, dimnames(ensemble_coord)[[2]]), ncol=2)
colnames(singleton_idx_pairs) <- paste0("i", 1:2)

# calculate g values for singletons
singleton_r_array <- ke::coord_array_to_r_array(aperm(ensemble_coord, c(2,1,3)), singleton_atom_pairs)
singleton_d_array <- ke::r_array_to_d_array(singleton_r_array)
singleton_g_mat <- sapply(grouping_list, function(grouping) ke::d_array_to_g(singleton_d_array, grouping))
colnames(singleton_g_mat) <- names(grouping_list)

# create data frame with singleton data
singleton_data <- data.frame(
	singleton_atom_pairs,
	singleton_idx_pairs,
	singleton_g_mat,
	s2=group_s2[singleton_pair_idx],
	rm3=group_rm3_mean[singleton_pair_idx],
	rm6=group_rm6_mean[singleton_pair_idx],
	group_r_range[singleton_pair_idx,],
	rsd=group_r_sd[singleton_pair_idx]
)

# write out data
write.csv(singleton_data, "singleton_data.csv")
