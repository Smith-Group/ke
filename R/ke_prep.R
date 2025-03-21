#### Permutations ####

#' Find permutations to coordinates of XH3 groups
#'
#' @param atomids character vector of atom IDs (PDB columns 13-27)
#'
#' @return A character matrix with columns corresponding to a different permutation. 
#'    Columns correspond to equivalent atoms whose coordinates should be swapped.
#'
#' @examples
#' pdb2lum <- read_ensemble("https://files.rcsb.org/download/2LUM.pdb", proton_only=TRUE)
#' find_methyl_permutations(dimnames(pdb2lum)[[2]])
find_methyl_permutations <- function(atomids) {

	pos_3 <- gregexpr("3", substr(atomids, 1, 4))
	pos_3 <- lapply(seq_along(pos_3), function(i) cbind(i, pos_3[[i]])[,pos_3[[i]] != -1,drop=FALSE])
	pos_3 <- pos_3[sapply(pos_3, length) > 0]
	pos_3 <- do.call(rbind, pos_3)
	
	atom3 <- atomids[pos_3[,1]]
	
	atom1 <- sapply(seq_len(nrow(pos_3)), function(i) { 
		x <- atomids[pos_3[i,1]]
		substr(x, pos_3[i,2], pos_3[i,2]) <- "1"
		x
	})
	
	atom2 <- sapply(seq_len(nrow(pos_3)), function(i) { 
		x <- atomids[pos_3[i,1]]
		substr(x, pos_3[i,2], pos_3[i,2]) <- "2"
		x
	})
	methyl_idx <- atom1 %in% atomids & atom2 %in% atomids
	
	methyl_array <- rbind(atom1, atom2, atom3, atom2, atom3, atom1, atom3, atom1, atom2)[,methyl_idx]
	dim(methyl_array) <- c(3, 3, sum(methyl_idx))
	dimnames(methyl_array) <- list(NULL, paste("methyl", 1:3, sep=""), NULL)
	
	methyl_list <- lapply(seq_len(dim(methyl_array)[3]), function(i) methyl_array[,,i])
	
	for (i in seq_along(methyl_list)) {
		colnames(methyl_list[[i]]) <- paste(trimws(substr(methyl_list[[i]][1,], 10, 14)), trimws(substr(methyl_list[[i]][1,], 1, 4)), sep=" ")
	}
	
	methyl_list
}

#' Find permutations to coordinates of phenylalanine and tyrosine rings
#'
#' @param atomids character vector of atom IDs (PDB columns 13-27)
#'
#' @return A character matrix with columns corresponding to a different permutation. 
#'    Columns correspond to equivalent atoms whose coordinates should be swapped.
#'
#' @examples
#' pdb2lum <- read_ensemble("https://files.rcsb.org/download/2LUM.pdb", proton_only=TRUE)
#' find_aromatic_permutations(dimnames(pdb2lum)[[2]])
find_aromatic_permutations <- function(atomids) {

	tyr_phe_hd1 <- grep("( HD1 PHE| HD1 TYR)", atomids, value=TRUE)
	tyr_phe_he1 <- sub(" HD1", " HE1", tyr_phe_hd1)
	tyr_phe_hd2 <- sub(" HD1", " HD2", tyr_phe_hd1)
	tyr_phe_he2 <- sub(" HD1", " HE2", tyr_phe_hd1)
	
	aromatic_array <- rbind(tyr_phe_hd1, tyr_phe_he1, tyr_phe_hd2, tyr_phe_he2, tyr_phe_hd2, tyr_phe_he2, tyr_phe_hd1, tyr_phe_he1)
	dimnames(aromatic_array) <- NULL
	dim(aromatic_array) <- c(4, 2, length(tyr_phe_hd1))
	
	aromatic_list <- lapply(seq_len(dim(aromatic_array)[3]), function(i) aromatic_array[,,i])
	
	for (i in seq_along(aromatic_list)) {
		colnames(aromatic_list[[i]]) <- paste(trimws(substr(aromatic_list[[i]][1,], 10, 14)), c("A1", "A2"), sep=" ")
	}
	
	aromatic_list
}

#' Create list of matrices for doing two or three atom permutations
#'
#' @param atomids character vector of atom IDs (PDB columns 13-27)
#' @param atom_permutations list of permutations from `find_methyl_permutations()` and 
#'    `find_aromatic_permutations()`
#'
#' @examples
#' pdb2lum <- read_ensemble("https://files.rcsb.org/download/2LUM.pdb", proton_only=TRUE)
#' perm_methyl <- find_methyl_permutations(dimnames(pdb2lum)[[2]])
#' perm_aro <- find_aromatic_permutations(dimnames(pdb2lum)[[2]])
#' perm_list <- make_atom_perm_list(dimnames(pdb2lum)[[2]], c(perm_methyl, perm_aro))
#' perm_list[[3]][unlist(lapply(perm_methyl, "[", , 1)),]
#' perm_list[[2]][unlist(lapply(perm_aro, "[", , 1)),]
make_atom_perm_list <- function(atomids, atom_permutations) {

	perm_lengths <- sapply(atom_permutations, ncol)
	
	perm_mat_list <- vector("list", max(perm_lengths))
	
	for (i in seq_along(perm_mat_list)) {
	
		perm_mat_list[[i]] <- matrix(seq_along(atomids), nrow=length(atomids), ncol=i, dimnames=list(atomids, NULL))
		atom_perm_idx <- which(perm_lengths == i)
		for (j in atom_perm_idx) {
			perm_mat_list[[i]][atom_permutations[[j]][,1],] <- match(atom_permutations[[j]], atomids)
		}
	}
	
	perm_mat_list
}

#' Generate matrix of possible permuted assignments to individual atoms in a group
#'
#' @param permutation single permutation matrix
#'
#' This is useful for calculating interactions between a permutation group and other atoms
#'
#' @examples
#' pdb2lum <- read_ensemble("https://files.rcsb.org/download/2LUM.pdb", proton_only=TRUE)
#' perm_methyl <- find_methyl_permutations(dimnames(pdb2lum)[[2]])
#' perm_aro <- find_aromatic_permutations(dimnames(pdb2lum)[[2]])
#' unique_atom_map(perm_methyl[[1]])
#' unique_atom_map(perm_aro[[1]])
unique_atom_map <- function(permutation) {

	#print(permutation)
	
	sorted_permutation <- t(apply(permutation, 1, sort))
	#print(sorted_permutation)
	
	unique_sorted_permutation <- unique(sorted_permutation)
	#print(unique_sorted_permutation)
	
	permutation_idx <- match(apply(sorted_permutation, 1, paste, collapse=""), apply(unique_sorted_permutation, 1, paste, collapse=""))
	#print(permutation_idx)

	stopifnot(tapply(permutation_idx, permutation_idx, length) == ncol(permutation))

	rownames(unique_sorted_permutation) <- unique_sorted_permutation[,1]

	unique_sorted_permutation
}

#' Generate matrix of possible permuted assignments of atom pairs within a group
#'
#' @param permutation single permutation matrix
#'
#' This is useful for calculating interactions within a permutation group
#' @examples
#' pdb2lum <- read_ensemble("https://files.rcsb.org/download/2LUM.pdb", proton_only=TRUE)
#' perm_methyl <- find_methyl_permutations(dimnames(pdb2lum)[[2]])
#' perm_aro <- find_aromatic_permutations(dimnames(pdb2lum)[[2]])
#' unique_atom_pair_map(perm_methyl[[1]])
#' unique_atom_pair_map(perm_aro[[1]])
unique_atom_pair_map <- function(permutation) {

	#print(permutation)
	
	unique_atoms <- sort(unique(as.vector(permutation)))
	#print(unique_atoms)
	
	unique_atom_pairs <- utils::combn(unique_atoms, 2)
	#print(unique_atom_pairs)
	
	match_mat <- apply(unique_atom_pairs, 2, function(x) {
		match_mat <- apply(permutation, 2, function(y) sort(match(x, y)))
		match_mat[,order(match_mat[1,], match_mat[2,])]
	})
	
	match_char <- apply(match_mat, 2, paste, collapse=" ")
	
	unique_match_idx <- match(match_char, unique(match_char))
	
	lapply(seq_len(max(unique_match_idx)), function(idx) {
		set_atom_pairs <- unique_atom_pairs[,unique_match_idx == idx,drop=FALSE]
		rownames(set_atom_pairs) <- set_atom_pairs[,1]
		set_atom_pairs
	})
}

#### Transition Matrix Creation ####

#' Create transition rate matrix with single rate
#'
#' @param k desired eigenvalue for the rate matrix
#' @param n_names character vector with row and columns
#'
#' @examples
#' rate_mat <- rate_mat_simple(4, c("a", "b", "c", "d"))
#' rate_mat <- rate_mat_diag(rate_mat)
#' eigen(rate_mat)
rate_mat_simple <- function(k, n_names) {

	rate_mat <- matrix(k/(length(n_names)), nrow=length(n_names), ncol=length(n_names), dimnames=list(n_names, n_names))
	diag(rate_mat) <- NA
	
	rate_mat
}

#' Create hierarchical transition rate matrix
#'
#' @param group_vec integer vector giving group membership of each state
#' @param k_intra desired eigenvalue for intra-group transitions
#' @param k_inter desired eigenvalue for inter-group transitions
#'
#' @examples
#' rate_mat <- rate_mat_intra_inter(c(0,0,1,1), 10, 4)
#' rate_mat <- rate_mat_diag(rate_mat)
#' eigen(rate_mat)
rate_mat_intra_inter <- function(group_vec, k_intra, k_inter) {

	rate_mat <- matrix(k_inter/length(group_vec), nrow=length(group_vec), ncol=length(group_vec), dimnames=list(names(group_vec), names(group_vec)))
	
	group_idx <- lapply(unique(group_vec), function(x) which(group_vec == x))
	
	for (idx in group_idx) {
		#print((length(group_vec)-length(idx))/length(group_vec))
		rate_mat[idx,idx] <- (k_intra-k_inter*(length(group_vec)-length(idx))/length(group_vec))/length(idx)
	}
	
	diag(rate_mat) <- NA
	
	rate_mat
}

#' Calculate Kronecker product of two transition rate matrices
#'
#' @param trans_rate_a first transition rate matrix
#' @param trans_rate_b second transition rate matrix
#'
#' @examples
#' rate_mat_fast <- rate_mat_simple(10, c("f1", "f2"))
#' rate_mat_slow <- rate_mat_simple(4, c("s1", "s2"))
#' rate_mat <- rate_mat_kronecker(rate_mat_fast, rate_mat_slow)
#' rate_mat <- rate_mat_diag(rate_mat)
#' eigen(rate_mat)
rate_mat_kronecker <- function(rate_mat_a, rate_mat_b) {

	a_names <- rep(colnames(rate_mat_a), ncol(rate_mat_b))
	b_names <- rep(colnames(rate_mat_b), each=ncol(rate_mat_a))

	new_names <- paste(a_names, b_names, sep="_")
	
	diag(rate_mat_a) <- 0
	diag(rate_mat_b) <- 0
	
	#rate_mat_new <- kronecker(rate_mat_b/prod(dim(rate_mat_a)[1]), rate_mat_a/prod(dim(rate_mat_b)[1]), pmin)
	rate_mat_new <- kronecker(diag(nrow(rate_mat_b)), rate_mat_a) + kronecker(rate_mat_b, diag(nrow(rate_mat_a)))
	
	dimnames(rate_mat_new) <- list(new_names, new_names)
	diag(rate_mat_new) <- NA
	
	rate_mat_new
}

#' Update transition rate diagonal using rates
#'
#' @param trans_rate_mat square transition rate matrix
rate_mat_diag <- function(rate_mat) {

	for (i in seq_len(nrow(rate_mat))) {
		stopifnot(all(rate_mat[i,-i] >= 0))
		rate_mat[i,i] <- -sum(rate_mat[i,-i])
	}
	
	rate_mat
}

#### Transition Rate Analysis ####

#' Get counts enabling calculation of eigenvalues from individual rates
#'
#' @param trans_rate 
#' @param all_permutations 
#' @param eps_factor 
#'
#' @examples
#' rate_mat <- rate_mat_intra_inter(c(0,0,1,1), 10, 4)
#' get_rate_count_mat(rate_mat, NULL)
get_rate_count_mat <- function(trans_rate, all_permutations, eps_factor=0.5) {

	base_trans_rates <- -eigen(rate_mat_diag(trans_rate))$values[-1]
	
	base_trans_rate_groups <- get_rate_groups(base_trans_rates, eps_factor)
	
	base_rate_means <- sapply(base_trans_rate_groups, function(idx) mean(base_trans_rates[idx]))
	if (length(base_rate_means) == 0) {
		base_rate_means <- numeric()
	}
	
	all_perm_mat <- cbind(as.numeric(names(all_permutations)), sapply(all_permutations, ncol))
	
	unique_permutations <- all_permutations[!duplicated(all_perm_mat)]
	
	#print(unique_permutations)
	
	stopifnot(!duplicated(names(unique_permutations)))
	
	permutation_rates <- as.numeric(names(unique_permutations))
	
	#print(c(rep(list(0:1), length(base_rate_means)), rep(list(0:2), length(permutation_rates))))
	
	base_rate_counts <- as.matrix(expand.grid(rep(list(0:1), length(base_rate_means))))
	if (length(base_rate_counts)) {
		base_rate_rates <- colSums(t(base_rate_counts)*base_rate_means)
		base_rate_comp <- apply(abs(outer(base_rate_rates, base_trans_rates, "-")), 1, min)/min(base_trans_rates)
		base_rate_counts <- base_rate_counts[base_rate_rates == 0 | base_rate_comp < eps_factor,,drop=FALSE]
	}
	
	permutation_rate_counts <- as.matrix(expand.grid(rep(list(0:2), length(permutation_rates))))
	permutation_rate_counts <- permutation_rate_counts[rowSums(permutation_rate_counts) <= 2,,drop=FALSE]
	
	if (nrow(base_rate_counts)) {
		
		if (nrow(permutation_rate_counts)) {
		
			rate_count_exp <- as.matrix(expand.grid(seq_len(nrow(base_rate_counts)), seq_len(nrow(permutation_rate_counts))))
	
			rate_count_mat <- unname(t(cbind(base_rate_counts[rate_count_exp[,1],], permutation_rate_counts[rate_count_exp[,2],])))
		
		} else {
		
			rate_count_mat <- unname(t(cbind(base_rate_counts)))
		}
	
	} else {
	
		rate_count_mat <- unname(t(cbind(permutation_rate_counts)))
	}
	
	#print(rate_count_mat)
	
	#rate_count_mat <- unname(t(as.matrix(expand.grid(c(
	#	rep(list(0:1), length(base_rate_means)), 
	#	rep(list(0:2), length(permutation_rates))
	#)))))
	#rate_count_mat <- rate_count_mat[,colSums(rate_count_mat[-seq_along(base_rate_means),]) <= 2]
	
	#print(rate_count_mat)
	
	rownames(rate_count_mat) <- c(base_rate_means, permutation_rates)
	
	all_rates <- colSums(c(base_rate_means, permutation_rates)*rate_count_mat)
	
	colnames(rate_count_mat) <- all_rates
	
	rate_count_mat <- rate_count_mat[,order(all_rates),drop=FALSE]
	
	rate_count_mat
}

#' Find groups of numerically equivalent rates (degenerate eigenvalues)
#'
#' @param rates numeric vector of rates
#' @param eps_factor rates must differ by this amount times the lowest rate
#' @param eps_log10 log10(rates) must differ by this amount
#'
#' If set, eps_log10 is used and eps_factor is bypassed.
get_rate_groups <- function(rates, eps_factor=0.5, eps_log10=NULL) {

	if (length(rates) == 0) {
		return(list())
	}

	rate_ord <- order(rates)
	rates_ord <- rates[rate_ord]

	if (is.null(eps_log10)) {
		group_div <- which(diff(rates_ord) > min(rates_ord)*eps_factor)
	} else {
		rate_diffs <- diff(log10(rates_ord))
		group_div <- which(rate_diffs > eps_log10)
	}
	group_start <- c(1, group_div+1)
	group_end <- c(group_div, length(rates_ord))
	
	rate_groups <- vector("list", length(group_start))
	
	for (i in seq_along(group_start)) {
		rate_groups[[i]] <- rate_ord[seq(group_start[i], group_end[i])]
	}
	
	names(rate_groups) <- rates[sapply(rate_groups, "[", 1)]
	
	rate_groups
}

#' Determine subsets of states within transition rate matrices
#'
#' @trans_rate_eigen list with `values` and `vectors` as returned by `eigen()`
#' @rate_groups list with numerically equivalent rates organized into groups
#'
#' This only works with rate matrices returned by `rate_mat_simple()` and 
#' `rate_mat_intra_inter`. It does not work with `rate_mat_kronecker()`.
#'
#' @examples
#' rate_mat <- rate_mat_intra_inter(c(0,0,1,1), 10, 4)
#' rate_mat <- rate_mat_diag(rate_mat)
#' calc_subset_mat(eigen(rate_mat))
calc_subset_mat <- function(trans_rate_eigen, rate_groups=get_rate_groups(-trans_rate_eigen$values[-1])) {

	state_pop <- trans_rate_eigen$vectors[,1]^2
	
	trans_rate_prob <- tcrossprod(trans_rate_eigen$vectors[,1])*state_pop
	
	rate_eigenvectors <- trans_rate_eigen$vectors[,-1,drop=FALSE]

	subset_mat <- matrix(NA_real_, nrow=nrow(trans_rate_eigen$vectors), ncol=length(rate_groups), dimnames=list(NULL, names(rate_groups)))

	# Handle signed zero values (ran into this on Linux R with Intel compilers)
	acceptable_ratios <- round(c(-1e-16, 1e-16, 1))

	for (i in seq_along(rate_groups)) {
	
		idx <- rate_groups[[i]]
		unique_trans_rate_prob <- tcrossprod(rate_eigenvectors[,idx])*state_pop
		trans_rate_prob <- trans_rate_prob+unique_trans_rate_prob
		trans_rate_prob_ideal <- tcrossprod(sqrt(diag(trans_rate_prob)))
		trans_rate_prob_ratios <- round(trans_rate_prob/trans_rate_prob_ideal, 8)
		#print(trans_rate_prob_ratios)
		
		stopifnot(trans_rate_prob_ratios %in% acceptable_ratios)
		
		item_names <- apply(trans_rate_prob_ratios, 1, paste, collapse="")
		group_vector <- match(item_names, unique(item_names))
		group_list <- unname(tapply(seq_along(group_vector), group_vector, function(x) x, simplify=FALSE))
		
		for (j in seq_along(group_list)) {
			subset_mat[group_list[[j]],i] <- j
		}
	}
	
	subset_mat
}

#' Expand subset matrix consistent with applying Kronecker product to original rate matrix
#'
#' @param rate_subset_mat
#' @param permutation_counts
#'
#' @examples
#' rate_mat_fast <- rate_mat_simple(10, c("f1", "f2"))
#' rate_mat_fast <- rate_mat_diag(rate_mat_fast)
#' subset_mat <- calc_subset_mat(eigen(rate_mat_fast))
#' expand_subset_mat(subset_mat, list("4"=2))
expand_subset_mat <- function(rate_subset_mat, permutation_counts) {
	
	permutation_subset_mat <- as.matrix(expand.grid(lapply(permutation_counts, seq_len)))

	# swapped from ker to be compatible with array ordering having models on the outside
	#exp_grid <- expand.grid(seq_len(nrow(rate_subset_mat)), seq_len(nrow(permutation_subset_mat)))
	exp_grid <- expand.grid(seq_len(nrow(permutation_subset_mat)), seq_len(nrow(rate_subset_mat)))
	
	# swapped from ker to be compatible with array ordering having models on the outside
	#rate_subset_mat <- cbind(rate_subset_mat[exp_grid[,1],,drop=FALSE], permutation_subset_mat[exp_grid[,2],,drop=FALSE])
	rate_subset_mat <- cbind(rate_subset_mat[exp_grid[,2],,drop=FALSE], permutation_subset_mat[exp_grid[,1],,drop=FALSE])
	
	rownames(rate_subset_mat) <- NULL
	
	rate_subset_mat[,order(as.numeric(colnames(rate_subset_mat))),drop=FALSE]
}

#' Get `k` coefficients for creating linear combinations of `a` values from `g` values
#' 
#' @param trans_rate_eigen
#' @param all_rates 
#' @param subset_mat created by `calc_subset_mat()` or `expand_subset_mat()`
#' @param validate
#' @param eps_factor
#'
#' @examples
#' rate_mat <- rate_mat_intra_inter(c(0,0,1,1), 10, 4)
#' rate_mat <- rate_mat_diag(rate_mat)
#' get_eigen_groups(eigen(rate_mat))
get_eigen_groups <- function(trans_rate_eigen, all_rates=NULL, subset_mat=NULL, validate=FALSE, eps_factor=0.5) {
	
	state_pop <- trans_rate_eigen$vectors[,1]^2

	rates <- signif(-trans_rate_eigen$values[-1], 8)
	
	# get list of groupings based on structure of the eigenvectors
	rate_groups <- get_rate_groups(rates)
	
	# determine subset matrix from the eigenvectors if it isn't already provided
	if (is.null(subset_mat)) {	
		subset_mat <- calc_subset_mat(trans_rate_eigen, rate_groups)
	}
	
	# calculate a grid with different combinations of subsets
	subset_comb <- as.matrix(expand.grid(rep(list(c(FALSE,TRUE)), ncol(subset_mat))))
	if (ncol(subset_mat) == 0) {
		subset_comb <- matrix(nrow=1, ncol=0)
	}
	colnames(subset_comb) <- colnames(subset_mat)
	rownames(subset_comb) <- apply(subset_comb, 1, function(x) paste(which(x), collapse=","))
	if (nrow(subset_comb)) {
		rownames(subset_comb)[1] <- "0"
	}
	
	# vector of all unique rates, 
	unique_rates <- rates[unlist(sapply(rate_groups, "[", 1))]
	
	# vector of all unique rates, including 0 as the first element
	zero_unique_rates <- c(0, unique_rates)
	
	# determine the rates represented by each row of the subset grid
	subset_rates <- colSums(t(subset_comb)*as.numeric(colnames(subset_comb)))
	
	# determine mapping between subset grid rates and unique rates
	rate_subset_idx <- sapply(subset_rates, function(x) {
		abs_diffs <- abs(zero_unique_rates-x)
		#print(abs_diffs)
		min_idx <- which.min(abs_diffs)
		if (abs_diffs[min_idx] < min(unique_rates,1)*eps_factor) {
			min_idx
		} else {
			NA_integer_
		}
	})
	
	#print(subset_mat)
	#print(subset_rates)
	#print(zero_unique_rates)
	#print(rate_subset_idx)
	
	# remove subsets that don't map onto a unique rate
	subset_comb <- subset_comb[!is.na(rate_subset_idx),,drop=FALSE]
	rate_subset_idx <- rate_subset_idx[!is.na(rate_subset_idx)]
	
	# update rate grid to include implicit dependencies between grid rows
	for (i in seq_len(max(ncol(subset_comb)-1, 0))) {
		for (j in seq(i+1, ncol(subset_comb))) {
			if (!any(subset_comb[,i] & subset_comb[,j])) {
				subset_comb[,i] <- subset_comb[,i] | subset_comb[,j]
			}
		}
	}
	
	# determine coefficients that produce the eigenvector cross products from averaged subsets (the key algorithm)
	subset_coefficient_mat <- matrix(0L, nrow=nrow(subset_comb), ncol=nrow(subset_comb), dimnames=list(rownames(subset_comb), rownames(subset_comb)))
	subset_coefficient_mat[1,1] <- 1
	subset_precursor_mat <- matrix(FALSE, nrow=nrow(subset_comb), ncol=nrow(subset_comb), dimnames=list(rownames(subset_comb), rownames(subset_comb)))
	
	for (i in seq_len(ncol(subset_coefficient_mat))[-1]) {
		
		# precursors have one fewer true boolean
		precursor_num_true <- sum(subset_comb[i,])-1
		precursor_idx <- which(rowSums(subset_comb) == precursor_num_true & colSums(t(subset_comb) & subset_comb[i,]) == precursor_num_true)
		
		# update matrix with the history of precursors
		subset_precursor_mat[,i] <- apply(subset_precursor_mat[,precursor_idx,drop=FALSE], 1, any)
		subset_precursor_mat[precursor_idx,i] <- TRUE
		
		# set the coefficients equal to the current subset minus history of precursors
		subset_coefficient_mat[,i] <- -rowSums(subset_coefficient_mat[,subset_precursor_mat[,i],drop=FALSE])
		subset_coefficient_mat[i,i] <- 1
	}
	
	# collapse subset_coefficient_mat such that duplicate rates are summed
	subset_coefficient_mat <- t(apply(subset_coefficient_mat, 1, tapply, rate_subset_idx, sum))
	colnames(subset_coefficient_mat) <- zero_unique_rates
	
	# get a list of the subsets for later calculations
	subset_idx <- lapply(seq_len(nrow(subset_comb)), function(i) {
		row_chars <- apply(subset_mat[,subset_comb[i,],drop=FALSE], 1, paste, collapse=" ")
		unname(tapply(seq_along(row_chars), row_chars, "c", simplify=FALSE))
	})

	all_rate_idx <- if (is.null(all_rates)) NULL else sapply(zero_unique_rates, function(x) which.min(abs(all_rates-x)))
	stopifnot(!duplicated(all_rate_idx))
	
	if (validate) {
	
		# generate matrices representing the cross product of each subset
		subset_mats <- lapply(subset_idx, function(sset) {
			mat <- array(0, dim(trans_rate_eigen$vectors))
			for (idx in sset) {
				pop <- state_pop[idx]
				sqrt_norm_pop <- sqrt(pop/sum(pop))
				mat[idx,idx] <- sqrt_norm_pop %o% sqrt_norm_pop
			}
			mat/sum(mat)
		})
		
		rate_eigenvectors <- trans_rate_eigen$vectors[,-1,drop=FALSE]
	
		subset_mat_diff <- rep(NA_real_, length(unique_rates))
		
		for (i in seq_along(rate_groups)) {
		
			idx <- rate_groups[[i]]
			unique_trans_rate_prob <- tcrossprod(rate_eigenvectors[,idx])*state_pop
			subset_coefs <- subset_coefficient_mat[,i+1]
			subset_array <- simplify2array(lapply(which(subset_coefs != 0), function(i) subset_coefs[i]*subset_mats[[i]]))
			
			subset_mat_diff[i] <- sum(abs(unique_trans_rate_prob-rowSums(subset_array, dims=2)))
		}
		
		stopifnot(subset_mat_diff < 1e-6)
	}
	
	subset_group_invpop <- lapply(subset_idx, function(group_list) 1/sapply(group_list, function(idx) sum(state_pop[idx])))
	
	list(
		state_pop=state_pop, 
		unique_rates=unique_rates, 
		all_rates=all_rates, 
		all_rate_idx=all_rate_idx,
		subset_idx=subset_idx,
		subset_group_invpop=subset_group_invpop,
		subset_coef=subset_coefficient_mat,
		subset_mat=subset_mat
	)
}

#' Construct `k` coefficients in a recursive manner
#'
#' @param trans_rate transition rate matrix
#' @param parent_data previous data to modify assuming application of Kronecker product
#' @param all_rates vector of all rates
#' @param validate passed onto `get_eigen_groups()`
#'
#' @examples
#' rate_mat_fast <- rate_mat_simple(10, c("f1", "f2"))
#' rate_mat_slow <- rate_mat_simple(4, c("s1", "s2"))
#' rate_mat <- rate_mat_kronecker(rate_mat_fast, rate_mat_slow)
#' all_rates <- unique(as.numeric(colnames(get_rate_count_mat(rate_mat, NULL))))
#' rate_data_base <- get_rate_data(rate_mat_fast, rate_mat_slow)
#' rate_mat <- rate_mat_diag(rate_mat)
#' eigen(rate_mat)
get_rate_data <- function(trans_rate, parent_data=NULL, all_rates=NULL, validate=FALSE) {

	if (is.null(parent_data)) {
		data_name <- "base"
	} else {
		data_name <- paste(parent_data[["name"]], paste(1/trans_rate[!diag(nrow(trans_rate))], collapse=","), sep="/")
	}
	
	#cat(paste("Calculating ", data_name, sep=""), sep="\n")
	
	#print(system.time({
	
	subset_mat <- NULL
	
	if (!is.null(parent_data)) {
		permutation_counts <- list(nrow(trans_rate))
		names(permutation_counts) <- signif(-eigen(rate_mat_diag(trans_rate))$values[2])
		subset_mat <- expand_subset_mat(parent_data[["subset_mat"]], permutation_counts)
		#trans_rate <- rate_mat_kronecker(parent_data[["trans_rate"]], trans_rate)
		trans_rate <- rate_mat_kronecker(trans_rate, parent_data[["trans_rate"]])
	}
	
	trans_rate_eigen <- eigen(rate_mat_diag(trans_rate))
	
	eigen_groups <- get_eigen_groups(trans_rate_eigen, all_rates=all_rates, subset_mat=subset_mat, validate=validate)
	
	#}))
	
	rate_data <- c(
		list(
			name=data_name, 
			trans_rate=trans_rate#, 
			#eigen=trans_rate_eigen
		),
		eigen_groups
	)
	
	rate_data
}

#' Convert a list of subset groups into a matrix with group indices
#'
#' @param subset_idx list of lists of integers giving members of each group
#'
#' @return a matrix with row for each outer list and group assignment for each ensemble member
subset_idx_to_grouping_mat <- function(subset_idx) {

	group_mat <- matrix(NA_integer_, nrow = length(subset_idx), ncol = max(unlist(subset_idx)))
	
	for (i in seq_along(subset_idx)) {
		for (j in seq_along(subset_idx[[i]])) {
			group_mat[i,subset_idx[[i]][[j]]] <- j
		}	
	}
	
	group_mat
}

update_numeric_names <- function(trial_names, master_names) {

	trial_names_numeric <- as.numeric(trial_names)
	master_names_numeric <- as.numeric(master_names)
	
	idx <- sapply(trial_names_numeric, function(x) which.min(abs(x-master_names_numeric)))
	
	stopifnot(!duplicated(idx))
	
	master_names[idx]
}

#' Convert rate data into a format more easily serializable to set of tables
#'
#' @param rate_data list returned by `get_rate_data()`
#' @param rate_count_mat matrix returned by `get_rate_count_mat()` including methyl and
#'     aromatic rates
#' @param numeric vector whose names give the rates to be used
#'
#' @return a list with elements: `groupings`, `a_coef`, and `lambda_coef`
rate_data_to_mat_list <- function(rate_data, rate_count_mat, rates_named = NULL) {

	groupings <- subset_idx_to_grouping_mat(rate_data$subset_idx)
	a_coef <- rate_data$subset_coef
	# fix possible roundoff/numeric errors in eigenvalues
	colnames(a_coef) <- update_numeric_names(colnames(a_coef), colnames(rate_count_mat))
	lambda_coef <- rate_count_mat[ , colnames(a_coef), drop = FALSE]
	lambda_coef <- lambda_coef[apply(lambda_coef != 0, 1, any), , drop = FALSE]
	
	name_map <- names(rates_named)
	names(name_map) <- rates_named
	
	rownames(a_coef) <- NULL
	rownames(lambda_coef) <- unname(name_map[rownames(lambda_coef)])
	
	name_mat <- apply(lambda_coef, 2, paste0, rownames(lambda_coef))
	name_mat <- sub("0.+", "", name_mat)
	name_mat <- sub("^([0-9])", "+\\1", name_mat)
	name_mat <- sub("^[+]1", "+", name_mat)
	if (!is.matrix(name_mat)) {
		name_mat <- t(name_mat)
	}
	name_vec <- apply(name_mat, 2, function(x) paste(x[x != ""], collapse=""))
	name_vec <- sub("^[+]", "", name_vec)
	name_vec[1] <- "0"
	
	colnames(a_coef) <- unname(name_vec)
	colnames(lambda_coef) <- unname(name_vec)
	
	list(groupings = groupings, a_coef = a_coef, lambda_coef = lambda_coef)
}

#### Relaxation Matrix Setup ####

#' Create list of equivalent atoms
#'
#' @param atomids character vector of atom IDs (PDB columns 13-27)
#' @param atom_permutations list of permutations from `find_methyl_permutations()` and 
#'    `find_aromatic_permutations()`
#'
#' @return list of atomid vectors where all atomids in a given vector are equivalent
make_equiv_list <- function(atomids, permutations) {

	equiv_vec <- rep(NA_integer_, length(atomids))
	names(equiv_vec) <- atomids
	for (i in seq_along(permutations)) {
		perm_unique <- unique(t(apply(permutations[[i]], 1, sort)))
		for (j in seq_len(nrow(perm_unique))) {
			equiv_vec[perm_unique[j,]] <- i+j*length(permutations)
		}
	}
	equiv_vec[is.na(equiv_vec)] <- -seq_len(sum(is.na(equiv_vec)))
	equiv_vec[] <- match(equiv_vec, unique(equiv_vec))
	
	equiv_list <- lapply(unique(equiv_vec), function(i) names(equiv_vec)[equiv_vec == i])
	
	equiv_list
}

#' Update names of equivalent atoms
#'
#' @param equiv_list
#' @param restype logical indicating whether to prefix with one-letter residue type
#' @param sep character giving separator between residue and atom names
#' @param atom_format represent multiatoms with "Q" (`q`) or with regular expression (`re`) syntax
#'
#' @return `equiv_list` with updated names
equiv_list_name <- function(equiv_list, restype = TRUE, sep = ":", multiatom_format = c("q", "re")) {

	res <- sapply(lapply(lapply(equiv_list, substr, 11, 15), trimws), "[", 1)
	
	if (restype) {
	
		req_type <- aa3_to_aa1[sapply(lapply(lapply(equiv_list, substr, 6, 8), trimws), "[", 1)]
		res <- paste0(req_type, res)
	}
	
	atom_list <- lapply(lapply(equiv_list, substr, 1, 4), trimws)
	atom_start <- sapply(atom_list, "[", 1)
	atom_length <- sapply(atom_list, length)
	atom_start[atom_length > 1] <- sapply(atom_start[atom_length > 1], function(x) substr(x, 1, nchar(x)-1))
	atom_end <- rep("", length(atom_list))
	atom_end[atom_length > 1] <- sapply(atom_list[atom_length > 1], function(x) paste(substr(x, nchar(x[1]), nchar(x[1])), collapse=""))
	
	multiatom_format <- match.arg(multiatom_format)
	
	if (multiatom_format == "q") {
		atom <- atom_start
		substr(atom[atom_length > 1], 1, 1) <- "Q"
	} else if (multiatom_format == "re") {
		atom <- paste0(atom_start, ifelse(nchar(atom_end), paste0("[", atom_end, "]"), ""))
	}
	
	names(equiv_list) <- paste0(res, sep, atom)
	
	equiv_list
}

#' Encapsulate kinetic ensemble data into a list data structure
#'
#' @param coord_array 
#' @param base_rate_mat 
#' @param base_rates 
#' @param kc reciprocal of molecular tumbling time
#' @param kmethyl methyl rotation eigenvalue
#' @param karo phenylalanine/tyrosine flip eigenvalue
#' @param proton_mhz 
#' @param mix_times mixing times to calculate
make_ke_data <- function(coord_array, base_rate_mat, base_rates, kc, kmethyl = 1/1e-12, karo = 1/100e-6, proton_mhz, mix_times) {

	perm_methyl <- find_methyl_permutations(dimnames(coord_array)[[2]])
	perm_aro <- find_aromatic_permutations(dimnames(coord_array)[[2]])
	
	names(perm_methyl) <- rep(kmethyl, length(perm_methyl))
	names(perm_aro) <- rep(karo, length(perm_aro))
	
	permutations <- c(perm_methyl, perm_aro)
	
	equiv_list <- equiv_list_name(make_equiv_list(dimnames(coord_array)[[2]], permutations))
	
	list(
		coord_array = coord_array,
		base_rate_mat = base_rate_mat,
		rates = c(base_rates, kc=kc, kmethyl=kmethyl, karo=karo),
		equiv_list = equiv_list,
		permutations = permutations,
		proton_mhz = proton_mhz,
		mix_times = mix_times
	)
}

#' Create data structures for calculating sigma/rho from dipole-dipole interaction tensors
#'
#' @param ke_data list structure returned by `make_ke_data()`
#' @param equiv_pair_mat matrix with pairs of equivalent atom labels matching names of
#'    `ke_data$equiv_list`
#' @param perm_internal logical indicating whether to include atom pair internal to
#'    permutation groups
#' @param optional numeric vector of cross relaxation rates
make_spec_den_data <- function(ke_data, equiv_pair_mat, perm_internal = FALSE, sigma = NULL) {

	# create rate matrices
	base_rate_mat <- ke_data$base_rate_mat
	aro_rate_mat <- rate_mat_simple(ke_data$rates["karo"], paste0("a", 1:2))
	methyl_rate_mat <- rate_mat_simple(ke_data$rates["kmethyl"], paste0("m", 1:3))
	
	rate_data_1_1 <- get_rate_data(base_rate_mat, validate = TRUE)
	rate_data_1_2 <- get_rate_data(aro_rate_mat, rate_data_1_1, validate = TRUE)
	rate_data_1_3 <- get_rate_data(methyl_rate_mat, rate_data_1_1, validate = TRUE)
	rate_data_2_2 <- get_rate_data(aro_rate_mat, rate_data_1_2, validate = TRUE)
	rate_data_2_3 <- get_rate_data(aro_rate_mat, rate_data_1_3, validate = TRUE)
	rate_data_3_3 <- get_rate_data(methyl_rate_mat, rate_data_1_3, validate = TRUE)
	
	rate_data <- list(
		"1-1"=rate_data_1_1,
		"1-2"=rate_data_1_2,
		"1-3"=rate_data_1_3,
		"2-2"=rate_data_2_2,
		"2-3"=rate_data_2_3,
		"3-3"=rate_data_3_3
	)

	rate_count_mat <- get_rate_count_mat(base_rate_mat, ke_data$permutations)
	
	mat_list <- lapply(rate_data, rate_data_to_mat_list, rate_count_mat, ke_data$rates)
	
	equiv_length <- sapply(ke_data$equiv_list, length)
	
	equiv_perm <- rep(NA_integer_, length(ke_data$equiv_list))
	names(equiv_perm) <- names(ke_data$equiv_list)
	for (i in seq_along(ke_data$permutations)) {
		equiv_idx <- sapply(ke_data$equiv_list, function(x) any(ke_data$permutations[[i]][,1] %in% x))
		equiv_perm[equiv_idx] <- i
	}
	
	# remove equivalent pairs that are within the same permutation group
	bad_idx <- which(equiv_perm[equiv_pair_mat[,1]] == equiv_perm[equiv_pair_mat[,2]])
	if (length(bad_idx)) {
		equiv_pair_mat <- equiv_pair_mat[-bad_idx,,drop=FALSE]
		if (!is.null(sigma)) {
			sigma <- sigma[-bad_idx]
		}
	}
	
	# sort equivalent group pairs into order of the number of atoms each has
	equiv_pair_mat <- t(apply(equiv_pair_mat, 1, function(x) x[order(equiv_length[x])]))
	equiv_pair_mat <- matrix(equiv_pair_mat, ncol=2)
	
	unique_length_pairs <- unique(matrix(equiv_length[equiv_pair_mat], ncol=2))
	
	# add unique length pairs for internal permutations if necessary
	if (perm_internal && any(sapply(ke_data$permutations, ncol) == 2)) {
		unique_length_pairs <- unique(rbind(unique_length_pairs, c(1L, 2L)))
	}
	if (perm_internal && any(sapply(ke_data$permutations, ncol) == 3)) {
		unique_length_pairs <- unique(rbind(unique_length_pairs, c(1L, 3L)))
	}
	
	unique_length_pairs <- unique_length_pairs[order(unique_length_pairs[,1], unique_length_pairs[,2]),,drop=FALSE]

	spec_den_data_list <- list()
	
	for (i in seq_len(nrow(unique_length_pairs))) {
	
		pair_name <- paste(unique_length_pairs[i,], collapse="-")
		
		idx <- equiv_length[equiv_pair_mat[,1]] == unique_length_pairs[i,1] & equiv_length[equiv_pair_mat[,2]] == unique_length_pairs[i,2]
		
		atomid_list <- lapply(which(idx), function(i) {
			as.matrix(expand.grid(ke_data$equiv_list[[equiv_pair_mat[i,1]]], ke_data$equiv_list[[equiv_pair_mat[i,2]]], stringsAsFactors=FALSE))
		})
		
		if (perm_internal) {
		
			if (all(unique_length_pairs[i,] == c(1L,2L))) {
				aro_perm <- ke_data$permutations[sapply(ke_data$permutations, ncol) == 2]
				aro_atomid_list <- unname(lapply(lapply(aro_perm, t), unname))
				atomid_list <- c(atomid_list, lapply(aro_atomid_list, "[", , c(1,2)))
				atomid_list <- c(atomid_list, lapply(aro_atomid_list, "[", , c(1,3)))
				atomid_list <- c(atomid_list, lapply(aro_atomid_list, "[", , c(1,4)))
				atomid_list <- c(atomid_list, lapply(aro_atomid_list, "[", , c(2,4)))
			}
			if (all(unique_length_pairs[i,] == c(1L,3L))) {
				methyl_perm <- ke_data$permutations[sapply(ke_data$permutations, ncol) == 3]
				methyl_atomid_list <- unname(lapply(lapply(methyl_perm, t), unname))
				atomid_list <- c(atomid_list, lapply(methyl_atomid_list, "[", , 1:2))
			}
		}
		
		atomid_array <- aperm(simplify2array(atomid_list), c(3, 1, 2))
	
		atom_pair_mat <- atomid_array
		dim(atom_pair_mat) <- c(prod(dim(atomid_array)[1:2]), dim(atomid_array)[3])
		
		colnames(atom_pair_mat) <- c("atom1", "atom2")
		
		if (!is.null(sigma)) {
			atom_pair_mat <- data.frame(atom_pair_mat, sigma=NA_real_)
			atom_pair_mat[seq_along(which(idx)),"sigma"] <- sigma[which(idx)]
		}
		
		spec_den_data_list[[pair_name]] <- c(
			list(atom_pairs = atom_pair_mat),
			mat_list[[pair_name]]
		)
	}
	
	spec_den_data_list
}

#' Tests with toy example from Smith 2020 J Biomol NMR
test_toy <- function() {

	toy_r_mat <- matrix(c(
		0.848351683690084, -0.529433112659379, 0,
		0.966177888683851, 0.257876496444355, 0,
		0.966177888683851, -0.257876496444355, 0,
		0.848351683690084, 0.529433112659379, 0
	), ncol=3, byrow=TRUE, dimnames=list(paste0("M", 1:4), c("x","y","z")))

	toy_grouping_list <- list(
		list(c(1L, 2L, 3L, 4L)),
		list(c(1L, 2L), c(3L, 4L)),
		list(1L, 2L, 3L, 4L)
	)

	toy_r_array <- array(toy_r_mat, dim=c(1, dim(toy_r_mat)), dimnames=c(list(NULL), dimnames(toy_r_mat)))
	toy_d_array <- r_array_to_d_array(toy_r_array, gradient=TRUE)
	toy_g_list <- list(
		g1=d_array_to_g(toy_d_array, toy_grouping_list[[1]], gradient=TRUE),
		g2=d_array_to_g(toy_d_array, toy_grouping_list[[2]], gradient=TRUE),
		g3=d_array_to_g(toy_d_array, toy_grouping_list[[3]], gradient=TRUE)
	)

	# derivative checks
	stopifnot(abs(deriv_check(r_array_to_d_array, toy_r_array, dv=1e-8, vdims=3, gdims=4)) < 1e-6)
	stopifnot(abs(deriv_check(d_array_to_g, toy_d_array, dv=1e-8, vdims=2:3, gdims=2:3, grouping=toy_grouping_list[[1]])) < 1e-7)
	stopifnot(abs(deriv_check(d_array_to_g, toy_d_array, dv=1e-8, vdims=2:3, gdims=2:3, grouping=toy_grouping_list[[2]])) < 1e-8)
	stopifnot(abs(deriv_check(d_array_to_g, toy_d_array, dv=1e-8, vdims=2:3, gdims=2:3, grouping=toy_grouping_list[[3]])) < 1e-8)
}

#' Tests with random internuclear vectors of length 1
test_random <- function() {

	if (!"r_mat_check" %in% ls()) {
		r_mat_check <- matrix(rnorm(30), ncol=3)
		r_mat_check <- r_mat_check/sqrt(rowSums(r_mat_check^2))
	}

	# derivative checks
	stopifnot(abs(deriv_check(r_array_to_d_array, r_mat_check, dv=1e-8, vdims=2, gdims=3)) < 1e-6)
}

#' Tests with EROS3 ensemble subset (First two ensemble members and five atoms from M1)
test_eros3 <- function() {

	eros3_coord <- read_ensemble("https://files.rcsb.org/download/6V5D.pdb")

	eros3_sub_atom_pairs <- matrix(c(
		" HA  MET A   1 ", " HB2 MET A   1 ",
		" HA  MET A   1 ", " HB3 MET A   1 ",
		" HA  MET A   1 ", " HG2 MET A   1 ",
		" HA  MET A   1 ", " HG3 MET A   1 "
	), ncol=2, byrow=TRUE)

	eros3_grouping_list <- list(list(1L, 2L), list(c(1L, 2L)))

	eros3_sub_coord <- aperm(eros3_coord[,dimnames(eros3_coord)[[2]] %in% eros3_sub_atom_pairs,1:2], c(2,1,3))
	eros3_sub_r_array <- coord_array_to_r_array(eros3_sub_coord, eros3_sub_atom_pairs)
	eros3_sub_d_array <- r_array_to_d_array(eros3_sub_r_array, gradient=TRUE)
	eros3_sub_g_list <- list(
		g1=d_array_to_g(eros3_sub_d_array, eros3_grouping_list[[1]], gradient=TRUE),
		g2=d_array_to_g(eros3_sub_d_array, eros3_grouping_list[[2]], gradient=TRUE)
	)

	eros3_sub_g <- coord_array_to_g(eros3_sub_coord, eros3_sub_atom_pairs, eros3_grouping_list)

	# get g values using M1 twice
	eros3_sub_1_g <- coord_array_to_g(eros3_sub_coord[,,c(1,1)], eros3_sub_atom_pairs, eros3_grouping_list)

	eros3_sub_energy <- coord_array_to_g_energy(eros3_sub_coord, eros3_sub_atom_pairs, eros3_grouping_list, eros3_sub_1_g, 1, gradient=TRUE)

	# derivative checks
	stopifnot(abs(deriv_check(r_array_to_d_array, eros3_sub_r_array, dv=1e-8, vdims=3, gdims=4)) < 1e-6)
	stopifnot(abs(deriv_check(d_array_to_g, eros3_sub_d_array, dv=1e-8, vdims=2:3, gdims=2:3, grouping=eros3_grouping_list[[1]])) < 1e-8)
	stopifnot(abs(deriv_check(d_array_to_g, eros3_sub_d_array, dv=1e-8, vdims=2:3, gdims=2:3, grouping=eros3_grouping_list[[2]])) < 1e-8)
	stopifnot(abs(deriv_check(coord_array_to_g_energy, eros3_sub_coord, dv=1e-8, vdims=1:3, gdims=1:3, atom_pairs=eros3_sub_atom_pairs, grouping_list=eros3_grouping_list, g0=eros3_sub_1_g, k=1)) < 1e-11)
}

test_gb3 <- function() {

	pdb2lum <- read_ensemble("https://files.rcsb.org/download/2LUM.pdb", proton_only = TRUE)[,,1:3]
	
	base_rates <- c(kens=1/2e-9)
	
	base_rate_mat <- rate_mat_simple(base_rates, dimnames(pdb2lum)[[3]])
	
	ke_data <- make_ke_data(pdb2lum, base_rate_mat, base_rates, kc=1/4e-9, proton_mhz=700, mix_times=70e-3)
	
	ke_data$equiv_list <- equiv_list_name(ke_data$equiv_list, restype = FALSE)
	
	atom_pairs <- structure(c("1:HA", "1:H1", "1:H1", "30:QD", "3:QD", "4:QZ", "1:H1", "3:QD", "1:QE", "3:QD", "1:QE", "1:QE"), dim = c(6L, 2L))
	
	spec_den_data <- make_spec_den_data(ke_data, atom_pairs, sigma = rep(0, nrow(atom_pairs)))
	
	coord_array <- aperm(pdb2lum, c(2,1,3))
	
	lapply(spec_den_data, function(dipole_kinetic_data) {
	
		# calculate internuclear vectors (convert from Å^-3 to m^-3)
		r_array <- coord_array_to_r_array(coord_array*1e-10, dipole_kinetic_data[["atom_pairs"]][,1:2,drop=FALSE])
	
		# calculate dipole-dipole interaction tensors
		d_array <- r_array_to_d_array(r_array)
		stopifnot(abs(deriv_check(r_array_to_d_array, r_array, dv=1e-18, vdims=3, gdims=4)) < 1e-6)
	
		# calculate the factor by which the number of models should be expanded
		n_shift <- ncol(dipole_kinetic_data[["groupings"]])/dim(coord_array)[3]
	
		# shift tensor components from atom pairs into virtual models
		d_array_shifted <- array_shift(d_array, n_shift)
		
		# calculate matrix of g values
		g_matrix <- d_array_to_g_matrix(d_array_shifted, dipole_kinetic_data[["groupings"]])
		
		# calculate matrix of a values
		a_matrix <- g_matrix_to_a_matrix(g_matrix, dipole_kinetic_data[["a_coef"]])
		
		# calculate lambda eigenvalues
		lambda_vec <- -colSums(ke_data$rates[rownames(dipole_kinetic_data[["lambda_coef"]])]*dipole_kinetic_data[["lambda_coef"]])
		
		# update eigenvalues to account for molecular tumbling
		lambda_prime_vec <- lambda_vec - ke_data$rates["kc"]
		
		#print(min(attr(a_matrix_to_sigma(a_matrix, lambda_prime_vec, ke_data$proton_mhz, gradient=TRUE), "gradient")))
		#print(a_matrix)
		
		stopifnot(abs(deriv_check(a_matrix_to_sigma, a_matrix, dv=min(a_matrix)*1e-2, vdims=2, gdims=2, lambda_prime_vec=lambda_prime_vec, proton_mhz=ke_data$proton_mhz)) < 1e-68)
	})
}
