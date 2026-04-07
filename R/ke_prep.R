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
#'
#' @export
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
#'
#' @export
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
#'
#' @export
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
#'
#' @export
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
#'
#' @examples
#' pdb2lum <- read_ensemble("https://files.rcsb.org/download/2LUM.pdb", proton_only=TRUE)
#' perm_methyl <- find_methyl_permutations(dimnames(pdb2lum)[[2]])
#' perm_aro <- find_aromatic_permutations(dimnames(pdb2lum)[[2]])
#' unique_atom_pair_map(perm_methyl[[1]])
#' unique_atom_pair_map(perm_aro[[1]])
#'
#' @export
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
#'
#' @export
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
#'
#' @export
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
#' @param rate_mat_a first transition rate matrix
#' @param rate_mat_b second transition rate matrix
#'
#' @examples
#' rate_mat_fast <- rate_mat_simple(10, c("f1", "f2"))
#' rate_mat_slow <- rate_mat_simple(4, c("s1", "s2"))
#' rate_mat <- rate_mat_kronecker(rate_mat_fast, rate_mat_slow)
#' rate_mat <- rate_mat_diag(rate_mat)
#' eigen(rate_mat)
#'
#' @export
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
#' @param rate_mat square transition rate matrix
#'
#' @export
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
#' @param trans_rate transition rate matrix
#' @param all_permutations all possible permutations (i.e. methyl, aromatic) applied
#' @param eps_factor rates must differ by this amount times the lowest rate
#'
#' @examples
#' rate_mat <- rate_mat_intra_inter(c(0,0,1,1), 10, 4)
#' get_rate_count_mat(rate_mat, NULL)
#'
#' @export
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
#' @param trans_rate_eigen list with `values` and `vectors` as returned by `eigen()`
#' @param rate_groups list with numerically equivalent rates organized into groups
#'
#' This only works with rate matrices returned by `rate_mat_simple()` and 
#' `rate_mat_intra_inter`. It does not work with `rate_mat_kronecker()`.
#'
#' @examples
#' rate_mat <- rate_mat_intra_inter(c(0,0,1,1), 10, 4)
#' rate_mat <- rate_mat_diag(rate_mat)
#' calc_subset_mat(eigen(rate_mat))
#'
#' @export
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
#' @param rate_subset_mat matrix (ensemble members, rates) giving the subset each ensemble member belongs to
#' @param permutation_counts vector with counts of states and names giving the rate
#'
#' @examples
#' rate_mat_fast <- rate_mat_simple(10, c("f1", "f2"))
#' rate_mat_fast <- rate_mat_diag(rate_mat_fast)
#' subset_mat <- calc_subset_mat(eigen(rate_mat_fast))
#' expand_subset_mat(subset_mat, list("4"=2))
#'
#' @export
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
#' @param trans_rate_eigen list with `values` and `vectors` as returned by `eigen()`
#' @param all_rates numeric vector of all rates
#' @param subset_mat created by `calc_subset_mat()` or `expand_subset_mat()`
#' @param validate check to see whether result is correct using another method
#' @param eps_factor epsilon used for heuristically matching rates generated in two ways
#'
#' @examples
#' rate_mat <- rate_mat_intra_inter(c(0,0,1,1), 10, 4)
#' rate_mat <- rate_mat_diag(rate_mat)
#' get_eigen_groups(eigen(rate_mat))
#'
#' @export
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
#' rate_data_base <- get_rate_data(rate_mat_fast)
#' rate_data <- get_rate_data(rate_mat_slow, rate_data_base)
#'
#' @export
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
#' @param rates_named numeric vector whose names give the rates to be used
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
#' @param permutations list of permutations from `find_methyl_permutations()` and 
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
#' @param equiv_list list of atomid vectors where all atomids in a given vector are equivalent
#' @param restype logical indicating whether to prefix with one-letter residue type
#' @param sep character giving separator between residue and atom names
#' @param multiatom_format represent multiatoms with "Q" (`q`) or with regular expression (`re`) syntax
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
#' @param coord_array 3D array (xyz, atoms, models) with atomic coordinates
#' @param base_rate_mat rate matrix for transitions between ensemble members
#' @param base_rates named vector of rates
#' @param kc reciprocal of molecular tumbling time
#' @param kmethyl methyl rotation eigenvalue
#' @param karo phenylalanine/tyrosine flip eigenvalue
#' @param proton_mhz proton frequency in MHz
#' @param mix_times mixing times to calculate
#'
#' @export
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
#' @param sigma optional numeric vector of cross relaxation rates
#'
#' @export
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

#' Simplify symbolic spectral-density term names
#'
#' This helper collapses identical sum and difference labels such as
#' `wHpwH -> 2wH` and `wHmwH -> 0`.
#'
#' @param term_names Character vector of spectral-density term names
#'
#' @return Character vector of simplified term names
#'
#' @noRd
.simplify_term_names <- function(term_names) {
	term_names <- sub("^w([A-Za-z0-9]+)pw\\1$", "2w\\1", term_names)
	term_names <- sub("^w([A-Za-z0-9]+)mw\\1$", "0", term_names)
	term_names
}

#' Build a spectral-density term array from term metadata
#'
#' Repeated terms with the same `term_name` are merged by summing their
#' coefficients while requiring a common frequency. In other words, if several
#' contributions share the same spectral density label \eqn{t}, this helper
#' constructs
#' \deqn{c_t = \sum_k c_{tk}}
#' while keeping a single associated \eqn{\omega_t}. Terms with merged
#' coefficient \eqn{c_t = 0} are dropped from the returned array.
#'
#' @param n_pairs Integer number of rows
#' @param term_names Character vector of term names
#' @param coef Numeric vector of term coefficients
#' @param freq Numeric vector of term frequencies
#'
#' @return Numeric array `(pairs, terms, 2)`
#'
#' @noRd
.spec_den_terms_to_array <- function(n_pairs, term_names, coef, freq) {
	stopifnot(
		is.numeric(n_pairs),
		length(n_pairs) == 1,
		n_pairs >= 0,
		n_pairs == as.integer(n_pairs),
		is.character(term_names),
		is.numeric(coef),
		is.numeric(freq),
		length(term_names) == length(coef),
		length(term_names) == length(freq)
	)

	term_names <- .simplify_term_names(term_names)
	term_coef <- vapply(unique(term_names), function(term_name) {
		sum(coef[term_names == term_name])
	}, numeric(1))
	term_levels <- names(term_coef)[term_coef != 0]
	spec_den_term_array <- array(
		0,
		dim = c(n_pairs, length(term_levels), 2),
		dimnames = list(NULL, term_levels, c("coef", "freq"))
	)

	for (term_name in term_levels) {
		idx <- term_names == term_name
		term_freq <- unique(freq[idx])
		if (length(term_freq) != 1) {
			stop("Term `", term_name, "` was assigned multiple frequencies")
		}
		spec_den_term_array[, term_name, "coef"] <- term_coef[[term_name]]
		spec_den_term_array[, term_name, "freq"] <- term_freq
	}

	spec_den_term_array
}

#' Make spectral density term array for heteronuclear R1 relaxation
#'
#' This function encodes the coefficients and frequencies needed for
#' `a_matrix_to_relax()` to calculate the heteronuclear longitudinal
#' relaxation rate \eqn{R_1} of the observed spin \eqn{I}, including both
#' dipole-dipole (DD) and chemical shift anisotropy (CSA) contributions.
#'
#' The implemented expression is
#' \deqn{
#' R_1(I) =
#' \frac{1}{10} d_{IS}^2
#' \left[
#' J(\omega_I - \omega_S) + 3J(\omega_I) + 6J(\omega_I + \omega_S)
#' \right]
#' +
#' c_I^2 J(\omega_I),
#' }
#' where
#' \deqn{
#' d_{IS}^2 =
#' \left(
#' \frac{\mu_0}{4 \pi} \hbar \gamma_I \gamma_S r_{IS}^{-3}
#' \right)^2
#' }
#' and
#' \deqn{
#' c_I^2 =
#' \frac{2}{15}\omega_I^2 \Delta\sigma_I^2.
#' }
#' Here \eqn{\Delta\sigma_I} is the CSA anisotropy of the observed spin in
#' fractional units, obtained from `delta_sigma_ppm`.
#'
#' The dipolar part is obtained from Solomon's longitudinal relaxation
#' coefficient \eqn{p} for two unlike spins,
#' \deqn{p = w_0 + 2w_1 + w_2,}
#' given in Solomon (1955, Eq. 15), together with the unlike-spin transition
#' probabilities \eqn{w_0}, \eqn{w_1}, and \eqn{w_2} for a pure dipole-dipole
#' interaction in Solomon (1955, Eq. 35). Writing those transition probabilities
#' in spectral-density form gives the standard dipolar contribution
#' \deqn{
#' R_1^{DD}(I) = \frac{1}{10} d_{IS}^2
#' \left[
#' J(\omega_I - \omega_S) + 3J(\omega_I) + 6J(\omega_I + \omega_S)
#' \right].
#' }
#' In Abragam's notation, the same dipolar longitudinal relaxation expression
#' appears in Chapter VIII, Eq. (88) of \emph{The Principles of Nuclear
#' Magnetism} (Abragam, 1961), one notation layer upstream in terms of the
#' rank-specific spectral densities \eqn{J^{(0)}} \eqn{J^{(1)}}, and
#' \eqn{J^{(2)}}.
#'
#' The CSA term corresponds to Abragam's treatment of shielding-anisotropy
#' relaxation in the same chapter, where for an
#' axially symmetric CSA tensor (\eqn{\eta = 0}) his Eq. (141) can be written
#' in the present notation as
#' \deqn{R_1^{CSA} = c_I^2 J(\omega_I)}
#' with the identification
#' \deqn{c_I^2 = (2/15)\omega_I^2 \Delta\sigma_I^2.}
#' The present function then adds this CSA contribution to the dipolar result.
#'
#' This function uses a single spectral density function for both DD and CSA
#' terms. That is an effective approximation appropriate when the CSA and
#' dipolar interaction tensors are assumed to experience the same underlying
#' motional spectral density.
#'
#' The underlying transition-probability treatment follows Solomon (1955)
#' \doi{10.1103/PhysRev.99.559}, which in turn builds on the rapid-motion
#' treatment of Abragam and Pound (1953) \doi{10.1103/PhysRev.92.943} and the
#' dipolar correlation-function formalism of Bloembergen, Purcell, and Pound
#' (1948) \doi{10.1103/PhysRev.73.679}.
#'
#' @param n_pairs integer number of atom pairs
#' @param proton_mhz spectrometer proton field strength in MHz
#' @param nucleus_i character scalar giving the observed nucleus, such as
#'   `"15N"`
#' @param nucleus_s character scalar giving the coupled partner nucleus, such
#'   as `"1H"`
#' @param r_is_angstrom numeric scalar internuclear distance \eqn{r_{IS}} in
#'   Angstrom
#' @param delta_sigma_ppm numeric scalar CSA anisotropy of the observed spin in
#'   ppm
#'
#' @return Array with dimensions `(pairs, terms, components)`. The term
#'   dimension is named with simplified frequency labels such as
#'   `c("wHmwN", "wN", "wHpwN")`, and the component dimension is named
#'   `c("coef", "freq")`.
#'
#' @export
make_r1_spec_den_term_array <- function(n_pairs, proton_mhz, nucleus_i = "15N", nucleus_s = "1H", r_is_angstrom, delta_sigma_ppm) {
	stopifnot(
		is.numeric(n_pairs),
		length(n_pairs) == 1,
		n_pairs >= 0,
		n_pairs == as.integer(n_pairs),
		is.numeric(proton_mhz),
		length(proton_mhz) == 1,
		is.character(nucleus_i),
		length(nucleus_i) == 1,
		is.character(nucleus_s),
		length(nucleus_s) == 1,
		is.numeric(r_is_angstrom),
		length(r_is_angstrom) == 1,
		r_is_angstrom > 0,
		is.numeric(delta_sigma_ppm),
		length(delta_sigma_ppm) == 1
	)
	if (nucleus_i == nucleus_s) {
		stop("`make_r1_spec_den_term_array()` currently implements the heteronuclear case only")
	}

	d_sq <- .dipolar_prefactor_sq(nucleus_i, nucleus_s) / (r_is_angstrom * 1e-10)^6
	c_sq <- .csa_prefactor_sq(proton_mhz, nucleus_i, delta_sigma_ppm)
	omega_i <- .larmor_omega(proton_mhz, nucleus_i)
	omega_s <- .larmor_omega(proton_mhz, nucleus_s)
	label_i <- .nucleus_label(nucleus_i)
	label_s <- .nucleus_label(nucleus_s)

	.spec_den_terms_to_array(
		n_pairs = n_pairs,
		term_names = c(
			paste0("w", label_s, "mw", label_i),
			paste0("w", label_i),
			paste0("w", label_s, "pw", label_i),
			paste0("w", label_i)
		),
		coef = c(
			0.1 * d_sq,
			0.3 * d_sq,
			0.6 * d_sq,
			c_sq
		),
		freq = c(
			omega_s - omega_i,
			omega_i,
			omega_s + omega_i,
			omega_i
		)
	)
}

#' Make spectral density term array for heteronuclear R2 relaxation
#'
#' This function encodes the coefficients and frequencies needed for
#' `a_matrix_to_relax()` to calculate the heteronuclear transverse relaxation
#' rate \eqn{R_2} of the observed spin \eqn{I}, including both dipole-dipole
#' (DD) and chemical shift anisotropy (CSA) contributions.
#'
#' The implemented expression is
#' \deqn{
#' R_2(I) =
#' \frac{1}{20} d_{IS}^2
#' \left[
#' 4J(0) + J(\omega_I - \omega_S) + 3J(\omega_I) +
#' 6J(\omega_S) + 6J(\omega_I + \omega_S)
#' \right]
#' +
#' \frac{1}{6} c_I^2
#' \left[
#' 4J(0) + 3J(\omega_I)
#' \right],
#' }
#' where
#' \deqn{
#' d_{IS}^2 =
#' \left(
#' \frac{\mu_0}{4 \pi} \hbar \gamma_I \gamma_S r_{IS}^{-3}
#' \right)^2
#' }
#' and
#' \deqn{
#' c_I^2 =
#' \frac{2}{15}\omega_I^2 \Delta\sigma_I^2.
#' }
#' Here \eqn{\Delta\sigma_I} is the CSA anisotropy of the observed spin in
#' fractional units, obtained from `delta_sigma_ppm`.
#'
#' The dipolar part is obtained from Solomon's transverse relaxation
#' coefficient \eqn{\nu} for two unlike spins,
#' \deqn{\nu = w_0 + 2w_1' + w_2,}
#' given in Solomon (1955, Eq. 21), together with the unlike-spin transition
#' probabilities in Solomon (1955, Eq. 35). Rewriting those transition
#' probabilities in spectral-density form yields the standard dipolar
#' contribution
#' \deqn{
#' R_2^{DD}(I) = \frac{1}{20} d_{IS}^2
#' \left[
#' 4J(0) + J(\omega_I - \omega_S) + 3J(\omega_I) +
#' 6J(\omega_S) + 6J(\omega_I + \omega_S)
#' \right].
#' }
#' In Abragam's notation, the same dipolar transverse relaxation expression
#' appears in Chapter VIII, Eq. (89) of \emph{The Principles of Nuclear
#' Magnetism} (Abragam, 1961), again written in terms of the rank-specific
#' spectral densities \eqn{J^{(0)}} \eqn{J^{(1)}}, and \eqn{J^{(2)}}.
#'
#' The CSA term corresponds to Abragam's treatment of shielding-anisotropy
#' relaxation in the same chapter, where for an
#' axially symmetric CSA tensor (\eqn{\eta = 0}) his Eq. (142) can be written
#' in the present notation as
#' \deqn{
#' R_2^{CSA} = \frac{1}{6} c_I^2 \left[4J(0) + 3J(\omega_I)\right]
#' }
#' with the identification
#' \deqn{c_I^2 = (2/15)\omega_I^2 \Delta\sigma_I^2.}
#' The function then adds this CSA contribution to the dipolar result.
#'
#' This function uses a single spectral density function for both DD and CSA
#' terms. That is an effective approximation appropriate when the CSA and
#' dipolar interaction tensors are assumed to experience the same underlying
#' motional spectral density.
#'
#' The underlying transition-probability treatment follows Solomon (1955)
#' \doi{10.1103/PhysRev.99.559}, which in turn builds on the rapid-motion
#' treatment of Abragam and Pound (1953) \doi{10.1103/PhysRev.92.943} and the
#' dipolar correlation-function formalism of Bloembergen, Purcell, and Pound
#' (1948) \doi{10.1103/PhysRev.73.679}.
#'
#' @param n_pairs integer number of atom pairs
#' @param proton_mhz spectrometer proton field strength in MHz
#' @param nucleus_i character scalar giving the observed nucleus, such as
#'   `"15N"`
#' @param nucleus_s character scalar giving the coupled partner nucleus, such
#'   as `"1H"`
#' @param r_is_angstrom numeric scalar internuclear distance \eqn{r_{IS}} in
#'   Angstrom
#' @param delta_sigma_ppm numeric scalar CSA anisotropy of the observed spin in
#'   ppm
#'
#' @return Array with dimensions `(pairs, terms, components)`. The term
#'   dimension is named with simplified frequency labels such as
#'   `c("0", "wHmwN", "wN", "wH", "wHpwN")`, and the component dimension is
#'   named `c("coef", "freq")`.
#'
#' @export
make_r2_spec_den_term_array <- function(n_pairs, proton_mhz, nucleus_i = "15N", nucleus_s = "1H", r_is_angstrom, delta_sigma_ppm) {
	stopifnot(
		is.numeric(n_pairs),
		length(n_pairs) == 1,
		n_pairs >= 0,
		n_pairs == as.integer(n_pairs),
		is.numeric(proton_mhz),
		length(proton_mhz) == 1,
		is.character(nucleus_i),
		length(nucleus_i) == 1,
		is.character(nucleus_s),
		length(nucleus_s) == 1,
		is.numeric(r_is_angstrom),
		length(r_is_angstrom) == 1,
		r_is_angstrom > 0,
		is.numeric(delta_sigma_ppm),
		length(delta_sigma_ppm) == 1
	)
	if (nucleus_i == nucleus_s) {
		stop("`make_r2_spec_den_term_array()` currently implements the heteronuclear case only")
	}

	d_sq <- .dipolar_prefactor_sq(nucleus_i, nucleus_s) / (r_is_angstrom * 1e-10)^6
	c_sq <- .csa_prefactor_sq(proton_mhz, nucleus_i, delta_sigma_ppm)
	omega_i <- .larmor_omega(proton_mhz, nucleus_i)
	omega_s <- .larmor_omega(proton_mhz, nucleus_s)
	label_i <- .nucleus_label(nucleus_i)
	label_s <- .nucleus_label(nucleus_s)

	.spec_den_terms_to_array(
		n_pairs = n_pairs,
		term_names = c(
			"0",
			paste0("w", label_s, "mw", label_i),
			paste0("w", label_i),
			paste0("w", label_s),
			paste0("w", label_s, "pw", label_i),
			"0",
			paste0("w", label_i)
		),
		coef = c(
			0.2 * d_sq,
			0.05 * d_sq,
			0.15 * d_sq,
			0.3 * d_sq,
			0.3 * d_sq,
			(4 / 6) * c_sq,
			0.5 * c_sq
		),
		freq = c(
			0,
			omega_s - omega_i,
			omega_i,
			omega_s,
			omega_s + omega_i,
			0,
			omega_i
		)
	)
}

#' Make spectral density term array for dipolar sigma cross relaxation
#'
#' This function encodes the coefficients and frequencies needed for
#'    `a_matrix_to_relax()` to calculate the dipolar cross-relaxation
#'    rate sigma for a pair of nuclei.
#'
#' The sigma cross-relaxation rate is defined
#' \deqn{
#' \sigma_{IS} =
#' \frac{1}{10} d_{IS}^2
#' \left( 6 J(\omega_I + \omega_S) - J(\omega_I - \omega_S) \right)
#' }
#' with
#' \deqn{
#' d_{IS}^2 =
#' \left(
#' \frac{\mu_0}{4 \pi} \hbar \gamma_I \gamma_S r_{IS}^{-3}
#' \right)^2.
#' }
#'
#' The returned array stores the two spectral density terms:
#' \deqn{J(\omega_I - \omega_S) \text{ with coefficient } -\frac{1}{10} d_{IS}^2}
#' \deqn{J(\omega_I + \omega_S) \text{ with coefficient } \frac{3}{5} d_{IS}^2.}
#'
#' This dipolar cross-relaxation term is obtained from Solomon's longitudinal
#' cross coefficient \eqn{\sigma = w_2 - w_0} for two unlike spins in Solomon
#' (1955, Eq. 15), together with the unlike-spin transition probabilities
#' \eqn{w_0} and \eqn{w_2} in Solomon (1955, Eq. 35). The same unlike-spin 
#' longitudinal cross-relaxation term appears in Chapter VIII, Eq. (88) of 
#' \emph{The Principles of Nuclear Magnetism} (Abragam, 1961).
#'
#' The static magnetic field is inferred from `proton_mhz` using the proton
#' gyromagnetic ratio, and the Larmor frequencies of the requested nuclei are
#' then calculated from an internal lookup table of gyromagnetic ratios.
#' Supported nucleus names currently include `"1H"`, `"13C"`, `"15N"`,
#' `"19F"`, `"31P"`, and `"2H"`.
#'
#' If `r_is_angstrom` is `NA`, the returned coefficients omit the
#' \eqn{r_{IS}^{-6}} factor and are therefore suitable for the distance-
#' dependent representation in which internuclear distance is already encoded in
#' the dipolar interaction amplitudes. If `r_is_angstrom` is supplied, the
#' returned coefficients include the \eqn{r_{IS}^{-6}} factor explicitly and
#' are suitable for a unit-tensor representation.
#'
#' The underlying transition-probability treatment follows Solomon (1955)
#' \doi{10.1103/PhysRev.99.559}, which in turn builds on the rapid-motion
#' treatment of Abragam and Pound (1953) \doi{10.1103/PhysRev.92.943} and the
#' dipolar correlation-function formalism of Bloembergen, Purcell, and Pound
#' (1948) \doi{10.1103/PhysRev.73.679}.
#'
#' @param n_pairs integer number of atom pairs
#' @param proton_mhz spectrometer proton field strength in MHz
#' @param nucleus_i character scalar giving the first nucleus, such as `"1H"`
#' @param nucleus_s character scalar giving the second nucleus. Defaults to
#'   `nucleus_i`.
#' @param r_is_angstrom optional numeric scalar internuclear distance
#'   \eqn{r_{IS}} in Angstrom. If `NA`, distance dependence is omitted from the
#'   coefficients.
#'
#' @return Array with dimensions `(pairs, terms, components)`. The term
#'    dimension is named according to the simplified frequency labels implied
#'    by `nucleus_i` and `nucleus_s`, such as `c("0", "2wH")` for
#'    homonuclear proton sigma or `c("wHmwN", "wHpwN")` for heteronuclear
#'    proton-nitrogen sigma. The component dimension is named
#'    `c("coef", "freq")`.
#'
#' @export
make_sigma_spec_den_term_array <- function(n_pairs, proton_mhz, nucleus_i = "1H", nucleus_s = nucleus_i, r_is_angstrom = NA_real_) {
	stopifnot(
		is.numeric(n_pairs),
		length(n_pairs) == 1,
		n_pairs >= 0,
		n_pairs == as.integer(n_pairs),
		is.numeric(proton_mhz),
		length(proton_mhz) == 1,
		is.character(nucleus_i),
		length(nucleus_i) == 1,
		is.character(nucleus_s),
		length(nucleus_s) == 1,
		is.numeric(r_is_angstrom),
		length(r_is_angstrom) == 1
	)
	if (!is.na(r_is_angstrom) && r_is_angstrom <= 0) {
		stop("`r_is_angstrom` must be positive when supplied")
	}

	omega_i <- .larmor_omega(proton_mhz, nucleus_i)
	omega_s <- .larmor_omega(proton_mhz, nucleus_s)
	d_sq <- .dipolar_prefactor_sq(nucleus_i, nucleus_s)
	if (!is.na(r_is_angstrom)) {
		d_sq <- d_sq / (r_is_angstrom * 1e-10)^6
	}
	label_i <- .nucleus_label(nucleus_i)
	label_s <- .nucleus_label(nucleus_s)

	.spec_den_terms_to_array(
		n_pairs = n_pairs,
		term_names = c(
			paste0("w", label_i, "mw", label_s),
			paste0("w", label_i, "pw", label_s)
		),
		coef = c(-0.1 * d_sq, 0.6 * d_sq),
		freq = c(omega_i - omega_s, omega_i + omega_s)
	)
}

#' Convert steady-state NOE values to sigma cross-relaxation rates
#'
#' For observed nucleus \eqn{X}, the steady-state heteronuclear NOE is related
#' to the dipolar cross-relaxation rate \eqn{\sigma_{HX}} through
#' \deqn{
#' \mathrm{NOE} = 1 + \frac{\gamma_H}{\gamma_X} \frac{\sigma_{HX}}{R_1}.
#' }
#' This function rearranges that relation to evaluate
#' \deqn{
#' \sigma_{HX} = (\mathrm{NOE} - 1) R_1 \frac{\gamma_X}{\gamma_H}.
#' }
#'
#' When `dnoe` and `dr1` are supplied, the uncertainty in `sigma` is propagated
#' under the assumption that the NOE and \eqn{R_1} uncertainties are
#' independent:
#' \deqn{
#' d\sigma =
#' \sqrt{
#' \left[(\mathrm{NOE} - 1)\frac{\gamma_X}{\gamma_H} dR_1\right]^2 +
#' \left[R_1 \frac{\gamma_X}{\gamma_H} d\mathrm{NOE}\right]^2
#' }.
#' }
#'
#' The returned numeric matrix always contains a `sigma` column and includes
#' `sigma_err` when both `dnoe` and `dr1` are supplied.
#'
#' @param noe Numeric vector of steady-state heteronuclear NOE values
#' @param r1 Numeric vector of longitudinal relaxation rates corresponding to
#'   `noe`
#' @param dnoe Optional numeric vector of NOE uncertainties
#' @param dr1 Optional numeric vector of \eqn{R_1} uncertainties
#' @param nucleus_x Character scalar observed nucleus identifier such as
#'   `"15N"` or `"13C"`
#'
#' @return Numeric matrix with column `sigma`, and optional column `sigma_err`
#'   when both error vectors are supplied
#'
#' @export
noe_to_sigma <- function(noe, r1, dnoe = NULL, dr1 = NULL, nucleus_x = "15N") {
	stopifnot(
		is.numeric(noe),
		is.numeric(r1),
		length(noe) == length(r1),
		is.character(nucleus_x),
		length(nucleus_x) == 1
	)
	if (xor(is.null(dnoe), is.null(dr1))) {
		stop("`dnoe` and `dr1` must either both be supplied or both be `NULL`")
	}
	if (!is.null(dnoe)) {
		stopifnot(
			is.numeric(dnoe),
			is.numeric(dr1),
			length(dnoe) == length(noe),
			length(dr1) == length(r1)
		)
	}

	gamma_ratio <- .nucleus_gamma(nucleus_x) / .nucleus_gamma("1H")
	sigma <- (noe - 1) * r1 * gamma_ratio

	out <- matrix(sigma, ncol = 1, dimnames = list(NULL, "sigma"))
	if (!is.null(dnoe)) {
		sigma_err <- sqrt(
			(((noe - 1) * gamma_ratio * dr1) ^ 2) +
			((r1 * gamma_ratio * dnoe) ^ 2)
		)
		out <- cbind(out, sigma_err = sigma_err)
	}

	out
}

#' Normalize permutation-rate metadata for spectral-density relaxation input
#'
#' This helper validates a length-two `perm_rates` specification. Missing
#' values indicate no permutation on that side, while non-missing values give
#' the symbolic kinetic rate constants to associate with the inferred
#' permutation processes on the `atom1` and `atom2` sides.
#'
#' @param perm_rates numeric vector of length 2
#'
#' @return A list with elements `rate_values` and `rate_names`
#'
#' @noRd
.normalize_spec_den_perm_rates <- function(perm_rates) {
	stopifnot(
		is.numeric(perm_rates),
		length(perm_rates) == 2
	)

	rate_names <- names(perm_rates)
	if (is.null(rate_names)) {
		rate_names <- rep("", length(perm_rates))
	}
	rate_names[is.na(perm_rates)] <- ""

	bad_named <- !is.na(perm_rates) & rate_names == ""
	if (any(bad_named)) {
		stop("Non-missing `perm_rates` entries must be named with symbolic rate constants")
	}
	if (any(!is.na(perm_rates) & perm_rates <= 0)) {
		stop("Non-missing `perm_rates` entries must be positive rate constants")
	}

	rate_values <- perm_rates
	names(rate_values) <- c("atom1", "atom2")
	names(rate_names) <- c("atom1", "atom2")

	list(
		rate_values = rate_values,
		rate_names = rate_names
	)
}

#' Normalize unique symbolic permutation rates
#'
#' @param perm_info output from `.normalize_spec_den_perm_rates()`
#'
#' @return named numeric vector of unique symbolic permutation rates
#'
#' @noRd
.unique_spec_den_perm_rates <- function(perm_info) {
	side_rate_names <- perm_info[["rate_names"]][perm_info[["rate_names"]] != ""]
	side_rate_values <- perm_info[["rate_values"]][perm_info[["rate_names"]] != ""]

	if (!length(side_rate_names)) {
		return(numeric())
	}

	unique_rate_values <- tapply(unname(side_rate_values), side_rate_names, unique, simplify = FALSE)
	bad_idx <- lengths(unique_rate_values) != 1L
	if (any(bad_idx)) {
		stop(
			"Permutation rate constants with the same symbolic name must have identical numeric values: ",
			paste(names(unique_rate_values)[bad_idx], collapse = ", ")
		)
	}

	unlist(unique_rate_values, use.names = TRUE)
}

#' Detect internal permutation blocks from atom-name patterns
#'
#' This helper identifies the special within-permutation blocks used for
#' methyl and aromatic internal interactions. These blocks have size 2 or 3
#' and, by convention, the atom identifiers on each side differ pairwise by a
#' single character.
#'
#' @param block two-column character matrix or data frame containing one block
#'   of expanded atom-pair rows
#'
#' @return Logical scalar indicating whether `block` matches the internal-block
#'   heuristic
#'
#' @noRd
.is_internal_spec_den_block <- function(block) {
	block <- as.matrix(block)
	block_size <- nrow(block)
	if (!(block_size %in% c(2L, 3L)) || ncol(block) < 2) {
		return(FALSE)
	}

	atom1_levels <- unique(block[, 1])
	atom2_levels <- unique(block[, 2])
	if (length(atom1_levels) != block_size || length(atom2_levels) != block_size) {
		return(FALSE)
	}

	count_char_diffs <- function(x, y) {
		sx <- strsplit(x, "", fixed = TRUE)[[1]]
		sy <- strsplit(y, "", fixed = TRUE)[[1]]
		if (length(sx) != length(sy)) {
			return(Inf)
		}
		sum(sx != sy)
	}

	atom1_diff <- outer(atom1_levels, atom1_levels, Vectorize(count_char_diffs))
	atom2_diff <- outer(atom2_levels, atom2_levels, Vectorize(count_char_diffs))

	all(atom1_diff[upper.tri(atom1_diff)] == 1L) &&
		all(atom2_diff[upper.tri(atom2_diff)] == 1L)
}

#' Classify one expanded atom-pair block
#'
#' @param block two-column character matrix or data frame containing one block
#'   of expanded atom-pair rows
#'
#' @return A list with elements `type`, `multiplicity`, `order`, and
#'   `block_index_order`, or `NULL` if the block does not match a supported
#'   layout
#'
#' @noRd
.classify_spec_den_relax_block <- function(block) {
	block <- as.matrix(block)
	block_size <- nrow(block)
	atom1_levels <- unique(block[, 1])
	atom2_levels <- unique(block[, 2])
	multiplicity <- c(atom1 = length(atom1_levels), atom2 = length(atom2_levels))

	if (block_size %in% c(2L, 3L) && .is_internal_spec_den_block(block)) {
		return(list(
			type = "internal",
			multiplicity = c(atom1 = 1L, atom2 = block_size),
			order = "atom1_fastest",
			block_index_order = seq_len(block_size)
		))
	}

	if (prod(multiplicity) != block_size || any(!(multiplicity %in% c(1L, 2L, 3L)))) {
		return(NULL)
	}

	expected_atom1_fastest <- cbind(
		atom1 = rep(atom1_levels, times = length(atom2_levels)),
		atom2 = rep(atom2_levels, each = length(atom1_levels))
	)
	expected_atom2_fastest <- cbind(
		atom1 = rep(atom1_levels, each = length(atom2_levels)),
		atom2 = rep(atom2_levels, times = length(atom1_levels))
	)

	if (identical(unname(block), unname(expected_atom1_fastest))) {
		return(list(
			type = "cartesian",
			multiplicity = multiplicity,
			order = "atom1_fastest",
			block_index_order = seq_len(block_size)
		))
	}
	if (identical(unname(block), unname(expected_atom2_fastest))) {
		return(list(
			type = "cartesian",
			multiplicity = multiplicity,
			order = "atom2_fastest",
			block_index_order = match(
				paste(expected_atom1_fastest[, 1], expected_atom1_fastest[, 2], sep = "\r"),
				paste(expected_atom2_fastest[, 1], expected_atom2_fastest[, 2], sep = "\r")
			)
		))
	}

	NULL
}

#' Infer and validate block structure for expanded atom-pair input
#'
#' Expanded atom-pair tables are assumed to be written so that rows belonging
#' to the same relaxation-rate are separated by \eqn{N}, the number of
#' relaxation-rates. This helper tests candidate block sizes against the full
#' table, allowing each `4`, `6`, or `9` row Cartesian block to have its own
#' atom1-fastest or atom2-fastest ordering while also recognizing special
#' internal-permutation blocks of size `2` and `3`.
#'
#' @param atom_pairs data frame or matrix whose first two columns define the
#'   expanded atom pairs
#' @return A list with elements:
#'   \describe{
#'     \item{`n_relax_rates`}{Number of relaxation-rates represented by the
#'       expanded atom-pair table.}
#'     \item{`inferred_order`}{Character scalar giving the within-block
#'       ordering after canonicalization. Supported layouts are always mapped
#'       onto `"atom1_fastest"`.}
#'     \item{`row_order`}{Integer permutation that maps the input rows onto the
#'       canonical atom1-fastest expanded ordering.}
#'     \item{`inferred_multiplicity`}{Named integer vector giving the inferred
#'       multiplicities of the `atom1` and `atom2` sides for kinetic-model
#'       construction.}
#'   }
#'
#' @noRd
.infer_spec_den_relax_blocks <- function(atom_pairs) {
	atom_pair_mat <- as.matrix(atom_pairs[, 1:2, drop = FALSE])
	max_block_size <- min(nrow(atom_pair_mat), 9L)

	for (block_size in c(9L, 6L, 4L, 3L, 2L, 1L)) {
		if (block_size > max_block_size || nrow(atom_pair_mat) %% block_size != 0) {
			next
		}

		n_relax_rates <- nrow(atom_pair_mat) / block_size
		row_index_mat <- matrix(seq_len(nrow(atom_pair_mat)), nrow = n_relax_rates, ncol = block_size)
		block_info <- vector("list", n_relax_rates)
		valid_block_size <- TRUE
		for (i in seq_len(n_relax_rates)) {
			block_info[[i]] <- .classify_spec_den_relax_block(
				atom_pair_mat[row_index_mat[i, ], , drop = FALSE]
			)
			if (is.null(block_info[[i]])) {
				valid_block_size <- FALSE
				break
			}
		}
		if (!valid_block_size) {
			next
		}

		block_type <- vapply(block_info, `[[`, character(1), "type")
		block_order <- vapply(block_info, `[[`, character(1), "order")
		block_multiplicity <- do.call(rbind, lapply(block_info, `[[`, "multiplicity"))

		if (block_size %in% c(4L, 6L, 9L)) {
			if (any(block_type != "cartesian")) {
				next
			}
			if (!all(apply(block_multiplicity, 1, paste, collapse = "-") ==
				apply(block_multiplicity[1, , drop = FALSE], 1, paste, collapse = "-"))) {
				next
			}
			inferred_multiplicity <- block_info[[1]][["multiplicity"]]
		} else if (block_size %in% c(2L, 3L)) {
			if (all(block_type == "internal")) {
				inferred_multiplicity <- c(atom1 = 1L, atom2 = block_size)
			} else {
				cart_idx <- block_type == "cartesian"
				cart_multiplicity <- block_multiplicity[cart_idx, , drop = FALSE]
				if (!nrow(cart_multiplicity)) {
					next
				}
				if (!all(apply(cart_multiplicity, 1, paste, collapse = "-") ==
					apply(cart_multiplicity[1, , drop = FALSE], 1, paste, collapse = "-"))) {
					next
				}
				inferred_multiplicity <- block_info[[which(cart_idx)[1]]][["multiplicity"]]
				if (sum(inferred_multiplicity > 1L) != 1L) {
					next
				}
			}
		} else {
			if (!all(block_type == "cartesian") || !all(block_multiplicity == 1L)) {
				next
			}
			inferred_multiplicity <- c(atom1 = 1L, atom2 = 1L)
		}

		inferred_order <- "atom1_fastest"

		row_order <- unlist(
			lapply(seq_len(n_relax_rates), function(i) row_index_mat[i, block_info[[i]][["block_index_order"]]]),
			use.names = FALSE
		)

		return(list(
			n_relax_rates = n_relax_rates,
			inferred_order = inferred_order,
			row_order = row_order,
			inferred_multiplicity = inferred_multiplicity
		))
	}

	stop("Could not match the expanded atom-pair table to a supported block layout")
}

#' Build spectral-density relaxation data from expanded atom-pair input
#'
#' This function constructs the matrices that
#' implement a user-supplied kinetic model for a given number of ensemble
#' members. It converts expanded atom-pair input, together with the
#' corresponding `relax_data_list`, into the heterogeneous data structure
#' needed for spectral-density relaxation calculations.
#'
#' The function infers the permutation block layout, if any, from the atom-pair
#' row ordering, validates that all blocks follow the same Cartesian-product
#' structure, and constructs the internal-motion grouping and coefficient
#' matrices from `base_rate_mat` and `perm_rates`. The spectral-density term
#' arrays in `relax_data_list` are passed through unchanged.
#'
#' The resulting object is intended to work with any rate definition encoded in
#' `relax_data_list`, including term arrays generated by
#' [make_sigma_spec_den_term_array()], [make_r1_spec_den_term_array()],
#' [make_r2_spec_den_term_array()], or other functions that produce the same
#' `spec_den_term_array` structure.
#'
#' The current implementation assumes that the rows in `atom_pairs` follow the
#' expanded ordering convention where rows belonging to the same
#' relaxation-rate are separated by \eqn{N}, the number of relaxation-rates.
#' Reordering rows breaks this interpretation and causes validation to fail.
#'
#' @param atom_relax_data heterogeneous list with at least `atom_pairs`,
#'   `unit`, and `relax_data_list`
#' @param base_rate_mat transition-rate matrix for ensemble-member exchange
#' @param base_rates named numeric vector of symbolic base-process rates used to
#'   label `lambda_int_coef`
#' @param perm_rates length-two numeric vector describing permutation processes
#'   on the `atom1` and `atom2` sides. `NA` indicates no permutation, while a
#'   non-missing named value supplies the symbolic rate constant for the
#'   inferred permutation on that side.
#'
#' @return A list with elements `atom_pairs`, `unit`, `relax_data_list`,
#'   `groupings`, `a_int_coef`, `lambda_int_coef`, and `inferred_multiplicity`
#'
#' @export
make_spec_den_relax_data <- function(atom_relax_data, base_rate_mat, base_rates, perm_rates = c(NA_real_, NA_real_)) {
	stopifnot(
		is.list(atom_relax_data),
		"atom_pairs" %in% names(atom_relax_data),
		"unit" %in% names(atom_relax_data),
		"relax_data_list" %in% names(atom_relax_data),
		is.matrix(base_rate_mat),
		is.numeric(base_rates),
		!is.null(names(base_rates))
	)
	stopifnot(is.logical(atom_relax_data[["unit"]]), length(atom_relax_data[["unit"]]) == 1)

	perm_info <- .normalize_spec_den_perm_rates(perm_rates)
	block_info <- .infer_spec_den_relax_blocks(atom_relax_data[["atom_pairs"]])

	atom_pairs <- atom_relax_data[["atom_pairs"]][block_info[["row_order"]], , drop = FALSE]
	unit <- atom_relax_data[["unit"]]

	relax_data_list <- atom_relax_data[["relax_data_list"]]
	if (length(relax_data_list) == 0) {
		stop("`make_spec_den_relax_data()` requires `atom_relax_data[['relax_data_list']]`")
	}
	stopifnot(is.list(relax_data_list))
	for (rate_name in names(relax_data_list)) {
		rate_obj <- relax_data_list[[rate_name]]
		stopifnot(
			is.list(rate_obj),
			"value" %in% names(rate_obj),
			"spec_den_term_array" %in% names(rate_obj)
		)
		if (length(rate_obj[["value"]]) != block_info[["n_relax_rates"]]) {
			stop(
				"`relax_data_list[['", rate_name, "']][['value']]` must have length ",
				block_info[["n_relax_rates"]],
				", matching the number of relaxation-rates"
			)
		}
		if (dim(rate_obj[["spec_den_term_array"]])[1] != block_info[["n_relax_rates"]]) {
			stop(
				"`relax_data_list[['", rate_name, "']][['spec_den_term_array']]` must have ",
				block_info[["n_relax_rates"]],
				" rows in its first dimension, matching the number of relaxation-rates"
			)
		}
		if ("k" %in% names(rate_obj)) {
			if (!is.numeric(rate_obj[["k"]])) {
				stop("`relax_data_list[['", rate_name, "']][['k']]` must be numeric when supplied")
			}
			if (!length(rate_obj[["k"]]) %in% c(1, block_info[["n_relax_rates"]])) {
				stop(
					"`relax_data_list[['", rate_name, "']][['k']]` must have length 1 or ",
					block_info[["n_relax_rates"]],
					", matching the number of relaxation-rates"
				)
			}
		}
	}

	expected_has_perm <- block_info[["inferred_multiplicity"]] > 1L
	declared_has_perm <- !is.na(perm_info[["rate_values"]])
	if (!identical(unname(expected_has_perm), unname(declared_has_perm))) {
		stop(
			"Inferred permutation structure (",
			paste(block_info[["inferred_multiplicity"]], collapse = "-"),
			") does not match which sides have non-missing `perm_rates`"
		)
	}
	if (any(block_info[["inferred_multiplicity"]][declared_has_perm] == 1L)) {
		stop("Permutation rate constants may only be supplied for sides whose inferred multiplicity is greater than 1")
	}
	perm_rate_values_named <- .unique_spec_den_perm_rates(perm_info)

	make_perm_matrix <- function(side_label) {
		mult <- block_info[["inferred_multiplicity"]][[side_label]]
		rate_value <- perm_info[["rate_values"]][[side_label]]
		if (is.na(rate_value)) {
			return(NULL)
		}
		rate_mat_simple(rate_value, paste0(side_label, seq_len(mult)))
	}

	left_perm_mat <- make_perm_matrix("atom1")
	right_perm_mat <- make_perm_matrix("atom2")

	rate_data <- get_rate_data(base_rate_mat, validate = TRUE)
	if (!is.null(right_perm_mat)) {
		rate_data <- get_rate_data(right_perm_mat, rate_data, validate = TRUE)
	}
	if (!is.null(left_perm_mat)) {
		rate_data <- get_rate_data(left_perm_mat, rate_data, validate = TRUE)
	}

	all_permutations <- list()
	if (!is.na(perm_info[["rate_values"]][["atom1"]])) {
		all_permutations <- c(
			all_permutations,
			list(matrix(0, nrow = 1, ncol = block_info[["inferred_multiplicity"]][["atom1"]]))
		)
	}
	if (!is.na(perm_info[["rate_values"]][["atom2"]])) {
		all_permutations <- c(
			all_permutations,
			list(matrix(0, nrow = 1, ncol = block_info[["inferred_multiplicity"]][["atom2"]]))
		)
	}
	if (length(all_permutations)) {
		all_perm_names <- c(
			if (!is.na(perm_info[["rate_values"]][["atom1"]])) as.character(perm_info[["rate_values"]][["atom1"]]),
			if (!is.na(perm_info[["rate_values"]][["atom2"]])) as.character(perm_info[["rate_values"]][["atom2"]])
		)
		names(all_permutations) <- all_perm_names
	} else {
		all_permutations <- list()
	}

	rate_count_mat <- get_rate_count_mat(base_rate_mat, all_permutations)
	if (ncol(rate_count_mat) == 0) {
		mat_list <- list(
			groupings = matrix(1L, nrow = 1, ncol = 1),
			a_coef = matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "0")),
			lambda_coef = matrix(numeric(), nrow = 0, ncol = 1, dimnames = list(NULL, "0"))
		)
	} else {
		mat_list <- rate_data_to_mat_list(rate_data, rate_count_mat, c(base_rates, perm_rate_values_named))
	}

	list(
		atom_pairs = atom_pairs,
		unit = unit,
		relax_data_list = relax_data_list,
		groupings = mat_list[["groupings"]],
		a_int_coef = mat_list[["a_coef"]],
		lambda_int_coef = mat_list[["lambda_coef"]],
		inferred_multiplicity = block_info[["inferred_multiplicity"]]
	)
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
		r_mat_check <- matrix(stats::rnorm(30), ncol=3)
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

	eros3_sub_g <- coord_array_to_g_matrix(eros3_sub_coord, eros3_sub_atom_pairs, eros3_grouping_list)

	# get g values using M1 twice
	eros3_sub_1_g <- coord_array_to_g_matrix(eros3_sub_coord[,,c(1,1)], eros3_sub_atom_pairs, eros3_grouping_list)

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
