#' Convert PDB formatted atom lines into a matrix
#'
#' Generate a matrix containing coordinates in a set of PDB atom lines
#'
#' The column names in the resulting matrix correspond to columns 13 through 27 of the 
#' ATOM lines. The definitions of the PDB columns and corresponding colnames columns are
#' as follows:
#' \tabular{lllll}{
#'   PDB \tab colnames \tab Data Type \tab Field \tab Definition \cr
#'   13-16 \tab 1-4 \tab Atom \tab name \tab Atom name \cr
#'   17 \tab 5 \tab Character \tab altLoc \tab Alternate location indicator \cr
#'   18-20 \tab 6-8 \tab Residue name \tab resName \tab Residue name \cr
#'   22 \tab 10 \tab Character \tab chainID \tab Chain identifier \cr
#'   23-26 \tab 11-4 \tab Integer \tab resSeq \tab Residue sequence number \cr
#'   27 \tab 15 \tab AChar \tab iCode \tab Code for insertion of residues
#' }
#'
#' See \url{http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM}
#' for more information.
#'
#' @param atom_lines character vector of ATOM lines
#'
#' @return 3xN matrix with a column for every atom
#'
#' @export
atom_lines_to_coord <- function(atom_lines) {

	coord <- rbind(
		as.numeric(substr(atom_lines, 31, 38)),
		as.numeric(substr(atom_lines, 39, 46)),
		as.numeric(substr(atom_lines, 47, 54))
	)
	
	colnames(coord) <- substr(atom_lines, 13, 27)
	
	coord
}

#' Convert a matrix of coordinates into PDB formatted atom lines
#'
#' Make PDB ATOM lines from a 3xN matrix of coordinates
#'
#' The names of the columns are expected to correspond to columns 13 to 27 of the PDB
#' atom records.
#'
#' @param coord 3xN matrix with a column for every atom
#'
#' @return a character vector with the PDB ATOM lines
#'
#' @export
coord_to_atom_lines <- function(coord) {

	sprintf("ATOM  %5i %s   %8.3f%8.3f%8.3f  1.00  0.00", seq_len(ncol(coord)), colnames(coord), coord[1,], coord[2,], coord[3,])
}

#' Read an ensemble of PDB files
#'
#' @param pdb_files character vector paths to PDB format file(s)
#' @param model_idx vector with indices of models to read
#' @param proton_only logical indicating whether to only read hydrogen ATOM records
#'
#' @export
read_ensemble <- function(pdb_files, model_idx=NULL, proton_only=FALSE) {

	coord_list <- lapply(pdb_files, function(pdb_file) {
	
		pdb_lines <- readLines(pdb_file)
		model_line_idx <- grep("^MODEL", pdb_lines)
		
		if (length(model_line_idx)) {
			start_line_idx <- model_line_idx+1
			end_line_idx <- c(model_line_idx[-1]-1, length(pdb_lines))
			model_idx_list <- lapply(seq_along(start_line_idx), function(i) seq(start_line_idx[i], end_line_idx[i]))
		} else {
			model_idx_list <- list(seq_along(pdb_lines))
		}
		
		if (!is.null(model_idx)) {
			model_idx_list <- model_idx_list[model_idx]
		}
		
		coord_model_list <- lapply(model_idx_list, function(idx) {
			atom_lines <- grep("^ATOM  ..... ....[ A]", pdb_lines[idx], value=TRUE)
			coord <- atom_lines_to_coord(atom_lines)
			substr(colnames(coord), 5, 5) <- " "
			
			coord
		})
		
		coord_model_list
	})
	
	if (is.null(model_idx)) {
		model_idx <- seq_along(coord_list[[1]])
	}
	
	coord_list <- unlist(coord_list, recursive=FALSE)
	
	unique_atoms <- colnames(coord_list[[1]])
	
	for (i in seq_along(coord_list)[-1]) {
		unique_atoms <- intersect(unique_atoms, colnames(coord_list[[i]]))
	}
	
	#unique_atoms <- intersect(unique_atoms, c(exchangable_atomids(unique_atoms), n_atomids))
	
	if (proton_only) {
		unique_atoms <- grep("^(H|.H)", unique_atoms, value=TRUE)
	}
	
	for (i in seq_along(coord_list)) {
		coord_list[[i]] <- coord_list[[i]][,unique_atoms]
	}
	
	coord_array <- simplify2array(coord_list)
	
	if (dim(coord_array)[3] == length(pdb_files)) {
	
		dimnames(coord_array)[[3]] <- pdb_files
	
	} else {
	
		coord_array <- array(
			coord_array,
			dim=c(dim(coord_array)[1:2], dim(coord_array)[3]/length(pdb_files), length(pdb_files)),
			dimnames=c(dimnames(coord_array)[1:2], list(paste0("M", model_idx), pdb_files))
		)
	}
	
	coord_array
}

#' Check analytical derivatives via finite difference approximation
#'
#' @param func function to evaluate
#' @param value first argument to pass to `func`
#' @param dv delta used to increment elements of `value`
#' @param vdims vector of dimensions that must be incremented separately
#' @param gdims vector of dimensions in gradient corresponding to `vdims`
#' @param ... additional arguments passed to `func`
#'
#' @return difference between `gradient` attribute returned by `func` and finite difference approximation
#'
#' @export
deriv_check <- function(func, value, dv, vdims, gdims, ...) {

	result <- func(value, ...)
	
	# gradient calculated analytically to check
	gradient <- attr(result, "gradient")
	attr(result, "gradient") <- NULL
	
	# gradient calculated using finite difference
	gradient_fd <- gradient
	gradient_fd[] <- NA

	stopifnot(dim(value)[vdims] == dim(gradient)[gdims])

	idx_mat <- as.matrix(expand.grid(lapply(vdims, function(i) seq_len(dim(value)[i]))))
	v_list <- lapply(seq_along(dim(value)), function(i) seq_len(dim(value)[i]))
	g_list <- lapply(seq_along(dim(gradient)), function(i) seq_len(dim(gradient)[i]))

	for (i in seq_len(nrow(idx_mat))) {
	
		v_list[vdims] <- as.list(idx_mat[i,])
		v_idx <- as.matrix(expand.grid(v_list))
		g_list[gdims] <- as.list(idx_mat[i,])
		g_idx <- as.matrix(expand.grid(g_list))
		
		new_value <- value
		new_value[v_idx] <- new_value[v_idx] + dv
		
		new_result <- func(new_value, ..., gradient=FALSE)
		
		gradient_fd[g_idx] <- (new_result-result)/dv
	}
	
	gradient-gradient_fd
}

#' Calculate internuclear vectors from atomic coordinates
#'
#' @param coord_array 3D array (atoms, xyz, models) with atomic coordinates
#' @param atom_pairs matrix with each row having the names or indices of an atom pair (first dimension in `coord_array`)
#'
#' @return 3D array (pairs, models, xyz) with internuclear vectors. Atom pair names follow the format `resSeq:Atom-resSeq:Atom`.
#;
#' @export
coord_array_to_r_array <- function(coord_array, atom_pairs) {

	r_array <- aperm(coord_array[atom_pairs[,2],,]-coord_array[atom_pairs[,1],,], c(1,3,2))
	
	if (!is.null(dimnames(coord_array))) {
		atom_names <- sub(" *([^ ]+) *([A-Z]{3}) (.) +([^ ]+) *", "\\4:\\1", dimnames(coord_array)[[1]])
		names(atom_names) <- dimnames(coord_array)[[1]]
		dimnames(r_array)[[1]] <- paste(atom_names[atom_pairs[,1]], atom_names[atom_pairs[,2]], sep="-")
		dimnames(r_array)[[3]] <- c("x", "y", "z")
	}
	
	r_array
}

#' Calculate dipole-dipole interaction tensors from internuclear vectors
#'
#' @details
#' Given an internuclear vector with components (`x`, `y`, `z`) and distance (`r`), the vector form (`d`) of the dipole-dipole interaction tensor is defined:
#'
#' \eqn{d = 1/r^5 [ z^2 - 1/2 (x^2+y^2), \sqrt{3}/2 (x^2-y^2), \sqrt{3} xz, \sqrt{3} yz, \sqrt{3} xy) ]}
#'
#' Note that the coefficients of the last four elements are different than what was given in \href{https://doi.org/10.1007/s10858-019-00288-8}{Smith 2020}. Those published coefficients do not account for the 3/2 factor in the cartesian tensor.
#'
#' @param r_array 3D array (pairs, models, xyz) with internuclear vectors
#' @param gradient a logical value indicating whether to calculate the derivative
#'
#' @return 3D array (pairs, models, tensor elements) with interaction tensors
#'
#' The optional derivative is contained in the `"gradient"` attribute. It is a 4D array (pairs, models, tensor elements, xyz).
#'
#' @export
r_array_to_d_array <- function(r_array, gradient=FALSE) {

	input_dim <- dim(r_array)
	input_dimnames <- dimnames(r_array)
	
	# convert to matrix if necessary
	if (length(dim(r_array)) == 3) {
		dim(r_array) <- c(prod(dim(r_array)[1:2]), dim(r_array)[3])
	}

	x <- r_array[,1]
	y <- r_array[,2]
	z <- r_array[,3]
	sr3x <- 1.73205080756888*x #sqrt(3)
	sr3y <- 1.73205080756888*y #sqrt(3)
	xsq <- x^2
	ysq <- y^2
	zsq <- z^2
	xsq_ysq <- xsq+ysq
	r <- sqrt(xsq_ysq+zsq)
	da <- 1/r^5
	db <- cbind(
		zsq-0.5*(xsq_ysq),
		0.866025403784439*(xsq-ysq), # sqrt(3)/2
		sr3x*z,
		sr3y*z,
		sr3x*y
	)
	
	value <- da*db
	colnames(value) <- c("d1", "d2", "d3", "d4", "d5")
	
	# convert back to array if necessary
	if (length(input_dim) == 3) {
		dim(value) <- c(input_dim[1:2], 5)
		dimnames(value) <- c(input_dimnames[1:2], list(c("d1", "d2", "d3", "d4", "d5")))
	}
	
	if (gradient) {
	
		# it is numerically more efficient to calculate the derivative with the product rule
		r7inv <- 1/r^7
		#ddadxyz <- c(-5*x/r7, -5*y/r7, -5*z/r7)
		ddadxyz <- -5*r_array*r7inv
		sr3z <- 1.73205080756888*z
		zero <- numeric(nrow(r_array))
		ddbdxyz <- array(c(
			-x, sr3x, sr3z, zero, sr3y,
			-y, -sr3y, zero, sr3z, sr3x,
			2*z, zero, sr3x, sr3y, zero
		), dim=c(nrow(r_array), 5, 3))
		
		# calculate derivative using product rule: db*ddadxyz + da*ddbdxyz
		# da could be premultiplied during construction of ddbdxyz
		grad <- as.vector(db)*as.vector(ddadxyz[,c(1L,1L,1L,1L,1L,2L,2L,2L,2L,2L,3L,3L,3L,3L,3L)]) + da*ddbdxyz
		
		# convert to array if necessary
		if (length(input_dim) == 3) {
			dim(grad) <- c(input_dim[1:2], 5, input_dim[3])
			dimnames(grad) <- c(input_dimnames[1:2], list(c("d1", "d2", "d3", "d4", "d5")), input_dimnames[3])
		}
		
		attr(value, "gradient") <- grad
	}
	
	value
}

#' Calculate group norm squared from dipole-dipole interaction tensors
#'
#' @param d_array 3D array (pairs, models, tensor elements) with interaction tensors
#' @param grouping list of integer vectors giving groupings of models to average interaction tensors
#' @param gradient a logical value indicating whether to calculate the derivative
#'
#' @return vector with norm squared for each atom pair. 
#'
#' The optional derivative is contained in the `"gradient"` attribute. It is a 3D array (pairs, models, tensor elements).
#'
#' @export
d_array_to_g <- function(d_array, grouping, gradient=FALSE) {

	value <- numeric(dim(d_array)[1])
	
	if (gradient) {
	
		grad <- d_array
		grad[] <- 0
		attr(grad, "gradient") <- NULL
	}
	
	for (i in seq_along(grouping)) {
		
		# average dipole interaction tensor
		d_matrix <- 0
		
		# sum the dipole interaction tensors within each group
		for (j in seq_along(grouping[[i]])) {
			d_matrix <- d_matrix + d_array[,grouping[[i]][[j]],,drop=FALSE]
		}
		
		if (gradient) {
		
			# D[((dself + dother)/grouplength)^2*grouplength/totallength, {dself}]
			d_matrix_grad <- d_matrix*(2/dim(d_array)[2]/length(grouping[[i]]))
			# the gradients are the same for all members of the group because they contribute equally
			for (j in seq_along(grouping[[i]])) {
				grad[,grouping[[i]][[j]],] <- d_matrix_grad
			}
		}
		
		# calculate the mean by dividing by the length
		d_matrix <- d_matrix*(1/length(grouping[[i]]))
		
		# calculate self dot product (norm squared) and accumulate group's contribution to mean g
		value <- value + rowSums(d_matrix^2)*(length(grouping[[i]])/dim(d_array)[2])
	}
	
	if (gradient) {
	
		attr(value, "gradient") <- grad
	}
	
	value
}

#' Calculate restraint energy from group norm squared values
#'
#' @param g current group norm squared values
#' @param g0 target group norm squared values
#' @param k force constant
#' @param gradient a logical value indicating whether to calculate the derivative
#'
#' @return restraint energy calculated using \eqn{k*(g-g0)^2}
#'
#' @export
g_to_energy <- function (g, g0, k=1, gradient=FALSE) {

    expr1 <- g - g0
    value <- k * expr1^2
    if (gradient) {
		grad <- k * (2 * expr1)
		attr(value, "gradient") <- grad
    } else {
    	attr(value, "gradient") <- NULL
    }
    value
}

#' Calculate group norm squared values from atomic coordinates
#'
#' @param coord_array 3D array (atoms, xyz, models) with atomic coordinates
#' @param atom_pairs matrix with each row having the names or indices of an atom pair (first dimension in `coord_array`)
#' @param grouping_list list of lists of integer vectors giving groupings of models to average interaction tensors
#'
#' @return matrix (pairs, groupings) with group norm squared values
#'
#' @export
coord_array_to_g <- function(coord_array, atom_pairs, grouping_list) {

	r_array <- coord_array_to_r_array(coord_array, atom_pairs)
	
	d_array <- r_array_to_d_array(r_array, gradient=FALSE)
	
	g_list <- lapply(grouping_list, function(grouping) d_array_to_g(d_array, grouping, gradient=FALSE))
	
	simplify2array(g_list)
}

#' Calculate restraint energy from atomic coordinates
#'
#' @param coord_array 3D array (atoms, xyz, models) with atomic coordinates
#' @param atom_pairs matrix with each row having the names or indices of an atom pair (first dimension in `coord_array`)
#' @param grouping_list list of lists of integer vectors giving groupings of models to average interaction tensors
#' @param g0 target group norm squared values
#' @param k force constant
#' @param gradient a logical value indicating whether to calculate the derivative
#'
#' @return total restraint energy calculated using \eqn{k*(g-g0)^2}
#'
#' The optional derivative is contained in the `"gradient"` attribute. It is a 3D array (atoms, xyz, models).
#'
#' @export
coord_array_to_energy <- function(coord_array, atom_pairs, grouping_list, g0, k, gradient=FALSE) {

	# calculate internuclear vectors
	r_array <- coord_array_to_r_array(coord_array, atom_pairs)
	
	# calculate dipole-dipole interaction tensors
	d_array <- r_array_to_d_array(r_array, gradient=gradient)
	
	# calculate norm squared for different groupings of dipole-dipole interaction tensors
	g_list <- lapply(grouping_list, function(grouping) d_array_to_g(d_array, grouping, gradient=gradient))
	
	# convert the list above into a matrix (pairs, groupings)
	g_matrix <- simplify2array(g_list)
	
	# calculate energies from the norm squared values
	energy_matrix <- g_to_energy(g_matrix, g0, k, gradient=gradient)
	
	# return the sum of all the individual restraint energies
	value <- sum(energy_matrix)
	
	if (gradient) {
	
		# calculate de/dd = dg/dd * de/dg for all individual interaction tensor components
		d_energy_d_d_array <- 0
		# sum the contributions from the different norm squared values
		for (i in seq_along(g_list)) {
			d_energy_d_d_array <- d_energy_d_d_array + attr(g_list[[i]], "gradient")*attr(energy_matrix, "gradient")[,i]
		}
		
		# calculate de/dr = dd/dr * de/dd for each xyz component of the internuclear vectors
		# coerce d_energy_d_d_array to a vector so it gets repeated over dd/dr
		d_energy_d_r_array_all <- attr(d_array, "gradient")*as.vector(d_energy_d_d_array)
		# sum the individual interaction tensor component derivatives associated with x, y, and z
		d_energy_d_r_array <- apply(d_energy_d_r_array_all, c(1,2,4), sum)
		
		# initialize empty gradient with dimensions equal to input coordinates
		grad <- coord_array
		grad[] <- 0
		
		# propagate the internuclear vector derivatives back onto the atomic coordinates
		for (i in seq_len(dim(d_energy_d_r_array)[1])) {
			
			pair_grad <- t(d_energy_d_r_array[i,,])
			# first component is subtracted so derivative is too
			grad[atom_pairs[i,1],,] <- grad[atom_pairs[i,1],,] - pair_grad
			# second component is added so derivative is too
			grad[atom_pairs[i,2],,] <- grad[atom_pairs[i,2],,] + pair_grad
		}
		
		attr(value, "gradient") <- grad
	}
	
	value
}
