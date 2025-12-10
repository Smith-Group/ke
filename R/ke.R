#' Convert between one and three letter residue names 
#'
#' @docType data
#' @rdname residue_names
#' @examples
#' aa1_to_aa3["W"]
#' aa1_to_aa3[c("W", "F", "Y")]
#'
#' @export
aa1_to_aa3 <- c("TRP", "PHE", "TYR", "MET", "LEU", "ILE", "VAL", "ALA", "GLY", "SER", "THR", "ARG", "LYS", "HIS", "ASN", "GLN", "ASP", "GLU", "PRO", "CYS")

#' @docType data
#' @rdname residue_names
#' @examples
#' aa3_to_aa1["TRP"]
#' aa3_to_aa1[c("TRP", "PHE", "TYR")]
#'
#' @export
aa3_to_aa1 <- c("W", "F", "Y", "M", "L", "I", "V", "A", "G", "S", "T", "R", "K", "H", "N", "Q", "D", "E", "P", "C")

names(aa1_to_aa3) <- aa3_to_aa1
names(aa3_to_aa1) <- aa1_to_aa3

#' Convert between proton atom/residue name pairs
#'
#' @rdname atom_residue_names
#' @export
amber_to_pdb <- c(" H   ALA", " HA  ALA", " HB1 ALA", " HB2 ALA", " HB3 ALA", " H   ARG", " HA  ARG", " HB2 ARG", " HB3 ARG", " HD2 ARG", " HD3 ARG", " HE  ARG", " HG2 ARG", " HG3 ARG", "HH11 ARG", "HH12 ARG", "HH21 ARG", "HH22 ARG", " H   ASN", " HA  ASN", " HB2 ASN", " HB3 ASN", "HD21 ASN", "HD22 ASN", " H   ASP", " HA  ASP", " HB2 ASP", " HB3 ASP", " H   CYS", " HA  CYS", " HB2 CYS", " HB3 CYS", " HG  CYS", " H   GLN", " HA  GLN", " HB2 GLN", " HB3 GLN", " HG2 GLN", " HG3 GLN", "HE21 GLN", "HE22 GLN", " H   GLU", " HA  GLU", " HB2 GLU", " HB3 GLU", " HG2 GLU", " HG3 GLU", " H   GLY", " HA2 GLY", " HA3 GLY", " H   HIS", " HA  HIS", " HB2 HIS", " HB3 HIS", " HD1 HIS", " HD2 HIS", " HE1 HIS", " HE2 HIS", " H   ILE", " HA  ILE", " HB  ILE", "HD11 ILE", "HD12 ILE", "HD13 ILE", "HG12 ILE", "HG13 ILE", "HG21 ILE", "HG22 ILE", "HG23 ILE", " H   LEU", " HA  LEU", " HB2 LEU", " HB3 LEU", " HG  LEU", "HD11 LEU", "HD12 LEU", "HD13 LEU", "HD21 LEU", "HD22 LEU", "HD23 LEU", " H   LYS", " HA  LYS", " HB2 LYS", " HB3 LYS", " HD2 LYS", " HD3 LYS", " HE2 LYS", " HE3 LYS", " HG2 LYS", " HG3 LYS", " HZ1 LYS", " HZ2 LYS", " HZ3 LYS", " H1  MET", " H2  MET", " H3  MET", " HA  MET", " HB2 MET", " HB3 MET", " HE1 MET", " HE2 MET", " HE3 MET", " HG2 MET", " HG3 MET", " H   PHE", " HA  PHE", " HB2 PHE", " HB3 PHE", " HD1 PHE", " HD2 PHE", " HE1 PHE", " HE2 PHE", " HZ  PHE", " HA  PRO", " HB2 PRO", " HB3 PRO", " HD2 PRO", " HD3 PRO", " HG2 PRO", " HG3 PRO", " H   SER", " HA  SER", " HB2 SER", " HB3 SER", " HG  SER", " H   THR", " HA  THR", " HB  THR", " HG1 THR", "HG21 THR", "HG22 THR", "HG23 THR", " H   TRP", " HA  TRP", " HB2 TRP", " HB3 TRP", " HD1 TRP", " HE1 TRP", " HE3 TRP", " HZ2 TRP", " HZ3 TRP", " HH2 TRP", " H   TYR", " HA  TYR", " HB2 TYR", " HB3 TYR", " HD1 TYR", " HD2 TYR", " HE1 TYR", " HE2 TYR", " HH  TYR", " H   VAL", " HA  VAL", " HB  VAL", "HG11 VAL", "HG12 VAL", "HG13 VAL", "HG21 VAL", "HG22 VAL", "HG23 VAL")

#' @rdname atom_residue_names
#' @export
pdb_to_amber <- c(" H   ALA", " HA  ALA", " HB1 ALA", " HB2 ALA", " HB3 ALA", " H   ARG", " HA  ARG", " HB1 ARG", " HB2 ARG", " HD1 ARG", " HD2 ARG", " HE  ARG", " HG1 ARG", " HG2 ARG", "1HH1 ARG", "2HH1 ARG", "1HH2 ARG", "2HH2 ARG", " H   ASN", " HA  ASN", " HB1 ASN", " HB2 ASN", "1HD2 ASN", "2HD2 ASN", " H   ASP", " HA  ASP", " HB1 ASP", " HB2 ASP", " H   CYS", " HA  CYS", " HB1 CYS", " HB2 CYS", " HG  CYS", " H   GLN", " HA  GLN", " HB1 GLN", " HB2 GLN", " HG1 GLN", " HG2 GLN", "1HE2 GLN", "2HE2 GLN", " H   GLU", " HA  GLU", " HB1 GLU", " HB2 GLU", " HG1 GLU", " HG2 GLU", " H   GLY", " HA1 GLY", " HA2 GLY", " H   HIS", " HA  HIS", " HB1 HIS", " HB2 HIS", " HD1 HIS", " HD2 HIS", " HE1 HIS", " HE2 HIS", " H   ILE", " HA  ILE", " HB  ILE", " HD1 ILE", " HD2 ILE", " HD3 ILE", "1HG1 ILE", "2HG1 ILE", "1HG2 ILE", "2HG2 ILE", "3HG2 ILE", " H   LEU", " HA  LEU", " HB1 LEU", " HB2 LEU", " HG  LEU", "1HD1 LEU", "2HD1 LEU", "3HD1 LEU", "1HD2 LEU", "2HD2 LEU", "3HD2 LEU", " H   LYS", " HA  LYS", " HB1 LYS", " HB2 LYS", " HD1 LYS", " HD2 LYS", " HE1 LYS", " HE2 LYS", " HG1 LYS", " HG2 LYS", " HZ1 LYS", " HZ2 LYS", " HZ3 LYS", " H1  MET", " H2  MET", " H3  MET", " HA  MET", " HB1 MET", " HB2 MET", " HE1 MET", " HE2 MET", " HE3 MET", " HG1 MET", " HG2 MET", " H   PHE", " HA  PHE", " HB1 PHE", " HB2 PHE", " HD1 PHE", " HD2 PHE", " HE1 PHE", " HE2 PHE", " HZ  PHE", " HA  PRO", " HB1 PRO", " HB2 PRO", " HD1 PRO", " HD2 PRO", " HG1 PRO", " HG2 PRO", " H   SER", " HA  SER", " HB1 SER", " HB2 SER", " HG  SER", " H   THR", " HA  THR", " HB  THR", " HG1 THR", "1HG2 THR", "2HG2 THR", "3HG2 THR", " H   TRP", " HA  TRP", " HB1 TRP", " HB2 TRP", " HD1 TRP", " HE1 TRP", " HE3 TRP", " HZ2 TRP", " HZ3 TRP", " HH2 TRP", " H   TYR", " HA  TYR", " HB1 TYR", " HB2 TYR", " HD1 TYR", " HD2 TYR", " HE1 TYR", " HE2 TYR", " HH  TYR", " H   VAL", " HA  VAL", " HB  VAL", "1HG1 VAL", "2HG1 VAL", "3HG1 VAL", "1HG2 VAL", "2HG2 VAL", "3HG2 VAL")

names(amber_to_pdb)	<- pdb_to_amber
names(pdb_to_amber) <- amber_to_pdb

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

#' Get or set four-character atom names
#'
#' @param coord matrix or array with second dimension having columns 13 to 27 of the PDB ATOM records
#'
#' @rdname coord_atomnames
#' @export
coord_atomnames <- function(coord) {

	if (is.array(coord) && is.character(dimnames(coord)[[2]])) {
		substr(dimnames(coord)[[2]], 1, 4)
	} else if (is.character(coord)) {
		substr(coord, 1, 4)
	} else {
		stop("coord not named array or character")
	}
}

#' @param value character vector with four characters per element
#'
#' @rdname coord_atomnames
#' @export
`coord_atomnames<-` <- function(coord, value) {

	stopifnot(nchar(value) == 4)

	if (is.array(coord) && is.character(dimnames(coord)[[2]])) {
		substr(dimnames(coord)[[2]], 1, 4) <- value
	} else if (is.character(coord)) {
		substr(coord, 1, 4) <- value
	} else {
		stop("coord not named array or character")
	}
	
	coord
}

#' Get or set three-character residue names
#'
#' @param coord matrix or array with second dimension having columns 13 to 27 of the PDB res records
#'
#' @rdname coord_resnames
#' @export
coord_resnames <- function(coord) {

	substr(dimnames(coord)[[2]], 6, 8)
}

#' @param value character vector with three characters per element
#'
#' @rdname coord_resnames
#' @export
`coord_resnames<-` <- function(coord, value) {

	stopifnot(nchar(value) == 4)

	substr(dimnames(coord)[[2]], 6, 8) <- value
	
	coord
}

#' Get or set atom/residue name pairs
#'
#' @param coord matrix or array with second dimension having columns 13 to 27 of the PDB ATOM records
#'
#' @rdname coord_atomresnames
#' @export
coord_atomresnames <- function(coord) {

	atomresnames <- substr(dimnames(coord)[[2]], 1, 8)
	substr(atomresnames, 5, 5) <- " "
	
	atomresnames
}

#' @param value character vector with eight characters per element (1-4: atom name, 6-8: residue name)
#'
#' @rdname coord_atomresnames
#' @export
`coord_atomresnames<-` <- function(coord, value) {

	stopifnot(nchar(value) == 4)

	substr(dimnames(coord)[[2]], 1, 4) <- substr(value, 1, 4)
	substr(dimnames(coord)[[2]], 6, 8) <- substr(value, 6, 8)
	
	coord
}

#' Get or set residue sequence number
#'
#' @param coord matrix or array with second dimension having columns 13 to 27 of the PDB ATOM records
#'
#' @rdname coord_resseq
#' @export
coord_resseq <- function(coord) {

	if (is.array(coord) && is.character(dimnames(coord)[[2]])) {
		as.integer(substr(dimnames(coord)[[2]], 11, 14))
	} else if (is.character(coord)) {
		as.integer(substr(coord, 11, 14))
	} else {
		stop("coord not named array or character")
	}
}

#' @param value integer or character vector with no more than four digits
#'
#' @rdname coord_resseq
#' @export
`coord_resseq<-` <- function(coord, value) {

	value <- as.character(value)
	stopifnot(nchar(value) <= 4)

	if (is.array(coord) && is.character(dimnames(coord)[[2]])) {
		substr(dimnames(coord)[[2]], 11, 14) <- sprintf("%-4s", value)
	} else if (is.character(coord)) {
		substr(coord, 11, 14) <- sprintf("%-4s", value)
	} else {
		stop("coord not named array or character")
	}
	
	coord
}

#' Determine whether an atom is rapidly exchangeable at neutral pH
#'
#' @param coord matrix or array with second dimension having columns 13 to 27 of the PDB ATOM records
#' @param amber logical indicating whether atom names are from AMBER force field
#'
#' @export
coord_rapidly_exchangeable <- function(coord, amber=FALSE) {

	rapidly_exchangeable_re <- c(
		"^ HH  TYR", # tyrosine hydroxyl
		"^ HZ[123] LYS", # lysine amino group
		"^ HE  ARG", # arginine epsilon NH
		"^HH[12][12] ARG", # arginine eta NH2
		"^ HG  SER", # serine hydroxyl
		"^ HG1 THR", # threonine hydroxyl
		"^ HD1 HIS", # histidine delta amide hydrogen
		"^ HE2 HIS", # histidine epsilon amide hydrogen
		"^ H[123]  ...", # N-terminal amino group
		"^ HG  CYS" # cysteine hydrogen
	)
	
	if (amber) {
		rapidly_exchangeable_re <- c(
			"^ HH  TYR", # tyrosine hydroxyl
			"^ HZ[123] LYS", # lysine amino group
			"^ HE  ARG", # arginine epsilon NH
			"^[12]HH[12] ARG", # arginine eta NH2
			"^ HG  SER", # serine hydroxyl
			"^ HG1 THR", # threonine hydroxyl
			"^ HD1 HIS", # histidine delta amide hydrogen
			"^ HE2 HIS", # histidine epsilon amide hydrogen
			"^ H[123]  ..." # N-terminal amino group
		)
	}
	
	match_mat <- sapply(rapidly_exchangeable_re, grepl, dimnames(coord)[[2]])
	
	apply(match_mat, 1, any)
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
	
	#unique_atoms <- intersect(unique_atoms, c(exchangeable_atomids(unique_atoms), n_atomids))
	
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
		
		if (dim(coord_array)[4] == 1) {
			
			new_dim <- dim(coord_array)[-4]
			new_dimnames <- dimnames(coord_array)[-4]
			dim(coord_array) <- new_dim
			dimnames(coord_array) <- new_dimnames
		}
	}
	
	coord_array
}

#' Select groups of protons that are typically detectable with unique chemical shifts
#'
#' @param coord_mat 3xN matrix with a column for every atom
#' @param alpha_group logical indicating whether to group aliphatic alpha protons
#' @param beta_group logical indicating whether to group aliphatic beta protons
#' @param gamma_group logical indicating whether to group aliphatic gamma protons
#' @param delta_group logical indicating whether to group aliphatic delta protons
#' @param epsilon_group logical indicating whether to group aliphatic epsilon protons
#' @param methyl_group logical indicating whether to group methyl protons
#' @param aromatic_group logical indicating whether to group phenylalanine/tyrosine delta/epsilon protons
#' @param amine_group logical indicating whether to group asparagine/glutamine amine protons
#' @param exclude_exchangeable logical indicating whether rapidly exchangeable protons should be excluded
#' @param include character vector of exchangeable protons to be explicitly included
#' @param amber logical indicating if names are from Amber force field
#'
#' @return Vector with selected protons. The names give the group the proton belongs to. The name will be the heavy atom name if all of the protons bound to that heavy atom are in the same group. Otherwise, the name will be the proton name.
#'
#' @export
coord_proton_groups <- function(coord_mat, alpha_group=FALSE, beta_group=FALSE, gamma_group=FALSE, delta_group=FALSE, epsilon_group=FALSE, amine_group=FALSE, methyl_group=FALSE, aromatic_group=TRUE, exclude_exchangeable=TRUE, include=character(), amber=FALSE) {

	heavy_idx <- substr(colnames(coord_mat), 2, 2) %in% c("C","N","O","S")
	
	atom_dist <- as.matrix(stats::dist(t(coord_mat)))
	
	group_names <- colnames(coord_mat)[heavy_idx][apply(atom_dist[!heavy_idx,heavy_idx], 1, which.min)]
	names(group_names) <- colnames(coord_mat)[!heavy_idx]
	
	if (exclude_exchangeable) {
		exclude_idx <- which(coord_rapidly_exchangeable(coord_mat[,!heavy_idx], amber))
		exclude_idx <- setdiff(exclude_idx, which(names(group_names) %in% include))
		if (length(exclude_idx)) {
			group_names <- group_names[-exclude_idx]
		}
	}
	
	if (aromatic_group) {
		group_names <- sub("( CD| CE)[12] (PHE|TYR)", "\\1* \\2", group_names)
	}
	
	if (methyl_group) {
		group_names <- sub("( CD)[12] (LEU)", "\\1* \\2", group_names)
		group_names <- sub("( CG)[12] (VAL)", "\\1* \\2", group_names)
	}
	
	groups_with_nonequivalent_protons <- character()
	
	if (!alpha_group) {
		groups_with_nonequivalent_protons <- c(" CA  GLY")
	}
	
	if (!beta_group) {
		groups_with_nonequivalent_protons <- c(groups_with_nonequivalent_protons, " CB  ARG", " CB  ASN", " CB  ASP", " CB  CYS", " CB  GLN", " CB  GLU", " CB  HIS", " CB  LEU", " CB  LYS", " CB  MET", " CB  PHE", " CB  PRO", " CB  SER", " CB  TRP", " CB  TYR")
	}
	
	if (!gamma_group) {
		groups_with_nonequivalent_protons <- c(groups_with_nonequivalent_protons, " CG  ARG", " CG  GLN", " CG  GLU", " CG1 ILE", " CG  LYS", " CG  MET", " CG  PRO")
	}
	
	if (!delta_group) {
		groups_with_nonequivalent_protons <- c(groups_with_nonequivalent_protons, " CD  ARG", " CD  LYS", " CD  PRO")
	}
	
	if (!epsilon_group) {
		groups_with_nonequivalent_protons <- c(groups_with_nonequivalent_protons, " CE  LYS")
	}
	
	if (!amine_group) {
		groups_with_nonequivalent_protons <- c(groups_with_nonequivalent_protons, " ND2 ASN", " NE2 GLN")
	}
	
	nonequivalent_idx <- substr(group_names, 1, 8) %in% groups_with_nonequivalent_protons
	group_names[nonequivalent_idx] <- names(group_names)[nonequivalent_idx]
	
	group_names
}

#' Calculate bond distances between a set of atoms
#'
#' @param coord_mat 3xN matrix with a column for every atom
#' @param depth maximum bond distance to calculate
#'
#' @return NxN matrix with bond distances between atoms (or NA if not calculated)
#'
#' @export
bond_separation <- function(coord_mat, depth=4) {

	bond_sep_mat <- matrix(NA_integer_, nrow=ncol(coord_mat), ncol=ncol(coord_mat), dimnames=list(colnames(coord_mat), colnames(coord_mat)))
	diag(bond_sep_mat) <- 0

	heavy_idx <- substr(colnames(coord_mat), 2, 2) %in% c("C","N","O","S")
	
	atom_dist <- as.matrix(stats::dist(t(coord_mat)))

	heavy_names <- colnames(coord_mat)[heavy_idx][apply(atom_dist[!heavy_idx,heavy_idx], 1, which.min)]
	proton_names <- colnames(coord_mat)[!heavy_idx]

	bond_sep_mat[cbind(heavy_names,proton_names)] <- bond_sep_mat[cbind(proton_names,heavy_names)] <- 1

	heavy_bond_mat <- which(atom_dist[heavy_idx,heavy_idx] > 0 & atom_dist[heavy_idx,heavy_idx] < 2, arr.ind=TRUE)
	heavy_bond_mat <- matrix(colnames(coord_mat)[heavy_idx][heavy_bond_mat], ncol=2)

	bond_sep_mat[heavy_bond_mat] <- 1

	for (i in seq_len(depth-1)) {
		for (j in seq_len(nrow(bond_sep_mat))) {
			bond_idx <- which(bond_sep_mat[j,] == i)
			onebond_idx <- unique(which(bond_sep_mat[,bond_idx,drop=FALSE] == 1, arr.ind=TRUE)[,1])
			nbond_idx <- onebond_idx[is.na(bond_sep_mat[j,onebond_idx])]
			bond_sep_mat[j,nbond_idx] <- i+1
		}
	}

	bond_sep_mat
}

#' Collapse first two dimensions of an array
#'
#' @param arr two or three dimensional array with the first two dimensions having the same length
#' @param idx_list list of character or integer indices to accumulate using func
#' @param func function to apply
#'
#' @export
collapse_array <- function(arr, idx_list, func=`+`) {

	idx_length <- sapply(idx_list, length)
	
	source_idx <- lapply(seq_len(max(idx_length)), function(idx) sapply(idx_list, function(x) x[idx]))
	dest_idx <- lapply(seq_along(source_idx), function(idx) which(!is.na(source_idx[[idx]])))
	source_idx <- lapply(seq_along(source_idx), function(idx) source_idx[[idx]][dest_idx[[idx]]])
	
	stopifnot(length(dim(arr)) %in% c(2L, 3L))
	
	if (length(dim(arr)) == 2) {
		
		arr_rows_collapsed <- arr[source_idx[[1]],,drop=FALSE]
		for (i in seq_along(source_idx)[-1]) {
			arr_rows_collapsed[dest_idx[[i]],] <- func(arr_rows_collapsed[dest_idx[[i]],],arr[source_idx[[i]],])
		}
	
		arr_rows_cols_collapsed <- arr_rows_collapsed[,source_idx[[1]],drop=FALSE]
		for (i in seq_along(source_idx)[-1]) {
			arr_rows_cols_collapsed[,dest_idx[[i]]] <- func(arr_rows_cols_collapsed[,dest_idx[[i]]],arr_rows_collapsed[,source_idx[[i]]])
		}
	
	} else if (length(dim(arr)) == 3) {
		
		arr_rows_collapsed <- arr[source_idx[[1]],,,drop=FALSE]
		for (i in seq_along(source_idx)[-1]) {
			arr_rows_collapsed[dest_idx[[i]],,] <- func(arr_rows_collapsed[dest_idx[[i]],,],arr[source_idx[[i]],,])
		}
	
		arr_rows_cols_collapsed <- arr_rows_collapsed[,source_idx[[1]],,drop=FALSE]
		for (i in seq_along(source_idx)[-1]) {
			arr_rows_cols_collapsed[,dest_idx[[i]],] <- func(arr_rows_cols_collapsed[,dest_idx[[i]],],arr_rows_collapsed[,source_idx[[i]],])
		}
	}
	
	dimnames(arr_rows_cols_collapsed)[1:2] <- list(names(idx_list), names(idx_list))
	
	arr_rows_cols_collapsed
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

	result <- func(value, ..., gradient=TRUE)
	
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

	r_array <- aperm(coord_array[atom_pairs[,2],,,drop=FALSE]-coord_array[atom_pairs[,1],,,drop=FALSE], c(1,3,2))
	
	if (!is.null(dimnames(coord_array))) {
		atom_names <- sub(" *([^ ]+) *([A-Z]{3}) (.) +([^ ]+) *", "\\4:\\1", dimnames(coord_array)[[1]])
		names(atom_names) <- dimnames(coord_array)[[1]]
		dimnames(r_array)[[1]] <- paste(atom_names[atom_pairs[,1]], atom_names[atom_pairs[,2]], sep="-")
		dimnames(r_array)[[3]] <- c("x", "y", "z")
	}
	
	r_array
}

#' Back-propagate energy derivative from r array to coordinates
#'
#' @param d_energy_d_coord_array 3D array (atoms, xyz, models) to accumulate derivatives into
#' @param atom_pairs matrix with each row having the names or indices of an atom pair (first dimension in `coord_array`
#' @param d_energy_d_r_array 3D array (pairs, models, xyz)
#'
#' @return 3D array (pairs, xyz, models) with d_energy_d_coord_array
coord_array_to_r_array_backprop <- function(d_energy_d_coord_array, atom_pairs, d_energy_d_r_array) {

	#print(str(d_energy_d_coord_array))
	#print(str(d_energy_d_r_array))
	# propagate the internuclear vector derivatives back onto the atomic coordinates
	for (i in seq_len(dim(d_energy_d_r_array)[1])) {
		
		pair_grad <- t(d_energy_d_r_array[i,,])
		# first component is subtracted so derivative is too
		d_energy_d_coord_array[atom_pairs[i,1],,] <- d_energy_d_coord_array[atom_pairs[i,1],,] - pair_grad
		# second component is added so derivative is too
		d_energy_d_coord_array[atom_pairs[i,2],,] <- d_energy_d_coord_array[atom_pairs[i,2],,] + pair_grad
	}
	
	d_energy_d_coord_array
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

#' Back-propagate energy derivative from d array to r array
#'
#' @param d_d_array_d_r_array 4D array (pairs, models, tensor elements, xyz) from
#'    `gradient` attribute of `r_array_to_d_array()`
#' @param d_energy_d_d_array 3D array (pairs, models, tensor elements)
#'
#' @return 3D array (pairs, models, xyz) with d_energy_d_r_array
r_array_to_d_array_backprop <- function(d_d_array_d_r_array, d_energy_d_d_array) {

	# calculate de/dr = dd/dr * de/dd for each xyz component of the internuclear vectors
	# coerce d_energy_d_d_array to a vector so it gets repeated over dd/dr
	d_energy_d_r_array_all <- d_d_array_d_r_array*as.vector(d_energy_d_d_array)
	# sum the individual interaction tensor component derivatives associated with x, y, and z
	apply(d_energy_d_r_array_all, c(1,2,4), sum)
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
	
	if (is.integer(grouping)) {
		
		# convert integer vector grouping to list grouping if necessary
		grouping <- lapply(lapply(seq_len(max(grouping)), "==", grouping), which)
	}
	
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

#' Calculate group norm squared from dipole-dipole interaction tensors
#'
#' @param d_array 3D array (pairs, models, tensor elements) with interaction tensors
#' @param grouping_mat integer matrix (groups, models) giving groupings of models to average interaction tensors
#' @param gradient a logical value indicating whether to calculate the derivative
#'
#' @return vector with norm squared for each atom pair. 
#'
#' The optional derivative is contained in the `"gradient"` attribute. It is a 3D array (pairs, models, tensor elements).
#'
#' @export
d_array_to_g_matrix <- function(d_array, grouping_mat, gradient=FALSE) {

	g_list <- apply(grouping_mat, 1, d_array_to_g, d_array=d_array, gradient=gradient, simplify=FALSE)
	
	value <- simplify2array(g_list)
	
	if (!is.matrix(value)) {
		value <- matrix(value, nrow=dim(d_array)[1])
	}
	
	if (gradient) {
	
		attr(value, "gradient") <- simplify2array(lapply(g_list, attr, "gradient"))
	}
	
	value
}

#' Back-propagate energy derivative from g matrix to d array
#'
#' @param d_g_matrix_d_d_array 4D array (pairs, models, tensor components, groups) fom
#'    `gradient` attribute of `d_array_to_g_matrix()`
#' @param d_energy_d_g_matrix matrix (pairs, groups)
#'
#' @return 3D array (pairs, models, tensor elements) with d_energy_d_d_array
d_array_to_g_matrix_backprop <- function(d_g_matrix_d_d_array, d_energy_d_g_matrix) {

	# calculate de/dd = dg/dd * de/dg for all individual interaction tensor components
	d_energy_d_d_array <- 0
	# sum the contributions from the different norm squared values
	for (i in seq_len(dim(d_g_matrix_d_d_array)[4])) {
		d_energy_d_d_array <- d_energy_d_d_array + d_g_matrix_d_d_array[,,,i]*d_energy_d_g_matrix[,i]
	}
	
	d_energy_d_d_array
}

#' Calculate matrix of a values from a matrix of g values
#'
#' @param g_matrix matrix of g values with columns associated with different groupings
#' @param a_coef coefficients used for calculating a values
g_matrix_to_a_matrix <- function(g_matrix, a_coef) {

	a_matrix <- matrix(0, nrow=nrow(g_matrix), ncol=ncol(a_coef), dimnames=list(rownames(g_matrix), colnames(a_coef)))
	
	for (i in seq_len(ncol(a_coef))) {
		for (j in which(a_coef[,i] != 0)) {
			a_matrix[,i] <- a_matrix[,i] + a_coef[j,i]*g_matrix[,j]
		}
	}
	
	a_matrix
}

#' Back-propagate energy derivative from a matrix to g matrix
#'
#' @param a_coef coefficients used for calculating a values
#' @param d_energy_d_a_matrix matrix (pairs, eigenvalues)
#'
#' @return matrix (pairs, groupings) with d_energy_d_g_matrix
g_matrix_to_a_matrix_backprop <- function(a_coef, d_energy_d_a_matrix) {

	d_energy_d_g_matrix <- matrix(0, nrow=nrow(d_energy_d_a_matrix), ncol=nrow(a_coef), dimnames=list(rownames(d_energy_d_a_matrix), NULL))

	for (i in seq_len(nrow(a_coef))) {
		for (j in which(a_coef[i,] != 0)) {
			d_energy_d_g_matrix[,i] <- d_energy_d_g_matrix[,i] + d_energy_d_a_matrix[,j]/a_coef[i,j]
		}
	}
	
	d_energy_d_g_matrix
}

#' Calculate sigma from a matrix of a values
#'
#' @param a_matrix matrix of a values with columns associated with eigenvalues
#' @param lambda_prime_vec eigenvalues augmented with tumbling rate
#' @param proton_mhz spectrometer proton field strength in MHz
#' @param gradient a logical value indicating whether to calculate the derivative
a_matrix_to_sigma <- function(a_matrix, lambda_prime_vec, proton_mhz, gradient=FALSE) {

	# K^2
	# ((1.2566370614e-6 kilogram meters / (ampere^2 second^2))/(4*pi)*(1.054571726e-34 meter^2 kilograms / second)*(267.513e6 radian s^-1 T^-1)^2)^2
	# 5.69549944e-49 rad m^6 / s^2
	K_sq <- 5.69549944e-49
	omega0 <- proton_mhz*1e6*2*pi

	a_matrix <- t(a_matrix)
	
	J0 <- -2*colSums(a_matrix/lambda_prime_vec)
	#Jomega <- -2*colSums(a_matrix*lambda_prime_vec/(lambda_prime_vec^2 + (omega0)^2))
	J2omega <- -2*colSums(a_matrix*lambda_prime_vec/(lambda_prime_vec^2 + (2*omega0)^2))
	
	value <- 0.1*K_sq*(3*J2omega - 0.5*J0)
	
	if (gradient) {
	
		d_J0_d_a_matrix <- -2/lambda_prime_vec
		d_J2omega_d_a_matrix <- -2*lambda_prime_vec/(lambda_prime_vec^2 + (2*omega0)^2)
		
		d_sigma_d_a_matrix <- 0.1*K_sq*(3*d_J2omega_d_a_matrix - 0.5*d_J0_d_a_matrix)
		
		attr(value, "gradient") <- matrix(d_sigma_d_a_matrix, nrow=ncol(a_matrix), ncol=nrow(a_matrix), byrow=TRUE)
	}
	
	value
}

#' Back-propagate energy derivative from sigma to a matrix
#'
#' @param d_sigma_d_a_matrix matrix (pairs, eigenvalues) fom `gradient` attribute of
#'    `a_matrix_to_sigma()`
#' @param d_energy_d_sigma vector (pairs)
#'
#' @return matrix (pairs, eigenvalues) with d_energy_d_a_matrix
a_matrix_to_sigma_backprop <- function(d_sigma_d_a_matrix, d_energy_d_sigma) {

	d_sigma_d_a_matrix*d_energy_d_sigma
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
g_to_energy <- function(g, g0, k=1, gradient=FALSE) {

	attr(g, "gradient") <- NULL
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

#' Scale real or back calculated values towards zero by taking a power
#'
#' @param x numerical values to be rescaled
#' @param p numeric power to raise the values to, usually 1 or less
#'
#' @return \deqn{ ( ( |x| + 1 )^p - 1 ) \operatorname{sgn}(x) }
power_scale <- function(x, p) {
	( ( abs(x) + 1 )^p - 1 ) * sign(x)
}

#' Traditional squared loss function scaled by a power
#'
#' @param x current values
#' @param x0 target values
#' @param p numeric power to raise the values to, usually 1 or less
#' @param k force constant
#' @param gradient a logical value indicating whether to calculate the derivative
#'
#' @return \deqn{ k \left ( ( ( |x| + 1 )^p - 1 ) \operatorname{sgn}(x) - ( ( |x_0| + 1 )^p - 1 ) \operatorname{sgn}(x_0) \right )^2 }
#'
#' @examples
#' x <- seq(-10, 10, by=0.1)
#' loss_x <- sapply(x, power_scaled_loss, 2, 0.25)
#' loss_x_grad <- power_scaled_loss(x, 2, 0.25, gradient=TRUE)
#' par(mfrow=c(2, 1))
#' plot(x, loss_x, type="l", ylab="loss")
#' plot(x, attr(loss_x_grad, "gradient"), type="l", ylab="dloss/dx")
#' points(x[-1]-mean(diff(x))/2, diff(loss_x)/mean(diff(x)), type="l", col="blue")
#' abline(h = 0, col="gray")
#' legend("bottomright", legend=c("Analytical", "Finite Difference"), bty="n", lwd=1, col=c("black", "blue"))
#'
#' @export
power_scaled_loss <- function(x, x0, p=1, k=1, gradient=FALSE) {
	
	expr1 <- ( ( abs(x) + 1 )^p - 1 ) * sign(x) - ( ( abs(x0) + 1 )^p - 1 ) * sign(x0)
	value <- sum(k * expr1^2)
	if (gradient) {
		# D[(((Abs[x] + 1)^p - 1)*Sign[x] - ((Abs[x0] + 1)^p - 1)*Sign[x0])^2, x]
		grad <- k * (2 * expr1) * p * ( abs(x) + 1 )^(p-1) #* sign(x) * sign(x)
		attr(value, "gradient") <- grad
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
coord_array_to_g_matrix <- function(coord_array, atom_pairs, grouping_list) {

	r_array <- coord_array_to_r_array(coord_array, atom_pairs)
	
	d_array <- r_array_to_d_array(r_array, gradient=FALSE)
	
	g_list <- lapply(grouping_list, function(grouping) d_array_to_g(d_array, grouping, gradient=FALSE))
	
	simplify2array(g_list)
}

#' Calculate g value restraint energy from atomic coordinates
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
coord_array_to_g_energy <- function(coord_array, atom_pairs, grouping_list, g0, k, gradient=FALSE) {

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

#' Calculate g value restraint energy from atomic coordinates
#'
#' @param coord_array 3D array (atoms, xyz, models) with atomic coordinates
#' @param atom_pairs matrix with each row having the names or indices of an atom pair (first dimension in `coord_array`)
#' @param grouping_mat integer matrix (groups, models) giving groupings of models to average interaction tensors
#' @param g0 target group norm squared values
#' @param loss_func loss function to use to calculate energy
#' @param ... additional parameters passed to `loss_func`
#' @param gradient a logical value indicating whether to calculate the derivative
#'
#' @return total restraint energy calculated using `loss_func`
#'
#' The optional derivative is contained in the `"gradient"` attribute. It is a 3D array (atoms, xyz, models).
#'
#' @export
coord_array_to_g_energy_refactored <- function(coord_array, atom_pairs, grouping_mat, g0, loss_func = power_scaled_loss, ..., gradient=FALSE) {

	# calculate internuclear vectors
	r_array <- coord_array_to_r_array(coord_array, atom_pairs)
	
	# calculate dipole-dipole interaction tensors
	d_array <- r_array_to_d_array(r_array, gradient=gradient)
	
	# calculate norm squared for different groupings of dipole-dipole interaction tensors
	g_matrix <- d_array_to_g_matrix(d_array, grouping_mat, gradient=gradient)
	
	# calculate energies from the norm squared values
	energy_matrix <- loss_func(g_matrix, g0, ..., gradient=gradient)
	
	# return the sum of all the individual restraint energies
	value <- sum(energy_matrix)
	
	if (gradient) {
	
		# initialize empty gradient with dimensions equal to input coordinates
		d_energy_d_coord_array <- coord_array
		d_energy_d_coord_array[] <- 0

		d_energy_d_g_matrix <- attr(energy_matrix, "gradient")

		# back-propagate derivatives from g matrix to d array
		d_energy_d_d_array <- d_array_to_g_matrix_backprop(attr(g_matrix, "gradient"), d_energy_d_g_matrix)
		
		# back-propagate derivatives from d array to r array
		d_energy_d_r_array <- r_array_to_d_array_backprop(attr(d_array, "gradient"), d_energy_d_d_array)
		
		# accumulate back-propagated derivatives from r array to coord array gradient
		d_energy_d_coord_array <- coord_array_to_r_array_backprop(d_energy_d_coord_array, atom_pairs, d_energy_d_r_array)
		
		# set coordinate gradient to accumulated values
		attr(value, "gradient") <- d_energy_d_coord_array
	}
	
	value
}

#' Read data for calculating spectral density functions
#'
#' @param prefix_path to prefix of four CSV files
#'
#' @return a list with elements: `atom_pairs`, `groupings`, `a_coef`, and `lambda_coef`
#'
#' @export
read_spec_den_data <- function(prefix_path) {

	atom_pairs <- utils::read.csv(paste0(prefix_path, "_atom_pairs.csv"))
	groupings <- unname(as.matrix(utils::read.csv(paste0(prefix_path, "_groupings.csv"), header=FALSE, row.names=NULL)))
	a_coef <- as.matrix(utils::read.csv(paste0(prefix_path, "_a_coef.csv"), check.names=FALSE))
	lambda_coef <- as.matrix(utils::read.csv(paste0(prefix_path, "_lambda_coef.csv"), check.names=FALSE, row.names=1))
	
	list(
		atom_pairs = atom_pairs,
		groupings = groupings,
		a_coef = a_coef,
		lambda_coef = lambda_coef
	)
}

#' Shift from one array dimension to another
#'
#' @param a array whose dimensions should be shifted
#' @param n integer factor to divide `dfrom` dimension and multiply `dto` dimension
#' @param dnames list with dimension names as alternative to specifying `n`
#' @param dfrom integer with dimension to move from
#' @param dto integer with dimension to move to (should be one greater than `dfrom`)
array_shift <- function(a, n=NULL, dnames=NULL, dfrom=1, dto=2) {

	if (!is.null(n)) {
	
		if (n == 1) {
			return(a)
		}
	
		new_a <- a
	
		new_dim <- dim(a)
		new_dim[dfrom] <- dim(a)[dfrom]/n
		new_dim[dto] <- dim(a)[dto]*n
		dim(new_a) <- new_dim
		
		new_dimnames <- dimnames(a)
		if (!is.null(dimnames(a)[[dfrom]])) {
			new_dimnames[[dfrom]] <- dimnames(a)[[dfrom]][seq_len(length(dimnames(a)[[dfrom]])/n)]
		}
		if (!is.null(dimnames(a)[[dto]])) {
			new_dimnames[[dto]] <- paste(rep(dimnames(a)[[dto]], each=n), seq_len(n), sep="-")
		}
		dimnames(new_a) <- new_dimnames
	
	} else if (!is.null(dnames)) {
	
		new_a <- a
		dim(new_a) <- sapply(dnames, length)
		dimnames(new_a) <- dnames
	
	} else {
	
		stop()
	}
	
	new_a
}

#' Calculate cross relaxation rates from atomic coordinates
#'
#' @param coord_array 3D array (atoms, xyz, models) with atomic coordinates
#' @param rates named numeric vector with rates
#' @param spec_den_data_list list of data for calculating spectral density functions
#' @param proton_mhz spectrometer proton field strength in MHz
#'
#' @export
coord_array_to_sigma <- function(coord_array, rates, spec_den_data_list, proton_mhz) {

	lapply(spec_den_data_list, function(spec_den_data) {
	
		# calculate internuclear vectors (convert from Å^-3 to m^-3)
		r_array <- coord_array_to_r_array(coord_array*1e-10, spec_den_data[["atom_pairs"]][,1:2,drop=FALSE])
	
		# calculate dipole-dipole interaction tensors
		d_array <- r_array_to_d_array(r_array)
	
		# calculate the factor by which the number of models should be expanded
		n_shift <- ncol(spec_den_data[["groupings"]])/dim(coord_array)[3]
	
		# shift tensor components from atom pairs into virtual models
		d_array_shifted <- array_shift(d_array, n_shift)
		
		# calculate matrix of g values
		g_mat <- apply(spec_den_data[["groupings"]], 1, d_array_to_g, d_array=d_array_shifted)
		
		# calculate matrix of a values
		a_mat <- g_matrix_to_a_matrix(g_mat, spec_den_data[["a_coef"]])
		
		# calculate lambda eigenvalues
		lambda_vec <- -colSums(rates[rownames(spec_den_data[["lambda_coef"]])]*spec_den_data[["lambda_coef"]])
		
		# update eigenvalues to account for molecular tumbling
		lambda_prime_vec <- lambda_vec - rates["kc"]
		
		a_matrix_to_sigma(a_mat, lambda_prime_vec, proton_mhz)
	})
}

#' Calculate sigma restraint energy from atomic coordinates
#'
#' @param coord_array 3D array (atoms, xyz, models) with atomic coordinates
#' @param rates named numeric vector with rates
#' @param spec_den_data_list list of data for calculating spectral density functions
#' @param proton_mhz spectrometer proton field strength in MHz
#' @param loss_func loss function to use
#' @param ... additional parameters passed to `loss_func`
#' @param gradient a logical value indicating whether to calculate the derivative
#'
#' @return total restraint energy calculated using `loss_func`
#'
#' The optional derivative is contained in the `"gradient"` attribute. It is a 3D array 
#' (atoms, xyz, models).
#'
#' Testing with `deriv_check` showed a slight systematic (~0.2% underestimation) of the 
#' gradient between two methyl groups. Perhaps there's some unaccounted correlation in 
#' the purely additive derivative calculation?
#'
#' @export
coord_array_to_sigma_energy <- function(coord_array, rates, spec_den_data_list, proton_mhz, loss_func = power_scaled_loss, ..., gradient=FALSE) {

	# intermediate derivatives that need to be captured for back-propagation
	d_array_list <- g_matrix_list <- sigma_vector_list <- vector("list", length(spec_den_data_list))
	
	for (i in seq_along(spec_den_data_list)) {

		# calculate internuclear vectors
		r_array <- coord_array_to_r_array(coord_array*1e-10, spec_den_data_list[[i]][["atom_pairs"]][,1:2,drop=FALSE])
	
		# calculate dipole-dipole interaction tensors
		d_array_list[[i]] <- r_array_to_d_array(r_array, gradient=gradient)
	
		# calculate the factor by which the number of models should be expanded
		n_shift <- ncol(spec_den_data_list[[i]][["groupings"]])/dim(coord_array)[3]
	
		# shift tensor components from atom pairs into virtual models
		d_array_shifted <- array_shift(d_array_list[[i]], n_shift)
	
		# calculate norm squared for different groupings of dipole-dipole interaction tensors
		g_matrix_list[[i]] <- d_array_to_g_matrix(d_array_shifted, spec_den_data_list[[i]][["groupings"]], gradient=gradient)
		
		# calculate matrix of a values
		a_matrix <- g_matrix_to_a_matrix(g_matrix_list[[i]], spec_den_data_list[[i]][["a_coef"]])
		
		# calculate lambda eigenvalues
		lambda_vector <- -colSums(rates[rownames(spec_den_data_list[[i]][["lambda_coef"]])]*spec_den_data_list[[i]][["lambda_coef"]])
		
		# update eigenvalues to account for molecular tumbling
		lambda_prime_vector <- lambda_vector - rates["kc"]
		
		# calculate cross relaxation rates
		sigma_vector_list[[i]] <- a_matrix_to_sigma(a_matrix, lambda_prime_vector, proton_mhz, gradient=TRUE)
	}
	
	# combine sigmas into a single vector
	sigma_vector <- unlist(sigma_vector_list)

	# combine target sigmas into a single vector
	sigma0_vector <- unlist(lapply(spec_den_data_list, function(x) {
		idx <- seq_len(nrow(x$atom_pairs)/(ncol(x$groupings)/dim(coord_array)[[3]]))
		x$atom_pairs[idx,"sigma"]
	}))
	
	# calculate energies from the sigma vectors
	energy_vector <- loss_func(sigma_vector, sigma0_vector, ..., gradient=gradient)
	
	# return the sum of all the individual restraint energies
	value <- sum(energy_vector)
	
	if (gradient) {
		
		# calculate starting offsets within sigma gradients
		sigma_vector_list_lengths <- sapply(sigma_vector_list, length)
		start_idx <- cumsum(c(1, sigma_vector_list_lengths[-length(sigma_vector_list)]))
	
		# initialize empty gradient with dimensions equal to input coordinates
		d_energy_d_coord_array <- coord_array
		d_energy_d_coord_array[] <- 0

		for (i in seq_along(spec_den_data_list)) {

			# extract particular segment of energy gradient
			d_energy_d_sigma <- attr(energy_vector, "gradient")[seq(start_idx[i], length.out=length(sigma_vector_list[[i]]))]
			
			# back-propagate derivatives from sigma to a matrix
			d_energy_d_a_matrix <- a_matrix_to_sigma_backprop(attr(sigma_vector_list[[i]], "gradient"), d_energy_d_sigma)
			
			# back-propagate derivatives from a array to g matrix
			d_energy_d_g_matrix <- g_matrix_to_a_matrix_backprop(spec_den_data_list[[i]][["a_coef"]], d_energy_d_a_matrix)
			
			# back-propagate derivatives from g matrix to d array
			d_energy_d_d_array <- d_array_to_g_matrix_backprop(attr(g_matrix_list[[i]], "gradient"), d_energy_d_g_matrix)
			
			# because of the way the arrays are structured shifting back isn't necessary
			# calculate the factor by which the number of models should be expanded
			#n_shift <- ncol(spec_den_data_list[[i]][["groupings"]])/dim(coord_array)[3]
			
			# shift tensor components from virtual models back to atom pairs
			#if (n_shift > 1) {
			#	d_energy_d_d_array <- array_shift(d_energy_d_d_array, dnames=dimnames(d_array_list[[i]]))
			#}
		
			# back-propagate derivatives from d array to r array
			d_energy_d_r_array <- r_array_to_d_array_backprop(attr(d_array_list[[i]], "gradient"), d_energy_d_d_array)
		
			# accumulate back-propagated derivatives from r array into coord array gradient
			d_energy_d_coord_array <- coord_array_to_r_array_backprop(d_energy_d_coord_array, spec_den_data_list[[i]][["atom_pairs"]][,1:2,drop=FALSE], d_energy_d_r_array)
		}
		
		# set coordinate gradient to accumulated values
		attr(value, "gradient") <- d_energy_d_coord_array*1e-10
	}
	
	value
}
